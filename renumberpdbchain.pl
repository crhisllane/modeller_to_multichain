#!/usr/bin/perl   -wT

# Copyright 2010 Emma Rath. All rights reserved.
# Soon this program will be released under an open source software license such as GNU General Public License or 
# Creative Commons license for Free Software Foundation's GNU General Public License at creativecommons.org

# perl -T renumberpdbchain.pl -infile 1ZOY.pdb -add 200 -chain D -from 230 -to 357

# This perl program runs as a as a cgi script on a web server that displays web pages.
# This perl program also runs as a batch program on the command line.

# This program reads in a PDB 3D file, and for the given chain (column 22 in the file),
# renumbers all the residues (residue number is columns 23-26) of that chain, 
# for lines starting with ATOM or HETATM,
# by adding the given number to the residue number.
# The given number can be a negative number (so that a number is subtracted from each residue number).

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use Math::Complex;
use Math::Trig;
use CGI;
use Data::Dumper;
# $Data::Dumper::Purity = 1;
# print "var x = " . Dumper( $x ) . "\n"; # debug
# my $wait_for_key_press = <STDIN>; # debug

my $q = new CGI;
my $global;
my $input;
my $debug = '';

#===================================================================================================
#=================================              MAIN              ==================================
#===================================================================================================

init();

get_program_mode();

if ($global->{'program_mode'} eq 'cgi-display-form') {

	display_html_form();

} elsif ($global->{'program_mode'} eq 'cgi-output') {

	get_html_input();
	if ($global->{'got_error'} == 0) {
		build_and_output_results_for_html_output();
	} else {
		display_html_form();
	}

} else { # $global->{'program_mode'} eq 'command-line'

	get_command_line_input();
	if ($global->{'got_error'} == 0) {

		my $untainted_infile = untaint( $input->{'infile'} );
		if ($untainted_infile eq '') {
			$untainted_infile = "renumberpdbchain";
		}
		my $output_file = $untainted_infile . ".pdb";
		open( OUTFILE, ">$output_file") or
			die "Cannot open $output_file for writing\n";
		build_and_output_results_for_command_line();
		close OUTFILE;
	} else {
		output_command_line_help();
	}
}



#===================================================================================================
#	this is the main logic subroutine of this program for command-line mode.
#	given an input PDB file, process the input and produce the output.
#===================================================================================================

sub build_and_output_results_for_command_line {

	build_the_output();

	my @input_lines = @{$input->{'input_lines'}};
	foreach my $line (@input_lines) {
		print OUTFILE "$line\n";
	}
}



#===================================================================================================
#	this is the main logic subroutine of this program for cgi-output mode.
#	given an input PDB file, process the input and produce the output.
#===================================================================================================

sub build_and_output_results_for_html_output {

	build_the_output();

	if ($input->{'whereto'} eq 'S') { # write output to the screen

	print $q->header();
	print	#$q->start_html(-title=>$text_title),
		"<html><head>\n",
		'<style type="text/css">' . "\n",
		'<!--' . "\n",
		'.tt {' . "\n",
		'	font-size: 14px;' . "\n",
		'	font-family: courier, monospace;' . "\n",
		'	color: #000000;' . "\n",
		'}' . "\n",
		'-->' . "\n",
		'</style>' . "\n",
		"</head><body>\n";
		print "<span class='tt'>\n";
		my @input_lines = @{$input->{'input_lines'}};
		foreach my $line (@input_lines) {
			$line =~ s/ /&nbsp;/g;
			$line = "$line<br>\n";
			print $line;
		}
		print "</span></body></html>\n";

	} else { # $input->{'whereto'} eq 'F' # write output to a file to be downloaded

		print $q->header("Content-type: text/plain");
		my @input_lines = @{$input->{'input_lines'}};
		foreach my $line (@input_lines) {
			print "$line\n";
		}

	}

	if ($debug ne '') {
		print "$debug\n";
	}
}



#===================================================================================================
#	move the input PDB file contents to output, changing the residue numbering
#===================================================================================================

sub build_the_output {

	# ATOM   4726  CE2 TYR A 622      85.569  38.124 138.803  1.00115.05           C  
	# ATOM   4727  CZ  TYR A 622      86.737  38.576 139.390  1.00115.85           C  
	# ATOM   4728  OH  TYR A 622      86.852  39.904 139.730  1.00116.74           O  
	# ATOM   4729  OXT TYR A 622      84.283  32.238 138.506  1.00111.29           O  
	# TER    4730      TYR A 622                                                      
	# ATOM   4731  N   PRO B   9     111.971  33.330  97.856  1.00 92.26           N  
	# ATOM   4732  CA  PRO B   9     112.160  34.798  97.760  1.00 92.31           C  
	# ATOM   4733  C   PRO B   9     111.127  35.407  96.816  1.00 93.29           C  
	# ATOM   4734  O   PRO B   9     110.966  34.956  95.677  1.00 92.89           O  
	# HETATM    8  C4' FAD A 700      91.658  54.405 126.549  1.00 52.15           C
	# HETATM    9  O4' FAD A 700      91.040  54.647 125.317  1.00 56.02           O
	# HETATM   10  C3' FAD A 700      92.030  52.882 126.685  1.00 54.68           C
	# CONECT    1    3    4    5    7
	# CONECT    2    5   18   23   39
	# CONECT    3    1
	# CONECT    4    1
	# CONECT    5    1    2

	my @input_lines = @{$input->{'input_lines'}};

	my $add = $input->{'add'};

	for ( my $i = 0; $i < @input_lines; $i++ ) {
		my $line = $input_lines[$i];
		my $l = $i + 1;
		$line = trim($line);
		if ($line ne '') {
			my $change_this_line = 0;
			if ($input->{'chain'} ne '') {
				if (length($line) >= 22) {
					my $this_chain = trim(substr($line,21,1));
					if ($this_chain eq $input->{'chain'}) {
						$change_this_line = 1;
					}
				}
			} else {
				$change_this_line = 1;
			}
			if ($change_this_line == 1) {
				if (length($line) >= 6) {
					my $line_id = substr($line,0,6);
					if (($line_id eq 'ATOM  ') or ($line_id eq 'HETATM') or ($line_id eq 'TER   ')) {
						if (($input->{'from'} ne '') or ($input->{'to'} ne '')) {
							my $res_num = trim(substr($line,22,4));
							if ($input->{'from'} ne '') {
								if ($res_num < $input->{'from'}) {
									$change_this_line = 0;
								}
							}
							if ($input->{'to'} ne '') {
								if ($res_num > $input->{'to'}) {
									$change_this_line = 0;
								}
							}
						}
					}
				}
			}
			if ($change_this_line == 1) {
				if (length($line) >= 6) {
					my $line_id = substr($line,0,6);
					if (($line_id eq 'ATOM  ') or ($line_id eq 'HETATM') or ($line_id eq 'TER   ')) {
						my $bit1 = substr($line,0,22);
						my $bit2 = substr($line,26);
						my $old_res_num = trim(substr($line,22,4));
						my $new_res_num = abs($old_res_num + $add);
						my $display_res_num = sprintf("%4s", $new_res_num);
						my $new_line = $bit1 . $display_res_num . $bit2;
						$input_lines[$i] = $new_line;
					}
				}
			}
		}
	}

	$input->{'input_lines'} = \@input_lines;
}



#===================================================================================================
#	this program is being called to display the HTML input form
#===================================================================================================

sub display_html_form {

	html_header();

	my $request_uri = $ENV{'REQUEST_URI'};
	my @bits = split(/\//, $request_uri);
	my $pgm_name = $bits[$#bits];

	print "<br><br>\n";

	print "	<form action='$pgm_name' method='post' enctype='multipart/form-data'>\n",
		"	<table cellspacing=0 cellpadding=2 border=0>\n";

	if ($global->{'got_error'} == 1) {
print "		<tr>\n",
"			<td valign='top' align='left'>\n",
"				<p><font color='red'>" . $global->{'err_msg'} . "<br><b>Didn't process because of input error(s).</b><br><br></font></p>\n",
"			</td>\n",
"		</tr>\n";
	}

print "		<tr>\n",
"			<td valign='top' align='left'>\n",
"				<table cellspacing=0 cellpadding=0 border=0>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"Please provide a PDB file with the chain id in column 22.<br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"Upload the PDB file from your computer.<br>Please note that this PDB file is not validated, so make sure that it is a valid PDB file.<br>\n",
							"<input type='file' name='upload_file'><br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"<table cellspacing=10 cellpadding=0 border=0>\n",
							"<tr><td width='20'></td><td valign='top' align='left'>or</td></tr>\n",
							"</table>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"<table cellspacing=10 cellpadding=0 border=0>\n",
							"<tr><td width='20'></td>	<td valign='top' align='left'>Enter the PDB file contents here.<br>Please note that these PDB file contents are not validated, so make sure that they are valid PDB file contents.</td></tr>\n",
							"<tr><td></td>			<td valign='top' align='left'><textarea name='pdb_contents' rows='5' cols='50'></textarea></td></tr>\n",
							"<tr><td></td>			<td valign='top' align='left'>(Provide either a file or enter the contents here, not both.<br>If both provided, then contents here will be processed instead of uploaded file.)</td></tr>\n",
							"</table><br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"Add this number to each residue number in the chain : ",
							"<input type='text' name='add' size=3 value=''>",
							" (mandatory field, can be negative)<br><br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"Chain whose residues are to be renumbered : ",
							"<input type='text' name='chain' size=3 value=''>",
							" (optional field, all residues residues will be renumbered if left blank)<br><br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"Renumber residues from residue number : ",
							"<input type='text' name='from' size=3 value=''>",
							" (optional field)<br><br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"Renumber residues to residue number : ",
							"<input type='text' name='to' size=3 value=''>",
							" (optional field)<br><br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"<input type='radio' name='whereto' value='S' checked>Display results on screen<br>\n",
							"<input type='radio' name='whereto' value='F'>Write results to output file<br><br>\n",
"						</td>\n",
"					</tr>\n",

"					<tr><td height='15'></td></tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"<input type='submit' name='submit' value='Upload'>\n",
"						</td>\n",
"					</tr>\n",

"				</table>\n",
"			</td>\n",
"		</tr>\n";

	print "	</table>\n",
		"</form>\n";

	html_footer();
}



#===================================================================================================
#	in command line mode - get the command line options, read the input file.
#===================================================================================================

sub get_command_line_input {

	my $input_file;
	my $input_add;
	my $input_chain;
	my $input_from;
	my $input_to;
	my $input_h;
	my $input_help;

	GetOptions( "infile=s" => \$input_file, "add=s" => \$input_add, "chain=s" => \$input_chain, "from=s" => \$input_from, "to=s" => \$input_to, "h" => \$input_h, "help" => \$input_help );
	if (!defined $input_file) {
		die "Usage $0 -infile INPUTFILE -add ADDNUMBER -chain CHAINID -from FROMNUMBER -to TONUMBER\n";
	}

	if (defined $input_file) {
		$input->{'infile'} = $input_file;
	}
	if (defined $input_add) {
		$input->{'add'} = $input_add;
	}
	if (defined $input_chain) {
		$input->{'chain'} = $input_chain;
	}
	if (defined $input_from) {
		$input->{'from'} = $input_from;
	}
	if (defined $input_to) {
		$input->{'to'} = $input_to;
	}
	if (defined $input_h) {
		$input->{'h'} = $input_h;
	}
	if (defined $input_help) {
		$input->{'help'} = $input_help;
	}
	if ((defined $input_h) || (defined $input_help)) {
		$global->{'got_error'} = 1;
	}

	open INFILE, $input->{'infile'} or die $!;
	my @input_lines = <INFILE>;
	my @new_input_lines;
	foreach my $line (@input_lines) {
		chomp($line);
		push(@new_input_lines, $line);
	}
	$input->{'input_lines'} = \@new_input_lines;
}



#===================================================================================================
#	get the input from the html form, so can process the input
#===================================================================================================

sub get_html_input {

	my $params = $q->Vars;

	if (defined($params->{'add'})) {
		if ($params->{'add'} ne '') {
			my $input_add = trim( $params->{'add'} );
			if (is_valid_add_number($input_add) == 1) {
				$input->{'add'} = $input_add;
			} else {
				$global->{'got_error'} = 1;
				$global->{'err_msg'} .= $global->{'err_msg'} . "The number to add that you entered is not a valid integer.<br>\n";
			}
		}
	}

	if (defined($params->{'chain'})) {
		my $input_chain = trim( $params->{'chain'} );
		if (is_valid_chain($input_chain) == 1) {
			$input->{'chain'} = $input_chain;
		} else {
			$global->{'got_error'} = 1;
			$global->{'err_msg'} .= $global->{'err_msg'} . "The chain that you entered is invalid. It must only be 1 character long.<br>\n";
		}
	}

	if (defined($params->{'from'})) {
		my $input_from = trim( $params->{'from'} );
		if (is_valid_residue_number($input_from) == 1) {
			$input->{'from'} = $input_from;
		} else {
			$global->{'got_error'} = 1;
			$global->{'err_msg'} .= $global->{'err_msg'} . "The from-residue-number that you entered is invalid. It must be a positive integer.<br>\n";
		}
	}

	if (defined($params->{'to'})) {
		my $input_to = trim( $params->{'to'} );
		if (is_valid_residue_number($input_to) == 1) {
			$input->{'to'} = $input_to;
		} else {
			$global->{'got_error'} = 1;
			$global->{'err_msg'} .= $global->{'err_msg'} . "The to-residue-number that you entered is invalid. It must be a positive integer.<br>\n";
		}
	}

	if (defined($params->{'whereto'})) {
		if ($params->{'whereto'} ne '') {
			my $input_whereto = trim( $params->{'whereto'} );
			$input->{'whereto'} = uc($input_whereto);
			if (($input->{'whereto'} ne 'S') and ($input->{'whereto'} ne 'F')) {
				$input->{'whereto'} = 'S';
			}
		}
	}

	my @input_lines;
	$input->{'input_lines'} = \@input_lines;

	my $user_provided_pdb_contents = 0;
	if (defined($params->{'pdb_contents'})) {
		my $pdb_contents = trim( $params->{'pdb_contents'} );
		if ($pdb_contents ne '') {
			$user_provided_pdb_contents = 1;
			my @pdb_lines = split(/\n|\r/, $pdb_contents);
			foreach my $pdb_line (@pdb_lines) {
				$pdb_line = trim($pdb_line);
				if ($pdb_line ne '') {
					push( @input_lines, $pdb_line );
				}
			}
		}
	}

	if ( (defined($params->{'upload_file'})) && ($user_provided_pdb_contents == 0)) {
		my $upload_file = $q->param('upload_file');
		my $upload_filehandle = $q->upload('upload_file');
		if ($upload_file ne '') {
			if ($upload_filehandle eq '') {
				$global->{'got_error'} = 1;
				$global->{'err_msg'} .= "The upload didn't work. Probably your file is too big - bigger than 5MB.<br>";
			} else {

				# read the uploaded file

				my @uploaded_lines;

				my $user_data = '';
				while ( <$upload_filehandle> ) {
					my $input_line = $_;
					chomp($input_line);
					push( @input_lines, $input_line );
				}
			}
		}
	}

	if ($global->{'got_error'} == 0) {
		my $num_input_lines = $#{$input->{'input_lines'}};
		if ($num_input_lines < 1) {
			$global->{'got_error'} = 1;
			$global->{'err_msg'} .= "No input PDB lines were given. Please enter your PDB file contents or upload a PDB file.<br>";
		}
	}
}



#===================================================================================================
#	figure out whether this program is being called to display the HTML input form,
#	or whether the form has been filled and now the results need to be displayed,
#	or whether this program is being run in batch mode from the command line.
#===================================================================================================

sub get_program_mode {

	$global->{'program_mode'} = 'command-line';
	if (defined($ENV{'REQUEST_METHOD'})) {
		if ($ENV{'REQUEST_METHOD'} eq 'POST') {
			$global->{'program_mode'} = 'cgi-output';
		} else {
			$global->{'program_mode'} = 'cgi-display-form';
		}
	} else {
		$global->{'program_mode'} = 'command-line';
	}
}



#===================================================================================================
#	html code for the graphics part of the header of the html page
#===================================================================================================

sub html_code_for_graphics_header {

	my $title = shift;
	my $html_code = '';

	$html_code .= "<table width='100%' cellspacing='0' cellpadding='0' border='0'><tr><td align='left'>\n" .
		"<table height='160' cellspacing='0' cellpadding='0' border='0' align='left' bgcolor='#a4a5e9'>\n" .
		#"	style='background-image: url(\"bluegold.gif\");\n" .
		#"		background-attachment: fixed;\n" .
		#"		background-position: top right;\n" .
		#"		background-repeat: no-repeat;'>\n" .
		"	<tr>\n" .
		"		<td height='160' valign='center'>\n" .
		"			<table height='160' cellspacing='0' cellpadding='0' border='0'>\n" .
		"				<tr>\n" .
		"					<td width='20'></td>\n" .
		"					<td width='740' valign='center'>\n" .
		"						<font size='5' color='black'>$title</font>\n" .
		"					</td>\n" .
		#"					<td width='560' valign='center'>\n" .
		#"						<font size='5' color='black'>$title</font>\n" .
		#"					</td>\n" .
		#"					<td width='180' valign='center'>\n" .
		#"						<img src='greenquoll.gif' width='180' height='100' border='1'>\n" .
		#"					</td>\n" .
		"					<td width='240'>\n" .
		"						<img src='bluegold.gif' width='240' height='160' border='0'>\n" .
		"					</td>\n" .
		"				</tr>\n" .
		"			</table>\n" .
		"		</td>\n" .
		"	</tr>\n" .
		"</table>\n" .
		"</td></tr><tr><td align='left'>\n" .
		"<table height='20' cellspacing='0' cellpadding='0' border='0'><tr><td></td></tr></table>\n<br>";

	return $html_code;
}



#===================================================================================================
#	html code for the header of the html page
#===================================================================================================

sub html_header {

	my $title = 'RENUMBER RESIDUE NUMBERS IN PDB FILE';

	print $q->header();
	print	#$q->start_html(-title=>$text_title),
		"<html><head><title>$title</title>\n",
		'<style type="text/css">' . "\n",
		'<!--' . "\n",
		'.basestyle {' . "\n",
		'	color: #FFFFFF;' . "\n",
		'	font-size: 14px;' . "\n",
		'	font-family: Arial, Helvetica, Verdana, sans-serif;' . "\n",
		'}' . "\n",
		'h1 {' . "\n",
		'	font-size: 16px;' . "\n",
		'	font-weight: bold;' . "\n",
		'	font-family: Arial, Helvetica,Verdana, sans-serif;' . "\n",
		'}' . "\n",
		'body, td, p, input {' . "\n",
		'	font-size: 14px;' . "\n",
		'	font-family: Arial, Helvetica, Verdana, sans-serif;' . "\n",
		'	color: #333333;' . "\n",
		'}' . "\n",
		'.navlinks {' . "\n",
		'	font-size: 14px;' . "\n",
		'	font-family: Arial, Helvetica, Verdana, sans-serif;' . "\n",
		'	font-weight: bold;' . "\n",
		'	color: #FFFFFF;' . "\n",
		'}' . "\n",
		'.tt {' . "\n",
		'	font-size: 12px;' . "\n",
		'	font-family: courier, monospace;' . "\n",
		'	color: #000000;' . "\n",
		'}' . "\n",
		'-->' . "\n",
		'</style>' . "\n",
		"</head><body>\n";

		my $output_line = html_code_for_graphics_header($title);
		print $output_line;
}



#===================================================================================================
#	html code for the footer of the html page
#===================================================================================================

sub html_footer {

	print $debug . "<br>\n";

	print "<br><br>\n";
	print "<b>SOURCE CODE :</b> <a href='renumberpdbchain_pl.txt'>renumberpdbchain_pl.txt</a><br>\n";

	print "<br><br>\n";
	print "<b>COPYRIGHT :</b><br>\n";
	print "Copyright &copy; 2010 Emma Rath. All rights reserved.<br>\n";
	print "Soon this program will be released under an open source software license such as GNU General Public License or<br>\n";
	print "Creative Commons license for Free Software Foundation's GNU General Public License at creativecommons.org<br>\n";

	print "<br><br>\n";
	print "</td></tr></table>\n";
	print $q->end_html();
}



#===================================================================================================
#	initialise global and program variables before starting processing.
#===================================================================================================

sub init {

	my %global_hash;
	$global = \%global_hash;
	$global->{'got_error'} = 0;
	$global->{'err_msg'} = '';

	my %input_hash;
	$input = \%input_hash;
	$input->{'infile'} = '';
	$input->{'add'} = '';
	$input->{'chain'} = '';
	$input->{'from'} = '';
	$input->{'to'} = '';
	$input->{'whereto'} = '';
	$input->{'h'} = '';
	$input->{'help'} = '';
	my @input_lines;
	$input->{'input_lines'} = \@input_lines;
}



#===================================================================================================
sub is_valid_chain {

	my $chain = shift;
	my $is_ok = 0;
	if (length($chain) <= 1) {
		$is_ok = 1;
	}
	return $is_ok;
}



#===================================================================================================
sub is_valid_decimal {

	my $input_integer = shift;
	my $return = 1;

	if ($input_integer !~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/) {
		$return = 0;
	} else {
		if ($input_integer < 1) {
			$return = 0;
		}
	}

	return $return;
}



#===================================================================================================
sub is_valid_integer {

	my $input_integer = shift;
	my $return = 1;

	if ($input_integer !~ /^\d+$/) {
		$return = 0;
	} else {
		if ($input_integer < 1) {
			$return = 0;
		}
	}

	return $return;
}



#===================================================================================================
sub is_valid_residue_number {

	my $input_integer = shift;
	my $return = 1;

	if ($input_integer !~ /^\d+$/) {
		$return = 0;
	} else {
		if ($input_integer < 1) {
			$return = 0;
		}
	}

	return $return;
}



#===================================================================================================
sub is_valid_add_number {

	my $input_integer = shift;
	my $return = 1;

	if ($input_integer !~ /^-?\d+$/) {
		$return = 0;
	}

	return $return;
}



#===================================================================================================
#	in command line mode, output some text help
#===================================================================================================

sub output_command_line_help {

	my $command_line_help_text = '';

	$command_line_help_text .= "Here are some examples of calling this program in command line mode :\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          perl -T renumberpdbchain.pl -infile 1ZOY.pdb -add 200 -chain D -from 230 -to 357\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "This program reads in a PDB 3D file, and for the given chain (column 22 in the file),\n";
	$command_line_help_text .= "renumbers all the residues (residue number is columns 23-26) of that chain,\n";
	$command_line_help_text .= "for lines starting with ATOM or HETATM,\n";
	$command_line_help_text .= "by adding the given number to the residue number.\n";
	$command_line_help_text .= "The given number can be a negative number (so that a number is subtracted from each residue number).\n";
	$command_line_help_text .= "\n";

	print $command_line_help_text;
}



#===================================================================================================
#	a perl utility subroutine - remove leading and trailing blanks from a string
#===================================================================================================

sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}



#===================================================================================================
#	untaint file name
#===================================================================================================

sub untaint {
	my $path = shift;
	# user_error("path cannot contain metacharacters") if $path=~/[\n|<>&!;\'\"]/;

	# extract only legal alphanumerics from the filename
	# $path =~ m!(/[a-zA-Z/0-9._~\-]+)!;
	# $path =~ /([\w"-"".""_"]+)/;
	$path =~ /([a-zA-Z0-9._~\-]+)/;
	if (defined($1)) {
		$path = $1;
	} else {
		$path = '';
	}

	# user_error("path cannot contain relative directories") 
	#	if $path=~m!\.\.!;

	return $path;
}
