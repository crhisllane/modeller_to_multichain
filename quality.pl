#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Parallel::ForkManager;

my $usage="\nCommand Line:\n$0 [tabela com a primeira coluna sendo a query e a segunda contendo seu template] [diretorio dos scripts] [numero de processos] [caminho para o arquivo procheck.src do seu procheck] [caminho para o arquivo modpy.sh da seu modeller]\n\n na fatboy o caminho seria: /home/crhisllane/Doutorado/programas/bin/Modeller9.19/modeller-9.19/bin/modpy.sh\n\n"; 

my $tablefile = $ARGV[0]|| die "$usage";
my $dir = $ARGV[1]|| die "$usage";
my $maxProcess = $ARGV[2]|| die "$usage";
my $seuProcheck = $ARGV[3]|| die "$usage";
my $seuModeller = $ARGV[4]|| die "$usage";


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CORE TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $numMaxprocessor=`cat /proc/cpuinfo | grep -c 'processor'`;
if ($maxProcess >= $numMaxprocessor){
	die "out of processors\n";
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


open (TAB, "$tablefile");
my $line; my @querysTemp;
while($line=<TAB>){
	chomp($line);
	push(@querysTemp,$line)
}
close TAB;
$dir=~s/\/$//;
my $pm = Parallel::ForkManager->new($maxProcess);
foreach my $queryTemp (@querysTemp){
   	my $pid=$pm->start and next;

    my @query_templine = split ('\t', $queryTemp);
    if (@query_templine[0]!~m/^P/){
        my $query = @query_templine[0];
        my $template = @query_templine[1];
        my $templateLC = lc @query_templine[1];
        my $templateUC = uc @query_templine[1];
        chdir ("$query");
        my @pdbslist = glob ("*.pdb");
        foreach my $pdblist (@pdbslist){
            mkdir ("PROCHECK_$pdblist");
            chdir ("PROCHECK_$pdblist");
            system ("cp ../$pdblist .");
            system ("$seuProcheck $pdblist 2.0");
            system ("$seuModeller $dir/ga341Dope_InputFiles $pdblist > $pdblist.log");            

            chdir ("..");
        }
        chdir ("..");

    }

}
