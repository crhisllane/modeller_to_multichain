#!/usr/bin/perl
use Bio::SeqIO;
#my $seqs_in = Bio::SeqIO->new(-format => 'fasta', -file => "AnGamPim_LCWJ01005811_1_1029608_1031407_MC_6ESC_FIRST:A.ali");
my $chains = 3;
my $query = "AnGamPim_LCWJ01005811_1_1029608_1031407";
my $template = "6ESC";
my $dir = "../../Git_modeller_to_multichain/modeller_to_multichain";

package criando {

    sub idAli {

        my ($aliEach) = @_;
        open (ALI, "$aliEach");
        my @aliLines, my $line;
        while(my $line=<ALI>){
            chomp($line);
            push(@aliLines,$line)
        }
        close ALI;
        my @cabecalhoCompleto;

        foreach my $aliLine (@aliLines){
            if ($aliLine=~m/^\>/){
                #print ("CAB $aliLine\n");
                push (@cabecalhoCompleto,$aliLine);
            } elsif(($aliLine=~m/^struc/) || ($aliLine=~m/^sequence/)) {
                push (@cabecalhoCompleto,$aliLine);

            }
        }
        return @cabecalhoCompleto;

    }
    sub aliFinal {
        my ($aliEach) = @_;
        open (ALI, "$aliEach");
        my @aliLines, my $line;
        while(my $line=<ALI>){
            chomp($line);
            push(@aliLines,$line)
        }
        close ALI;
        my $controlAli = 0;
        open (TEMP1, ">>SEQ1");
        open (TEMP2, ">>SEQ2");
        foreach my $aliLine (@aliLines){
            if (($aliLine=~m/^\>/) && ($controlAli == 0)){
                print ("CAB1 $aliLine");
            }elsif (($aliLine=~m/^\>/) && ($controlAli != 0)){
                print ("CAB2 $aliLine");
                $controlAli ++;
            }elsif ($aliLine=~m/^struc/){
                print ("SUB1 $aliLine");
                $controlAli ++;
            }elsif ($aliLine=~m/^sequence/){
                print ("SUB2 $aliLine\n");
            }elsif (($aliLine!~m/^\>/) && ($aliLine!~m/^struc/) && ($aliLine!~m/^sequence/) && ($controlAli == 1)){
                print TEMP1 ("\n$aliLine");
                
            }elsif (($aliLine!~m/^\>/) && ($aliLine!~m/^struc/) && ($aliLine!~m/^sequence/) && ($controlAli == 2)){
                #print ("SEQ2 $aliLine\n");
                print TEMP2 ("\n$aliLine");


            }
        }
        print TEMP1 ("/"); 
        close TEMP1;
        print TEMP2 ("/");
        close TEMP2;

    }
}


my @controlChain = ("", "A".."Z");
my @allAlis = glob "$query\_MC\_$template*.ali";
#print ("$allAlis[0]\n");
my @cab = criando::idAli($allAlis[0]);
#print ("@cab[0]\n@cab[1]\n@cab[2]\n@cab[3]\n");


for (my $i=1; $i <= $chains; $i++){
    my $ali = "$query\_$template\_$controlChain[$i]\_$controlChain[$i]\.ali";
    #print ("\n\n\n::::$ali\n\n\n");
    criando::aliFinal("$ali");
}

open (ALICOM, ">completechains.ali");

#fazendo a seqTemp
system ("sed -i \'/^\$/d\' SEQ1");
system ("sed -i -e \'s/*//\' SEQ1");
system ("sed -i -e \'s/*//\' SEQ2");


print ALICOM ("@cab[0]\n@cab[1]\n");
open (TEMP1, "SEQ1");
while( my $line=<TEMP1>){
	chomp($line);
	print ALICOM ("$line");
}
close TEMP1;

#fazendo a seqQuery
print ALICOM ("\n\n@cab[2]\n@cab[3]\n");
open (TEMP2, "SEQ2"); 
while( my $line=<TEMP2>){
	chomp($line);
	print ALICOM ("$line");
}
close ALICOM;
system ("rm SEQ1 SEQ2");
print ("\n$dir/toBarras.sh\n");
system ("$dir/toBarras.sh completechains.ali");