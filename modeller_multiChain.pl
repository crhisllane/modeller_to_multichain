#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Parallel::ForkManager;

my $usage="\nCommand Line:\n$0 [arquivo com os fastas das querys] [tabela com a primeira coluna sendo a query e a segunda contendo seu template] [diretorio dos scripts] [numero de processos] [2 para dimeros 3 para trimeros] [caminho para o arquivo modpy.sh da seu modeller]\n\n na fatboy o caminho seria: /home/crhisllane/Doutorado/programas/bin/Modeller9.19/modeller-9.19/bin/modpy.sh\n\n"; 

my $queryFastaFile = $ARGV[0]|| die "$usage"; 
my $tablefile = $ARGV[1]|| die "$usage";
my $dir = $ARGV[2]|| die "$usage";
my $maxProcess = $ARGV[3]|| die "$usage";
my $chains = $ARGV[4]|| die "$usage";
my $seuModeller = $ARGV[5]|| die "$usage";
 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CORE TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $numMaxprocessor=`cat /proc/cpuinfo | grep -c 'processor'`;
if ($maxProcess >= $numMaxprocessor){
	die "out of processors\n";
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

package criando {
  
  sub dir {
    my ($n1) = @_;
    mkdir ("$n1") ;
  }

  sub file {
    my ($n1) = @_;
    open (TEMP, ">>temp.lst");
    print TEMP ("$n1\n");
    close TEMP;
  }

  sub fasta{
    my ($nid, $nseq) = @_;
    open (FILESEQ, ">>$nid\.fasta");
	print FILESEQ ">$nid\n$nseq\n";
    #print ">$nid\n$nseq\n";
	close FILESEQ;
  }

}
package teste {
    sub whichChain {
        my ($ni) = @_;
        if ($ni == 1){
            return "A";
            print ("chegou aqui");
        } elsif ($ni == 2){
            return "B";
        } elsif ($ni == 3){
            return "C";
        }
    }
}


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
    my $seqs_in = Bio::SeqIO->new(-format => 'fasta', -file => "$queryFastaFile");
    if (@query_templine[0]!~m/^P/){
        my $query = @query_templine[0];
        my $template = @query_templine[1];
        my $templateUC = lc @query_templine[1];
        print ("prot $query temp $template\n");
        criando::dir($query);
        chdir ("$query");
        criando::file($template);
        system("python $dir/pdb_download.py temp.lst");

		while(my $seq = $seqs_in->next_seq){
			my $id = $seq->display_id;
            my $sequence=$seq->seq;
            if ($id eq $query){
                criando::fasta($id, $sequence);
            }
		}
        
        #########alinhamentos
        for (my $i=1; $i <= $chains; $i++){
            my $nameChain = teste::whichChain($i);
            #print ("$seuModeller python $dir/align.py $templateUC  $templateUC\_$nameChain $templateUC\.pdb $query\.fasta FASTA $query $nameChain $nameChain");
            system ("$seuModeller python $dir/align.py $templateUC  $templateUC\_$nameChain $templateUC\.pdb $query\.fasta FASTA $query $nameChain $nameChain");
        }
        chdir ("..");
    }
    $pm->finish;
}
$pm->wait_all_children;


