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
    sub PDB{
        my ($diretorio, $pdbname) = @_;
        #system("python $diretorio/pdb_download.py temp.lst");
        system("wget https://files.rcsb.org/download/$pdbname\.pdb1.gz");
        system("gunzip $pdbname\.pdb1.gz"); 
        
    }
    sub PDBnewA{
        my ($numChain, $pdbname, $namChain) = @_;
        my $next = $numChain + 1;
        system("sed -n \'/MODEL        $numChain\/,\/MODEL        $next\/p\' $pdbname\.pdb1 > $pdbname\_$namChain.pdb");
        system("sed -i '/^MODEL        $next/d' $pdbname\_$namChain.pdb");
    }
    sub PDBnewC{
        my ($numChain, $pdbname, $namChain) = @_;
        my $next = $numChain + 1;
        system("sed -n \'/MODEL        $numChain\/,\/ENDMODEL\/p\' $pdbname\.pdb1 > $pdbname\_$namChain.pdb");
    }  

}
package work {
    sub whichChain {
        my ($ni, $template) = @_;
        if ($ni == 1){
            my $nameChain = "A";
            criando::PDBnewA($ni, $template, $nameChain);
            return "A";
        } elsif ($ni == 2){
            my $nameChain = "B";
            criando::PDBnewA($ni, $template, $nameChain);
            return "B";
        } elsif ($ni == 3){
            my $nameChain = "C";
            criando::PDBnewC($ni, $template, $nameChain);
            return "C";
        }
    }
    sub whichNumAA {
        my ($filename) = @_;
        my @allRow;
        if (open(my $fh, '<:encoding(UTF-8)', $filename)) {
            while (my $row = <$fh>) {
            chomp $row;
            push(@allRow,$row);
        }
        } else {
            warn "Não foi possível abrir o arquivo '$filename' $!";
        }
        my @IniAtom = grep {/ATOM/} @allRow; 
        return $IniAtom[1]; 
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
        my $templateLC = lc @query_templine[1];
        my $templateUC = uc @query_templine[1];

        print ("prot $query temp $template\n");
        criando::dir($query);
        chdir ("$query");
        criando::file($template);
        criando::PDB($dir, $templateUC);


		while(my $seq = $seqs_in->next_seq){
			my $id = $seq->display_id;
            my $sequence=$seq->seq;
            if ($id eq $query){
                criando::fasta($id, $sequence);
            }
		}
        
        #########alinhamentos
        for (my $i=1; $i <= $chains; $i++){
            my $nameChain = work::whichChain($i, $templateUC);
            #criando::PDBnew($i, $templateUC, $nameChain);
            my $fileName=$templateUC . "_" . $nameChain . ".pdb";
            my $NumAA = work::whichNumAA($fileName);
            my $integerSerial = substr $NumAA, 7, 4;
            my $chainID = substr $NumAA, 21, 1;
            my $resSeq = substr $NumAA, 23, 3;
            print "\n\n\n ####Primeira Num $NumAA\n >>> $integerSerial\t$chainID\t$resSeq\n\n";
            #print ("$seuModeller python $dir/align.py $templateLC  $templateLC\_$nameChain $templateLC\.pdb $query\.fasta FASTA $query $nameChain $nameChain");
            #system ("$seuModeller python $dir/align.py $templateLC  $templateLC\_$nameChain $templateLC\.pdb $query\.fasta FASTA $query $nameChain $nameChain");
        }
        chdir ("..");
    }
    $pm->finish;
}
$pm->wait_all_children;


