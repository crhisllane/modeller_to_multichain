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
    sub dividePDB{
        my ($filename, $templateUC, $newChain) = @_;
        print ("::::::: $filename, $templateUC, $newChain :::::::\n");
        my @allRow;
        open (FILECHAIN, "$filename");
        while(my $row=<FILECHAIN>){
	        chomp($row);
	        push(@allRow,$row)
        }    
        my $testeChain = "WW";
        foreach my $oneRow (@allRow){
            if ($oneRow=~m/ATOM/){
                my $integerSerial = substr $oneRow, 7, 4;
                my $chainID = substr $oneRow, 21, 1;
                my $resSeq = substr $oneRow, 23, 3;
                if ($testeChain eq $chainID){

                    print BYCHAIN "$oneRow\n";

                }elsif (($testeChain eq "WW" || $testeChain ne $chainID)){
                    if ($testeChain ne $chainID){
                        close BYCHAIN;
                    }
                    $testeChain = $chainID;
                    open (BYCHAIN, ">>$templateUC\_$newChain\_$chainID\.pdb");
                    print BYCHAIN "$oneRow\n";
                
                }
                
            }
        }
        close BYCHAIN;

    }
    sub PDBrenum{
        my ($templateUC, $newChain, $resSeq, $dirL) = @_;
        my @files = glob("$templateUC\_$newChain\_*");
        my $control;
        my $resnumb = $resSeq;
        foreach my $file (@files){                       
            if ($file=~m/$templateUC\_$newChain\_$newChain/){
                $control = "S";
                my $resnumb1 = $resSeq - 1;
                $resnumb = work::renumber($file, $resnumb1, $newChain, $dirL, $control);
                print ("\n\n:::::::PRIMEIRO resnum $resSeq e virou $resnumb :::::::\n\n");
            } 
        }
        foreach my $file (@files){
            if ($file!~m/$templateUC\_$newChain\_$newChain/) {
                $control = "N";
                my $firsline = `head -1 $file`;
                my $firstres = substr $firsline, 23, 3;
                my $firstresOK;
                if ($firstres != 1){
                    $firstresOK = $firstres - 1;
                    system ("$dirL/renamepdbchain.pl -infile $file -tochain $newChain");
                    system ("$dirL/renumberpdbchain.pl -infile $file\.pdb -add -$firstresOK -chain $newChain");           
                    system ("mv $file\.pdb\.pdb $file");
                    system ("rm $file\.pdb*");

                    $resnumb = work::renumber($file, $resnumb, $newChain, $dirL, $control);
                } else {

                    $resnumb = work::renumber($file, $resnumb, $newChain, $dirL, $control);
                }

            }            

        }
  

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
    sub renumber {
        my ($file, $resnumb, $newChain, $dirL, $control) =@_;
        my $resSeqN = $resnumb;
        system ("$dirL/renamepdbchain.pl -infile $file -tochain $newChain");
        if ($control eq "S"){
            system ("$dirL/renumberpdbchain.pl -infile $file\.pdb -add -$resSeqN -chain $newChain");
        } 
        if ($control eq "N"){
            system ("$dirL/renumberpdbchain.pl -infile $file\.pdb -add $resSeqN -chain $newChain");
        }
        my $lastLine = `tail -1 $file\.pdb\.pdb`;
        my $resOut = substr $lastLine, 23, 3;
        print (":::$lastLine - ANTIGO $resSeqN - NOVO $resOut:::\n");
        return $resOut;
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
            my $fileName=$templateUC . "_" . $nameChain . ".pdb";
            my $NumAA = work::whichNumAA($fileName);
            my $integerSerial = substr $NumAA, 7, 4;
            my $chainID = substr $NumAA, 21, 1;
            my $resSeq = substr $NumAA, 23, 3;
            print "\n\n\n ####Primeira Num $NumAA\n >>> $integerSerial\t$chainID\tmudarpara $nameChain \t$resSeq\n\n";
            criando::dividePDB($fileName, $templateUC, $nameChain);
            criando::PDBrenum($templateUC, $nameChain, $resSeq, $dir);
            #system ("$dir/renamepdbchain.pl -infile $fileName -tochain $nameChain");
            #system ("mv $fileName.pdb $fileName");
            #print ("$seuModeller python $dir/align.py $templateLC  $templateLC\_$nameChain $templateLC\.pdb $query\.fasta FASTA $query $nameChain $nameChain");
            #system ("$seuModeller python $dir/align.py $templateLC  $templateLC\_$nameChain $templateLC\.pdb $query\.fasta FASTA $query $nameChain $nameChain");
        }
        chdir ("..");
    }
    $pm->finish;
}
$pm->wait_all_children;



