from modeller import *
import sys

template=sys.argv[1]				#nome do template sem o ".pdb". ex: 1R26
codtemplate=sys.argv[2]				#nome do codigo, esse nome vai sair no resultado do alinhamento tambem. ex: 1r26_A	
pdb=sys.argv[3]					#arquivo pdb. ex: 1R26.pdb
query=sys.argv[4]				#arquivo da query. ex:LbrM.01.0300.fasta
formato=sys.argv[5]				#formato que esta o arquivo. ex: FASTA
codquery=sys.argv[6] 				#nome da query. ex:LbrM.01.0300
first=sys.argv[7]
last=sys.argv[8]

first='FIRST:' +first
last='LAST:' +last
fileali = codquery + "_" + template + ".ali"

env = environ()
aln = alignment(env)
mdl = model(env, file=template, model_segment=(first,last))
aln.append_model(mdl, align_codes=codtemplate, atom_files=pdb)
aln.append(file=query, alignment_format=formato, align_codes=codquery
)
aln.align2d()
aln.write(file=fileali, alignment_format='PIR') 	#precisa corrigir essa parte
aln.write(file='templatequery.ali', alignment_format='PIR') 	#precisa corrigir essa parte
#aln.write(file='templatequery.pap', alignment_format='PAP')  
