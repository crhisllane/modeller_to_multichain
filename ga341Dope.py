from modeller import *
from modeller.scripts import complete_pdb
import sys

#programa para gerar o valor de ga341 e dope a partir de arquivos.pdb gerados pela modelagem


arqpdb=sys.argv[1]

env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read a model previously generated by Modeller's automodel class
mdl = complete_pdb(env, arqpdb)

# Set template-model sequence identity. (Not needed in this case, since
# this is written by Modeller into the .pdb file.)
mdl.seq_id = 37.037

score = mdl.assess_ga341()


env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read a model previously generated by Modeller's automodel class
mdl = complete_pdb(env, arqpdb)

zscore = mdl.assess_normalized_dope()