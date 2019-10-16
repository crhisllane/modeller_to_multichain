# Comparative modeling by the automodel class - editado por Crhisllane para inserir o processo de otmizacao e restricoes das cadeias
#
# Demonstrates how to build multi-chain models, and symmetry restraints
#
from modeller import *
from modeller.automodel import *    
import sys

fileali=sys.argv[1]			#arquivo gerado pelo alinhamento.ali ex:TvLDH-1bdmA.ali
template=sys.argv[2]			#nome do arquivo pdb ex:1bdmA
query=sys.argv[3]			#nome da query ex:TvLDH

log.verbose()
env = environ()
env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)

# Override the 'special_restraints' and 'user_after_single_model' methods:
class MyModel(automodel):
    def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = selection(self.chains['A']).only_atom_types('CA')
        s2 = selection(self.chains['B']).only_atom_types('CA')
        s3 = selection(self.chains['C']).only_atom_types('CA')
        self.restraints.symmetry.append(symmetry(s1, s2, 1.0))
        self.restraints.symmetry.append(symmetry(s2, s3, 1.0))
    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)

env = environ()
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Be sure to use 'MyModel' rather than 'automodel' here!
a = MyModel(env,
            alnfile  = fileali ,     # alignment filename
            knowns   = template,              # codes of the templates
            sequence = query)              # code of the target

a.starting_model= 1                # index of the first model
a.ending_model  = 5                # index of the last model
                                   # (determines how many models to calculate)

# Very thorough Variable Target Function Method (VTFM) optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 300

# Thorough MD optimization:
a.md_level = refine.slow

# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
a.repeat_optimization = 2
a.max_molpdf = 1e6

a.make()

# Get clusters
a.cluster(cluster_cut=1.00)
# END OF MODEL CONSTRUCTION


# PRINT RESULTS
# Open a file
fo = open("model-single-opt.out", "w")

# Get a list of all successfully built models from a.outputs
ok_models = filter(lambda x: x['failure'] is None, a.outputs)

# Printing out a summary of all successfully generated models
print >> fo, '\n>> Summary of successfully produced model'
fields = [x for x in ok_models[0].keys() if x.endswith(' score')]
fields.sort()
fields = ['molpdf'] + fields
header = '%-25s ' % 'Filename' + " ".join(['%14s' % x for x in fields])
print >> fo, header
print >> fo, '-' * len(header)
for mdl in ok_models:
    text = '%-25s' % mdl['name']
    for field in fields:
	if isinstance(mdl[field], (tuple, list)):
	    text = text + ' %14.5f' % mdl[field][0]
	else:
	    text = text + ' %14.5f' % mdl[field]
    print >> fo, text
print >> fo, ''

# Printing top model results
print >> fo, '>> Top model results:'
  
# Rank models by molpdf score
key = 'molpdf'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
# Get top model - molpdf
m = ok_models[0]
print "Top model_molpdf: %s (molpdf %.3f)" % (m['name'], m[key])
print >> fo, 'molpdf: ', m[key], '(file: ', m['name'], ')'

# Rank models by DOPE score
key = 'DOPE score'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
# Get top model - DOPE
m = ok_models[0]
print "Top model_DOPE: %s (DOPE score %.3f)" % (m['name'], m[key])
print >> fo, 'DOPE score: ', m[key], '(file: ', m['name'], ')'

# Rank models by normalized DOPE score
key = 'GA341 score'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
# Get top model - normalized DOPE
m = ok_models[0]
print "Top model_GA341: %s (GA341 score %.3f)" % (m['name'], m[key][0])
print >> fo, 'GA341 score: ', m[key][0], '(file: ', m['name'], ')'

# Rank models by normalized DOPE score
key = 'Normalized DOPE score'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
# Get top model - normalized DOPE
m = ok_models[0]
print "Top model_nDOPE (z): %s (Normalized DOPE score %.3f)" % (m['name'], m[key])
print >> fo, 'Normalized DOPE score: ', m[key], '(file: ', m['name'], ')'

# Read a model previously generated by Modeller's automodel class
mdl = complete_pdb(env, './cluster.opt')

# Select all atoms in the first chain
atmsel = selection(mdl)

score = atmsel.assess_dope()
zscore = mdl.assess_normalized_dope()
score2 = mdl.assess_ga341()

# Printing assess results
print >> fo, '\n>> Cluster results:'

fo2 = open("cluster.opt", "r")
lines = [ i.rstrip() for i in fo2.readlines()]
# 3rd line
print >> fo, lines[1], '(molpdf)'

print >> fo, 'DOPE score: ', score
print >> fo, 'GA341 score: ', score2[0]
print >> fo, 'Normalized DOPE score: ', zscore

# Close opened file
fo.close()
#END OF PRINT RESULTS
