import sys

from modeller import *
from modeller.scripts import complete_pdb

PDBToComplete = sys.argv[1]

env = environ()
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')

m = complete_pdb(env, PDBToComplete, transfer_res_num="true")
m.write(file=PDBToComplete+"_completed")
