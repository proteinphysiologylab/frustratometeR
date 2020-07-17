#!/usr/bin/python3

###estoy corriendo script desde un folder que scripts, que tiene el nombre de la secuencia
###tengo pdb en mismo folder que scripts
###tengo uniprot en folder /scripts/UniprotsLahey

from modeller import *
from sys import argv

molde=argv[1]
secuencia=argv[2]
chain=argv[3]

env = environ()
aln = alignment(env)
mdl=()
if chain != "NULL":
    mdl = model(env, file=molde+'.pdb', model_segment=('FIRST:'+chain, 'LAST:'+chain))
else:
    mdl = model(env, file=molde+'.pdb')
        

aln.append_model(mdl, align_codes=molde, atom_files=molde+'.pdb')
aln.append(file=secuencia+'.ali', align_codes=secuencia)
aln.align2d()
aln.write(file=secuencia+'-'+molde+'.ali', alignment_format='PIR')
aln.write(file=secuencia+'-'+molde+'.pap', alignment_format='PAP')
