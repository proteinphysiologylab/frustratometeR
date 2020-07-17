#!/usr/bin/python3

###estoy corriendo script desde un folder que scripts, que tiene el nombre de la secuencia
###alineamiento esta en mismo folder de corrida


from modeller import *
from modeller.automodel import *
from sys import argv
#from modeller import soap_protein_od


molde=argv[1]
secuencia=argv[2]

env = environ()
a = automodel(env, alnfile=secuencia+'-'+molde+'.ali',
              knowns=molde, sequence=secuencia,
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 1
a.make()
