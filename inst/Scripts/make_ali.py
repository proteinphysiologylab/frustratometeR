#!/usr/bin/python3

###makes .PIR file from .fas

from sys import argv

secuencia=argv[1]

fasta=open(''+secuencia+'.fa', 'r')
salida=open(''+secuencia+'.ali', 'a')
salida.write('>P1;'+str(secuencia)+'\nsequence:'+str(secuencia)+':::::::0.00: 0.00')
for linea in fasta:
	if linea[0]!='>':
		salida.write('\n'+str(linea.rstrip('\n')))
salida.write('*')

fasta.close()
salida.close()


