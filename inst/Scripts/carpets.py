import sys, os
import os.path as path

fname=sys.argv[1]+"/FrustrationData"

if path.exists(fname):
	print 'Exist\n'
else:
	mkdirFD="cd "+sys.argv[1]+"; mkdir FrustrationData"
	os.system(mkdirFD)

fname=sys.argv[1]+"/Images"

if path.exists(fname):
		print 'Exist\n'
else:
	mkdirI="cd "+sys.argv[1]+"; mkdir Images"
	os.system(mkdirI)

fname=sys.argv[1]+"/VisualizationScripts"

if path.exists(fname):
		print 'Exist\n'
else:
	mkdirVS="cd "+sys.argv[1]+"; mkdir VisualizationScripts"	
	os.system(mkdirVS)

mv="mv "+sys.argv[1]+"/"+sys.argv[2]+".pdb_"+sys.argv[3]+"_5adens "+sys.argv[1]+"/FrustrationData/"+sys.argv[2]+".pdb_"+sys.argv[3]+"_5adens"
os.system(mv)

mv="mv "+sys.argv[1]+"/"+sys.argv[2]+".pdb_"+sys.argv[3]+" "+sys.argv[1]+"/FrustrationData/"+sys.argv[2]+".pdb_"+sys.argv[3]
os.system(mv)

mv="mv "+sys.argv[1]+"/"+sys.argv[2]+".pdb_"+sys.argv[3]+".pml"+" "+sys.argv[1]+"/VisualizationScripts/"+sys.argv[2]+".pdb_"+sys.argv[3]+".pml"
os.system(mv)

mv="cp "+sys.argv[1]+"/"+sys.argv[2]+".pdb"+" "+sys.argv[1]+"/VisualizationScripts/"+sys.argv[2]+".pdb"
os.system(mv)

mv="mv "+sys.argv[1]+"/"+sys.argv[2]+"_"+sys.argv[3]+".tcl"+" "+sys.argv[1]+"/VisualizationScripts/"+sys.argv[2]+"_"+sys.argv[3]+".tcl"
os.system(mv)

mv="mv "+sys.argv[1]+"/draw_links.py"+" "+sys.argv[1]+"/VisualizationScripts/draw_links.py"
os.system(mv)

mv="mv "+sys.argv[1]+"/"+sys.argv[2]+"_"+sys.argv[3]+".jml"+" "+sys.argv[1]+"/VisualizationScripts/"+sys.argv[2]+"_"+sys.argv[3]+".jml"
os.system(mv)

rm="cd "+sys.argv[1]+"/; rm *.in *.coord *vps *.txt *.seq *auxiliar *.dat *data* *.lammpstrj *.lammps *log *12 *.tcl"
os.system(rm)
