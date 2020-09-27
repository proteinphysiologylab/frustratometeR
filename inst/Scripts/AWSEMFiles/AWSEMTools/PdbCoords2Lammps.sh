#!/bin/bash

pdb_file=$1
output_file=$2

echo Pdb file: $pdb_file
echo Output file: $output_file

python3 $3/AWSEMFiles/AWSEMTools/PDBToCoordinates.py $pdb_file $output_file".coord"
python3 $3/AWSEMFiles/AWSEMTools/CoordinatesToWorkLammpsDataFile.py $3 $output_file".coord" "data."$output_file -b
