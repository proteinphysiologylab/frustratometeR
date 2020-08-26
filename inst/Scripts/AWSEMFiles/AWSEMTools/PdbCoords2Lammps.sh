#!/bin/bash

export s=$HOME/opt/script

if [[ ! $# -eq 2 ]]
then
        echo
        echo "> $0 pdb_id project_name"
        echo
        exit
fi

pdb_file=$1
output_file=$2

echo $pdb_file
echo $output_file

python3 $s/PDBToCoordinates.py $pdb_file $output_file".coord"
python3 $s/CoordinatesToWorkLammpsDataFile.py $output_file".coord" "data."$output_file -b
