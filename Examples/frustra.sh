#! /bin/bash
echo "FrustratometeR"

docker run -v $1:/pdb -it proteinphysiologylab/frustratometer env KEY_MODELLER=$4 R -e "library(frustratometeR)" -e "dir_frustration(PdbsDir = '/pdb/', Mode = '$2', ResultsDir = '/pdb/')"