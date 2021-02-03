#! /bin/bash
echo "FrustratometeR"

if [ "configurational" = "$2" ]; then
    docker run -v $1:/pdb -it cesarleonetti/frustratometer:version1 R -e "library(frustratometeR)" -e "dir_frustration(PdbsDir = '/pdb/', Mode = '$2', ResultsDir = '/pdb/')"
    exit
fi

if [ "mutational" = "$1" ]; then
    docker run -v $1:/pdb -it cesarleonetti/frustratometer:version1 R -e "library(frustratometeR)" -e "dir_frustration(PdbsDir = '/pdb/', Mode = '$2', ResultsDir = '/pdb/')"
    exit
fi

if [ "singleresidue" = "$1" ]; then
    docker run -v $1:/pdb -it cesarleonetti/frustratometer:version1 R -e "library(frustratometeR)" -e "dir_frustration(PdbsDir = '/pdb/', Mode = '$2', ResultsDir = '/pdb/')"
    exit
fi
