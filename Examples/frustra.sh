#! /bin/bash

docker run -v $1:/pdb --rm -it proteinphysiologylab/frustratometer:latest /bin/bash -c "sh script1.sh inicio $2 $3"
cd $1
sudo chown -R $(whoami)