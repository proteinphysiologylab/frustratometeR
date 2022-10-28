#!/bin/bash
#FROM ubuntu:22.04

#MAINTAINER Freiberger Maria Ines
#For R 

apt-get update
apt install -y r-base r-base-core r-recommended r-base-dev

#For frustraR package installation

apt update -y
apt install wget mc r-base python3 python3-pip pymol libmagick++-dev libcurl4-openssl-dev libssl-dev libgit2-dev -y
apt install libcurl4-gnutls-dev libxml2-dev -y
python3 -m pip install numpy biopython leidenalg
apt install software-properties-common -y

Rscript -e "update.packages(ask=FALSE, checkBuilt=TRUE)"
Rscript -e "install.packages('usethis', dependencies = TRUE)"
Rscript -e "install.packages('devtools', dependencies = TRUE)"
Rscript -e "install.packages(c('ggplot2','ggrepel','dplyr','igraph','FactoMineR','Hmisc','msa','magick','leiden'))"
Rscript -e "devtools::install_github('proteinphysiologylab/frustratometeR')"

