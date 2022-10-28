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
apt install libharfbuzz-dev libfribidi-dev
apt install libudunits2-dev
apt install libgeos-dev
apt-get install libgdal-dev
apt install cmake
apt-get install cargo
apt-get install -y libpoppler-cpp-dev
apt-get install -y libavfilter-dev
apt install libtesseract-dev libleptonica-dev
apt-get install tesseract-ocr-eng

Rscript -e "update.packages(ask=FALSE, checkBuilt=TRUE)"
Rscript -e "install.packages('usethis', dependencies = TRUE)"
Rscript -e "install.packages('devtools', dependencies = TRUE)"
Rscript -e "install.packages('ggplot2', dependencies = TRUE)"
Rscript -e "install.packages('ggrepel', dependencies = TRUE)"
Rscript -e "install.packages('dplyr', dependencies = TRUE)"
Rscript -e "install.packages('igraph', dependencies = TRUE)"
Rscript -e "install.packages('FactoMineR', dependencies = TRUE)"
Rscript -e "install.packages('Hmisc', dependencies = TRUE)"
Rscript -e "install.packages('magick', dependencies = TRUE)"
Rscript -e "install.packages('leiden', dependencies = TRUE)"
Rscript -e "options(timeout=9999999)"
Rscript -e "devtools::install_github('proteinphysiologylab/frustratometeR')"


