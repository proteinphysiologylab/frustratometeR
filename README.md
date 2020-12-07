# FrustratometeR
An R package for calculating energetic local frustration in protein structures. Additionally to the functionalities present at the web server, frustratometeR offers the possibility to analyse frustration across molecular dynamics simulations and to assess the effect of aminoacid variants.

# Installation 

`invisible(lapply(c("usethis", "devtools"), library, character.only = TRUE))`

`devtools::install_github("proteinphysiologylab/frustratometeR")`

# Dependencies

**R packages (R version >= 3.6.3 required. Compatible with R 4.0.X)**

`install.packages(c("usethis", "devtools"))`

**Python 3**

* numpy 
`python3 -m pip install numpy`

* biopython
`python3 -m pip install biopython`

* leiden(optional)
`python3 -m pip install leidenalg`

**others**

* pymol
`sudo apt install pymol`

* magick
`sudo apt-get install -y libmagick++-dev`

* perl
`sudo apt-get install perl`

**modeller**
(you need to install a licenced version of modeller in order to use some of the frustratometeR functionalities, https://salilab.org/modeller/)
modeller needs to be available via python3. To check if this is working execute python3 and try to do "form modeller import * " after installing modeller. We suggest to install modeller using conda environments. 

# Minimum code to calculate frustration in a protein
`library(frustratometeR)`

`Pdb=calculate_frustration(PdbID = "1n0r",Chain = "A",  ResultsDir = "/Users/parra/Desktop/" )`

`view_frustration_pymol(Pdb)`

## **You can find an example of how to use the package at:**

https://github.com/proteinphysiologylab/frustratometeR/tree/master/Examples

## **You can also find useful examples in our wiki!!:**

https://github.com/proteinphysiologylab/frustratometeR/wiki

