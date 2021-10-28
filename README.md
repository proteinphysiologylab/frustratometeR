# frustratometeR
An R package for calculating energetic local frustration in protein structures. Additionally to the functionalities present at the web server, frustratometeR offers the possibility to analyse frustration across molecular dynamics simulations and to assess the effect of aminoacid variants.

# Installation 

`invisible(lapply(c("usethis", "devtools"), library, character.only = TRUE))`

`devtools::install_github("proteinphysiologylab/frustratometeR")`

# Dependencies

** If you are using MacOs, FrustratometeR is compatible with systems using MacOS version 10.15 or higher.

** If you are using UBUNTU, FrustratometeR is compatible with Ubuntu 16 or higher.


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
modeller needs to be available via python3. To check if this is working execute python3 and try to do "from modeller import * " after installing modeller. We suggest to install modeller using conda environments. 

if modeller does not work from RStudio, you can use reticulate to tell R which conda environment to use. Execute these lines before you execute the frustration functions.

`library(reticulate)`

`use_python("/url/bin/python3")`

`use_condaenv("condaenvironment") `

`Sys.setenv(RETICULATE_PYTHON = "/url/bin/python3")`

`reticulate::py_config()`


# Minimum code to calculate frustration in a protein
`library(frustratometeR)`

`Pdb_conf <- calculate_frustration(PdbID = "1n0r",Chain = "A",  ResultsDir = "/Home/Desktop/")`

`view_frustration_pymol(Pdb_conf)`

## **You can find an example of how to use the package at:**

https://github.com/proteinphysiologylab/frustratometeR/tree/master/Examples

## **You can also find useful examples in our wiki!!:**

https://github.com/proteinphysiologylab/frustratometeR/wiki

## **You can find the data used in the Wiki examples at:**

https://github.com/proteinphysiologylab/frustratometer-data.git
