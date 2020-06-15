# frustratometeR
An R package for calculating energetic local frustration in protein structures. Additionally to the functionalities present at the web server, frustratometeR offers the possibility to analyse frustration across molecular dynamics simulations and to assess the effect of aminoacid variants.

# installation 

`library(usethis)`

`library(devtools)`

`library(bio3d)`

`library(ggplot2)`

`library(reshape2)`

`library(magick)`

`install_github("proteinphysiologylab/frustratometeR")`

# dependencies

**R packages**
`install.packages(c("bio3d", "usethis", "devtools", "ggplot2", "reshape2", "magick"))`

**Python 3**

numpy 

`python3 -m pip install numpy`

biopython

`python3 -m pip install biopython`


**others**

pymol

`sudo apt install pymol`

magick

'sudo apt-get install -y libmagick++-dev'

perl

'sudo apt-get install perl'

modeller (you need to install a licenced version of modeller in order to use some of the frustratometeR functionalities)

https://salilab.org/modeller/

**You can find an example of how to use the package at:**
