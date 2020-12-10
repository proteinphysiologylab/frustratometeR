library(devtools)
library(bio3d)
library(usethis)
library(devtools)
library(ggplot2)
library(reshape2)
library(magick)

install_github("proteinphysiologylab/frustratometeR")

library(frustratometeR)

PdbFile="/home/atilio/FrustratometeR/Examples/1n0r.pdb"

ResultsDir="/home/atilio/Escritorio/"
Modes="configurational"

# Calculate frustration for a givenfile
Pdb=calculate_frustration(PdbFile=PdbFile, Mode = "configurational", ResultsDir = ResultsDir,Graphics = TRUE)

# Calculate frustration for a structure to be downloaded from the PDB
Pdb=calculate_frustration(PdbID="1vkx", Mode = "configurational", ResultsDir = ResultsDir,Graphics = FALSE)

#Calculate frustration for a particular chain in a structure to be downloaded from the PDB
Pdb=calculate_frustration(PdbID="1ikn", Chain="D", Mode = "configurational", ResultsDir = ResultsDir)


plot_5Andens(Pdb, chain=NULL)
plot_5Adens_proportions(Pdb, chain=NULL)
plot_contact_map(Pdb, chain="A")

view_frustration_pymol(Pdb)
