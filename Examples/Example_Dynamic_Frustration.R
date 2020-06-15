library(devtools)
library(bio3d)
library(usethis)
library(devtools)
library(ggplot2)
library(reshape2)
library(magick)

install_github("proteinphysiologylab/frustratometeR")

library(frustratometeR)

PdbDir="/home/atilio/FrustratometeR/inst/Dynamic_Pdbs/"

ResultsDir="/home/atilio/Escritorio/Result/"

# Calculate frustration of the frames along the dynamics in PdbDir. With Gifs. Without analyzing specific residue
dynamic_frustration(PdbDir = PdbDir, Modes = "configurational", ResultsDir = ResultsDir, GIFs = TRUE)

# Calculate frustration of the frames along the dynamics in PdbDir. With Gifs. Analyzing residue Resno
dynamic_frustration(PdbDir = PdbDir, Modes = "mutational", ResultsDir = ResultsDir,Resno=32, GIFs = FALSE)

# Calculate frustration of frames throughout the dynamics in PdbDir, calculated for OrderList. No gifs. Analyzing residue Resno
OrderList=c("pdb2.pdb","pdb14.pdb","pdb26.pdb")
dynamic_frustration(PdbDir = PdbDir, OrderList=OrderList, Modes = "singleresidue", ResultsDir = ResultsDir,Resno=32, GIFs = FALSE)


# In general cases, the Resno Chain must be specified since the algorithm checks that said residue exists.
