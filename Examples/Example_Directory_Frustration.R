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

# Calculate frustration of Pdbs contained in PdbDir
dir_frustration(PdbDir = PdbDir,Modes = "configurational",ResultsDir = ResultsDir)

# Calculate frustration of the Pdbs contained in PdbDir. With electrostatics = -1 and seqdist = 3
dir_frustration(PdbDir = PdbDir,Modes = "mutational",Electrostatics_K = -1, seqdist = 3, ResultsDir = ResultsDir)

# Calculate frustration of the Pdbs contained in PdbDir and indicated in OrderList
OrderList=c("pdb2.pdb","pdb14.pdb","pdb26.pdb")
dir_frustration(PdbDir = PdbDir, OrderList = OrderList, Modes = "configurational", ResultsDir = ResultsDir)
