library(devtools)
library(bio3d)
library(usethis)
library(devtools)
library(ggplot2)
library(reshape2)
library(magick)

install_github("proteinphysiologylab/frustratometeR")

library(frustratometeR)

PdbPath="/home/atilio/FrustratometeR/Examples/1n0r.pdb"

JobDir="/home/atilio/Escritorio/Result/"

# Calculate frustration for the 20 AA in the Resno position
mutate_res(PdbPath = PdbPath, Resno = 91, Chain = "A", JobDir = JobDir)
#DataFile is the full path of the resulting mutate_res table, stored in Jobdir / FrustrationData
DataFile = "/home/atilio/Escritorio/Result/FrustrationData/configurational_Res91.txt"
#Graph the frustration index of each AA variant for each existing contact
plot_mutate_res(PdbPath = PdbPath, DataFile = DataFile, Resno = 91, Chain = "A", Modes = "configurational", ResultDir = JobDir)


# Calculate frustration for the 20 AA in the Resno position. In singleresidue mode.
mutate_res(PdbPath = PdbPath, Modes="singleresidue", Resno = 91, Chain = "A", JobDir = JobDir)
mutate_res(PdbPath = PdbPath, Modes="singleresidue", Resno = 92, Chain = "A", JobDir = JobDir)
mutate_res(PdbPath = PdbPath, Modes="singleresidue", Resno = 93, Chain = "A", JobDir = JobDir)
#DataDir is the full path of the directory where the resulting tables of mutate_res are located in singleresidue mode
DataDir = "/home/atilio/Escritorio/Result/FrustrationData/"
#Graph the variation of the Frustration Index with respect to native AA. There can be n tables in DataDir of singleresidue results, the algorithm is able to graph them
plot_delta_frus(PdbPath = PdbPath, DataDir = DataDir, Chain = "A", ResultDir = JobDir)
