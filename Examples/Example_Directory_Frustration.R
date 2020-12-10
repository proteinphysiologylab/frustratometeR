library(frustratometeR)

PdbsDir <- "/your/path/here/"
ResultsDir <- "/your/path/here/Results/"

# Calculate frustration of Pdbs contained in PdbDir
dir_frustration(PdbsDir = PdbsDir, Mode = "configurational", ResultsDir = ResultsDir)

# Calculate frustration of the Pdbs contained in PdbDir. With electrostatics = -1 and seqdist = 3
dir_frustration(PdbsDir = PdbsDir, Mode = "mutational", Electrostatics_K = -1, seqdist = 12, ResultsDir = ResultsDir)

# Calculate frustration of the Pdbs contained in PdbDir and indicated in OrderList
# The listed structures are found in: https://github.com/proteinphysiologylab/frustratometeR/tree/master/Examples/Dynamic_Pdbs
OrderList <- c("pdb2.pdb", "pdb14.pdb", "pdb26.pdb")
dir_frustration(PdbsDir = PdbsDir, OrderList = OrderList, Mode = "configurational", ResultsDir = ResultsDir)
