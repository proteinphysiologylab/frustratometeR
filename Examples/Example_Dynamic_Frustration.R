#Example adapted for use with the dynamics frames found in https://github.com/proteinphysiologylab/frustratometeR/tree/master/Examples/Dynamic_Pdbs

library(frustratometeR)

PdbsDir <- "/your/path/here/"
ResultsDir <- "/your/path/here/Results/"

# With singleresidue mode ----

# Calculate frustration of the frames along the dynamics in PdbsDir. With Gifs. 
Dynamic_sing <- dynamic_frustration(PdbsDir = PdbsDir, ResultsDir = ResultsDir, 
                                    GIFs = TRUE, Mode = "singleresidue")

# Analyzing residue Resno
Dynamic_sing <- dynamic_res(Dynamic = Dynamic_sing, Resno = 32, Chain = "A")

# Visualization
plot_dynamic_res(Dynamic = Dynamic_sing, Resno = 32, Chain = "A")

# With configurational or mutational modes ----

# Calculate frustration of frames throughout the dynamics in PdbDir, calculated for OrderList. No gifs. Analyzing residue Resno
OrderList <- c("pdb2.pdb", "pdb14.pdb", "pdb26.pdb")
Dynamic_conf <- dynamic_frustration(PdbsDir = PdbsDir, OrderList = OrderList, 
                                    ResultsDir = ResultsDir, GIFs = FALSE)

# Analyzing residue Resno
Dynamic_conf <- dynamic_res(Dynamic = Dynamic_conf, Resno = 32, Chain = 'A', Graphics = TRUE)

#Visualisations
plot_dynamic_res(Dynamic = Dynamic_conf, Resno = 32, Chain = "A") 
plot_dynamic_res_5Adens_proportion(Dynamic = Dynamic_conf, Resno = 32, Chain = "A")
