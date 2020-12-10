#Example adapted for the dynamics found in https://github.com/proteinphysiologylab/frustratometer-data/tree/main/IkappaBalfa_dynamics

#Load frustratometeR
library(frustratometeR)

# Defining the order of the structures
OrderList <- c()
for(i in 1:2040){ OrderList <- c(OrderList, paste("612_", i, ".pdb", sep = "")) }

#Calculation of local energy frustration for all frames of dynamics
Dynamic_sing <- dynamic_frustration(PdbsDir = "/home/Dynamic/pdbs/", OrderList = OrderList, ResultsDir = "/home/Dynamic/results/", Mode = "singleresidue")

#Detect dynamic clusters
Dynamic_sing <- detect_dynamic_clusters(Dynamic = Dynamic_sing)

#Visualisations
plot_variable_res_filter(Dynamic = Dynamic_sing)
plot_dynamic_clusters_graph(Dynamic = Dynamic_sing)
plot_res_dynamics(Dynamic = Dynamic_sing, Resno = 85, Chain = "A")
plot_res_dynamics(Dynamic = Dynamic_sing, Resno = 167, Chain = "A")