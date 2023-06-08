#to mutate for a single amino acid of interest
#ResToMut=c('L','V') where L in the native residue and V is the mutation
library(frustratometeR)

Pdb_conf <- calculate_frustration(PdbFile = "XXXX.pdb", Mode = "mutational", ResultsDir ="FrustrationResults/", Chain="A")
Pdb_conf <- mutate_res_only(Pdb = Pdb_conf, Resno = 6, Chain = "A",ResToMut=c('L','V'))
plot_mutate_res(Pdb = Pdb_conf, Resno = 6, Chain = "A", Save = TRUE)
