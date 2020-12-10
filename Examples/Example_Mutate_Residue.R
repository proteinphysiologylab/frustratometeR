library(frustratometeR)

# With frustration index at the contact level (configurational or mutational) ----

# Obtain Pdb frustration object
Pdb_conf <- calculate_frustration(PdbID = "1nfi", Chain = "F")

# Calculate local frustration for 20 amino acid variants for residue 178 in chain A
Pdb_conf <- mutate_res(Pdb = Pdb_conf, Resno = 178, Chain = "F")

# View results
plot_mutate_res(Pdb = Pdb_conf, Resno = 178, Chain = "F")

# With frustration index at the single-residue level ----

# Obtain Pdb frustration object
Pdb_sing <- calculate_frustration(PdbID = "1nfi", Chain = "F", Mode = "singleresidue")

# Calculate local frustration for 20 amino acid variants for residue 178 in chain A
Pdb_sing <- mutate_res(Pdb = Pdb_sing, Resno = 178, Chain = "F")

# View results
plot_delta_frus(Pdb = Pdb_sing, Resno = 178, Chain = "F")