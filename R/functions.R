#get_os----
#' @title Get OS
#'
#' @description Gets the operating system on which FrustratometeR is running
#'
#' @return Character string indicating the operating system
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
#replace_Expr----
#' @title Expression replacement in a file
#'
#' @description Search and replace pattern by replacement in the file lines.
#'
#' @param pattern Pattern string to replace.
#' @param replacement The character string to be replaced.
#' @param file Full path of the file where the pattern will be replaced.
#' 
#' @return File with pattern replaced by replacement.
replace_Expr <- function(pattern, replacement, file){

  f <- file(file, open = "r")
  document <- readLines(f)
  close(f)
  document <- gsub(pattern, replacement, document)
  writeLines(document, con = file)
}
#XAdens----
#' @title X Adens.
#'
#' @description Calculate the proportion of each type of contact (neutral, highly minimally frustrated) of each residue Pdb within a sphere of radius ratio
#' Is performed by obtaining information from the "tertiary_frustration.dat" file (intermediate processing).
#' 
#' @param Pdb Pdb frustration object.
#' @param ratio Sphere radius. Type: numeric.
#' 
#' @return File (Mode)_5adens.
XAdens <- function(Pdb, ratio = 5){
  
  options(stringsAsFactors = F)
  CA_xyz <- atom.select(Pdb, elety = "CA", value = T)
  CA_x <- c()
  CA_y <- c()
  CA_z <- c()
  for(i in seq(1, length(CA_xyz$xyz), 3)){
    CA_x <- c(CA_x, CA_xyz$xyz[i])
    CA_y <- c(CA_y, CA_xyz$xyz[i + 1])
    CA_z <- c(CA_z, CA_xyz$xyz[i + 2])
  }
  Conts_coords <- read.table(paste(Pdb$JobDir, "tertiary_frustration.dat", sep = ""))[, c(1, 2, 5, 6, 7, 8, 9, 10, 19)]
  Positions <- Pdb$equivalences[, 3]
  ResChain <- Pdb$equivalences[, 1]
  
  vps <- data.frame(Conts_coords[, 1], Conts_coords[, 2],
                    (Conts_coords[, 6] + Conts_coords[, 3]) / 2.0,
                    (Conts_coords[, 7] + Conts_coords[, 4]) / 2.0,
                    (Conts_coords[, 8] + Conts_coords[, 5]) / 2.0,
                    Conts_coords[, 9])
  
  write.table(paste(Pdb$JobDir, Pdb$PdbBase, ".pdb.vps", sep = ""), col.names = F, row.names = F, x = vps)
  
  if(!file.exists(paste(Pdb$JobDir, Pdb$PdbBase, ".pdb.", Pdb$Mode, "_5adens", sep = "")))  
      file.create(paste(Pdb$JobDir, Pdb$PdbBase, ".pdb.", Pdb$Mode, "_5adens", sep = ""))
  
  write("Res ChainRes Total nHighlyFrst nNeutrallyFrst nMinimallyFrst relHighlyFrustrated relNeutralFrustrated relMinimallyFrustrated", 
        file = paste(Pdb$JobDir, Pdb$PdbBase, ".pdb_", Pdb$Mode, "_5adens", sep = ""), append = TRUE)
  for(i in 1:length(CA_x)){
    total_density <- nrow(vps[sqrt((CA_x[i] - vps[, 3]) ^ 2 + (CA_y[i] - vps[, 4]) ^ 2 + (CA_z[i] - vps[, 5]) ^ 2 ) < ratio, ])
    highly_frustrated <- nrow(vps[sqrt((CA_x[i] - vps[, 3]) ^ 2 + (CA_y[i] - vps[, 4]) ^ 2 + (CA_z[i] - vps[, 5]) ^ 2 ) < ratio & vps[, 6] <= (-1), ])
    neutral_frustrated <- nrow(vps[sqrt((CA_x[i] - vps[, 3]) ^ 2 + (CA_y[i] - vps[, 4]) ^ 2 + (CA_z[i] - vps[, 5]) ^ 2 ) < ratio & (vps[, 6] > (-1) & vps[, 6] < (0.78)), ])
    minimally_frustrated <- nrow(vps[sqrt((CA_x[i] - vps[, 3]) ^ 2 + (CA_y[i] - vps[, 4]) ^ 2 + (CA_z[i] - vps[, 5]) ^ 2 ) < ratio & vps[, 6] >= 0.78, ])
    relHighlyFrustratedDensity <- 0
    relNeutralFrustratedDensity <- 0
    relMinimallyFrustratedDensity <- 0
    
    if(total_density > 0)
    {
      relHighlyFrustratedDensity = highly_frustrated / total_density
      relNeutralFrustratedDensity = neutral_frustrated / total_density
      relMinimallyFrustratedDensity = minimally_frustrated / total_density
    }
    write(paste(Positions[i], ResChain[i], total_density, highly_frustrated, neutral_frustrated, minimally_frustrated,
                relHighlyFrustratedDensity, relNeutralFrustratedDensity, relMinimallyFrustratedDensity), 
                file = paste(Pdb$JobDir, Pdb$PdbBase, ".pdb_", Pdb$Mode, "_5adens", sep = ""), append = TRUE)
  }
}
#get_frustration----
#' @title Get frustration data.
#'
#' @description Returns the frustration of all Pdb residues, of a specific Chain or residue (Resno)
#' By default the complete Pdb frustration table is obtained and returned.
#' 
#' @param Pdb Pdb frustration object obtained by calculate_frustration().
#' @param Resno Specific residue in Pdb. Type: numeric. Default: NULL.
#' @param Chain Specific chain in Pdb. Type: character. Default: NULL.
#' 
#' @return Frustration table. Type: data frame.
#' 
#' @export
get_frustration <- function(Pdb, Resno = NULL, Chain = NULL){
  
  if(!is.null(Chain))
    if(!(Chain %in% unique(Pdb$atom$chain))) stop("Chain not exist!")
  if(!is.null(Resno))
    if(!(Resno %in% unique(Pdb$atom$resno))) stop("Resno not exist!")
      
  frustraTable <- read.table(paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb_", Pdb$Mode, sep = ""), header = T)
  if(Pdb$Mode == "singleresidue" ){
    if(!is.null(Chain)) frustraTable <- frustraTable[frustraTable$ChainRes == Chain, ]
    if(!is.null(Resno)) frustraTable <- frustraTable[frustraTable$Res == Resno, ]
  }
  else if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational"){
    if(!is.null(Chain) & !is.null(Resno)) frustraTable <- frustraTable[(frustraTable$Res1 == Resno & frustraTable$ChainRes1 == Chain) |
                                                                         (frustraTable$Res2 == Resno & frustraTable$ChainRes2 == Chain), ]
    else if(!is.null(Chain))  frustraTable <- frustraTable[frustraTable$ChainRes1 == Chain | frustraTable$ChainRes2 == Chain, ]
    else if(!is.null(Resno))  frustraTable <- frustraTable[frustraTable$Res1 == Resno | frustraTable$Res2 == Resno, ]
  }

  print(paste("Frustration obtained!. Frustration index ", Pdb$Mode, sep = ""))
  
  return(frustraTable)
}
#get_frustration_dynamic----
#' @title Get frustration dynamic.
#'
#' @description Obtains and returns the table of frustration of a specific chain and residue in a complete dynamic or in the indicated frames.
#' The frustration of a specific residue must have previously been calculated using dynamic_res().
#' 
#' @param Dynamic Dynamic frustration object obtained by dynamic_frustration().
#' @param Resno Specific residue. Type: numeric.
#' @param Chain Specific chain. Type: character.
#' @param Frames Specific frames. Type: vector.
#' 
#' @return Frustration table. Type: data frame.
#' 
#' @export
get_frustration_dynamic <- function(Dynamic, Resno, Chain, Frames = NULL){
  
  if(!is.null(Dynamic$ResiduesDynamic[[Chain]])){
    if(is.null(Dynamic$ResiduesDynamic[[Chain]][[paste("Res_", Resno, sep = "")]]))
      stop(paste("No analysis to ", Resno, " residue from ", Chain, " chain", sep = ""))
  }
  else  stop(paste("No analysis ", Chain, " chain", sep = ""))
  
  frustraTable <- read.table(Dynamic$ResiduesDynamic[[Chain]][[paste("Res_", Resno, sep = "")]], header = T)
  if(!is.null(Frames)){
    frustraTable <- frustraTable[Frames, ]
  }
  
  print(paste("Frustration of residue ", Resno, " from chain ", Chain," obtained!", sep = ""))
  
  return(frustraTable)
}
#get_clusters----
#' @title Get clusters information
#'
#' @description Obtain information about the clusters obtained from detect_dynamic_clusters(), 
#' the name and number of the residues belonging to each cluster indicated in Clusters.
#' 
#' @param Dynamic Dynamic Frustration Object.
#' @param Clusters Indicates the clusters, for example, c(1, 2, 3), clusters 1, 2 and 3. Default: "all".
#' 
#' @return Data frame containing name, number and cluster of each residue belonging to Clusters
#' 
#' @export
get_clusters <- function(Dynamic, Clusters = "all"){
  
  if(is.null(Dynamic$Clusters[["Graph"]]))
    stop("Cluster detection failed, run detect_dynamic_clusters()")
  
  clusterData <- c()
  clusterData <- cbind(substr(V(Dynamic$Clusters$Graph)$name, 0, 3),
                       substr(V(Dynamic$Clusters$Graph)$name, 5, length(V(Dynamic$Clusters$Graph)$name)),
                       Dynamic$Clusters$LeidenClusters$cluster)
  clusterData <- as.data.frame(clusterData)
  colnames(clusterData) <- c("AA", "Res", "Cluster")
  clusterData$Res <- as.numeric(clusterData$Res)
  clusterData$Cluster <- as.numeric(clusterData$Cluster)
  
  if(Clusters[1] != "all"){
    clusterData <- clusterData[clusterData$Cluster %in% Clusters, ]
  }
  
  return(clusterData)
}

#pdb_equivalences----
#' @title Pdb Equivalences.
#'
#' @description Internal function that produces auxiliary files in the execution of the script, numerical equivalences between the Pdb and its sequence.
#'
#' @param Pdb Pdb frustration object.
#' 
#' @return Pdb equivalences in matrix and file "pdb_equivalences.txt".
pdb_equivalences <- function(Pdb){
  
  existing_res <- unique(cbind(Pdb$atom$chain[which(Pdb$atom$type == "ATOM")], 
                               Pdb$atom$resno[which(Pdb$atom$type == "ATOM")], 
                               Pdb$atom$resid[which(Pdb$atom$type == "ATOM")]))
  equivalences <- cbind(existing_res[, 1],
                        1:length(existing_res[, 1]),
                        existing_res[, 2],
                        existing_res[, 3])
  write.table(equivalences, file = paste(Pdb$JobDir, "/", Pdb$PdbBase, ".pdb_equivalences.txt", sep = ""),
              quote = F, col.names = F, row.names = F, sep = "\t")
  
  return(equivalences)
}
#check_backbone_complete----
#' @title Check Backbone Complete.
#'
#' @description Checks the backbone of a given protein structure to be processed by the frustratometeR pipeline.
#'
#' @param Pdb Pdb frustration object.
#' 
#' @return Flag indicating if the backbone should be completed.
check_backbone_complete <- function(Pdb){
  # Select backbone atoms
  backpdb <- atom.select(Pdb, elety = c("N", "CA", "C", "O", "CB"), value = TRUE)
  # Only for those atoms with ATOM coords
  atom_res <- backpdb$atom[backpdb$atom$type == "ATOM", ]
  # Check there are no missing atoms in the backbone
  Complete <- length(which(atom_res$elety == "CA")) == dim(Pdb$equivalences)[1] & 
              length(which(atom_res$elety == "O")) == dim(Pdb$equivalences)[1] & 
              (length(which(atom_res$elety == "CB")) + length(which(Pdb$equivalences[, 4] == "GLY"))) == dim(Pdb$equivalences)[1] &
              length(which(atom_res$elety == "N")) == dim(Pdb$equivalences)[1]

  return(Complete)
}
#complete_backbone----
#' @title Complete Backbone.
#'
#' @description Completes the backbone of a given protein structure to be processed by the frustratometeR pipeline.
#' 
#' @param Pdb Pdb Frustration object.
#' 
#' @return Pdb Frustration object with backbone completed and flag indicated if the backbone was completed.
complete_backbone <- function(Pdb){
  Completed <- FALSE
  if(!check_backbone_complete(Pdb))
  {
      system(paste("python3 ", Pdb$scriptsDir, "/MissingAtoms.py ",  Pdb$JobDir, Pdb$PdbBase, ".pdb", sep=""))
      system(paste("mv ", Pdb$JobDir, Pdb$PdbBase, ".pdb_completed", " ", Pdb$JobDir, Pdb$PdbBase, ".pdb", sep=""))
      Completed<-TRUE
  }
  return(Completed)
}
#calculate_frustration----
#' @title Calculate Frustration.
#'
#' @description Calculates local energetic frustration for a protein structure.
#'
#' @param PdbFile Full path of the file containing the protein structure. Default: NULL.
#' @param PdbID PdbID of the protein structure. Default: NULL.
#' @param Chain Chain of the protein structure. Type: character. Default: NULL.
#' @param Mode Local frustration index to be calculated (configurational, mutational, singleresidue). Type: string. Default: configurational.
#' @param Electrostatics_K K constant to use in the electrostatics Mode. Type: numeric. Default: NULL (no electrostatics is considered).
#' @param SeqDist Sequence at which contacts are considered to interact (3 or 12). Type: numeric. Default: 12.
#' @param Graphics The corresponding graphics are made (TRUE or FALSE). Type: bool. Default: TRUE.
#' @param Visualization Make visualizations, including pymol. Type: bool. Default: TRUE.
#' @param ResultsDir Path to the folder where results will be stored. If not specified, it will be stored in the directory returned by tempdir(). Default: NULL.
#' 
#' @return Pdb frustration object.
#'
#' @importFrom bio3d pdbsplit basename.pdb read.pdb write.pdb get.pdb
#' @export
calculate_frustration <- function(PdbFile = NULL, PdbID = NULL, Chain = NULL, Electrostatics_K = NULL, SeqDist = 12,
                                  Mode = "configurational", Graphics = TRUE, Visualization = TRUE, ResultsDir = NULL){
  
  if(is.null(ResultsDir)) 
    ResultsDir <- paste(tempdir(), "/", sep = "")
  tempfolder <- tempdir()
  #system(paste0("rm -r ", tempfolder, "/*"))
  
  #Chain
  if(is.null(Chain))  boolsplit = F
  else  boolsplit = T
  
  if(is.null(PdbFile)){
    #Get URL and download
    cat("-----------------------------Download files-----------------------------\n")
    pdbURL <- get.pdb(id = PdbID, split = boolsplit, URLonly = TRUE)
    system(paste("wget --no-check-certificate -P ", tempfolder, pdbURL, " -q --progress=bar:force:noscroll --show-progress", sep = ' '))
    
    if(!boolsplit){
      PdbFile = paste(tempfolder, "/", PdbID, ".pdb", sep = "")
    }else{
      #Chain splitting
      pdbsplit(pdb.files = paste(tempfolder, "/", PdbID, ".pdb", sep = ""), path = file.path(tempfolder, "/split_chain"))
      #Nonexistent Chain exception
      if(!file.exists(paste(tempfolder, "/split_chain/", PdbID, "_", Chain, ".pdb", sep = ""))) stop("Nonexistent Chain")
      PdbFile = paste(tempfolder, "/split_chain/", PdbID, "_", Chain, ".pdb", sep = "")
    }
  }
  else{
    if(boolsplit){
      pdbsplit(pdb.files = PdbFile, path = file.path(tempfolder, "/split_chain"))
      #Nonexistent Chain exception
      if(!file.exists(paste(tempfolder, "/split_chain/", basename.pdb(PdbFile), "_", Chain, ".pdb", sep = ""))){ stop("Nonexistent Chain") }
      PdbFile = paste(tempfolder, "/split_chain/", basename.pdb(PdbFile), "_", Chain, ".pdb", sep = "")
    }
  }
  
  PdbBase <- basename.pdb(PdbFile)
  JobDir <- paste(ResultsDir, PdbBase, ".done/", sep = "")
  
  #Creates JobDir
  if(!dir.exists(JobDir))  dir.create(JobDir)
  system(paste("cp ", PdbFile, " ", JobDir, sep = ""))
  
  #Set job directory
  setwd(JobDir)
  
  PdbFile <- paste(JobDir, PdbBase, ".pdb", sep = "")

  cat("-----------------------------Filtering-----------------------------\n")
  # Here we filter out DNA chains
  system(paste("awk '{if($0!~/DT/ && $0!~/DG/ && $0!~/DA/ && $0!~/DC/ && $4!=",'"',"T",'"'," && $4!=",'"',"G",'"'," && $4!=",'"',"A",'"'," && $4!=",'"',"C",'"',") print }' ", PdbFile , " > ", JobDir, "/aux", sep=""))
  system(paste("mv ", JobDir, "/aux ", PdbFile, sep = ""))

  # We read the pdb file and ignore alternative conformations, aminoacids and only read ATOM lines
  Pdb <- read.pdb(PdbFile, ATOM.only = T, rm.alt = T, rm.insert = T)

  #Fix chain NA
  Pdb$atom$chain[is.na(Pdb$atom$chain[])] <- "A"

  #MSE to MET filter
  Pdb$atom$type[which(Pdb$atom$resid == "MSE")] <- "ATOM"
  Pdb$atom$resid[which(Pdb$atom$resid == "MSE")] <- "MET"

  Pdb[["PdbBase"]] <- PdbBase
  Pdb[["JobDir"]] <- JobDir
  Pdb[["scriptsDir"]] <- paste(find.package("frustratometeR"), "/Scripts", sep = "")
  Pdb[["Mode"]] <- Mode
  Pdb[["PdbPath"]] <- paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb", sep = "")

  # Save equivalences
  equivalences = pdb_equivalences(Pdb)
  Pdb[["equivalences"]] <- equivalences

  write.pdb(Pdb, paste(JobDir, PdbBase, ".pdb", sep = ""))

  if(complete_backbone(Pdb)){
  	Pdb <- read.pdb(paste(JobDir, PdbBase, ".pdb", sep = ""))
  	Pdb <- atom.select(Pdb, type = "ATOM", value = TRUE)
  	write.pdb(Pdb, paste(JobDir, PdbBase, ".pdb", sep = ""))

  	Pdb[["PdbBase"]] <- PdbBase
  	Pdb[["JobDir"]] <- JobDir
  	Pdb[["scriptsDir"]] <- paste(find.package("frustratometeR"), "/Scripts", sep = "")
  	Pdb[["Mode"]] <- Mode
  	Pdb[["PdbPath"]] <- paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb", sep = "")
  	Pdb[["equivalences"]] <- equivalences
  }

  cat("-----------------------------Preparing files-----------------------------\n")
  system(paste("sh ", Pdb$scriptsDir, "/AWSEMFiles/AWSEMTools/PdbCoords2Lammps.sh ", Pdb$PdbBase, " ", Pdb$PdbBase, " ", Pdb$scriptsDir, sep = ""))
  system(paste("cp ", Pdb$scriptsDir, "/AWSEMFiles/*.dat* ", Pdb$JobDir, sep = ""))

  cat("-----------------------------Setting options-----------------------------\n")
  replace_Expr("run\t\t10000", "run\t\t0", paste(Pdb$JobDir, Pdb$PdbBase, ".in", sep = ""))
  replace_Expr("mutational", Pdb$Mode, paste(Pdb$JobDir, "fix_backbone_coeff.data", sep = ""))
  
  if(!is.null(Electrostatics_K))
  {
    cat("-----------------------------Setting electrostatics-----------------------------\n")
    replace_Expr("\\[DebyeHuckel\\]-", "\\[DebyeHuckel\\]", paste(Pdb$JobDir, "fix_backbone_coeff.data", sep = ""))
    replace_Expr("4.15 4.15 4.15", paste(Electrostatics_K, Electrostatics_K, Electrostatics_K, sep = " "),
                 paste(Pdb$JobDir, "fix_backbone_coeff.data", sep = ""))
    print("Setting electrostatics...")
    system(paste("python3 ", Pdb$scriptsDir, "/Pdb2Gro.py ", Pdb$PdbBase, ".pdb ", Pdb$PdbBase, ".pdb.gro; perl ",
                 Pdb$scriptsDir, "/GenerateChargeFile.pl ", Pdb$PdbBase, ".pdb.gro > ", JobDir, "charge_on_residues.dat", sep = ""))
  }

  cat("-----------------------------Calculating-----------------------------\n")
  OperativeSystem <- get_os()
  if(OperativeSystem == "linux"){
     system(paste("cp ", Pdb$scriptsDir, "/lmp_serial_", SeqDist, "_Linux ", Pdb$JobDir, "; chmod +x lmp_serial_", SeqDist,
                 "_Linux ; ./lmp_serial_", SeqDist, "_Linux < ", Pdb$PdbBase, ".in", sep = ""))
  }
  else if(OperativeSystem == "osx"){
    system(paste("cp ", Pdb$scriptsDir, "/lmp_serial_", SeqDist, "_MacOS ", Pdb$JobDir, "; chmod +x lmp_serial_", SeqDist,
                 "_MacOS ; ./lmp_serial_", SeqDist, "_MacOS < ", Pdb$PdbBase, ".in", sep = ""))
  }
 
  system(paste("perl ", Pdb$scriptsDir, "/RenumFiles.pl ", Pdb$PdbBase, " ", Pdb$JobDir, " ", Pdb$Mode, sep = "" ))
  
  if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational"){
    XAdens(Pdb)
  }
  
  cat("-----------------------------Reorganization-----------------------------\n")
  Frustration <- paste(Pdb$JobDir, "FrustrationData", sep = "")
  if (!dir.exists(Frustration))  dir.create(Frustration)
  system(paste("mv *.pdb_", Pdb$Mode, " ", Frustration, sep = ""))
  if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational"){
    system(paste("mv *_", Pdb$Mode, "_5adens ", Frustration, sep = ""))
  }
  system(paste("mv *.pdb ", Frustration, sep = ""))
  
  if(Graphics & Mode != "singleresidue"){
    cat("-----------------------------Images-----------------------------\n")
    Images <- paste(Pdb$JobDir, "Images", sep = "")
    if (!dir.exists(Images))  dir.create(Images)
    plot_5Andens(Pdb, Chain = Chain, Save = T)
    plot_5Adens_proportions(Pdb, Chain = Chain, Save = T)
    plot_contact_map(Pdb, Chain = Chain, Save = T)
  }
  
  if(Visualization & Mode != "singleresidue"){
    cat("-----------------------------Visualizations-----------------------------\n")
    system(paste("perl ", Pdb$scriptsDir, "/GenerateVisualizations.pl ", Pdb$PdbBase, "_", Pdb$Mode,
                 ".pdb_auxiliar ", Pdb$PdbBase, " ", gsub(".$", "", Pdb$JobDir), " ", Pdb$Mode, sep = ""))
    VisualizationDir <- paste(Pdb$JobDir, "VisualizationScrips", sep = "")
    if (!dir.exists(VisualizationDir))  dir.create(VisualizationDir)
    system(paste("cp ", Frustration, "/", Pdb$PdbBase, ".pdb ", VisualizationDir, "/", Pdb$PdbBase, ".pdb", sep = ""))
    system(paste("mv *_", Pdb$Mode, ".pml ", VisualizationDir, sep = ""))
    system(paste("mv *_", Pdb$Mode, ".tcl ", VisualizationDir, sep = ""))
    system(paste("mv *_", Pdb$Mode, ".jml ", VisualizationDir, sep = ""))
    system(paste("cp ", Pdb$scriptsDir, "/draw_links.py ", VisualizationDir, sep = ""))
  }
  
  cat("\n\n****Storage information****\n")
  cat(paste("The frustration data was stored in ", Frustration, "\n", sep = ""))
  if(Graphics & Mode != "singleresidue")
    cat(paste("Graphics are stored in ", Images, "\n", sep = ""))
  if(Visualization & Mode != "singleresidue")
    cat(paste("Visualizations are stored in ", VisualizationDir, "\n", sep = ""))
  
  #Unnecessary files are removed
  system(paste("ls -F1 > output;", sep = ""))
  files <- read.table(paste(JobDir, "output", sep = ""), header = F)
  files[,1] <- as.character(files[, 1])
  for(i in seq(1, dim(files)[1])){
  	finalCharacter <- substr(files[i, 1], nchar(files[i, 1]), nchar(files[i, 1]))
  	if(finalCharacter != '/' & finalCharacter != '.'){
  		system(paste("rm -f ", files[i, 1], sep = ""));
  	}
  }
  setwd(tempfolder)
  
  #We delete temporary files
  if(!is.null(Chain))
  {
  	system(paste("rm -f -R ", tempfolder, "/split_chain", sep = "" ))
  }
  system(paste("rm -f ", tempfolder, "/", PdbID, ".pdb", sep = "" ))

  return(Pdb)
}
#dir_frustration----
#' @title Calculate Directory Frustration.
#'
#' @description Calculate local energy frustration for all protein structures in one directory.
#'
#' @param PdbsDir Directory containing all protein structures. The full path to the file is needed.
#' @param OrderList Orderer list of PDB files to calculate frustration. If it is NULL, frustration is calculated for all PDBs. Type: vector. Default: NULL.
#' @param Chain Chain of the protein structure. Type: character. Default: NULL.
#' @param Electrostatics_K K constant to use in the electrostatics Mode. Type: numeric. Default: NULL (no electrostatics is considered).
#' @param Mode Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational.
#' @param SeqDist Sequence at which contacts are considered to interact (3 or 12). Type: numeric. Default: 12.
#' @param Graphics The corresponding graphics are made (TRUE or FALSE). Type: bool. Default: TRUE.
#' @param Visualization Make visualizations, including pymol. Type: bool. Default: TRUE.
#' @param ResultsDir Path to the folder where results will be stored.
#'
#' @return Files resulting from the frustration of each structure in ResultsDir.
#'
#' @export
dir_frustration <- function(PdbsDir, OrderList = NULL, Chain = NULL, Electrostatics_K = NULL, SeqDist = 12,
                            Mode = "configurational", Graphics = TRUE, Visualization = TRUE, ResultsDir){
  
  setwd(PdbsDir)
  
  CalculationEnabled = TRUE
  
  if(file.exists(paste(ResultsDir, "Modes.log", sep = ""))){
    Modes <- read.table(paste(ResultsDir, "Modes.log", sep = ""), header = FALSE)
    for(i in seq(1, length(Modes[, 1]))){
      if(Modes[i, 1] == Mode)
        CalculationEnabled = FALSE
    }
  }
  
  if(CalculationEnabled){
    if(is.null(OrderList)){
      
      OrderList <- as.vector(list.files(path = PdbsDir))
      listPdb <- data.frame(OrderList)
      
      for(index in seq(1, dim(listPdb)[1])){
        
        if(Mode == "configurational"){
          Pdb = calculate_frustration(PdbFile = paste(PdbsDir, listPdb[index, 1], sep = ""), Chain = Chain, Electrostatics_K = Electrostatics_K, SeqDist = SeqDist,
                                      Mode = "configurational", Graphics = Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
        }else if(Mode == "mutational"){
          Pdb = calculate_frustration(PdbFile = paste(PdbsDir, listPdb[index, 1], sep = ""), Chain = Chain, Electrostatics_K = Electrostatics_K, SeqDist = SeqDist,
                                      Mode = "mutational", Graphics = Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
        }else if(Mode == "singleresidue"){
          Pdb = calculate_frustration(PdbFile = paste(PdbsDir, listPdb[index, 1], sep = ""), Chain = Chain, Electrostatics_K = Electrostatics_K, SeqDist = SeqDist,
                                      Mode = "singleresidue" , Graphics = Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
        }
      }
    }else{
      
      for(index in seq(1, length(OrderList))){
        if(Mode == "configurational"){
          Pdb = calculate_frustration(PdbFile = paste(PdbsDir, OrderList[index], sep = ""), Chain = Chain, Electrostatics_K = Electrostatics_K, SeqDist = SeqDist,
                                      Mode = "configurational", Graphics = Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
        }else if(Mode == "mutational"){
          Pdb = calculate_frustration(PdbFile = paste(PdbsDir, OrderList[index], sep = ""), Chain = Chain, Electrostatics_K = Electrostatics_K, SeqDist = SeqDist, 
                                      Mode = "mutational", Graphics = Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
        }else if(Mode == "singleresidue"){
          Pdb = calculate_frustration(PdbFile = paste(PdbsDir, OrderList[index], sep = ""), Chain = Chain, Electrostatics_K = Electrostatics_K, SeqDist = SeqDist,
                                      Mode = "singleresidue", Graphics = Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
        }
      }
    }
    write(Mode, file = paste(ResultsDir, "Modes.log", sep = ""), append = TRUE)
    
    cat("\n\n****Storage information****\n")
    cat(paste("Frustration data for all Pdb's directory ", PdbsDir, " are stored in ", ResultsDir, "\n", sep = ""))
  }
}
#dynamic_frustration----
#' @title Dynamic Frustration.
#'
#' @description Calculates local energetic frustration for a dynamic.
#'
#' @param PdbsDir Directory containing all protein structures. The full path to the file is needed.
#' @param OrderList Orderer list of PDB files to calculate frustration. If it is NULL, frustration is calculated for all PDBs. Type: vector. Default: NULL.
#' @param Chain Chain of the protein structure. Type: character. Default: NULL.
#' @param Electrostatics_K K constant to use in the electrostatics Mode. Type: numeric. Default: NULL (no electrostatics is considered).
#' @param SeqDist Sequence at which contacts are considered to interact (3 or 12). Type: numeric. Default: 12.
#' @param Mode Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational.
#' @param GIFs If it is TRUE, the contact map gifs and 5 adens proportion of all the frames of the dynamic will be stored, otherwise they will not be stored. Type: bool. Default: FALSE
#' @param ResultsDir Path to the folder where results will be stored. If not specified, it will be stored in the directory returned by tempdir(). Default: NULL.
#'
#' @return Dynamic frustration object
#'
#' @export
dynamic_frustration <- function(PdbsDir, OrderList = NULL, Chain = NULL, Electrostatics_K = NULL, SeqDist = 12, Mode = "configurational", GIFs = FALSE, ResultsDir = NULL){
  
  cat("-----------------------------Object Dynamic Frustration-----------------------------\n")
  if(is.null(ResultsDir)) 
    ResultsDir <- paste(tempdir(), "/", sep = "")
  if(is.null(OrderList)) 
    OrderList <- as.vector(list.files(path = PdbsDir))
  Dynamic <- c()
  Dynamic$Mode <- Mode
  Dynamic$Chain <- Chain
  Dynamic$PdbsDir <- PdbsDir
  Dynamic$ResultsDir <- ResultsDir
  Dynamic$Electrostatics_K <- Electrostatics_K
  Dynamic$OrderList <- OrderList
  Dynamic$SeqDist <- SeqDist
  
  #Pdbs indicated in OrderList found in PdbsDir are frustrated
  cat("-----------------------------Calculating Dynamic Frustration-----------------------------\n")
  dir_frustration(PdbsDir = PdbsDir, OrderList = OrderList, Chain = Chain, Electrostatics_K = Electrostatics_K,
                  SeqDist = SeqDist, Mode = Mode, ResultsDir = ResultsDir)
  
  cat("\n\n****Storage information****\n")
  cat(paste("The frustration of the full dynamic is stored in ", ResultsDir, "\n", sep = ""))
  
  #Gifs
  if(GIFs){
    if (Mode == "configurational" | Mode == "mutational"){
      gif_5adens_proportions(Dynamic)
      gif_contact_map(Dynamic)
    }
  }
  
  return(Dynamic)
}
#dynamic_res----
#' @title Dynamic Res.
#'
#' @description Obtain the local frustration of a specific residue in the dynamics.
#'
#' @param Dynamic Dynamic frustration object.
#' @param Resno Residue specific analyzed. Type: numeric.
#' @param Chain Chain of specific residue. Type: character.
#' @param Graphics If it is TRUE, the graphs corresponding to the residual frustration in the dynamics are made and stored according to the frustration index used. Type: bool. Default: TRUE.
#' 
#' @return Dynamic frustration object adding Mutations attribute for the residue of the indicated chain.
#'
#' @export
dynamic_res <- function(Dynamic, Resno, Chain, Graphics = TRUE){

  PlotsDir <- paste(Dynamic$ResultsDir, "Dynamic_plots_res_", Resno, "_", Chain, sep = "")
  ResultFile <- paste(PlotsDir, "/", Dynamic$Mode, "_Res_", Resno, "_", Chain, sep = "")
  
  setwd(Dynamic$ResultsDir)
  
  if(!dir.exists(PlotsDir)) dir.create(PlotsDir)
  if(file.exists(ResultFile)) file.remove(ResultFile)
  
  cat("-----------------------------Getting frustrated data-----------------------------\n")
  if(Dynamic$Mode == "configurational" | Dynamic$Mode == "mutational"){
    write("Res ChainRes Total nHighlyFrst nNeutrallyFrst nMinimallyFrst relHighlyFrustrated relNeutralFrustrated relMinimallyFrustrated",
          file = ResultFile)
    for(i in 1:length(Dynamic$OrderList)){
      DataRes <- read.table(paste(Dynamic$ResultsDir, basename.pdb(Dynamic$OrderList[i]), ".done/FrustrationData/", basename.pdb(Dynamic$OrderList[i]), ".pdb_", Dynamic$Mode, "_5adens", sep = ""),
                            header = T)
      DataRes <- DataRes[DataRes$Res == Resno & DataRes$ChainRes == Chain, ]
      write.table(DataRes, file = ResultFile, append = T, col.names = F, row.names = F, quote = F)
    }
  }
  else{
    write("Res ChainRes DensityRes AA NativeEnergy DecoyEnergy SDEnergy FrstIndex", file = ResultFile)
    for(i in seq(1, length(Dynamic$OrderList))){
      DataRes <- read.table(paste(Dynamic$ResultsDir, basename.pdb(Dynamic$OrderList[i]), ".done/FrustrationData/", basename.pdb(Dynamic$OrderList[i]), ".pdb_singleresidue", sep = ""),
                            header = T)
      DataRes <- DataRes[DataRes$Res == Resno & DataRes$ChainRes == Chain, ]
      write.table(DataRes, file = ResultFile, append = T, col.names = F, row.names = F, quote = F)
    }
  }
  Dynamic$ResiduesDynamic[[Chain]][[paste("Res_", Resno, sep = "")]] <- ResultFile
  
  cat("\n\n****Storage information****\n")
  cat(paste("The frustration of the residue is stored in ", ResultFile, "\n", sep = ""))
  #Graphics
  if(Graphics){
    plot_dynamic_res(Dynamic = Dynamic, Resno = Resno, Chain = Chain, Save = T)
    if(Dynamic$Mode != "singleresidue")
      plot_dynamic_res_5Adens_proportion(Dynamic = Dynamic, Resno = Resno, Chain = Chain, Save = T)
  }
  
  return(Dynamic)
}
#mutate_res----
#' @title Frustration of mutated residue
#'
#' @description Calculate the local energy frustration for each of the 20 residual variants in the Resno position and Chain chain.
#' Use the frustration index indicated in the Pdb object.
#' 
#' @param Pdb Pdb frustration object.
#' @param Chain Chain of the residue to be mutated. Type: character. Default: NULL.
#' @param Resno Resno of the residue to be mutated. Type: numeric. Default: NULL.
#' @param Split Split that you are going to calculate frustration. If it is TRUE specific string, if it is FALSE full complex. Default: TRUE.
#' @param Method Method indicates the method to use to perform the mutation (Threading or Modeller). Default: Threading
#' 
#' @return Returns Pdb frustration object with corresponding Mutation attribute
#' 
#' @importFrom msa seqbind seqaln
#' @importFrom bio3d write.pdb atom.select get.seq read.fasta
#' 
#' @export
mutate_res <- function(Pdb, Resno, Chain, Split = TRUE, Method = "Threading"){
  
  #Exceptions
  if(length(atom.select(Pdb, resno = Resno, chain = Chain, elety = "CA")$atom) == 0)
    stop("Resno of chain not exist")
  else if(is.na(atom.select(Pdb, resno = Resno, elety = "CA", value = TRUE)$atom$chain[1]))
    Chain = "A"
  if(Split == FALSE & Method == "Modeller")  stop("Complex Modeling not available")
  
  #Output file
  FrustraMutFile <- paste(Pdb$JobDir, "MutationsData/", Pdb$Mode, "_Res", Resno, "_", Method, "_", Chain, ".txt", sep = "")
  if(!dir.exists(paste(Pdb$JobDir, "MutationsData", sep = "")))  dir.create(paste(Pdb$JobDir, "MutationsData", sep = ""))
  setwd(paste(Pdb$JobDir, "MutationsData", sep = ""))
  if(file.exists(FrustraMutFile)){
    file.remove(FrustraMutFile)
    file.create(FrustraMutFile)
  }
  if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational")
    write("Res1 Res2 ChainRes1 ChainRes2 AA1 AA2 FrstIndex FrstState", file = FrustraMutFile, append = TRUE)
  else if(Pdb$Mode == "singleresidue")  write("Res ChainRes AA FrstIndex", file = FrustraMutFile, append = TRUE)  
  
  #Data types, files
  colClasses <- c()
  if(Pdb$Mode == "singleresidue")
    colClasses <- c("integer", "character","numeric","character",
                    "numeric", "numeric", "numeric", "numeric")
  else if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational")
    colClasses <- c("integer", "integer", "character", "character",
                    "numeric", "numeric", "character", "character",
                    "numeric", "numeric", "numeric", "numeric",
                    "character", "character")
  
  if(Method == "Threading"){
    
    AAvector <- c('LEU', 'ASP', 'ILE', 'ASN', 'THR', 'VAL', 'ALA', 'GLY', 'GLU', 'ARG',
                  'LYS', 'HIS', 'GLN', 'SER', 'PRO', 'PHE', 'TYR', 'MET', 'TRP', 'CYS')
    Glycine <- FALSE
    
    #Indices of all the atoms of the residue to mutate
    indexTotales <- atom.select(Pdb, chain = Chain, resno = Resno)
    #If it were glycine
    indexBackboneGly <- atom.select(Pdb, chain = Chain, resno = Resno, elety = c("N", "CA", "C", "O"))
    
    #Check if it is glycine
    if(Pdb$atom$resid[indexTotales$atom[1]] == "GLY") Glycine = TRUE
    else Glycine = FALSE
    
    #Backbone indices
    if(Glycine) indexBackbone <- atom.select(Pdb, chain = Chain, resno = Resno, elety = c("N", "CA", "C", "O"))
    else indexBackbone <- atom.select(Pdb, chain = Chain, resno = Resno, elety = c("N", "CA", "C", "O", "CB"))
    
    rename <-c()
    #It is mutated by 20 AA
    for(AA in AAvector){
      
      cat("\n-----------------------------Getting variant ", AA, "-----------------------------\n")
      
      PdbMut <- Pdb
      #If AA is equal to the residue it is not necessary to mutate and neither is it glycine
      if(AA != Pdb$atom$resid[indexTotales$atom[1]] & !Glycine){
        #if the residue to be inserted is not glycine, insert backbone with CB
        if(AA != "GLY"){
          diffAtom <- setdiff(indexTotales$atom, indexBackbone$atom)
          diffxyz <- setdiff(indexTotales$xyz, indexBackbone$xyz)
        }
        #if the residue to be inserted is glycine, a backbone without CB is inserted
        else{
          diffAtom <- setdiff(indexTotales$atom, indexBackboneGly$atom)
          diffxyz <- setdiff(indexTotales$xyz, indexBackboneGly$xyz)
        }
        #If the previous subtraction is not empty, the corresponding atoms and coordinates are removed
        if(length(diffAtom) > 0)  PdbMut$atom <- PdbMut$atom[-diffAtom, ]
        if(length(diffxyz) > 0) PdbMut$xyz <- PdbMut$xyz[-diffxyz]
      }
      
      #Residues are renamed
      rename <- atom.select(PdbMut, chain = Chain, resno = Resno)
      PdbMut$atom$resid[rename$atom] <- AA
      
      #Muted PDB is saved
      if(Split == TRUE) write.pdb(PdbMut, paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, ".pdb", sep = ""))
      else  write.pdb(PdbMut, paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""))
      
      cat("----------------------------Calculating frustration-----------------------------\n")
      if(Split == TRUE) calculate_frustration( PdbFile = paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, ".pdb", sep = ""),
                                               Mode = Pdb$Mode, ResultsDir = Pdb$JobDir, Graphics = F, Visualization = F, Chain = Chain)
      else  calculate_frustration(PdbFile = paste(Pdb$JobDir, "/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""), 
                                  Mode = Pdb$Mode, ResultsDir = Pdb$JobDir, Graphics = FALSE, Visualization = F)
      
      cat("----------------------------Storing-----------------------------\n")
      if(Pdb$Mode == "singleresidue"){
        system(paste("mv ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/FrustrationData/", Pdb$PdbBase, "_", 
                     Resno, "_", AA, "_", Chain, ".pdb_singleresidue ", Pdb$JobDir, "MutationsData/", sep = ""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_singleresidue", sep = ""),
                                   header = T, colClasses = colClasses)
        frustraTable <- frustraTable[frustraTable$ChainRes == Chain & frustraTable$Res == Resno, c(1, 2, 4, 8)]
        write(paste(frustraTable[, 1], frustraTable[, 2], frustraTable[, 3], frustraTable[, 4]), file = FrustraMutFile, append = TRUE)
      }else if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational"){
        system(paste("mv ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/FrustrationData/", Pdb$PdbBase, "_", 
                     Resno, "_", AA, "_", Chain, ".pdb_", Pdb$Mode, " ", Pdb$JobDir, "MutationsData/", sep = ""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_", Pdb$Mode, sep = ""),
                                   header = T, colClasses = colClasses)
        frustraTable <- frustraTable[frustraTable$ChainRes1 == Chain & frustraTable$Res1 == Resno |
                                     frustraTable$ChainRes2 == Chain & frustraTable$Res2 == Resno,
                                     c(1, 2, 3, 4, 7, 8, 12, 14)]
        write(paste(frustraTable[, 1], frustraTable[, 2], frustraTable[, 3], frustraTable[, 4],
                    frustraTable[, 5], frustraTable[, 6], frustraTable[, 7], frustraTable[, 8]),
                    file = FrustraMutFile, append = TRUE)
      }
      
      #Unnecessary files are removed
      file.remove(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_", Pdb$Mode, sep = ""))
      system(paste("rm -R ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/", sep = ""))
      system(paste("cd ", Pdb$JobDir, " ; rm *pdb", sep = ""))
    }
  }
  else if(Method == "Modeller"){
    
    if(!requireNamespace("msa", quietly = TRUE)){
      stop("Please install msa package from Bioconductor to continue")
    }
    else library(msa)
    
    AAvector <- c('L', 'D', 'I', 'N', 'T', 'V', 'A', 'G', 'E', 'R', 'K', 'H', 'Q', 'S', 'P', 'F', 'Y', 'M', 'W', 'C');
    
    cat("-----------------------------Getting sequence-----------------------------\n")
    sequence <- get.seq(ids = c(Pdb$PdbBase), db = "PDB", outfile = paste(Pdb$JobDir, "seqs.fasta", sep = ""))
    fasta <- read.fasta(paste(Pdb$JobDir, "seqs.fasta", sep = ""))
    
    Seq <- c()
    rowname <- c()
    for(i in 1:length(fasta$id)){
      Seq <- rbind(Seq, fasta$ali[i, ])
      rowname <- c(rowname, paste(substr(fasta$id[i], nchar(fasta$id[i]), nchar(fasta$id[i]))))
    }
    
    rownames(Seq) <- toupper(rowname)
    cat("-----------------------------Equivalences-----------------------------\n")
    SeqPdb <- as.data.frame(pdbseq(Pdb, atom.select(Pdb, chain = Chain, elety = "CA")))
    SeqPdb <- cbind(SeqPdb, as.numeric(row.names(SeqPdb)), 1:length(SeqPdb[, 1]))
    colnames(SeqPdb) <- c("AA", "resno", "index")

    aln <- seqbind(Seq[Chain, 1:table(Seq[Chain,] == "-")[1]], SeqPdb$AA)
    align <- seqaln(aln, outfile = tempfile(), exefile = "msa")
    SeqGap <- align$ali["seq2", ]
    
    SeqGap <- cbind(SeqGap, rep(0, length(SeqGap)), 1:length(SeqGap))
    
    j <- 1
    for(i in 1:length(SeqGap[,1])){
      if(SeqGap[i, 1] != "-" && SeqGap[i, 1] == SeqPdb$AA[j]){
        SeqGap[i, 2] <- SeqPdb$resno[j]
        j <- j + 1
      }
    }
    pos <- as.numeric(SeqGap[SeqGap[, 2] == Resno, 3])
    
    system(paste("cp ", Pdb$scriptsDir, "/align2d.py ", Pdb$JobDir, sep = ""))
    system(paste("cp ", Pdb$scriptsDir, "/make_ali.py ", Pdb$JobDir, sep = ""))
    system(paste("cp ", Pdb$scriptsDir, "/model-single.py ", Pdb$JobDir, sep = ""))
    system(paste("cp ", Pdb$PdbPath," ", Pdb$JobDir, sep = ""))

    for (AA in AAvector){
      cat("\n-----------------------------Getting variant ", AA, "-----------------------------\n")
      if(Split){
        write(">Modelo", file = paste(Pdb$JobDir, "Modelo.fa", sep = ""))
        SeqMut = Seq[Chain, 1:table(Seq[Chain,] == "-")[1]]
        SeqMut[pos] <- AA
        write(paste(as.vector(SeqMut), collapse = ""), file = paste(Pdb$JobDir, "Modelo.fa", sep = ""), append = TRUE)
      }
      
      cat("-----------------------------Aligning-----------------------------\n")
      system(paste("cd ", Pdb$JobDir, " ;python3 make_ali.py Modelo", sep = ""))
      if(Split) system(paste("cd ", Pdb$JobDir, " ;python3 align2d.py ", Pdb$PdbBase, " Modelo ", Chain, sep = ""))
  
      cat("-----------------------------Modeling-----------------------------\n")
      system(paste("cd ", Pdb$JobDir, " ;python3 model-single.py ", Pdb$PdbBase, " Modelo", sep = ""))
      system(paste("cd ", Pdb$JobDir, " ;mv Modelo.B99990001.pdb ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""))
      system(paste("cd ", Pdb$JobDir, " ;rm *D00000001 *ini *rsr *sch *V99990001 *ali *pap *fa"))
      
      cat("----------------------------Calculating frustration-----------------------------\n")
      calculate_frustration(PdbFile = paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""),
                            Mode = Pdb$Mode, ResultsDir = Pdb$JobDir, Graphics = F, Visualization = F)
      
      cat("----------------------------Storing-----------------------------\n")
      if(Pdb$Mode == "singleresidue"){
        system(paste("mv ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/FrustrationData/", Pdb$PdbBase,
                     "_", Resno, "_", AA, "_", Chain, ".pdb_singleresidue ", Pdb$JobDir, "MutationsData/", sep = ""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_singleresidue", sep = ""),
                                   header = T, colClasses = colClasses)
        frustraTable <- frustraTable[frustraTable$Res == pos, c(1, 2, 4, 8)]
        write(paste(frustraTable[, 1], frustraTable[, 2], frustraTable[, 3], frustraTable[, 4]), file = FrustraMutFile, append = TRUE)
      }else if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational"){
        system(paste("mv ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/FrustrationData/", Pdb$PdbBase,
                     "_", Resno, "_", AA, "_", Chain, ".pdb_", Pdb$Mode, " ", Pdb$JobDir, "MutationsData/", sep = ""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_", Pdb$Mode, sep = ""),
                                   header = T, colClasses = colClasses)
        frustraTable <- frustraTable[frustraTable$Res1 == pos |
                                     frustraTable$Res2 == pos,
                                     c(1, 2, 3, 4, 7, 8, 12, 14)]
        write(paste(frustraTable[, 1], frustraTable[, 2], frustraTable[, 3], frustraTable[, 4],
                    frustraTable[, 5], frustraTable[, 6], frustraTable[, 7], frustraTable[, 8]),
              file = FrustraMutFile, append = TRUE)
      }
      
      #Unnecessary files are removed
      file.remove(paste(Pdb$JobDir,"MutationsData/",Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_", Pdb$Mode, sep=""))
      system(paste("rm -R ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/", sep=""))
      system(paste("cd ", Pdb$JobDir, " ; rm ",Pdb$PdbBase,"_*", sep = ""))
      
    }
    system(paste("cd ", Pdb$JobDir, " ; rm *pdb seqs.fasta *py", sep = ""))
    
    cat("----------------------------Renumbering-----------------------------\n")
    if(Pdb$Mode == "singleresidue"){
      data <- read.table(paste(Pdb$JobDir, "MutationsData/singleresidue_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""), header = T)
      data[, 2] <- Chain
      for(i in 1:length(data[, 1])){
        data[i, 1] <- SeqGap[as.numeric(data[i,1]),2]
      }
      write.table(paste(Pdb$JobDir, "MutationsData/singleresidue_Res", Resno, "_", Method, "_", Chain, ".txt", sep = ""), col.names = T, row.names = F)
    }else if(Pdb$Mode == "configurational"){
      data <- read.table(paste(Pdb$JobDir, "MutationsData/configurational_Res", Resno, "_", Method, "_", Chain, ".txt", sep = ""), header = T)
      data[, c(3, 4)] <- Chain
      for(i in 1:length(data[,1])){
        data[i, 1] <- SeqGap[as.numeric(data[i, 1]), 2]
        data[i, 2] <- SeqGap[as.numeric(data[i, 2]), 2]
      }
      write.table(data, paste(Pdb$JobDir, "MutationsData/configurational_Res", Resno, "_", Method, "_", Chain, ".txt", sep = ""), col.names = T, row.names = F, quote = F)
    }else if(Pdb$Mode == "mutational"){
      data <- read.table(paste(Pdb$JobDir, "MutationsData/mutational_Res", Resno, "_", Method, "_", Chain, ".txt", sep = ""), header = T)
      data[, c(3, 4)] <- Chain
      for (i in 1:length(data[,1])) {
        data[i, 1] <- SeqGap[as.numeric(data[i, 1]), 2]
        data[i, 2] <- SeqGap[as.numeric(data[i, 2]), 2]
      }
      write.table(data, paste(Pdb$JobDir, "MutationsData/mutational_Res", Resno, "_", Method, "_", Chain , ".txt", sep = ""), col.names = T, row.names = F, quote = F)
    }
  }
  
  Pdb$Mutations[[Method]][[paste("Res_", Resno, "_", Chain, sep = "")]] <- data.frame(Method = Method ,Res = Resno, 
                                                                                               Chain = Chain, File = paste(Pdb$JobDir, "MutationsData/", Pdb$Mode, "_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""))

  cat("\n\n****Storage information****\n")
  cat(paste("The frustration of the residue is stored in ", Pdb$JobDir, "MutationsData/mutational_Res", Resno, "_", Method, "_", Chain , ".txt\n", sep = ""))
  
  return(Pdb)
}
#detect_dynamic_clusters----
#' @title Detect Dynamic clusters
#'
#' @description Detects residue modules with similar single-residue frustration dynamics.
#' It filters out the residuals with variable dynamics, for this, it adjusts a loess 
#' model with span = LoessSpan and calculates the dynamic range of frustration and the mean of single-residue frustration.
#' It is left with the residuals with a dynamic frustration range greater than the quantile defined by MinFrstRange and with a mean Mean <(-FiltMean) or Mean> FiltMean.
#' Performs an analysis of main components and keeps Ncp components, to compute the correlation(CorrType) between them and keep the residues that have greater correlation MinCorr and p-value> 0.05.
#' An undirected graph is generated and Leiden clustering is applied with LeidenResol resolution.
#'
#' @param Dynamic Dynamic Frustration Object.
#' @param LoessSpan Parameter α > 0 that controls the degree of smoothing of the loess() function of model fit. Type: numeric. Default: 0.05.
#' @param MinFrstRange Frustration dynamic range filter threshold. 0 <= MinFrstRange <= 1. Type: numeric. Default: 0.8.
#' @param FiltMean Frustration Mean Filter Threshold. FiltMean >= 0. Type: numeric. Default: 0.5.
#' @param Ncp Number of principal components to be used in PCA(). Ncp >= 1. Type: numeric. Default: 10.
#' @param MinCorr Correlation filter threshold. 0 <= MinCorr <= 1. Type: numeric. Default: 0.8.
#' @param LeidenResol Parameter that defines the coarseness of the cluster. LeidenResol > 0. Type: numeric. Default: 1.5.
#' @param CorrType Type of correlation index to compute. Values: "pearson" or "spearman". Type: character. Default: "pearson".
#' 
#' @return Dynamic Frustration Object and its Clusters attribute.
#' 
#' @importFrom FactoMineR PCA
#' @importFrom Hmisc rcorr
#' @importFrom bio3d basename.pdb
#' @importFrom igraph graph_from_adjacency_matrix delete_vertices
#' @importFrom leiden leiden
#' 
#' @export
detect_dynamic_clusters <- function(Dynamic = Dynamic, LoessSpan = 0.05, MinFrstRange = 0.8, FiltMean = 0.5, Ncp = 10, MinCorr = 0.8, LeidenResol = 1.5, CorrType = "pearson"){
  
  if(Dynamic$Mode != "singleresidue")
    stop("This functionality is only available for the singleresidue index, run dynamic_frustration() with Mode = 'singleresidue'")
  if(!(CorrType %in% c("pearson", "spearman")))
    stop("Correlation type(CorrType) indicated isn't available or no exist, indicate 'pearson' or 'spearman'")
  
  libraries <- c("leiden", "dplyr", "igraph", "FactoMineR", "bio3d", "Hmisc", "RColorBrewer")
  for(library in libraries){
    if(!requireNamespace(library, quietly = TRUE)){
      cat(paste("Please install ", library," package to continue", sep = ""))
      if(library == "leiden")
          cat("To import 'leiden' is the next module must first be installed with: 'python3 -m pip install leidenalg'")
    }
    else library(library, character.only = T)
  }
  
  #Loading residues and resno
  ini <- read.table(paste(Dynamic$ResultsDir, basename.pdb(Dynamic$OrderList[1]),".done/FrustrationData/",
                          basename.pdb(Dynamic$OrderList[1]), ".pdb_singleresidue", sep = ""), header = T)
  residues <- ini$AA
  resnos <- ini$Res
  rm(ini)
  
  #Loading data
  frustraData <- c()
  read <- c()
  for(i in 1:length(Dynamic$OrderList)){
    read <- read.table(paste(Dynamic$ResultsDir, basename.pdb(Dynamic$OrderList[i]), ".done/FrustrationData/",
                             basename.pdb(Dynamic$OrderList[i]), ".pdb_singleresidue", sep = ""), header = T)
    frustraData <- cbind(frustraData, read$FrstIndex)
  }
  frustraData <- as.data.frame(frustraData)
  colnames(frustraData) <- paste("frame_", 1:length(Dynamic$OrderList), sep = "")
  rownames(frustraData) <- paste(aa123(residues), "_", resnos, sep = "")
  
  #Model fitting and filter by difference and mean
  diferences <- c()
  means <- c()
  fitted <- c()
  res <- c()
  for (i in 1:length(residues)){
    res <- as.data.frame(cbind(t(frustraData[i, ]), 1:ncol(frustraData)))
    colnames(res) <- c("Frustration","Frames")
    modelo <- loess(Frustration ~ Frames, data = res, span = LoessSpan, degree = 1, family="gaussian")
    fitted <- as.data.frame(cbind(fitted, modelo$fitted))
    diferences <- c(diferences, max(modelo$fitted) - min(modelo$fitted))
    means <- c(means, mean(modelo$fitted))
  }
  
  estadistics <- data.frame(Diferences = diferences, Means = means)
  frustraData <- frustraData[(estadistics$Diferences > quantile(estadistics$Diferences, probs = MinFrstRange)) &
                      (estadistics$Means < -FiltMean | estadistics$Means > FiltMean) , ]
  
  #Principal component analysis
  pca <- PCA(frustraData, ncp = Ncp, graph = F)
  
  Cor <- as.matrix(rcorr(t(pca$ind$coord[,]), type = CorrType))
  Cor[[1]][lower.tri(Cor[[1]], diag = T)] <- 0
  
  Cor[[1]][!(Cor[[1]] < (-MinCorr)) & !(Cor[[1]] > (MinCorr)) | Cor[[3]] > 0.05] <- 0
  net <- graph_from_adjacency_matrix(adjmatrix = Cor[[1]], diag = F, mode = "undirected", weighted = T)
  
  #Leiden Clustering
  cluster <- data.frame(cluster = leiden(net, resolution_parameter = LeidenResol))
  cluster <- as.data.frame(cluster[-as.vector(V(net)[igraph::degree(net) == 0]), ])
  net <- delete_vertices(net, V(net)[igraph::degree(net) == 0])
  
  Dynamic$Clusters[["Graph"]] <- net
  Dynamic$Clusters[["LeidenClusters"]] <- cluster
  Dynamic$Clusters[["LoessSpan"]] <- LoessSpan
  Dynamic$Clusters[["MinFrstRange"]] <- MinFrstRange
  Dynamic$Clusters[["FiltMean"]] <- FiltMean
  Dynamic$Clusters[["Ncp"]] <- Ncp
  Dynamic$Clusters[["MinCorr"]] <- MinCorr
  Dynamic$Clusters[["LeidenResol"]] <- LeidenResol
  Dynamic$Clusters[["Fitted"]] <- fitted
  Dynamic$Clusters[["Means"]] <- means
  Dynamic$Clusters[["Diferences"]] <- diferences
  Dynamic$Clusters[["CorrType"]] <- CorrType
  
  return(Dynamic)
}
