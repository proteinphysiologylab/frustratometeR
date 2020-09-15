#' Replace expression
#'
#' Searches for a string in a file and replaces
#'
#' @param pattern Pattern string to replace
#' @param replacement Character string to be replaced
#' @param file Full path of the file where it will be replaced
#' @return

replace_Expr <- function(pattern = NULL, replacement = NULL, file = NULL){
  if(is.null(pattern)) stop("pattern parameter not indicated")
  if(is.null(replacement)) stop("replacement parameter not indicated")
  if(is.null(file)) stop("file parameter not indicated")
  
  f <- file(file, open = "r")
  document <- readLines(f)
  close(f)
  document <- gsub(pattern, replacement, document)
  writeLines(document, con = file)
}

#' Replace expression
#'
#' Calculate the proportion by type of contact in a sphere of radius X
#'
#' @param Pdb Pdb object
#' @param ratio Radius of the sphere for calculating proportions
#' @return

XAdens <- function(Pdb = NULL, ratio = 5){
  if(is.null(Pdb)) stop("Pdb parameter not indicated")
  
  options(stringsAsFactors = F)
  CA_xyz <- atom.select(Pdb,elety = "CA",value = T)
  CA_x <- c()
  CA_y <- c()
  CA_z <- c()
  for( i in seq(1,length(CA_xyz$xyz),3)) {
    CA_x <- c(CA_x, CA_xyz$xyz[i])
    CA_y <- c(CA_y, CA_xyz$xyz[i+1])
    CA_z <- c(CA_z, CA_xyz$xyz[i+2])
  }
  Conts_coords <- read.table(paste(Pdb$JobDir, "tertiary_frustration.dat",sep=""))[,c(1,2,5,6,7,8,9,10,19)]
  Positions <- Pdb$equivalences[,3]
  ResChain <- Pdb$equivalences[,1]
  
  vps <- data.frame(Conts_coords[,1],Conts_coords[,2],
                    (Conts_coords[,6]+Conts_coords[,3])/2.0,
                    (Conts_coords[,7]+Conts_coords[,4])/2.0,
                    (Conts_coords[,8]+Conts_coords[,5])/2.0,
                    Conts_coords[,9])
  
  write.table(paste(Pdb$JobDir,Pdb$PdbBase,".pdb.vps",sep=""),col.names=F,row.names=F,x = vps)
  
  if (!file.exists(paste(Pdb$JobDir,Pdb$PdbBase,".pdb.",Pdb$mode,"_5adens",sep="")))  
      file.create(paste(Pdb$JobDir,Pdb$PdbBase,".pdb.",Pdb$mode,"_5adens",sep=""))
  
  write("Res ChainRes Total nHighlyFrst nNeutrallyFrst nMinimallyFrst relHighlyFrustrated relNeutralFrustrated relMinimallyFrustrated", 
        file=paste(Pdb$JobDir,Pdb$PdbBase,".pdb_",Pdb$mode,"_5adens",sep=""), append=TRUE)
  for (i in 1:length(CA_x)) {
    total_density <- nrow(vps[sqrt((CA_x[i] - vps[,3])^2 + (CA_y[i] - vps[,4])^2 + (CA_z[i] - vps[,5])^2 )<ratio,])
    highly_frustrated <- nrow(vps[sqrt((CA_x[i] - vps[,3])^2 + (CA_y[i] - vps[,4])^2 + (CA_z[i] - vps[,5])^2 )<ratio & vps[,6]<=(-1),])
    neutral_frustrated <- nrow(vps[sqrt((CA_x[i] - vps[,3])^2 + (CA_y[i] - vps[,4])^2 + (CA_z[i] - vps[,5])^2 )<ratio & (vps[,6]>(-1)&vps[,6]<(0.78)),])
    minimally_frustrated <- nrow(vps[sqrt((CA_x[i] - vps[,3])^2 + (CA_y[i] - vps[,4])^2 + (CA_z[i] - vps[,5])^2 )<ratio & vps[,6]>=0.78,])
    relHighlyFrustratedDensity <- 0
    relNeutralFrustratedDensity <- 0
    relMinimallyFrustratedDensity <- 0
    
    if(total_density>0)
    {
      relHighlyFrustratedDensity=highly_frustrated/total_density
      relNeutralFrustratedDensity=neutral_frustrated/total_density
      relMinimallyFrustratedDensity=minimally_frustrated/total_density
    }
    write(paste(Positions[i],ResChain[i],total_density,highly_frustrated,neutral_frustrated,minimally_frustrated,
                relHighlyFrustratedDensity,relNeutralFrustratedDensity,relMinimallyFrustratedDensity), 
                file=paste(Pdb$JobDir,Pdb$PdbBase,".pdb_",Pdb$mode,"_5adens",sep=""), append=TRUE)
  }
}

#' Get frustration data
#'
#' Returns a data frame corresponding to the frustration data previously calculated by calculate_frustration
#'
#' @param Pdb Pdb object
#' @param Chain ******************************
#' @param Resno *****************************
#' @return Frustration table

get_frustration <- function(Pdb = NULL, Chain = NULL, Resno = NULL){
  
  if(is.null(Pdb)) stop("PDB parameter not indicated")
  if(!is.null(Chain))
    if(!(Chain %in% unique(Pdb$atom$chain))) stop("Chain not exist!")
  if(!is.null(Resno))
    if(!(Resno %in% unique(Pdb$atom$resno))) stop("Resno not exist!")
      
  frustraTable <- read.table(paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb_", Pdb$mode, sep = ""), header = T)
  if(Pdb$mode == "singleresidue" ){
    if(!is.null(Chain)) frustraTable <- frustraTable[frustraTable$ChainRes == Chain,]
    if(!is.null(Resno)) frustraTable <- frustraTable[frustraTable$Res == Resno,]
  }
  else if(Pdb$mode == "configurational" | Pdb$mode == "mutational"){
    if(!is.null(Chain) & !is.null(Resno)) frustraTable <- frustraTable[(frustraTable$Res1 == Resno & frustraTable$ChainRes1 == Chain) |
                                                                         (frustraTable$Res2 == Resno & frustraTable$ChainRes2 == Chain), ]
    else if(!is.null(Chain))  frustraTable <- frustraTable[frustraTable$ChainRes1 == Chain | frustraTable$ChainRes2 == Chain, ]
    else if(!is.null(Resno))  frustraTable <- frustraTable[frustraTable$Res1 == Resno | frustraTable$Res2 == Resno, ]
  }
  
  print(paste("Frustration obtained!. Frustration index ",Pdb$mode,sep = ""))
  
  return(frustraTable)
}

#' Pdb Equivalences
#'
#' Internal function that produces auxiliar files to run the frustratometer pipeline.
#'
#' @param Pdb Frustration object
#' @return Pdb equivalences

pdb_equivalences <- function(Pdb){
  existing_res <- unique(cbind(Pdb$atom$chain[which(Pdb$atom$type == "ATOM")], Pdb$atom$resno[which(Pdb$atom$type == "ATOM")], Pdb$atom$resid[which(Pdb$atom$type == "ATOM")]))
  equivalences <- cbind(existing_res[,1], 1:length(existing_res[,1]), existing_res[,2],existing_res[,3])
  write.table(equivalences, file = paste(Pdb$JobDir, "/", Pdb$PdbBase, ".pdb_equivalences.txt", sep = ""), quote = F, col.names = F, row.names = F, sep = "\t")
  return(equivalences)
}

#' Check Backbone Complete
#'
#' Checks the backbone of a given protein structure to be processed by the frustratometer pipeline
#'
#' @param Pdb Frustration object
#' @return Flag indicating if the backbone should be completed

check_backbone_complete <- function(Pdb){
  # Select backbone atoms
  backpdb <- atom.select(Pdb, elety=c("N", "CA", "C", "O", "CB"), value=TRUE)
  # Only for those atoms with ATOM coords
  atom_res <- backpdb$atom[backpdb$atom$type=="ATOM",]
  # Check there are no missing atoms in the backbone
  Complete <- length(which(atom_res$elety=="CA"))==dim(Pdb$equivalences)[1] & length(which(atom_res$elety=="O"))==dim(Pdb$equivalences)[1] & (length(which(atom_res$elety=="CB")) + length(which(Pdb$equivalences[,4]=="GLY")))==dim(Pdb$equivalences)[1] & length(which(atom_res$elety=="N"))==dim(Pdb$equivalences)[1]

  return(Complete)
}

#' Complete Backbone
#'
#' Completes the backbone of a given protein structure to be processed by the frustratometer pipeline
#'
#' @param Pdb Frustration object
#' @return Pdb Frustration object with backbone completed and flag indicated if the backbone was completed

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

#' Calculate Frustration
#'
#' Calculates local energetic frustration for a protein structure
#'
#' @param PdbFile File containing the protein structure. The full path to the file is needed
#' @param PdbID PdbID of the protein structure.
#' @param Chain Chain of the protein structure.
#' @param Mode Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#' @param Electrostatics_K K constant to use in the electrostatics mode. Default: NULL (no electrostatics is considered).
#' @param seqdist Sequence at which contacts are considered to interact (3 or 12). Default: 12
#' @param Graphics The corresponding graphics are made (TRUE or FALSE). Default: TRUE
#' @param Visualization Make visualizations, including pymol. Default: TRUE
#' @param ResultsDir Path to the folder where results will be stored.
#' @return Pdb Frustration object
#'
#' @export

calculate_frustration <- function(PdbFile = NULL, PdbID = NULL, Chain = NULL, Electrostatics_K = NULL, seqdist = 12, Mode = "configurational", Graphics = TRUE, Visualization = TRUE, ResultsDir = NULL){
  tempfolder <- tempdir()

  #Chain
  if(is.null(Chain))  boolsplit = F
  else  boolsplit = T
  
  if(is.null(PdbFile))
  {
    #Get URL and download
    print("-----------------------------Download files-----------------------------")
    pdbURL <- get.pdb(id = PdbID, split = boolsplit, URLonly = TRUE)
    system(paste("wget --no-check-certificate -P ", tempfolder, pdbURL, " -q --progress=bar:force:noscroll --show-progress", sep = ' '))
    
    if(!boolsplit)
    {
      PdbFile = paste(tempfolder, "/", PdbID, ".pdb", sep = "")
    }else{
      #Chain splitting
      pdbsplit(pdb.files = paste(tempfolder, "/", PdbID, ".pdb", sep = ""), path = file.path(tempfolder, "/split_chain"))
      #Nonexistent Chain exception
      if(!file.exists(paste(tempfolder, "/split_chain/", PdbID, "_", Chain, ".pdb", sep = ""))){ stop("Nonexistent Chain") }
      PdbFile = paste(tempfolder, "/split_chain/", PdbID, "_", Chain, ".pdb", sep = "")
    }
  }else{
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

  #**************FILTERS**************
  # Here we filter out DNA chains
  system(paste("awk '{if($0!~/DT/ && $0!~/DG/ && $0!~/DA/ && $0!~/DC/ && $4!=",'"',"T",'"'," && $4!=",'"',"G",'"'," && $4!=",'"',"A",'"'," && $4!=",'"',"C",'"',") print }' ", PdbFile , " > ", JobDir, "/aux", sep=""))
  system(paste("mv ", JobDir, "/aux ", PdbFile, sep = ""))

  # We read the pdb file and ignore alternative conformations, aminoacids and only read ATOM lines
  Pdb <- read.pdb(PdbFile, ATOM.only = T, rm.alt = T, rm.insert = T)

  #Fix chain NA
  Pdb$atom$chain[is.na(Pdb$atom$chain[])] <- "A"

  #MSE to MET filter
  Pdb$atom$type[ which( Pdb$atom$resid == "MSE" ) ] <- "ATOM"
  Pdb$atom$resid[ which( Pdb$atom$resid == "MSE" ) ] <- "MET"

  Pdb[["PdbBase"]] <- PdbBase
  Pdb[["JobDir"]] <- JobDir
  Pdb[["scriptsDir"]] <- paste(find.package("frustratometeR"), "/Scripts", sep = "")
  Pdb[["mode"]] <- Mode
  Pdb[["PdbPath"]] <- paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb", sep = "")

  # Save equivalences
  equivalences = pdb_equivalences(Pdb)
  Pdb[["equivalences"]] <- equivalences

  write.pdb(Pdb, paste(JobDir, PdbBase, ".pdb", sep = ""))

  if(complete_backbone(Pdb)){
  	Pdb <- read.pdb(paste(JobDir, PdbBase, ".pdb", sep = ""))
  	Pdb<-atom.select(Pdb, type = "ATOM", value = TRUE)
  	write.pdb(Pdb, paste(JobDir, PdbBase, ".pdb", sep = ""))

  	Pdb[["PdbBase"]] <- PdbBase
  	Pdb[["JobDir"]] <- JobDir
  	Pdb[["scriptsDir"]] <- paste(find.package("frustratometeR"), "/Scripts", sep = "")
  	Pdb[["mode"]] <- Mode
  	Pdb[["PdbPath"]] <- paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb", sep = "")
  	Pdb[["equivalences"]] <- equivalences
  }

  print("-----------------------------Preparing files-----------------------------")
  
  #Prepare the PDB file to get awsem input files, create the workdir and move neccessary files to it.
  system(paste("sh ", Pdb$scriptsDir, "/AWSEMFiles/AWSEMTools/PdbCoords2Lammps.sh ", Pdb$PdbBase, " ", Pdb$PdbBase, " ", Pdb$scriptsDir, sep = ""))
  system(paste("cp ", Pdb$scriptsDir, "/AWSEMFiles/*.dat* ", Pdb$JobDir, sep = ""))

  print("-----------------------------Setting options-----------------------------")

  replace_Expr("run\t\t10000", "run\t\t0", paste(Pdb$JobDir, Pdb$PdbBase, ".in", sep = ""))
  replace_Expr("mutational", Pdb$mode, paste(Pdb$JobDir, "fix_backbone_coeff.data", sep = ""))
  
  if(!is.null(Electrostatics_K))
  {
    print("Setting electrostatics...")
    replace_Expr("\\[DebyeHuckel\\]-", "\\[DebyeHuckel\\]", paste(Pdb$JobDir, "fix_backbone_coeff.data", sep = ""))
    replace_Expr("4.15 4.15 4.15", paste(Electrostatics_K, Electrostatics_K, Electrostatics_K, sep=" "), paste(Pdb$JobDir, "fix_backbone_coeff.data", sep = ""))
    print("Setting electrostatics...")
    system(paste("python3 ", Pdb$scriptsDir, "/Pdb2Gro.py ", Pdb$PdbBase, ".pdb ", Pdb$PdbBase, ".pdb.gro; perl ", Pdb$scriptsDir, "/GenerateChargeFile.pl ", Pdb$PdbBase, ".pdb.gro > ", JobDir, "charge_on_residues.dat", sep = ""))
  }

  print("-----------------------------Calculating-----------------------------")

  
  system(paste("cp ", Pdb$scriptsDir, "/lmp_serial_", seqdist, " ", Pdb$JobDir, "; chmod +x lmp_serial_", seqdist,
               "; ./lmp_serial_", seqdist, " < ", Pdb$PdbBase, ".in", sep = ""))

  print("-----------------------------RenumFiles-----------------------------")
  system(paste("perl ", Pdb$scriptsDir, "/RenumFiles.pl ", Pdb$PdbBase, " ", Pdb$JobDir, " ", Pdb$mode, sep = "" ))
  
  if(Pdb$mode == "configurational" | Pdb$mode == "mutational")
  {
    print("-----------------------------5Adens-----------------------------")
    XAdens(Pdb)
  }
  
  #Directory reorganization and results
  print("-----------------------------Reorganization-----------------------------")
  
  #Frustration Data
  Frustration <- paste(Pdb$JobDir, "FrustrationData", sep = "")
  if (!dir.exists(Frustration))  dir.create(Frustration)
  system(paste("mv *.pdb_", Pdb$mode, " ", Frustration, sep = ""))
  if(Pdb$mode == "configurational" | Pdb$mode == "mutational"){
    system(paste("mv *_", Pdb$mode, "_5adens ", Frustration, sep = ""))
  }
  system(paste("mv *.pdb ", Frustration, sep = ""))
  
  #Images
  if(Graphics & Mode != "singleresidue"){
    Images <- paste(Pdb$JobDir, "Images", sep = "")
    if (!dir.exists(Images))  dir.create(Images)
    plot_5Andens(Pdb, Chain = Chain, Show = F)
    plot_5Adens_proportions(Pdb, Chain = Chain, Show = F)
    plot_contact_map(Pdb, Chain = Chain, Show = F)
  }
  
  #Visualization
  if(Visualization & Mode != "singleresidue"){
    system(paste("perl ", Pdb$scriptsDir, "/GenerateVisualizations.pl ", Pdb$PdbBase, "_", Pdb$mode, ".pdb_auxiliar ", Pdb$PdbBase, " ", gsub(".$", "", Pdb$JobDir), " ", Pdb$mode, sep = ""))
    VisualizationDir <- paste(Pdb$JobDir, "VisualizationScrips", sep = "")
    if (!dir.exists(VisualizationDir))  dir.create(VisualizationDir)
    system(paste("cp ", Frustration, "/", Pdb$PdbBase, ".pdb ", VisualizationDir, "/", Pdb$PdbBase, ".pdb", sep = ""))
    system(paste("mv *_", Pdb$mode, ".pml ", VisualizationDir, sep = ""))
    system(paste("mv *_", Pdb$mode, ".tcl ", VisualizationDir, sep = ""))
    system(paste("mv *_", Pdb$mode, ".jml ", VisualizationDir, sep = ""))
    system(paste("cp ", Pdb$scriptsDir, "/draw_links.py ", VisualizationDir, sep = ""))
  }
  
  #Unnecessary files are removed
  system(paste("ls -F1 > output;", sep = ""))
  files <- read.table(paste(JobDir, "output", sep = ""), header = F)
  files[,1] <- as.character(files[,1])
  for (i in seq(1,dim(files)[1])) {
  	finalCharacter <- substr(files[i,1], nchar(files[i,1]), nchar(files[i,1]))
  	if(finalCharacter != '/' & finalCharacter != '.'){
  		system(paste("rm -f ", files[i,1], sep = ""));
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

#' Calculate Directory Frustration
#'
#' Calculate local energy frustration for all protein structures in one directory
#'
#' @param PdbDir File containing all proteins structures. The full path to the file is needed
#' @param OrderList Orderer list of PDB files to calculate frustration. If it is NULL, frustration is calculated for all PDBs. Default: NULL
#' @param Electrostatics_K K constant to use in the electrostatics mode. Default: NULL (no electrostatics is considered).
#' @param Mode Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#' @param seqdist Sequence at which contacts are considered to interact (3 or 12). Default: 12
#' @param Graphics The corresponding graphics are made (TRUE or FALSE). Default: TRUE
#' @param Visualization Make visualizations, including pymol. Default: TRUE
#' @param ResultsDir Path to the folder where results will be stored.
#'
#' @export

dir_frustration <- function(PdbDir=NULL, OrderList=NULL, Electrostatics_K=NULL, seqdist=12, Mode="configurational", Graphics=TRUE, Visualization=TRUE, ResultsDir=NULL){
	if(is.null(PdbDir)) stop("PdbDir not specified")
	if(is.null(ResultsDir)) stop("ResultsDir not specified")
  setwd(PdbDir)
  
	CalculationEnabled = TRUE

	if (file.exists(paste(ResultsDir, "modes.log", sep=""))){
		modes <- read.table(paste(ResultsDir, "modes.log", sep=""), header=FALSE)
		for (i in seq(1, length(modes[,1]))) {
			if (modes[i,1] == Mode)
				CalculationEnabled =FALSE
		}
	}

	if(CalculationEnabled){

		if(is.null(OrderList)){

			system(paste("ls *.pdb > listPdb"))
			listPdb<-read.table(paste(PdbDir, "listPdb", sep=""), header=FALSE)
			system(paste("rm listPdb"))

			for (index in seq(1,dim(listPdb)[1])) {

    			if( Mode == "configurational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,listPdb[index,1],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Mode = "configurational" , Graphics=Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
    			}else if( Mode == "mutational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,listPdb[index,1],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Mode = "mutational" , Graphics=Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
    			}else if( Mode == "singleresidue"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,listPdb[index,1],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Mode = "singleresidue" , Graphics=Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
    			}
    		}
		}else{

  			for (index in seq(1,length(OrderList))) {

    			if( Mode == "configurational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,OrderList[index],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Mode = "configurational" , Graphics=Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
    			}else if( Mode == "mutational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,OrderList[index],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Mode = "mutational" , Graphics=Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
    			}else if( Mode == "singleresidue"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,OrderList[index],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Mode = "singleresidue" , Graphics=Graphics, Visualization = Visualization, ResultsDir = ResultsDir)
    			}
    		}
		}

  		write(Mode, file=paste(ResultsDir,"modes.log",sep=""), append=TRUE)
	}

}

#' Dynamic Frustration
#'
#' Calculates local energetic frustration for a dynamic
#'
#' @param PdbDir File containing all proteins structures. The full path to the file is needed
#' @param OrderList Orderer list of PDB files to calculate frustration. If it is NULL, frustration is calculated for all PDBs. Default: NULL
#' @param Electrostatics_K K constant to use in the electrostatics mode. Default: NULL (no electrostatics is considered).
#' @param seqdist Sequence at which contacts are considered to interact (3 or 12). Default: 12
#' @param Mode Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#' @param Resno Resno of residue to analyze. If it is NULL, not analyze. Default: NULL
#' @param Chain Chain of residue to analyze
#' @param GIFs if it is TRUE, gifs of contact maps and proportion 5 adens are made for all the frames. Default: FALSE
#' @param ResultsDir Path to the folder where results will be stored.
#'
#' @export

dynamic_frustration <- function(PdbDir=NULL, OrderList=NULL, Electrostatics_K=NULL, seqdist=12, Mode="configurational", Resno=NULL, Chain=NULL, GIFs=FALSE, ResultsDir=NULL){

  if(is.null(OrderList)){
    system(paste("cd ",PdbDir,"; ls *.pdb > listPdb;",sep=""))
    OrderList<-read.table(paste(PdbDir,"listPdb",sep=""),header=FALSE)
    OrderList<-as.vector(OrderList[,1])
    system(paste("cd ",PdbDir,"; rm listPdb;",sep=""))
  }
	if(is.null(PdbDir)) stop("PdbDir not specified")
	if(is.null(ResultsDir)) stop("ResultsDir not specified")
	if(!is.null(Resno)){
		Pdb<-read.pdb(paste(PdbDir,OrderList[1],sep=""))
	  if(!is.null(Chain)){
	    if(length(atom.select(Pdb,resno=Resno,chain=Chain,elety="CA")$atom)==0)
	      stop("Resno of chain not exist")
	  }else if(length(atom.select(Pdb,resno=Resno,elety="CA")$atom)!=0 & is.na(atom.select(Pdb,resno=Resno,elety="CA",value=TRUE)$atom$chain)){
		  Chain="A"
		}
		rm(Pdb)
	}

	#Pdbs indicated in OrderList found in PdbDir are frustrated
	dir_frustration(PdbDir=PdbDir,OrderList=OrderList, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Mode=Mode, ResultsDir=ResultsDir)

	#It is checked if you want to graph on a particular residue
	if(!is.null(Resno)){

		PlotsDir<-paste(ResultsDir,"Dynamic_plots_res_",Resno,sep="")

		setwd(ResultsDir)
		system(paste("mkdir -p ",PlotsDir,sep=""))

		ResultFile <- paste(PlotsDir,"/",Mode,"_Res",Resno,sep="")

		if(file.exists(ResultFile)) system(paste("rm ",ResultFile,sep=""))


		if (Mode == "configurational" | Mode == "mutational"){

			write("MaximallyFrst	NeutrallyFrst	MinimallyFrst",file=ResultFile)
			for (i in seq(1,length(OrderList))) {
				system(paste("cd ",ResultsDir,basename.pdb(OrderList[i]),".done/FrustrationData;","awk '{if ($1 ==", Resno," && $2 ==",'"',Chain,'"',"){print $7,$8,$9 >> ",'"',ResultFile,'"',";}}' ",basename.pdb(OrderList[i]),".pdb_",Mode,"_5adens;",sep=""))
			}

		}
		else {
			write("FrstIndex",file=ResultFile)
			for (i in seq(1,length(OrderList))) {
				system(paste("cd ",ResultsDir,basename.pdb(OrderList[i]),".done/FrustrationData;","awk '{if ($1 ==", Resno," && $2 ==",'"',Chain,'"',"){print $8 >> ",'"',ResultFile,'"',";}}' ",basename.pdb(OrderList[i]),".pdb_singleresidue;",sep=""))
			}
		}
		#Residue dynamics in particular
		plot_dynamic_res(DynDir=PlotsDir,Resno=Resno,Mode=Mode)
	}

	#Gifs
	if(GIFs){

		if (Mode == "configurational" | Mode == "mutational"){
			gif_5adens_proportions(PdbDir=ResultsDir, OrderList=OrderList,Mode=Mode)
			gif_contact_map(PdbDir=ResultsDir, OrderList=OrderList,Mode=Mode)
		}
	}

	#FrustraMovie
	frustra_movie(PdbDir=ResultsDir,OrderList=OrderList,Mode=Mode)
}


#' Frustration of mutated residue
#'
#' Calculate the local energy frustration for each of the 20 residual variants in the Resno position and Chain chain.
#' Use the frustration index indicated in the Pdb object.
#' 
#' @param Pdb Pdb object
#' @param Chain Chain of the residue to be mutated. Default: NULL
#' @param Resno Resno of the residue to be mutated. Default: NULL
#' @param Split Split that you are going to calculate frustration. If it is TRUE specific string, if it is FALSE full complex. Default: TRUE
#' @param Method Method indicates the method to use to perform the mutation (Threading or Modeller). Default: Threading
#' @return Returns Pdb object with corresponding Mutation attribute
#' @export

mutate_res<-function(Pdb = NULL, Chain = NULL, Resno = NULL, Split = TRUE, Method = "Threading"){
  
  options(stringsAsFactors = FALSE)
  
  #Exceptions
  if(is.null(Pdb)) stop("Pdb not specified")
  if(!is.null(Resno)){
    if(!is.null(Chain)){
      if(length(atom.select(Pdb, resno = Resno, chain = Chain, elety = "CA")$atom) == 0)
        stop("Resno of chain not exist")
    }else if(length(atom.select(Pdb, resno = Resno, elety = "CA")$atom) != 0 & is.na(atom.select(Pdb, resno = Resno, elety = "CA", value = TRUE)$atom$chain)){
      Chain = "A"
    }
  }
  if(Split == FALSE & Method == "Modeller")  stop("Complex modeling not available")
  
  FrustraMutFile <- paste(Pdb$JobDir,"MutationsData/",Pdb$mode,"_Res",Resno,"_",Method,"_", Chain,".txt",sep="")
  if(!dir.exists(paste(Pdb$JobDir,"MutationsData",sep="")))  dir.create(paste(Pdb$JobDir,"MutationsData",sep=""))
  setwd(paste(Pdb$JobDir,"MutationsData",sep=""))
  if (!file.exists(FrustraMutFile)) file.create(FrustraMutFile)
  
  if(Pdb$mode == "configurational" | Pdb$mode == "mutational")
    write("Res1 Res2 ChainRes1 ChainRes2 AA1 AA2 FrstIndex FrstState",file=FrustraMutFile, append=TRUE)
  else if(Pdb$mode == "singleresidue")  write("Res ChainRes AA FrstIndex",file=FrustraMutFile, append=TRUE)  
  
  
  if(Method == "Threading"){
    
    AAvector <- c('LEU','ASP','ILE','ASN','THR','VAL','ALA','GLY','GLU','ARG','LYS','HIS','GLN','SER','PRO','PHE','TYR','MET','TRP','CYS')
    #AAvector<-c("ASP")
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
    
    #It is mutated by 20 AA
    for (AA in AAvector) {
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
        if(length(diffAtom))  PdbMut$atom <- PdbMut$atom[-diffAtom,]
        if(length(diffxyz)) PdbMut$xyz <- PdbMut$xyz[-diffxyz]
      }
      
      #Residues are renamed
      if(AA == "GLY") PdbMut$atom$resid[indexBackboneGly$atom] <- AA
      else PdbMut$atom$resid[indexBackbone$atom] <- AA
      
      #Muted PDB is saved
      if(Split == TRUE) write.pdb(PdbMut, paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, ".pdb", sep = ""))
      else  write.pdb(PdbMut, paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""))
      
      #Gets frustrated
      if(Split==TRUE) calculate_frustration( PdbFile = paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, ".pdb", sep = ""), Mode = Pdb$mode, ResultsDir = Pdb$JobDir, Graphics = F, Visualization = F,Chain = Chain)
      else  calculate_frustration(PdbFile = paste(Pdb$JobDir, "/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""), Mode = Pdb$mode, ResultsDir = Pdb$JobDir, Graphics = FALSE, Visualization = F)
      
      if(Pdb$mode == "singleresidue"){
        system(paste("mv ",Pdb$JobDir,Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".done/FrustrationData/",Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".pdb_singleresidue ",Pdb$JobDir,"MutationsData/",sep=""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_singleresidue", sep=""), header = T)
        frustraTable <- frustraTable[frustraTable$ChainRes == Chain & frustraTable$Res == Resno, c(1,2,4,8)]
        write(paste(frustraTable[,1], frustraTable[,2], frustraTable[,3], frustraTable[,4]), file = FrustraMutFile, append = TRUE)
      }else if(Pdb$mode == "configurational" | Pdb$mode == "mutational"){
        system(paste("mv ",Pdb$JobDir,Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".done/FrustrationData/",Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".pdb_",Pdb$mode," ",Pdb$JobDir,"MutationsData/",sep=""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".pdb_",Pdb$mode,sep=""), header = T)
        frustraTable <- frustraTable[frustraTable$ChainRes1 == Chain & frustraTable$Res1 == Resno |
                                     frustraTable$ChainRes2 == Chain & frustraTable$Res2 == Resno,
                                     c(1,2,3,4,7,8,12,14)]
        write(paste(frustraTable[,1], frustraTable[,2], frustraTable[,3], frustraTable[,4],
                    frustraTable[,5], frustraTable[,6], frustraTable[,7], frustraTable[,8]),
                    file = FrustraMutFile, append = TRUE)
      }
      
      #Unnecessary files are removed
      file.remove(paste(Pdb$JobDir,"MutationsData/",Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_", Pdb$mode, sep=""))
      system(paste("rm -R ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/", sep=""))
      system(paste("cd ", Pdb$JobDir, " ; rm *pdb", sep = ""))
    }
  }
  else if(Method=="Modeller"){
    
    if (!requireNamespace("msa", quietly = TRUE)){
      if(!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")} 
      BiocManager::install("msa")
    }
    
    AAvector <- c('L','D','I','N','T','V','A','G','E','R','K','H','Q','S','P','F','Y','M','W','C');
    
    sequence <- get.seq(ids = c(Pdb$PdbBase), db = "PDB", outfile = paste(Pdb$JobDir, "seqs.fasta", sep = ""))
    fasta <- read.fasta(paste(Pdb$JobDir, "seqs.fasta", sep = ""))
    
    Seq <- c()
    rowname <- c()
    for (i in 1:length(fasta$id)) {
      Seq <- rbind(Seq, fasta$ali[i,])
      rowname <- c(rowname,paste(substr(fasta$id[i], nchar(fasta$id[i]), nchar(fasta$id[i]))))
    }
    
    rownames(Seq) <- toupper(rowname)
    
    SeqPdb <- as.data.frame(pdbseq(Pdb, atom.select(Pdb, chain = Chain, elety = "CA")))
    SeqPdb <- cbind(SeqPdb, as.numeric(row.names(SeqPdb)), 1:length(SeqPdb[,1]))
    colnames(SeqPdb) <- c("AA", "resno", "index")

    aln <- seqbind(Seq[Chain, 1:table(Seq[Chain,] == "-")[1]], SeqPdb$AA)
    align <- seqaln(aln, outfile = tempfile(), exefile = "msa")
    SeqGap <- align$ali["seq2",]
    
    SeqGap <- cbind(SeqGap, rep(0, length(SeqGap)), 1:length(SeqGap))
    
    j <- 1
    for (i in 1:length(SeqGap[,1])) {
      if(SeqGap[i,1] != "-" && SeqGap[i,1] == SeqPdb$AA[j]){
        SeqGap[i,2] <- SeqPdb$resno[j]
        j <- j+1
      }
    }
    pos <- as.numeric(SeqGap[SeqGap[,2] == Resno, 3])
    
    system(paste("cp ", Pdb$scriptsDir, "/align2d.py ", Pdb$JobDir, sep=""))
    system(paste("cp ", Pdb$scriptsDir, "/make_ali.py ", Pdb$JobDir, sep=""))
    system(paste("cp ", Pdb$scriptsDir, "/model-single.py ", Pdb$JobDir, sep=""))
    system(paste("cp ", Pdb$PdbPath," ", Pdb$JobDir, sep=""))

    for (AA in AAvector){
      if(Split){
        write(">modelo", file = paste(Pdb$JobDir, "modelo.fa", sep = ""))
        SeqMut = Seq[Chain, 1:table(Seq[Chain,] == "-")[1]]
        SeqMut[pos] <- AA
        write(paste(as.vector(SeqMut), collapse = ""),file = paste(Pdb$JobDir, "modelo.fa", sep = ""), append = TRUE)
      }
      
      system(paste("cd ", Pdb$JobDir, " ;python3 make_ali.py modelo", sep = ""))
      if(Split) system(paste("cd ", Pdb$JobDir, " ;python3 align2d.py ", Pdb$PdbBase, " modelo ", Chain, sep = ""))
  
      system(paste("cd ", Pdb$JobDir, " ;python3 model-single.py ", Pdb$PdbBase, " modelo", sep = ""))
      system(paste("cd ", Pdb$JobDir, " ;mv modelo.B99990001.pdb ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""))
      system(paste("cd ", Pdb$JobDir, " ;rm *D00000001 *ini *rsr *sch *V99990001 *ali *pap *fa"))
      
      calculate_frustration(PdbFile = paste(Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb", sep = ""), Mode = Pdb$mode, ResultsDir = Pdb$JobDir, Graphics = F, Visualization = F)
      
      if(Pdb$mode == "singleresidue"){
        system(paste("mv ",Pdb$JobDir,Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".done/FrustrationData/",Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".pdb_singleresidue ",Pdb$JobDir,"MutationsData/",sep=""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_singleresidue", sep=""), header = T)
        frustraTable <- frustraTable[frustraTable$ChainRes == Chain & frustraTable$Res == Resno, c(1,2,4,8)]
        write(paste(frustraTable[,1], frustraTable[,2], frustraTable[,3], frustraTable[,4]), file = FrustraMutFile, append = TRUE)
      }else if(Pdb$mode == "configurational" | Pdb$mode == "mutational"){
        system(paste("mv ",Pdb$JobDir,Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".done/FrustrationData/",Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".pdb_",Pdb$mode," ",Pdb$JobDir,"MutationsData/",sep=""))
        frustraTable <- read.table(paste(Pdb$JobDir, "MutationsData/", Pdb$PdbBase,"_",Resno,"_",AA,"_",Chain,".pdb_",Pdb$mode,sep=""), header = T)
        frustraTable <- frustraTable[frustraTable$ChainRes1 == Chain & frustraTable$Res1 == Resno |
                                       frustraTable$ChainRes2 == Chain & frustraTable$Res2 == Resno,
                                     c(1,2,3,4,7,8,12,14)]
        write(paste(frustraTable[,1], frustraTable[,2], frustraTable[,3], frustraTable[,4],
                    frustraTable[,5], frustraTable[,6], frustraTable[,7], frustraTable[,8]),
              file = FrustraMutFile, append = TRUE)
      }
      
      #Unnecessary files are removed
      file.remove(paste(Pdb$JobDir,"MutationsData/",Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".pdb_", Pdb$mode, sep=""))
      system(paste("rm -R ", Pdb$JobDir, Pdb$PdbBase, "_", Resno, "_", AA, "_", Chain, ".done/", sep=""))
      system(paste("cd ", Pdb$JobDir, " ; rm ",Pdb$PdbBase,"_*", sep = ""))
      
    }
    system(paste("cd ", Pdb$JobDir, " ; rm *pdb seqs.fasta *py", sep = ""))
    
    #Renumber residues
    if(Pdb$mode == "singleresidue"){
      data <- read.table(paste(Pdb$JobDir, "MutationsData/singleresidue_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""), header = T)
      data[,2] <- Chain
      for (i in 1:length(data[,1])) {
        data[i,1] <- SeqGap[as.numeric(data[i,1]),2]
      }
      write.table(paste(Pdb$JobDir, "MutationsData/singleresidue_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""), col.names = T,row.names = F)
    }else if(Pdb$mode == "configurational"){
      data <- read.table(paste(Pdb$JobDir, "MutationsData/configurational_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""), header = T)
      data[,c(3,4)] <- Chain
      for (i in 1:length(data[,1])) {
        data[i,1] <- SeqGap[as.numeric(data[i,1]),2]
        data[i,2] <- SeqGap[as.numeric(data[i,2]),2]
      }
      write.table(data, paste(Pdb$JobDir, "MutationsData/configurational_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""), col.names = T, row.names = F, quote = F)
    }else if(Pdb$mode == "mutational"){
      data <- read.table(paste(Pdb$JobDir, "MutationsData/mutational_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""), header = T)
      data[,c(3,4)] <- Chain
      for (i in 1:length(data[,1])) {
        data[i,1] <- SeqGap[as.numeric(data[i,1]),2]
        data[i,2] <- SeqGap[as.numeric(data[i,2]),2]
      }
      write.table(data, paste(Pdb$JobDir, "MutationsData/mutational_Res", Resno, "_", Method,"_",Chain ,".txt", sep = ""), col.names = T, row.names = F, quote = F)
    }
  }
  
  Pdb$Mutations[[Method]][[paste("Res_", Resno, "_", Chain, sep = "")]] <- data.frame(Method = Method ,Res = Resno, 
                                                                                               Chain = Chain, File = paste(Pdb$JobDir, "MutationsData/", Pdb$mode, "_Res", Resno, "_", Method,"_",Chain, ".txt", sep = ""))

  return(Pdb)
}