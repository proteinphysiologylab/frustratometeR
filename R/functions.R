#' Pdb Equivalences
#'
#' Internal function that produces auxiliar files to run the frustratometer pipeline.
#'
#' @param Pdb Frustration object
#' @return Pdb equivalences

pdb_equivalences <- function(Pdb)
{
  existing_res<-unique( cbind( Pdb$atom$chain[ which( Pdb$atom$type=="ATOM") ], Pdb$atom$resno[ which( Pdb$atom$type=="ATOM" ) ], Pdb$atom$resid[ which( Pdb$atom$type=="ATOM" ) ] ) )
  equivalences <- cbind(existing_res[,1], 1:length(existing_res[,1]), existing_res[,2],existing_res[,3])
  write.table(equivalences, file = paste(Pdb$JobDir, "/", Pdb$PdbBase, ".pdb_equivalences.txt", sep=""), quote = F, col.names = F, row.names = F, sep="\t")
  return(equivalences)
}

#' Check Backbone Complete
#'
#' Checks the backbone of a given protein structure to be processed by the frustratometer pipeline
#'
#' @param Pdb Frustration object
#' @return Flag indicating if the backbone should be completed

check_backbone_complete <- function(Pdb)
{
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

complete_backbone <- function(Pdb)
{
  Completed<-FALSE
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
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default=configurational
#' @param Electrostatics_K K constant to use in the electrostatics mode. Default: NULL (no electrostatics is considered).
#' @param seqdist Sequence at which contacts are considered to interact (3 or 12). Default: 12
#' @param Graphics The corresponding graphics are made (TRUE or FALSE). Default: TRUE
#' @param ResultsDir Path to the folder where results will be stored.
#' @return Pdb Frustration object
#'
#' @export

calculate_frustration <- function(PdbFile=NULL, PdbID=NULL, Chain=NULL, Electrostatics_K=NULL, seqdist=12, Modes="configurational",Graphics=TRUE, ResultsDir)
{

  tempfolder <-tempdir()

  #Chain
  if(is.null(Chain))
  {
    boolsplit=F
  }else{
    boolsplit=T
  }

  if(is.null(PdbFile))
  {

    #Get URL and download
    print("-----------------------------Download files-----------------------------")
    pdbURL<-get.pdb( id = PdbID, split = boolsplit, URLonly = TRUE)
    system(paste("wget --no-check-certificate -P",tempfolder,pdbURL,sep = ' '))

    if(!boolsplit)
    {
    	PdbFile=paste(tempfolder,"/",PdbID,".pdb", sep="")

    }else{

      	#Chain splitting
    	pdbsplit(pdb.files = paste(tempfolder,"/",PdbID,".pdb",sep =""), path = file.path(tempfolder, "/split_chain"))

    	#Nonexistent Chain exception
    	if(!file.exists(paste(tempfolder, "/split_chain/",PdbID, "_", Chain, ".pdb", sep=""))){stop("Nonexistent Chain")}

        PdbFile=paste(tempfolder,"/split_chain/",PdbID, "_", Chain, ".pdb", sep="")
    }
  }
  else{

  	if(boolsplit){
  		pdbsplit(pdb.files = PdbFile, path = file.path(tempfolder, "/split_chain"))

  		#Nonexistent Chain exception
    	if(!file.exists(paste(tempfolder,"/split_chain/",basename.pdb(PdbFile), "_", Chain, ".pdb", sep=""))){stop("Nonexistent Chain")}

  		PdbFile=paste(tempfolder,"/split_chain/",basename.pdb(PdbFile), "_", Chain, ".pdb", sep="")
  	}
  }

  PdbBase <- basename.pdb(PdbFile)
  JobDir=paste(ResultsDir, PdbBase, ".done/", sep="")

  #Creates JobDir
  system(paste("mkdir -p ", JobDir, sep=""))
  system(paste("cp ", PdbFile, " ", JobDir, sep=""))

  PdbFile = paste(JobDir,PdbBase,".pdb",sep="")

  # Here we filter out DNA chains
  system(paste("awk '{if($0!~/DT/ && $0!~/DG/ && $0!~/DA/ && $0!~/DC/ && $4!=",'"',"T",'"'," && $4!=",'"',"G",'"'," && $4!=",'"',"A",'"'," && $4!=",'"',"C",'"',") print }' ", PdbFile ," > ", JobDir, "/aux", sep=""))
  system(paste("mv ", JobDir, "/aux ", PdbFile, sep=""))

  # We read the pdb file and ignore alternative conformations, aminoacids and only read ATOM lines
  Pdb <- read.pdb(PdbFile, ATOM.only=T, rm.alt = T, rm.insert = T)

  #Fix chain NA
  Pdb$atom$chain[is.na(Pdb$atom$chain[])] <- "A"

  #MSE to MET filter
  Pdb$atom$type[ which( Pdb$atom$resid == "MSE" ) ] <- "ATOM"
  Pdb$atom$resid[ which( Pdb$atom$resid == "MSE" ) ] <- "MET"

  Pdb[["PdbBase"]] <- PdbBase
  Pdb[["JobDir"]] <- JobDir
  Pdb[["scriptsDir"]] <- paste(find.package("frustratometeR"), "/Scripts", sep="")
  Pdb[["mode"]] <- Modes

  # Save equivalences
  equivalences = pdb_equivalences(Pdb)
  Pdb[["equivalences"]] <- equivalences

  write.pdb(Pdb,paste(JobDir,PdbBase,".pdb",sep=""))

  if(complete_backbone(Pdb)){
  	Pdb <- read.pdb(paste(JobDir,PdbBase,".pdb",sep=""))
  	Pdb<-atom.select(Pdb,type="ATOM",value=TRUE)
  	write.pdb(Pdb,paste(JobDir,PdbBase,".pdb",sep=""))

  	Pdb[["PdbBase"]] <- PdbBase
  	Pdb[["JobDir"]] <- JobDir
  	Pdb[["scriptsDir"]] <- paste(find.package("frustratometeR"), "/Scripts", sep="")
  	Pdb[["mode"]] <- Modes
  	Pdb[["equivalences"]] <- equivalences
  }

  print("-----------------------------Preparing files-----------------------------")
  #Prepare the PDB file to get awsem input files, create the workdir and move neccessary files to it.
  system(paste("cd ", Pdb$JobDir, "; pwd; ", Pdb$scriptsDir, "/AWSEMFiles/AWSEMTools/PdbCoords2Lammps.sh ", Pdb$PdbBase, " ", Pdb$PdbBase, " ", Pdb$scriptsDir, sep=""))
  system(paste("cp ", Pdb$scriptsDir, "/AWSEMFiles/*.dat* ", Pdb$JobDir, sep=""))

  print("-----------------------------Setting options-----------------------------")

  system(paste("cd ", Pdb$JobDir, "; sed -i  's/run		10000/run		0/g' ", Pdb$PdbBase, ".in; sed -i 's/mutational/", Pdb$mode, "/g' fix_backbone_coeff.data",  sep=""))

  if(!is.null(Electrostatics_K))
  {
    print("Setting electrostatics...")
    system(paste("cd ", Pdb$JobDir, "; sed -i 's/\\[DebyeHuckel\\]-/\\[DebyeHuckel\\]/g' fix_backbone_coeff.data; sed -i 's/4.15 4.15 4.15/", Electrostatics_K, " ", Electrostatics_K, " ", Electrostatics_K, "/g' fix_backbone_coeff.data;", sep=""))
    print("Setting electrostatics...")
    system(paste("cd ", Pdb$JobDir, "; python3 ", Pdb$scriptsDir, "/Pdb2Gro.py ", Pdb$PdbBase, ".pdb ", Pdb$PdbBase, ".pdb.gro; perl ", Pdb$scriptsDir, "/GenerateChargeFile.pl ", Pdb$PdbBase, ".pdb.gro > ", JobDir, "charge_on_residues.dat", sep=""))
  }

  print("-----------------------------Calculating-----------------------------")

  system(paste("cp ", Pdb$scriptsDir, "/lmp_serial_", seqdist, " ", Pdb$JobDir, "; cd ", Pdb$JobDir, "; ./lmp_serial_", seqdist, " < ", Pdb$PdbBase, ".in", sep=""))

  if(Pdb$mode == "configurational" | Pdb$mode == "mutational")
  {
    system(paste("perl ", Pdb$scriptsDir, "/5Adens.pl ", Pdb$PdbBase, ".pdb ", gsub(".$", "", Pdb$JobDir), " ", Pdb$mode, sep=""))
  }
  system(paste("perl ", Pdb$scriptsDir, "/RenumFiles.pl ", Pdb$PdbBase, " ", Pdb$JobDir, " ", Pdb$mode, sep="" ))
  system(paste("perl ", Pdb$scriptsDir, "/GenerateVisualizations.pl ", Pdb$PdbBase, "_", Pdb$mode, ".pdb_auxiliar ", Pdb$PdbBase, " ", gsub(".$", "", Pdb$JobDir), " ", Pdb$mode, sep=""))


  #Directory reorganization and results
  print("-----------------------------Reorganization-----------------------------")

  Frustration = paste(Pdb$JobDir,"FrustrationData",sep="")
  Images = paste(Pdb$JobDir,"Images",sep="")
  Visualization = paste(Pdb$JobDir,"VisualizationScrips",sep="")

  system(paste("mkdir -p ",Frustration,";","mkdir -p ",Images,";","mkdir -p ",Visualization,";",sep=""))

  system(paste("cd ",Pdb$JobDir,";","mv ","*.pdb_",Pdb$mode," ",Frustration,";",sep=""))
  if(Pdb$mode == "configurational" | Pdb$mode == "mutational"){
  	system(paste("cd ",Pdb$JobDir,";","mv ","*_",Pdb$mode,"_5adens ",Frustration,";",sep=""))
  }

  system(paste("cd ",Pdb$JobDir,";","mv ","*.pdb ",Visualization,sep=""))
  system(paste("cd ",Pdb$JobDir,";","mv ","*_",Pdb$mode,".pml ",Visualization,sep=""))
  system(paste("cd ",Pdb$JobDir,";","mv ","*_",Pdb$mode,".tcl ",Visualization,sep=""))
  system(paste("cp ", Pdb$scriptsDir, "/draw_links.py ",Visualization, sep=""))

  if(Graphics&Modes!="singleresidue"){
  	plot_5Andens(Pdb,Chain=Chain)
  	plot_5Adens_proportions(Pdb,Chain=Chain)
  	plot_contact_map(Pdb,Chain=Chain)
  }

  #Unnecessary files are removed
  system(paste("cd ",Pdb$JobDir,";ls -F1 > output;",sep=""));
  files <- read.table(paste(JobDir,"output",sep=""),header=F)
  files[,1] <- as.character(files[,1])
  for (i in seq(1,dim(files)[1])) {
  	finalCharacter<-substr(files[i,1],nchar(files[i,1]),nchar(files[i,1]))
  	if(finalCharacter!='/'&finalCharacter!='.'){
  		system(paste("cd ",Pdb$JobDir,";rm -f ",files[i,1],";",sep=""));
  	}
  }

  #We delete temporary files
  if(!is.null(Chain))
  {
  	system(paste("rm -f -R ",tempfolder,"/split_chain", sep="" ))
  }
  system(paste("rm -f ",tempfolder,"/",PdbID,".pdb", sep="" ))


  return(Pdb)
}

#' Calculate Directory Frustration
#'
#' Calculate local energy frustration for all protein structures in one directory
#'
#' @param PdbDir File containing all proteins structures. The full path to the file is needed
#' @param OrderList Orderer list of PDB files to calculate frustration. If it is NULL, frustration is calculated for all PDBs. Default: NULL
#' @param Electrostatics_K K constant to use in the electrostatics mode. Default: NULL (no electrostatics is considered).
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#' @param seqdist Sequence at which contacts are considered to interact (3 or 12). Default: 12
#' @param Graphics The corresponding graphics are made (TRUE or FALSE). Default: TRUE
#' @param ResultsDir Path to the folder where results will be stored.
#'
#' @export

dir_frustration <- function(PdbDir=NULL,OrderList=NULL, Electrostatics_K=NULL, seqdist=12, Modes="configurational", Graphics=TRUE, ResultsDir=NULL)
{
	if(is.null(PdbDir)) stop("PdbDir not specified")
	if(is.null(ResultsDir)) stop("ResultsDir not specified")

	CalculationEnabled=TRUE

	if (file.exists(paste(ResultsDir,"modes.log",sep=""))){
		modes<-read.table(paste(ResultsDir,"modes.log",sep=""),header=FALSE)
		for (i in seq(1,length(modes[,1]))) {
			if (modes[i,1] == Modes)
				CalculationEnabled=FALSE
		}
	}

	if(CalculationEnabled){

		if(is.null(OrderList)){

			system(paste("cd ",PdbDir,"; ls *.pdb > listPdb;",sep=""))
			listPdb<-read.table(paste(PdbDir,"listPdb",sep=""),header=FALSE)
			system(paste("cd ",PdbDir,"; rm listPdb;",sep=""))

			for (index in seq(1,dim(listPdb)[1])) {

    			if( Modes == "configurational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,listPdb[index,1],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Modes = "configurational" , Graphics=Graphics,ResultsDir = ResultsDir)
    			}else if( Modes == "mutational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,listPdb[index,1],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Modes = "mutational" , Graphics=Graphics, ResultsDir = ResultsDir)
    			}else if( Modes == "singleresidue"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,listPdb[index,1],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Modes = "singleresidue" , Graphics=Graphics, ResultsDir = ResultsDir)
    			}
    		}
		}else{

  			for (index in seq(1,length(OrderList))) {

    			if( Modes == "configurational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,OrderList[index],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Modes = "configurational" , Graphics=Graphics, ResultsDir = ResultsDir)
    			}else if( Modes == "mutational"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,OrderList[index],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Modes = "mutational" , Graphics=Graphics, ResultsDir = ResultsDir)
    			}else if( Modes == "singleresidue"){
    		  		Pdb = calculate_frustration(PdbFile = paste(PdbDir,OrderList[index],sep=""), Chain=NULL, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Modes = "singleresidue" , Graphics=Graphics, ResultsDir = ResultsDir)
    			}
    		}
		}

  		write(Modes, file=paste(ResultsDir,"modes.log",sep=""), append=TRUE)

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
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#' @param Resno Resno of residue to analyze. If it is NULL, not analyze. Default: NULL
#' @param Chain Chain of residue to analyze
#' @param GIFs if it is TRUE, gifs of contact maps and proportion 5 adens are made for all the frames. Default: FALSE
#' @param ResultsDir Path to the folder where results will be stored.
#'
#' @export

dynamic_frustration <- function(PdbDir=NULL, OrderList=NULL, Electrostatics_K=NULL, seqdist=12, Modes="configurational", Resno=NULL, Chain=NULL, GIFs=FALSE, ResultsDir=NULL)
{

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
	dir_frustration(PdbDir=PdbDir,OrderList=OrderList, Electrostatics_K=Electrostatics_K, seqdist=seqdist, Modes=Modes, ResultsDir=ResultsDir)

	#It is checked if you want to graph on a particular residue
	if(!is.null(Resno)){

		PlotsDir<-paste(ResultsDir,"Dynamic_plots_res_",Resno,sep="")

		setwd(ResultsDir)
		system(paste("mkdir -p ",PlotsDir,sep=""))

		ResultFile <- paste(PlotsDir,"/",Modes,"_Res",Resno,sep="")

		if(file.exists(ResultFile)) system(paste("rm ",ResultFile,sep=""))


		if (Modes == "configurational" | Modes == "mutational"){

			write("MaximallyFrst	NeutrallyFrst	MinimallyFrst",file=ResultFile)
			for (i in seq(1,length(OrderList))) {
				system(paste("cd ",ResultsDir,basename.pdb(OrderList[i]),".done/FrustrationData;","awk '{if ($1 ==", Resno," && $2 ==",'"',Chain,'"',"){print $7,$8,$9 >> ",'"',ResultFile,'"',";}}' ",basename.pdb(OrderList[i]),".pdb_",Modes,"_5adens;",sep=""))
			}

		}
		else {
			write("FrstIndex",file=ResultFile)
			for (i in seq(1,length(OrderList))) {
				system(paste("cd ",ResultsDir,basename.pdb(OrderList[i]),".done/FrustrationData;","awk '{if ($1 ==", Resno," && $2 ==",'"',Chain,'"',"){print $8 >> ",'"',ResultFile,'"',";}}' ",basename.pdb(OrderList[i]),".pdb_singleresidue;",sep=""))
			}
		}
		#Residue dynamics in particular
		plot_dynamic_res(DynDir=PlotsDir,Resno=Resno,Modes=Modes)
	}

	#Gifs
	if(GIFs){

		if (Modes == "configurational" | Modes == "mutational"){
			gif_5adens_proportions(PdbDir=ResultsDir, OrderList=OrderList,Modes=Modes)
			gif_contact_map(PdbDir=ResultsDir, OrderList=OrderList,Modes=Modes)
		}
	}

	#FrustraMovie
	frustra_movie(PdbDir=ResultsDir,OrderList=OrderList,Modes=Modes)
}


#' Frustration of mutated residue
#'
#' Calculates local energetic frustration for residue mutation
#'
#' @param PdbPath Pdb file path
#' @param JobDir Full directory path where the results will be stored
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#' @param Chain Chain of residue to analyze. Default: NULL
#' @param Resno Resno of residue to analyze. If it is NULL, not analyze. Default: NULL
#' @return Table of frustration results in file
#'
#' @export

mutate_res<-function(PdbPath=NULL,JobDir=NULL,Modes="configurational",Chain=NULL,Resno=NULL){

	#Read Pdb
	Pdb<-read.pdb(PdbPath)

	#Exceptions
	if(is.null(PdbPath)) stop("PdbPath not specified")
	if(is.null(JobDir)) stop("JobDir not specified")
	if(!is.null(Resno)){
	  Pdbcompr<-read.pdb(PdbPath)
	  if(!is.null(Chain)){
	    if(length(atom.select(Pdbcompr,resno=Resno,chain=Chain,elety="CA")$atom)==0)
	      stop("Resno of chain not exist")
	  }else if(length(atom.select(Pdbcompr,resno=Resno,elety="CA")$atom)!=0 & is.na(atom.select(Pdbcompr,resno=Resno,elety="CA",value=TRUE)$atom$chain)){
	    Chain="A"
	  }
	  rm(Pdbcompr)
	}

	AAvector<-c('LEU','ASP','ILE','ASN','THR','VAL','ALA','GLY','GLU','ARG','LYS','HIS','GLN','SER','PRO','PHE','TYR','MET','TRP','CYS')

	Glycine<-FALSE

	#Indices of all the atoms of the residue to mutate
	indexTotales<-atom.select(Pdb, chain=Chain, resno=Resno)
	#If it were glycine
	indexBackboneGly<-atom.select(Pdb, chain=Chain, resno=Resno,elety=c("N", "CA", "C", "O"))

	#Check if it is glycine
	if(Pdb$atom$resid[indexTotales$atom[1]]=="GLY") Glycine=TRUE
	else Glycine=FALSE

	#Backbone indices
	if(Glycine) indexBackbone<-atom.select(Pdb, chain=Chain, resno=Resno,elety=c("N", "CA", "C", "O"))
	else indexBackbone<-atom.select(Pdb, chain=Chain, resno=Resno,elety=c("N", "CA", "C", "O", "CB"))


	#It is mutated by 20 AA
	for (AA in AAvector) {
	  PdbMut<-Pdb
	  #If AA is equal to the residue it is not necessary to mutate and neither is it glycine
	  if(AA!=Pdb$atom$resid[indexTotales$atom[1]]&!Glycine){
	    #if the residue to be inserted is not glycine, insert backbone with CB
	    if(AA!="GLY"){
	      diffAtom<-setdiff(indexTotales$atom,indexBackbone$atom)
	      diffxyz<-setdiff(indexTotales$xyz,indexBackbone$xyz)
	    }
	    #if the residue to be inserted is glycine, a backbone without CB is inserted
	    else {
	      diffAtom<-setdiff(indexTotales$atom,indexBackboneGly$atom)
	      diffxyz<-setdiff(indexTotales$xyz,indexBackboneGly$xyz)
	    }
	    #If the previous subtraction is not empty, the corresponding atoms and coordinates are removed
		if(length(diffAtom))  PdbMut$atom<-PdbMut$atom[-diffAtom,]
		if(length(diffxyz)) PdbMut$xyz<-PdbMut$xyz[-diffxyz]
	}

		#Residues are renamed
		if(AA=="GLY") PdbMut$atom$resid[indexBackboneGly$atom]<-AA
		else PdbMut$atom$resid[indexBackbone$atom]<-AA

		#Muted PDB is saved
		write.pdb(PdbMut,paste(JobDir,"/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".pdb",sep=""))

		#Gets frustrated
		calculate_frustration(PdbFile = paste(JobDir,"/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".pdb",sep=""),Modes = Modes,ResultsDir = paste(JobDir,"/",sep=""),Graphics = FALSE)

		system(paste("mkdir -p ",JobDir,"/FrustrationData",sep=""))
		if(Modes=="singleresidue"){
			system(paste("mv ",JobDir,"/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".done/FrustrationData/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".pdb_singleresidue ",JobDir,"/FrustrationData/",sep=""))
			system(paste("cd ",JobDir,"/FrustrationData/ ; ","awk '{if($2==",'"',Chain,'"',"&&$1==",Resno,") print $1,$2,$4,$8 >> ",'"',"singleresidue_Res",Resno,".txt",'"'," }' *.pdb_singleresidue",sep=""))
			system(paste("cd ",JobDir,"/FrustrationData/ ; rm *pdb_singleresidue",sep = ""))
		}else if(Modes=="configurational"){
			system(paste("mv ",JobDir,"/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".done/FrustrationData/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".pdb_configurational ",JobDir,"/FrustrationData/",sep=""))
			system(paste("cd ",JobDir,"/FrustrationData/ ; ","awk '{if(($3==",'"',Chain,'"',"&&$1==",Resno,")||($4==",'"',Chain,'"',"&&$2==",Resno,")) print $1,$2,$7,$8,$12,$14 >> ",'"',"configurational_Res",Resno,".txt",'"'," }' *.pdb_configurational",sep=""))
			system(paste("cd ",JobDir,"/FrustrationData/ ; rm *pdb_configurational",sep = ""))
		}else if(Modes=="mutational"){
			system(paste("mv ",JobDir,"/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".done/FrustrationData/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".pdb_mutational ",JobDir,"/FrustrationData/",sep=""))
			system(paste("cd ",JobDir,"/FrustrationData/ ; ","awk '{if(($3==",'"',Chain,'"',"&&$1==",Resno,")||($4==",'"',Chain,'"',"&&$2==",Resno,")) print $1,$2,$7,$8,$12,$14 >> ",'"',"mutational_Res",Resno,".txt",'"'," }' *.pdb_mutational",sep=""))
			system(paste("cd ",JobDir,"/FrustrationData/ ; rm *pdb_mutational",sep = ""))
		}

		#Unnecessary files are removed
		system(paste("rm -R ",JobDir,"/",basename.pdb(PdbPath),"_",Resno,"_",AA,"_",Chain,".done/",sep=""))
		system(paste("cd ",JobDir," ; rm *pdb",sep = ""))
	}

}
