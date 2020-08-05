#' Plot 5Adens
#'
#' Generates plot to analyze the density of contacts around a sphere of 5 Armstrongs, centered in the C-alfa atom from the residue. The different classes of contacts based on the mutational frustration index are counted in absolute terms.
#'
#' @param Pdb Frustration object
#' @param Chain Chain of residue to analyze. Default: NULL
#' @export
#'
plot_5Andens <- function(Pdb, Chain=NULL)
{
  JobID=Pdb$PdbBase;
  Dir=Pdb$JobDir;
  mode=Pdb$mode;

  AdensTable=read.table(file=paste(Dir,"FrustrationData/", JobID,".pdb_", Pdb$mode, "_5adens", sep=""),stringsAsFactors = FALSE)
  AdensTable <- as.data.frame(AdensTable)
  colnames(AdensTable) <- c("Positions","Chains","Total","MaximallyFrst","NeutrallyFrst","MinimallyFrst")

  PositionsTotal=seq(from=1, to=length( AdensTable$Positions ), by=1)
  AdensTable <- cbind( AdensTable , PositionsTotal )

  if(is.null(Chain))
  {
    Maximum=max(c(max(AdensTable$MaximallyFrst),max(AdensTable$MinimallyFrst),max(AdensTable$NeutrallyFrst), max(AdensTable$Total) ))

    Graphic <- ggplot() + geom_line( aes(x = AdensTable$PositionsTotal, y = AdensTable$MaximallyFrst , colour = "a")  )
    Graphic <- Graphic + geom_line( aes( x = AdensTable$PositionsTotal, y = AdensTable$NeutrallyFrst , colour = "b" ) )
    Graphic <- Graphic + geom_line( aes( x = AdensTable$PositionsTotal, y = AdensTable$MinimallyFrst , colour =  "c" ) )
    Graphic <- Graphic + geom_line( aes( x = AdensTable$PositionsTotal, y = AdensTable$Total , colour = "d" ) )
    Graphic <- Graphic + ylab( "Local frustration density (5A sphere)" ) + xlab( "Position" )
    Graphic <- Graphic + scale_color_manual(name = "",labels = c( "Highly frustrated" , "Neutral" , "Minimally frustrated" , "Total" ) , values = c("red", "gray", "green", "black"))
    Graphic <- Graphic + scale_y_continuous( breaks = seq(0,Maximum,5) ) + scale_x_continuous( breaks = seq( 1 , length(AdensTable$Positions), trunc((length(AdensTable$Positions)-1)/10)) )
    Graphic <- Graphic + theme_classic()

    ggsave( filename=paste(Dir, "Images/", JobID, "_", mode, ".png_5Adens", ".png", sep=""), plot=Graphic,width=10, height= 6)
  }else{
    AdensTable<-AdensTable[AdensTable$Chains==Chain,]

    Maximum=max(c(max(AdensTable$MaximallyFrst),max(AdensTable$MinimallyFrst),max(AdensTable$NeutrallyFrst), max(AdensTable$Total) ))

    Graphic <- ggplot() + geom_line( aes( x = AdensTable$PositionsTotal,y = AdensTable$MaximallyFrst , colour = "a" )  )
    Graphic <- Graphic + geom_line( aes( x = AdensTable$PositionsTotal,y = AdensTable$NeutrallyFrst , colour = "b" ) )
    Graphic <- Graphic + geom_line( aes( x = AdensTable$PositionsTotal,y = AdensTable$MinimallyFrst , colour =  "c" ) )
    Graphic <- Graphic + geom_line( aes( x = AdensTable$PositionsTotal,y = AdensTable$Total , colour = "d" ) )
    Graphic <- Graphic + ylab( "Local frustration density (5A sphere)" ) + xlab( "Position" )
    Graphic <- Graphic + scale_color_manual(name = "", labels = c( "Highly frustrated" , "Neutral" , "Minimally frustrated" , "Total" ) , values = c("red", "gray", "green", "black"))
    Graphic <- Graphic + scale_y_continuous( breaks = seq(0,Maximum,5) ) + scale_x_continuous( breaks = seq( 1 , length(AdensTable$Positions), trunc((length(AdensTable$Positions)-1)/10) ))
    Graphic <- Graphic + theme_classic()

    ggsave( filename=paste(Dir, "Images/", JobID, "_", mode, "_5Adens__chain", Chain, ".png", sep=""), plot=Graphic,width=10, height= 6)
  }
}

#' Plot 5Adens proportions
#'
#' Generates plot to analyze the density of contacts around a sphere of 5 Armstrongs, centered in the C-alfa atom from the residue. The different classes of contacts based on the mutational frustration index are counted in relative terms.
#'
#' @param Pdb Frustration object
#' @param Chain Chain of residue to analyze. Default: NULL
#' @export
#'
plot_5Adens_proportions <- function(Pdb, Chain=NULL)
{
  JobID=Pdb$PdbBase;
  Dir=Pdb$JobDir;
  mode=Pdb$mode;

  AdensTable=read.table(file=paste(Dir,"FrustrationData/", JobID,".pdb_", Pdb$mode, "_5adens", sep=""),fill=T)
  Positions=as.numeric(AdensTable[,1])
  Chains=AdensTable[,2]
  MaximallyFrst= as.numeric(AdensTable[,4])
  NeutrallyFrst=as.numeric(AdensTable[,5])
  MinimallyFrst=as.numeric(AdensTable[,6])
  Total=as.numeric(AdensTable[,3])

  PositionsTotal=seq(from=1, to=length(Positions), by=1)
  MinimallyFrst=as.numeric(AdensTable[,9])
  NeutrallyFrst=as.numeric(AdensTable[,8])
  MaximallyFrst=as.numeric(AdensTable[,7])

  if(is.null(Chain))
  {
    png(filename = paste(Dir, "Images/", JobID, "_", mode, "_5Adens_around.png", sep=""), width = 540, height = 420)
    par(mar=c(5,5,3,2), xpd=TRUE)
    par(mgp=c(2,0.5,0))
    barplot(rbind(MinimallyFrst, NeutrallyFrst, MaximallyFrst), col=c("green", "gray", "red"), axis.lty=1, xlab="Position", ylab="Density arround 5A sphere (%)", cex.lab=1.5, border = NA, space = c(0), xaxt = "n")
    AuxAxis=seq(from=0, to=length(PositionsTotal))
    axis(side = 1, at = AuxAxis , labels = (AuxAxis+min(PositionsTotal)),  tick = FALSE)
    box(lwd=2)
    legend(x="top",inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated"), pch=c(15, 15, 15), col=c("green", "gray", "red"), horiz = T, bty = "n")
    par(mar=c(5,1,5,5))
    mtext(text = Pdb$PdbBase,line=4,side=4,cex=1.5)

    dev.off()
  }
  else{
    png(filename = paste(Dir, "Images/", JobID, "_", mode, "_5Adens_around_chain", Chain, ".png", sep=""), width = 540, height = 420)
    par(mar=c(5,5,3,2), xpd=TRUE)
    par(mgp=c(2,0.5,0))
    barplot(rbind(MinimallyFrst[which(Chains==Chain)], NeutrallyFrst[which(Chains==Chain)], MaximallyFrst[which(Chains==Chain)]), col=c("green", "gray", "red"), axis.lty=1, xlab="Position", ylab="Density arround 5A sphere (%)", cex.lab=1.5, border = NA, space = c(0), xaxt = "n")
    AuxAxis=seq(from=0, to=length(Positions[which(Chains==Chain)]))
    axis(side = 1, at = AuxAxis , labels = (AuxAxis+min(Positions[which(Chains==Chain)])),  tick = FALSE)
    legend(x="top",inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated"), pch=c(15, 15, 15), col=c("green", "gray", "red"), horiz = T, bty = "n")
    box(lwd=2)
    par(mar=c(5,1,5,5))
    mtext(text = Pdb$PdbBase,line=4,side=4,cex=1.5)
    dev.off()
  }
}

#' Plot frustration contact map
#'
#' Generates contact map plot to visualize the frustration values assigned to each contact
#'
#' @param Pdb Frustration object
#' @param Chain Chain of residue to analyze. Default: NULL
#' @export
#'
plot_contact_map <-function(Pdb, Chain=NULL)
{
  JobID=Pdb$PdbBase;
  Dir=Pdb$JobDir;
  mode=Pdb$mode;

  AdensTable=read.table(file=paste(Dir,"FrustrationData/", JobID,".pdb_", Pdb$mode, "_5adens", sep=""), fill=T)

  if(!is.null(Chain)){
    AdensTable<-AdensTable[AdensTable[,2]==Chain,]
  }

  Positions=as.numeric(AdensTable[,1])

  Chains=AdensTable[,2]
  MaximallyFrst= as.numeric(AdensTable[,4])
  NeutrallyFrst=as.numeric(AdensTable[,5])
  MinimallyFrst=as.numeric(AdensTable[,6])
  Total=as.numeric(AdensTable[,3])

  PositionsTotal=seq(from=1, to=length(Positions), by=1)

  datos<-read.table(file=paste(Dir,"FrustrationData/", JobID, ".pdb_", Pdb$mode ,sep=""),stringsAsFactors = F)

  if(!is.null(Chain)){
    datos<-datos[datos$V3==Chain,]
  }

  chains<-sort(unique(c(datos$V3,datos$V4)))
  positions<-matrix(ncol=3,nrow=length(chains))
  auxPosVec<-c()
  for(i in seq_along(chains)){
    positions[i,1:2]<-range(c(datos$V1[which(datos$V3==chains[i])],datos$V2[which(datos$V4==chains[i])]))
    positions[i,3]<-positions[i,2]-positions[i,1]+1
    auxPosVec<-c(auxPosVec,positions[i,1]:positions[i,2])
  }

  datos$pos1 <- NA
  datos$pos2 <- NA
  for(i in seq_along(chains)){
    if(i==1){bias <- 0}else{bias <- sum(positions[1:(i-1),3])}
    idx <- which(datos$V3==chains[i] )
    datos$pos1[idx] <- datos$V1[idx]-positions[i,1]+bias+1
    idx <- which(datos$V4==chains[i])
    datos$pos2[idx] <- datos$V2[idx]-positions[i,1]+bias+1
  }

  posNEW<-matrix(ncol=3,nrow=length(chains))
  for(i in seq_along(chains)){
    posNEW[i,1:2]<-range(c(datos$pos1[which(datos$V3==chains[i])],datos$pos2[which(datos$V4==chains[i])]))
    posNEW[i,3]<-posNEW[i,2]-posNEW[i,1]+1
  }

  total.positions<-sum(apply(positions,1,function(x){x[2]-x[1]+1}))
  matrz <- matrix(NA,ncol=total.positions,nrow=total.positions)
  for(i in 1:nrow(datos)){
    if(datos$V13[i]=="short") matrz[datos$pos1[i],datos$pos2[i]]<-20
    else if(datos$V13[i]=="long") matrz[datos$pos1[i],datos$pos2[i]]<-30
    else matrz[datos$pos1[i],datos$pos2[i]]<-40
    matrz[datos$pos2[i],datos$pos1[i]]<-datos$V12[i]
  }
  #PARA HACERLO TRIANGULAR SUPERIOR
  #matrz[upper.tri(matrz,diag = T)] <- 0

  longData<-melt(matrz)
  # Graphic<-ggplot() + geom_tile(data=longData[longData$value<20& !is.na(longData$value),],aes(x = Var2, y = Var1, fill=value))+
  # scale_fill_gradient2(low="red",mid="grey",high="green",limits=c(-4,4),breaks=c(-4,-3,-2,-1,0,1,2,3,4),labels=c("-4","-3","-2","-1","0","1","2","3","4"),
  # guide=guide_colourbar(title=paste("Local ",Pdb$mode," Frustration Index"),barwidth=2,barheight=20,title.position="right",title.hjust=0.5)) +
  # geom_tile(data=longData[longData$value>=20& !is.na(longData$value),],aes(x = Var2, y = Var1, color=value))+
  # scale_color_discrete("",breaks=c(20,30,40),labels=c("Short","Long","Water mediated"),values = c("yellow","violet","black"))+
  # labs(x="Residue i", y="Residue j") + ggtitle(paste("Contact map ",Pdb$PdbBase,sep=""))+
  #   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
  #                      axis.text.y=element_text(size=9),
  #                      plot.title=element_text(size=11,hjust=0.5),
  #                      legend.title=element_text(angle=-90))

  Graphic<-ggplot() + geom_tile(data=longData[longData$value<20& !is.na(longData$value),],aes(x = Var2, y = Var1, fill=value))+
    scale_fill_gradient2(low="red",mid="grey",high="green",limits=c(-4,4),breaks=c(-4,-3,-2,-1,0,1,2,3,4),labels=c("-4","-3","-2","-1","0","1","2","3","4"),
                         guide=guide_colourbar(barwidth=2,barheight=20,title=paste("Local ",Pdb$mode," Frustration Index"),title.position="right",title.hjust=0.5))+
    labs(x="Residue i", y="Residue j") + ggtitle(paste("Contact map ",Pdb$PdbBase,sep=""))+
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11,hjust=0.5),
                         legend.title=element_text(angle = -90))+
    geom_point(data=longData[longData$value>=20& !is.na(longData$value),],aes(x = Var2, y = Var1, color=as.factor(value)),shape=15,size=1)+
    scale_color_manual("Contact distance",breaks=c(20,30,40),labels=c("Short","Long","Water mediated"),values = c("black","gray","cyan"),guide=guide_legend( title.theme = element_text(angle = 360)))


  if(length(chains)>1){
    breaks <- round(seq(1,total.positions,total.positions/15))
    labels <- as.character(auxPosVec[breaks])
    Graphic <- Graphic + geom_vline(xintercept=cumsum(posNEW[,3]),color="gray",linetype="longdash")+
    geom_hline(yintercept=cumsum(posNEW[,3]),color="gray",linetype="longdash")+
    annotate("text", x = round(apply(posNEW[,1:2],1,mean)), y = total.positions+20, label = unique(chains),parse = TRUE,color="gray")+
    annotate("text", x = total.positions+20, y = round(apply(posNEW[,1:2],1,mean)), label = unique(chains),parse = TRUE,color="gray")+
    scale_x_continuous(breaks=breaks,labels=labels) + scale_y_continuous(breaks=breaks,labels=labels)
  }

  ggsave(paste(Dir, "Images/", JobID, "_", mode, "_map.png", sep=""),plot=Graphic)
}

#' View Frustration Pymol
#'
#' Generates a pymol session to observe the frustration patterns on top of the pdb protein structure
#'
#' @param Pdb Frustration object
#' @export
#'
view_frustration_pymol <- function(pdb)
{
  system(paste("cd ", Pdb$JobDir,"VisualizationScrips/",  "; pymol ", Pdb$PdbBase, ".pdb_", Pdb$mode, ".pml", sep=""))
}


#' View Dynamic Frustration of Resid
#'
#' Generate graphs corresponding to the frustration of the residue along the dynamics
#'
#' @param DynDir Full path of the table resulting from the execution of dynamic_frustration for the corresponding residue
#' @param Resno Resno of residue to analyze.
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#' @export
#'

plot_dynamic_res <- function(DynDir=NULL, Resno=NULL, Modes="configurational"){

  if(is.null(DynDir)) stop("DynDir not specified")
  if(is.null(Resno)) stop("Resno not specified")

  FrustrationResults <- read.table( paste(DynDir,"/",Modes,"_Res",Resno,sep=""), header=T)
  FrustrationResults <- as.data.frame(FrustrationResults)

  if(Modes == "configurational" | Modes == "mutational"){

    FrustrationResults <- cbind(seq(1, dim(FrustrationResults)[1]),FrustrationResults)
    colnames(FrustrationResults)<-c("Frame","MaximallyFrst","NeutrallyFrst","MinimallyFrst")
  }else {

    FrustrationResults <- cbind(seq(1, dim(FrustrationResults)[1]),FrustrationResults[1])
    FrustrationResults[ which(FrustrationResults[,2]>=0.58) ,3] <- "Minimally frustrated"
    FrustrationResults[ which( (FrustrationResults[,2]<0.58) & (FrustrationResults[,2]>-1.0) )  ,3] <- "Neutral"
    FrustrationResults[ which(FrustrationResults[,2]<=(-1.0)) ,3] <- "Highly frustrated"
    colnames(FrustrationResults)<-c("Frame","IndexFrst","Type")
    FrustrationResults$Type<-factor(x=FrustrationResults$Type,levels=c("Minimally frustrated","Neutral","Highly frustrated"))
  }


  if(Modes == "configurational" | Modes == "mutational"){

    #************GRAFICA 5ADENS POR FRAME**************

    Maximum=max(c(max(FrustrationResults$MaximallyFrst),max(FrustrationResults$MinimallyFrst),max(FrustrationResults$NeutrallyFrst)))

    adens5Frame <- ggplot() + geom_line( aes( x = FrustrationResults$Frame , y = FrustrationResults$MaximallyFrst , colour = "1" )  )
    adens5Frame <- adens5Frame + geom_line( aes( x = FrustrationResults$Frame , y = FrustrationResults$NeutrallyFrst , colour = "2" ) )
    adens5Frame <- adens5Frame + geom_line( aes( x = FrustrationResults$Frame , y = FrustrationResults$MinimallyFrst , colour =  "3" ) )
    adens5Frame <- adens5Frame + ylab( "Local frustration density (5A sphere)" ) + xlab( "Frame" ) + ggtitle(paste("Frustration of resid",Resno, Modes))
    adens5Frame <- adens5Frame + scale_colour_manual( name="",labels = c("Highly frustrated","Neutral","Minimally frustrated" ) , values = c( "red" ,"gray","green"))
    adens5Frame <- adens5Frame + scale_y_continuous(limits=c(0,1), breaks = seq(0.0,1.0,0.2),labels =as.character(seq(0.0,1.0,0.2)))
    adens5Frame <- adens5Frame + scale_x_continuous( breaks = seq( 1 , length(FrustrationResults$Frame),trunc(length(FrustrationResults$Frame)*0.05) ))
    adens5Frame <- adens5Frame + theme_classic() + theme(plot.title=element_text(size=11,hjust=0.5),axis.text.x = element_text(angle = 90))

    ggsave(paste(DynDir,"/dynamic5adens_",Modes,"_Res",Resno,".png", sep=""),plot=adens5Frame,width=10, height= 6)

    #************GRAFICA HISTOGRAMA 5ADENS**************

    png(filename = paste(DynDir,"/dynamicHist5adens_",Modes,"_Res",Resno,".png", sep=""), width = 540, height = 420)
    par(mar=c(5,5,3,2), xpd=TRUE)
    par(mgp=c(2,0.5,0))
    barplot(rbind(FrustrationResults$MinimallyFrst, FrustrationResults$NeutrallyFrst, FrustrationResults$MaximallyFrst), col=c("green", "gray", "red"), axis.lty=1, xlab="Frame", ylab="Density arround 5A sphere (%)", cex.lab=1.5, border = NA, space = c(0), xaxt = "n")
    AuxAxis=seq(from=0, to=dim(FrustrationResults)[1])
    axis(side = 1, at = AuxAxis , labels = AuxAxis,  tick = FALSE)
    box(lwd=2)
    legend(x="top",inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated"), pch=c(15, 15, 15), col=c("green", "gray", "red"), horiz = T, bty = "n")
    dev.off()

  }else{

    #************GRAFICA INDEX FRUSTRATION POR FRAME**************

    IndexFrusFrame <- ggplot( data = FrustrationResults,aes( x = Frame , y = IndexFrst,colour=Type )) + geom_point()
    IndexFrusFrame <- IndexFrusFrame + ylab( "Index Frustration" ) + xlab( "Frame" ) + ggtitle(paste("Index Frustration of resid",Resno, "Single residue"))
    cols <- c("Minimally frustrated" = "green", "Neutral" = "gray", "Highly frustrated" = "red")
    IndexFrusFrame <- IndexFrusFrame + scale_colour_manual( name="",values=cols)+ scale_y_reverse(limits=c(4,-4), breaks=seq(4.0, -4.0, -0.5),labels=as.character(seq(4.0, -4.0, -0.5)))
    IndexFrusFrame <- IndexFrusFrame + scale_x_continuous( breaks = seq( 1 , length(FrustrationResults$Frame),trunc(length(FrustrationResults$Frame)*0.05) ) )
    IndexFrusFrame <- IndexFrusFrame + theme_classic() + theme(plot.title=element_text(size=11,hjust=0.5),axis.text.x = element_text(angle = 90))

    ggsave(paste(DynDir,"/dynamic_IndexFrustration_",Modes,"_Res",Resno,".png", sep=""),plot=IndexFrusFrame,width=10, height= 6)
  }

}

#' Contact map gif
#'
#' Generate Gif of contact map of the structure under study along the dynamics
#'
#' @param PdbDir File containing all proteins structures. The full path to the file is needed
#' @param OrderList Orderer list of PDB files to calculate frustration.
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#'
#' @export
#'
#'
gif_contact_map<-function(PdbDir=NULL, OrderList=NULL,Modes="configurational"){

  if(is.null(PdbDir)) stop("PdbDir not specified")
  if(is.null(OrderList)) stop("OrderList not specified")

  paths<-c()
  for (i in seq(1,length(OrderList))) {
    paths[i]<-paste(PdbDir,basename.pdb(OrderList[i]),".done/Images/",basename.pdb(OrderList[i]),"_",Modes,"_map.png",sep="")
  }
  images<-list()
  for (i in seq(1,length(OrderList))) {
    images[[i]]<-image_read(paths[i])
  }
  images <- image_join(images)
  animation <- image_animate(images, fps = 2, optimize = TRUE)
  image_write(animation,paste(PdbDir,"contactMap_",Modes,".gif",sep=""))

}

#' 5Adens proportions gif
#'
#' Generate Gif of 5Adens proportions plot of the structure under study along the dynamics
#'
#' @param PdbDir File containing all proteins structures. The full path to the file is needed
#' @param OrderList Orderer list of PDB files to calculate frustration.
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#'
#' @export
#'
gif_5adens_proportions<-function(PdbDir=NULL, OrderList=NULL,Modes="configurational"){

  if(is.null(PdbDir)) stop("PdbDir not specified")
  if(is.null(OrderList)) stop("OrderList not specified")

  paths<-c()
  for (i in seq(1,length(OrderList))) {
    paths[i]<-paste(PdbDir,basename.pdb(OrderList[i]),".done/Images/",basename.pdb(OrderList[i]),"_",Modes,"_5Adens_around.png",sep="")
  }
  images<-list()
  for (i in seq(1,length(OrderList))) {
    images[[i]]<-image_read(paths[i])
  }

  images = image_join(images)

  animation <- image_animate(images, fps = 2,optimize = TRUE)
  image_write(image = animation,path = paste(PdbDir,"5Adens_proportions_",Modes,".gif",sep=""))

}

#' Frustration movie in Pymol
#'
#' It generates visualization in pymol of the frames of the indicated dynamics and script for its execution in pymol of generation of the film
#'
#' @param PdbDir File containing all proteins structures. The full path to the file is needed
#' @param OrderList Orderer list of PDB files to calculate frustration.
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default: configurational
#'
#' @export
#'
frustra_movie<-function(PdbDir=NULL,OrderList=NULL,Modes="configurational"){

  if(is.null(PdbDir)) stop("PdbDir not specified")
  if(is.null(OrderList)) stop("OrderList not specified")

  ResultDir<-paste(PdbDir,"FrustraMovie",sep="")
  system(paste("mkdir -p ",ResultDir,sep=""))

  fileChargePml<-paste(ResultDir,"/representations_",Modes,".pml",sep="")
  i<-1
  for (Pdb in OrderList) {
    write(paste("load ",PdbDir,basename.pdb(Pdb),".done/VisualizationScrips/",basename.pdb(Pdb),".pdb_",Modes,".pml",sep=""),file=fileChargePml,append=TRUE)
    write(paste("group frame",i,",*",basename.pdb(Pdb),sep=""),file=fileChargePml,append=TRUE)
    i<-i+1
  }

  fileMovie<-paste(ResultDir,"/GenerateMovie_",Modes,".pml",sep="")

  write(paste("mset 1x",length(OrderList),sep=""),file=fileMovie,append=TRUE)
  for (Nframe in 1:length(OrderList)) {
    write("hide all",file=fileMovie,append=TRUE)
    write(paste("show cgo, frame",Nframe,sep=""),file=fileMovie,append=TRUE)
    write(paste("show dashes, frame",Nframe,sep=""),file=fileMovie,append=TRUE)
    write(paste("show cartoon, frame",Nframe,sep=""),file=fileMovie,append=TRUE)
    write(paste("scene F",Nframe,",store",sep=""),file=fileMovie,append=TRUE)
    write(paste("mview store,",Nframe,",scene=F",Nframe,sep=""),file=fileMovie,append=TRUE)
  }
  write("show all",file=fileMovie,append=TRUE)
  write("show cgo,all",file=fileMovie,append=TRUE)
  write("show dashes,all",file=fileMovie,append=TRUE)
  write("show cartoon,all",file=fileMovie,append=TRUE)
}

#' Delta frustration of mutated residue
#'
#' Generate a graph of the singleresidue frustration difference for the mutation of the n mutated residues by using mutate_res
#'
#' @param PdbPath Full path of the Pdb to analyze
#' @param DataDir Full directory path where the singleresidue tables generated by mutate_res are located for the n residues analyzed
#' @param ResultDir Full path of the directory where the resulting graphic will be stored
#' @param Chain Chain of residue to analyze. Default: NULL
#'
#' @export
#'

plot_delta_frus<-function(PdbPath=NULL,DataDir=NULL,ResultDir=NULL,Chain=NULL){


  if(is.null(PdbPath)) stop("PdbPath not specified")
  if(is.null(DataDir)) stop("DataDir not specified")
  if(is.null(ResultDir)) stop("ResultDir not specified")

  Pdb<-read.pdb(PdbPath)

  system(paste("cd ",DataDir,"; cat singleresidue_Res*.txt > Data.txt;",sep=""))
  DataFrus=read.table(paste(DataDir,"/Data.txt",sep=""), header=F, stringsAsFactors=F)
  system(paste("cd ",DataDir,"; rm Data.txt;"))
  DataFrus<-as.data.frame(DataFrus)

  #Native residues
  Natives<-data.frame(stringsAsFactors = FALSE)
  Natives<-cbind(Pdb$atom[atom.select(Pdb,resno =unique(DataFrus$V1),elety="CA",chain=Chain)$atom,"resid"],Pdb$atom[atom.select(Pdb,resno =unique(DataFrus$V1),elety="CA",chain=Chain)$atom,"resno"])
  Natives[,1]<-aa321(Natives[,1])

  DataFrus<-cbind(DataFrus,1:length(DataFrus[,1]))
  DataFrus[DataFrus$V4>0.55,5]<-"green"
  DataFrus[DataFrus$V4<=-1.0,5]<-"red"
  DataFrus[DataFrus$V4>=-1.0&DataFrus$V4<=0.55,5]<-"gray"
  for (i in 1:dim(Natives)[1]) {
    DataFrus[Natives[i,2]==DataFrus$V1 & DataFrus$V2==Chain & DataFrus$V3==Natives[i,1] ,5]<-"blue"
  }

  DataFrus[,5]<-as.factor(DataFrus[,5])

  DataFrus<-cbind(DataFrus,seq(1,length(DataFrus[,1])))
  colnames(DataFrus)<-c("Res", "ChainRes", "AA", "FrstIndex","State Of Frustration","DeltaFrus")

  #DeltaFrus
  NativesFrus<-DataFrus[DataFrus$`State Of Frustration`=="blue",]
  DataMuts<-DataFrus[!DataFrus$`State Of Frustration`=="blue",]

  delta<-rep(0.0,length(NativesFrus$FrstIndex))-NativesFrus$FrstIndex
  diff<-data.frame(NativesFrus$Res,NativesFrus$ChainRes,delta)
  colnames(diff)<-c("Res","ChainRes","Diff")

  NativesFrus$DeltaFrus<-NativesFrus$FrstIndex+diff$Diff

  for (i in 1:length(DataMuts[,1])) {
    DataMuts[i,6]<-NativesFrus[NativesFrus$Res==DataMuts$Res[i]&NativesFrus$ChainRes==DataMuts$ChainRes[i],6]-(DataMuts[i,4]+diff[diff$Res==DataMuts$Res[i]&diff$ChainRes==DataMuts$ChainRes[i],3])
  }
  DataFrus<-rbind(DataMuts,NativesFrus)

  x1<-DataFrus$Res[1]-3
  x2<-DataFrus$Res[1]+3
  y1<--4
  y2<-4
  if(min(DataFrus$DeltaFrus)<(-4))  y1<-min(DataFrus$DeltaFrus)
  if(max(DataFrus$DeltaFrus)>(4))  y2<-max(DataFrus$DeltaFrus)

  Residues<-unique(DataFrus$Res)
  Residues<-Residues[order(Residues[])]
  Residues<-cbind(Residues,1:length(Residues))
  Residues<-as.data.frame(Residues)
  colnames(Residues)<-c("Res","Index")

  DataFrus$Res<-Residues[match(DataFrus$Res,Residues$Res),2]

  Resid<-atom.select(Pdb,resno=Residues$Res,elety="CA",chain=Chain,value=T)$atom$resid

  Graphic<-ggplot(data=DataFrus)+geom_point(aes(x=Res,y=DeltaFrus,color=`State Of Frustration`),shape=DataFrus$AA,size=3)
  Graphic<-Graphic+ xlab("Residue Position") + ylab("Delta Frustration") +  scale_x_continuous(breaks=Residues$Index,labels=paste(Resid,Residues$Res))
  Graphic<-Graphic+ scale_y_continuous(breaks=seq(y1,y2,0.5),labels=as.character(seq(y1,y2,0.5)),limits = c(y1,y2))
  Graphic<-Graphic+ scale_color_manual("Frustration",breaks=levels(DataFrus$`State Of Frustration`),labels=c("Native state","Neutral" , "Minimally frustrated","Highly frustrated" ),values=c("blue","gray","green","red"))
  Graphic<-Graphic+ theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))
  ggsave(plot=Graphic,paste(ResultDir,"/Delta_frus.png",sep=""),width=10, height= 6)

}

#' Frustration of mutated residue
#'
#' Generate a graph of the configurational or mutational frustration difference for the mutation of the n mutated residues by using mutate_res
#'
#' @param PdbPath Full path of the Pdb to analyze
#' @param DataFile Full path of the table resulting from the execution of mutate_res to parse
#' @param ResultDir Full path of the directory where the resulting graphic will be stored
#' @param Chain Chain of residue to analyze. Default: NULL
#' @param Modes Local frustration index to be calculated (configurational, mutational). Default: configurational
#' @param Resno Resno of residue to analyze.
#'
#' @export
#'

plot_mutate_res<-function(PdbPath=NULL,DataFile=NULL,ResultDir=NULL,Chain=NULL,Modes="configurational",Resno=NULL){

  if(is.null(PdbPath)) stop("PdbPath not specified")
  if(is.null(DataFile)) stop("DataFile not specified")
  if(is.null(ResultDir)) stop("ResultDir not specified")
  if(is.null(Chain)) stop("Chain not specified")
  if(is.null(Resno)) stop("Resno not specified")

  Pdb<-read.pdb(PdbPath)

  DataFrus=read.table(DataFile, header=F,stringsAsFactors=F)
  DataFrus<-as.data.frame(DataFrus)
  DataFrus<-cbind(DataFrus,seq(1,length(DataFrus[,1])))
  colnames(DataFrus)<-c("Res1","Res2","Chain1","Chain2","AA1","AA2","FrstIndex","FrstState","Color")
  DataFrus$FrstState<-as.factor(DataFrus$FrstState)

  Native<-atom.select(Pdb,resno=Resno,chain=Chain,value=TRUE)
  Native<-aa321(unique(Native$atom$resid))

  DataFrus[DataFrus$Res2==Resno,c(2,4,5,6)]<-c(DataFrus$Res1[DataFrus$Res2==Resno],DataFrus$Chain1[DataFrus$Res2==Resno],DataFrus$AA2[DataFrus$Res2==Resno],DataFrus$AA1[DataFrus$Res2==Resno])
  #Agregue por las cadenas
  DataFrus$Chain1[DataFrus$Chain1!=Chain]<-Chain
  DataFrus$Res1[DataFrus$Res1!=Resno]<-Resno

  DataFrus$Color[DataFrus$FrstState=="neutral"]<-"gray"
  DataFrus$Color[DataFrus$FrstState=="highly"]<-"red"
  DataFrus$Color[DataFrus$FrstState=="minimally"]<-"green"
  DataFrus$Color[DataFrus$AA1==Native]<-"blue"
  #Agregar color orange de las glicinas aca
  #DataFrus$Color[DataFrus$AA1=='G']<-"orange"

  DataFrus$Color<-as.factor(DataFrus$Color)

  Contacts<-unique(DataFrus[,c("Res2","Chain2")])
  Contacts<-Contacts[order(Contacts[,1]),]
  Contacts<-cbind(Contacts,1:length(Contacts[,1]))
  Contacts<-as.data.frame(Contacts)
  colnames(Contacts)<-c("Res","Chain","Index")

  Contacts$Res<-as.numeric(Contacts$Res)
  Resid<-atom.select(Pdb,resno=Contacts$Res,elety="CA",chain=Chain,value=T)$atom$resid

  DataFrus$Res2<-Contacts[match(DataFrus$Res2,Contacts$Res),3]

  DataFrus<-rbind(DataFrus[DataFrus$Color!="orange",],DataFrus[DataFrus$Color=="orange",])
  DataFrus<-rbind(DataFrus[DataFrus$Color!="blue",],DataFrus[DataFrus$Color=="blue",])

  y1<-(-4)
  y2<-4
  Graphic<-ggplot(data=DataFrus)+geom_point(aes(x=Res2,y=FrstIndex,color=Color),shape=DataFrus$AA1,size=3)
  Graphic<-Graphic+ xlab("Contact residue") + ylab("Frustration Index")
  Graphic<-Graphic+ scale_x_continuous(breaks=Contacts$Index,labels=paste(Resid,Contacts$Res,Contacts$Chain))+ scale_y_continuous(breaks=seq(y1,y2,0.5),labels=as.character(seq(y1,y2,0.5)),limits = c(y1,y2))
  Graphic<-Graphic+ geom_hline(yintercept=c(0.78,-1),color="gray",linetype="longdash")
  Graphic<-Graphic+ ggtitle(paste("Contact Frustration ",Modes," of residue ",DataFrus$Res1[1],sep=""))+ theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))
  Graphic<-Graphic+ scale_color_manual("",breaks=c("green","gray","red","blue","orange"),labels=c("Minimally frustrated","Neutral","Highly frustrated","Native","Glycine"),values=c("green","gray","red","blue","orange"))


  ggsave(plot=Graphic,paste(ResultDir,"/",Modes,"_",DataFrus$Res1[1],".png",sep=""),width=10, height= 6)

}
