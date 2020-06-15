#!/usr/bin/env Rscript

args= commandArgs(trailingOnly=T)

JobID=args[1];
Dir=args[2];
mode=args[3];

AdensTable=read.table(file=paste(Dir,"/", JobID,".pdb_", mode, "_5adens", sep=""))


Positions=as.numeric(AdensTable[,1])
Chains=AdensTable[,2]
MaximallyFrst= as.numeric(AdensTable[,4])
NeutrallyFrst=as.numeric(AdensTable[,5])
MinimallyFrst=as.numeric(AdensTable[,6])
Total=as.numeric(AdensTable[,3])

PositionsTotal=seq(from=1, to=length(Positions), by=1)

#Primer Grafico
png(filename = paste(Dir, "/", JobID, "_", mode, ".png_5Adens", ".png", sep=""), width = 540, height = 420)
Maximum=max(c(max(MaximallyFrst),max(MinimallyFrst),max(NeutrallyFrst), max(Total) ))
par(mar=c(5,5,3,2), xpd=TRUE)
plot(PositionsTotal, Total, type ="l", col="black", lwd=1, ylab="Local frustration density (5A sphere)", xlab="Position", ylim=(c(0, Maximum)), cex.axis=1.5, cex.lab=1.5)
lines(PositionsTotal, NeutrallyFrst, type ="l", col="gray", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
lines(PositionsTotal, MaximallyFrst, type ="l", col="red", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
lines(PositionsTotal, MinimallyFrst, type ="l", col="green", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
legend(x="top",  inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated", "total"), pch=c(15, 15, 15, 15), col=c("green", "gray", "red", "black"), horiz = T, x.intersp=0.5, xjust=0,yjust=0, bty = "n")
box(lwd=2)
dev.off()

#Primer grafico por chains
UChains=as.character(unique(Chains))
for(i in UChains)
{
  png(filename = paste(Dir, "/", JobID, "_", mode, "_5Adens_chain", i, ".png", sep=""), width = 540, height = 420)
  Maximum=max(c(max(MaximallyFrst[which(Chains==i)]),max(MinimallyFrst[which(Chains==i)]),max(NeutrallyFrst[which(Chains==i)]), max(Total[which(Chains==i)]) ))
  par(mar=c(5,5,3,2), xpd=TRUE)
  plot(Positions[which(Chains==i)], Total[which(Chains==i)], type ="l", col="black", lwd=1, ylab="Local frustration density (5A sphere)", xlab="Position", ylim=(c(0, Maximum)), cex.axis=1.5, cex.lab=1.5)
  lines(Positions[which(Chains==i)], NeutrallyFrst[which(Chains==i)], type ="l", col="gray", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
  lines(Positions[which(Chains==i)], MaximallyFrst[which(Chains==i)], type ="l", col="red", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
  lines(Positions[which(Chains==i)], MinimallyFrst[which(Chains==i)], type ="l", col="green", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
  legend(x="top",  inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated", "total"), pch=c(15, 15, 15, 15), col=c("green", "gray", "red", "black"), horiz = T, x.intersp=0.5, xjust=0,yjust=0, bty = "n")
  box(lwd=2)
  dev.off()
}

#Segundo Grafico total
MinimallyFrst=as.numeric(AdensTable[,9])
NeutrallyFrst=as.numeric(AdensTable[,8])
MaximallyFrst=as.numeric(AdensTable[,7])

png(filename = paste(Dir, "/", JobID, "_", mode, "_5Adens_around.png", sep=""), width = 540, height = 420)
par(mar=c(5,5,3,2), xpd=TRUE)
par(mgp=c(2,0.5,0))
barplot(rbind(MinimallyFrst, NeutrallyFrst, MaximallyFrst), col=c("green", "gray", "red"), axis.lty=1, xlab="Position", ylab="Density arround 5A sphere (%)", cex.lab=1.5, border = NA, space = c(0), xaxt = "n")
AuxAxis=seq(from=0, to=length(PositionsTotal))
axis(side = 1, at = AuxAxis , labels = (AuxAxis+min(PositionsTotal)),  tick = FALSE)
box(lwd=2)
legend(x="top",inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated"), pch=c(15, 15, 15), col=c("green", "gray", "red"), horiz = T, bty = "n")
dev.off()

#Segundo grafico por chain
for(i in UChains)
{
  png(filename = paste(Dir, "/", JobID, "_", mode, "_5Adens_around_chain", i, ".png", sep=""), width = 540, height = 420)
  par(mar=c(5,5,3,2), xpd=TRUE)
  par(mgp=c(2,0.5,0))
  barplot(rbind(MinimallyFrst[which(Chains==i)], NeutrallyFrst[which(Chains==i)], MaximallyFrst[which(Chains==i)]), col=c("green", "gray", "red"), axis.lty=1, xlab="Position", ylab="Density arround 5A sphere (%)", cex.lab=1.5, border = NA, space = c(0), xaxt = "n")
  AuxAxis=seq(from=0, to=length(Positions[which(Chains==i)]))
  axis(side = 1, at = AuxAxis , labels = (AuxAxis+min(Positions[which(Chains==i)])),  tick = FALSE)
  box(lwd=2)
  legend(x="top",inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated"), pch=c(15, 15, 15), col=c("green", "gray", "red"), horiz = T, bty = "n")
  dev.off()
}

# Grafico Contact Map gral

datos<-read.table(file=paste(Dir,"/", JobID, "_", mode ,sep=""),stringsAsFactors = F)

chains<-sort(unique(c(datos$V3,datos$V4)))
positions<-matrix(ncol=3,nrow=length(chains))
auxPosVec<-c()
for(i in seq_along(chains)){
  positions[i,1:2]<-range(c(datos$V1[which(datos$V3==chains[i])],datos$V2[which(datos$V4==chains[i])]))
  positions[i,3]<-positions[i,2]-positions[i,1]+1
  auxPosVec<-c(auxPosVec,positions[i,1]:positions[i,2])
}

#Renumero posiciones
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
  matrz[datos$pos1[i],datos$pos2[i]]<-datos$V12[i]
  matrz[datos$pos2[i],datos$pos1[i]]<-datos$V12[i]
  
}

png(filename = paste(Dir, "/", JobID, "_", mode, "_map.png", sep=""),width = 600,height = 540)
cotaM <- 4
n.breaks=50
vecvalores <- seq(-cotaM,cotaM,length.out = 2*n.breaks)
# vecvalores <- c(vecvalores[1:n.breaks],0,vecvalores[(n.breaks+1):(2*n.breaks)])
colores.paleta<-colorRampPalette(c("red","grey","green"))
colores<-colores.paleta(length(vecvalores)-1)
par(mar=c(5,5,5,5),cex.lab=1.5,cex.axis=1.5)
layout(mat = matrix(1:2,ncol=2),widths = c(5,1))
image(x=1:total.positions,y=1:total.positions,matrz,axes=F,xlab="Residue i",ylab="Residue j",col=colores,asp=1,breaks = vecvalores)

if(length(chains)>1){
  abline(v=cumsum(posNEW[,3])+.5,lty=3,col="darkgrey",lwd=1.5)
  abline(h=cumsum(posNEW[,3])+.5,lty=3,col="darkgrey",lwd=1.5)
  axis(side=3,at=apply(posNEW[,1:2],1,mean),tck=0,las=1,lwd=2,labels = chains,tick = 0)
  axis(side=4,at=apply(posNEW[,1:2],1,mean),tck=0,las=1,lwd=2,labels = chains,tick = 0)
  mtext(text = "Chain",line=3,side=3,cex=1.5)
  mtext(text = "Chain",line=3,side=4,cex=1.5)
}
box(lwd=2)

ticksPos <- axTicks(1)[which(axTicks(1)!=0)]
axis(side=1,at=ticksPos,labels = auxPosVec[ticksPos],tcl=.5,lwd=2)
axis(side=2,at=ticksPos,labels = auxPosVec[ticksPos],tcl=.5,lwd=2)
axis(side=3,at=ticksPos,labels = NA,tcl=.5,lwd=2)
axis(side=4,at=ticksPos,labels = NA,tcl=.5,lwd=2)

par(mar=c(5,1,5,5))
image(y=vecvalores,z=matrix(vecvalores,nrow=1),col=colores,axes=F,ylab="")
axis(side=4,at=seq(-3,3,by=1),tcl=0.3,las=1,lwd=2)
axis(side=4,at=c(-4,4),tcl=0,las=1,lwd=2)
mtext(text = "Local Configurational Frustration Index",line=3,side=4,cex=1.5)
axis(side=2,at=seq(-3,3,by=1),tcl=0.3,labels = NA,lwd=2)

box(lwd=2)

dev.off()

