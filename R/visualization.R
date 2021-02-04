#plot_5Andens----
#' @title Plot 5Adens.
#'
#' @description Generates plot to analyze the density of contacts around a sphere of 5 Armstrongs, centered in the C-alfa atom from the residue.
#'  The different classes of contacts based on the mutational frustration index are counted in absolute terms.
#'
#' @param Pdb Pdb Frustration object.
#' @param Chain Chain of residue. Type: character. Default: NULL.
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#' 
#' @import ggplot2
#' 
#' @export
plot_5Andens <- function(Pdb, Chain = NULL, Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  if(!is.null(Chain)){
    if(length(Chain) > 1)
      stop("You must enter only one Chain!")
    if(!all(Chain %in% unique(Pdb$atom$chain))) 
      stop(paste("The Chain ", Chain, " doesn't exist! The Chains are: ",
                                                          paste(unique(Pdb$atom$chain), collapse = " ", sep = ""), sep = ""))
  }
  if(!dir.exists(paste(Pdb$JobDir, "/Images", sep = "")))  
    dir.create(paste(Pdb$JobDir, "/Images", sep = ""))
  
  colClasses = c("integer", "character", "integer", 
                 "integer", "integer", "integer",
                 "numeric", "numeric", "numeric")
  AdensTable <- read.table(file = paste(Pdb$JobDir,"FrustrationData/", Pdb$PdbBase, ".pdb_", Pdb$Mode, "_5adens", sep = ""), header = T, colClasses = colClasses)
  AdensTable <- as.data.frame(AdensTable)
  colnames(AdensTable) <- c("Positions", "Chains", "Total", "MaximallyFrst", "NeutrallyFrst", "MinimallyFrst")
  PositionsTotal <- seq(from = 1, to = length(AdensTable$Positions), by = 1)
  AdensTable <- cbind(AdensTable, PositionsTotal)
  
  if(is.null(Chain)){
    
    Maximum = max(c(max(AdensTable$MaximallyFrst), max(AdensTable$MinimallyFrst),
                    max(AdensTable$NeutrallyFrst), max(AdensTable$Total)))
    
    Graphic <- ggplot() + geom_line(aes(x = AdensTable$PositionsTotal, y = AdensTable$MaximallyFrst , colour = "a"))
    Graphic <- Graphic + geom_line(aes(x = AdensTable$PositionsTotal, y = AdensTable$NeutrallyFrst , colour = "b" ))
    Graphic <- Graphic + geom_line(aes(x = AdensTable$PositionsTotal, y = AdensTable$MinimallyFrst , colour =  "c" ))
    Graphic <- Graphic + geom_line(aes(x = AdensTable$PositionsTotal, y = AdensTable$Total , colour = "d"))
    Graphic <- Graphic + ggtitle(paste("Density arround 5A sphere (%) in", Pdb$PdbBase))
    Graphic <- Graphic + ylab("Local frustration density (5A sphere)") + xlab("Position")
    Graphic <- Graphic + scale_color_manual(name = "", labels = c("Highly frustrated", "Neutral", "Minimally frustrated", "Total"),
                                            values = c("red", "gray", "green", "black"))
    Graphic <- Graphic + scale_y_continuous(breaks = seq(0, Maximum, 5)) + scale_x_continuous(breaks = seq(1 , length(AdensTable$PositionsTotal), trunc((length(AdensTable$PositionsTotal) - 1) / 10)))
    Graphic <- Graphic + theme(plot.title = element_text(size = 11, hjust = 0.5), panel.background = element_blank())
    
    if(Save){
      ggsave(filename = paste(Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, ".png_5Adens", ".png", sep = ""), plot = Graphic, width = 10, height = 6)
      cat(paste("5Adens plot is stored in ", Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, ".png_5Adens", ".png\n", sep = ""))
    }
    
  }else{
    AdensTable <- AdensTable[AdensTable$Chains == Chain, ]
    
    Maximum = max(c(max(AdensTable$MaximallyFrst), max(AdensTable$MinimallyFrst),
                    max(AdensTable$NeutrallyFrst), max(AdensTable$Total)))
    
    Graphic <- ggplot() + geom_line(aes(x = AdensTable$Positions, y = AdensTable$MaximallyFrst , colour = "a"))
    Graphic <- Graphic + geom_line(aes(x = AdensTable$Positions, y = AdensTable$NeutrallyFrst , colour = "b"))
    Graphic <- Graphic + geom_line(aes(x = AdensTable$Positions, y = AdensTable$MinimallyFrst , colour =  "c"))
    Graphic <- Graphic + geom_line(aes(x = AdensTable$Positions, y = AdensTable$Total , colour = "d"))
    Graphic <- Graphic + ggtitle(paste("Density arround 5A sphere (%) in", Pdb$PdbBase, "chain", Chain))
    Graphic <- Graphic + ylab("Local frustration density (5A sphere)") + xlab("Position")
    Graphic <- Graphic + scale_color_manual(name = "", labels = c("Highly frustrated", "Neutral", "Minimally frustrated", "Total"),
                                            values = c("red", "gray", "green", "black"))
    Graphic <- Graphic + scale_y_continuous(breaks = seq(0, Maximum, 5)) + scale_x_continuous(breaks = seq(AdensTable$Positions[1] , AdensTable$Positions[length(AdensTable$Positions)], trunc((length(AdensTable$Positions) - 1) / 10)))
    Graphic <- Graphic + theme(plot.title = element_text(size = 11, hjust = 0.5), panel.background = element_blank())
    
    if(Save){
      ggsave(filename = paste(Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_5Adens__chain", Chain, ".png", sep = ""), plot = Graphic, width = 10, height = 6)
      cat(paste("5Adens plot is stored in ", Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_5Adens__chain", Chain, ".png\n", sep = ""))
    }
  }
  
  return(Graphic)
}

#plot_5Adens_proportions----
#' @title Plot 5Adens proportions.
#'
#' @description Generates plot to analyze the density of contacts around a sphere of 5 Armstrongs, centered in the C-alfa atom from the residue. 
#' The different classes of contacts based on the mutational frustration index are counted in relative terms.
#'
#' @param Pdb Pdb Frustration object.
#' @param Chain Chain of residue. Type: character. Default: NULL.
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#' 
#' @import ggplot2
#' 
#' @export
plot_5Adens_proportions <- function(Pdb, Chain = NULL, Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  if(!is.null(Chain)){
    if(length(Chain) > 1)
      stop("You must enter only one Chain!")
    if(!all(Chain %in% unique(Pdb$atom$chain))) 
      stop(paste("The Chain ", Chain, " doesn't exist! The Chains are: ",
                 paste(unique(Pdb$atom$chain), collapse = " ", sep = ""), sep = ""))
  }
  
  if(!dir.exists(paste(Pdb$JobDir, "/Images", sep = "")))  
    dir.create(paste(Pdb$JobDir, "/Images", sep = ""))
  
  colClasses = c("integer", "character", "integer", 
                 "integer", "integer", "integer",
                 "numeric", "numeric", "numeric")
  AdensTable = read.table(file = paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb_", Pdb$Mode, "_5adens", sep = ""),
                          fill = T, header = T, colClasses = colClasses)
  AdensTable <- as.data.frame(AdensTable)
  colnames(AdensTable) <- c("Positions", "Chains", "Total", "MaximallyFrst", "NeutrallyFrst", "MinimallyFrst")
  
  if(!is.null(Chain)) AdensTable <- AdensTable[AdensTable$Chains == Chain, ]
  
  MinimallyFrst = as.numeric(AdensTable[, 9])
  NeutrallyFrst = as.numeric(AdensTable[, 8])
  MaximallyFrst = as.numeric(AdensTable[, 7])
  
  FrustrationData <- as.data.frame(rbind(MaximallyFrst, NeutrallyFrst, MinimallyFrst))
  FrustrationData$row <- seq_len(nrow(FrustrationData))
  FrustrationData <- melt(FrustrationData, id.vars = "row")
  FrustrationData$variable <- apply(as.matrix(FrustrationData$variable), c(1, 2),
                                    function(x){return(as.numeric(substr(start = 2, stop = nchar(x),x = x)))})
  FrustrationData$row <- as.factor(FrustrationData$row)
  
  Graphic <- ggplot(FrustrationData, aes(x = variable, y = value, fill = row), color = row) + geom_bar(stat = "identity", width = 1) 
  Graphic <- Graphic + scale_x_continuous(breaks = seq(1 , length(AdensTable[, 1]), by = trunc((length(AdensTable[, 1]) - 1) / 10)))
  Graphic <- Graphic + xlab("Positions") + ylab("Density arround 5A sphere (%)")
  if(is.null(Chain))
    Graphic <- Graphic + ggtitle(paste("Density arround 5A sphere (%) in", Pdb$PdbBase))
  else
    Graphic <- Graphic + ggtitle(paste("Density arround 5A sphere (%) in", Pdb$PdbBase, "chain", Chain))
  Graphic <- Graphic + scale_fill_manual(name = "",labels = c("Highly frustrated", "Neutral", "Minimally frustrated"), values = c("red", "gray", "green"))
  Graphic <- Graphic + scale_color_manual(name = "",values = c("red", "gray", "green"))
  Graphic <- Graphic + theme(plot.title = element_text(size = 11, hjust = 0.5), panel.background = element_blank())
  
  if(Save)
    if(is.null(Chain)){
      ggsave(filename = paste(Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_5Adens_around.png", sep = ""), plot = Graphic, width = 10, height = 6)
      cat(paste("5Adens proportion plot is stored in ", Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_5Adens_around.png\n", sep = ""))
    }else{
      ggsave(filename = paste(Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_5Adens_around_chain", Chain, ".png", sep = ""), plot = Graphic, width = 10, height = 6)
      cat(paste("5Adens proportion plot is stored in ", Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_5Adens_around_chain", Chain, ".png\n", sep = ""))
    }
  
  return(Graphic)
}
#plot_contact_map----
#' @title Plot frustration contact map.
#'
#' @description Generates contact map plot to visualize the frustration values assigned to each contact.
#'
#' @param Pdb Pdb Frustration object.
#' @param Chain Chain of residue. Type: character. Default: NULL.
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'   
#' @export
plot_contact_map <- function(Pdb, Chain = NULL, Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  if(!is.null(Chain)){
    if(length(Chain) > 1)
      stop("You must enter only one Chain!")
    if(!all(Chain %in% unique(Pdb$atom$chain))) 
      stop(paste("The Chain ", Chain, " doesn't exist! The Chains are: ",
                 paste(unique(Pdb$atom$chain), collapse = " ", sep = ""), sep = ""))
  }
  if(!dir.exists(paste(Pdb$JobDir, "/Images", sep = ""))) 
    dir.create(paste(Pdb$JobDir, "/Images", sep = ""))
  
  colClasses = c("integer", "character", "integer", 
                 "integer", "integer", "integer",
                 "numeric", "numeric", "numeric")
  AdensTable <- read.table(file = paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb_", Pdb$Mode, "_5adens", sep = ""),
                           fill = T, header = T, colClasses = colClasses)
  
  if(!is.null(Chain)){
    AdensTable <- AdensTable[AdensTable[, 2] == Chain, ]
  }
  
  Positions <- AdensTable[, 1]
  Chains <- AdensTable[, 2]
  MaximallyFrst <- AdensTable[, 4]
  NeutrallyFrst <- AdensTable[, 5]
  MinimallyFrst <- AdensTable[, 6]
  Total <- AdensTable[, 3]
  PositionsTotal <- seq(from = 1, to = length(Positions), by = 1)
  
  colClasses = c("integer", "integer", "character", 
                 "character", "numeric", "numeric",
                 "character", "character", "numeric",
                 "numeric", "numeric", "numeric",
                 "character", "character")
  datos <- read.table(file = paste(Pdb$JobDir, "FrustrationData/", Pdb$PdbBase, ".pdb_", Pdb$Mode, sep = ""),
                      header = T, colClasses = colClasses)
  
  if(!is.null(Chain)){
    datos <- datos[datos$ChainRes1 == Chain, ]
  }
  
  chains <- sort(unique(c(datos$ChainRes1, datos$ChainRes2)))
  positions <- matrix(ncol = 3, nrow = length(chains))
  auxPosVec <- c()
  for(i in seq_along(chains)){
    positions[i, 1:2] <- range(c(datos$Res1[which(datos$ChainRes1 == chains[i])],
                                datos$Res2[which(datos$ChainRes2 == chains[i])]))
    positions[i, 3] <- positions[i,2] - positions[i,1] + 1
    auxPosVec <- c(auxPosVec, positions[i, 1]:positions[i, 2])
  }
  
  datos$pos1 <- NA
  datos$pos2 <- NA
  for(i in seq_along(chains)){
    if(i == 1){bias <- 0} else{bias <- sum(positions[1:(i - 1), 3])}
    idx <- which(datos$ChainRes1 == chains[i])
    datos$pos1[idx] <- datos$Res1[idx] - positions[i, 1] + bias + 1
    idx <- which(datos$ChainRes2 == chains[i])
    datos$pos2[idx] <- datos$Res2[idx] - positions[i, 1] + bias + 1
  }
  
  posNEW <- matrix(ncol = 3, nrow = length(chains))
  for(i in seq_along(chains)){
    posNEW[i, 1:2] <- range(c(datos$pos1[which(datos$ChainRes1 == chains[i])],
                             datos$pos2[which(datos$ChainRes2 == chains[i])]))
    posNEW[i, 3] <- posNEW[i, 2] - posNEW[i, 1] + 1
  }
  
  total.positions <- sum(apply(positions, 1, function(x){x[2] - x[1] + 1}))
  matrz <- matrix(NA, ncol = total.positions, nrow = total.positions)
  for(i in 1:nrow(datos)){
    matrz[datos$pos2[i], datos$pos1[i]] <- datos$FrstIndex[i]
    matrz[datos$pos1[i], datos$pos2[i]] <- datos$FrstIndex[i]
  }
  #PARA HACERLO TRIANGULAR SUPERIOR
  matrz[upper.tri(matrz, diag = T)] <- 0
  
  longData <- melt(matrz)
  longData <- longData[longData$value != 0,]
  Graphic <- ggplot(longData, aes(x = Var2, y = Var1)) + geom_tile(aes(fill = value), na.rm = TRUE) 
  Graphic <- Graphic + scale_fill_gradient2(low = "red", mid = "grey", high = "green", limits = c(-4, 4), 
                                            breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
                                            guide = guide_colourbar(title = paste("Local ", Pdb$Mode, " Frustration Index"),
                                            barwidth = 2, barheight = 20, title.position = "right", title.hjust = 0.5)) 
  Graphic <- Graphic + labs(x = "Residue i", y = "Residue j") + ggtitle(paste("Contact map ", Pdb$PdbBase, sep = ""))
  Graphic <- Graphic + theme_bw() + theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
                                          axis.text.y = element_text(size = 9),
                                          plot.title = element_text(size = 11, hjust = 0.5),
                                          legend.title = element_text(angle = -90))
  
  if(length(chains) > 1){
    breaks <- round(seq(1, total.positions, total.positions / 15))
    labels <- as.character(auxPosVec[breaks])
    Graphic <- Graphic + geom_vline(xintercept = cumsum(posNEW[, 3]), color = "gray", linetype = "longdash") +
                         geom_hline(yintercept = cumsum(posNEW[, 3]), color = "gray", linetype = "longdash") +
                         annotate("text", x = round(apply(posNEW[, 1:2], 1, mean)), y = total.positions + 20, label = unique(chains), parse = TRUE, color = "gray")+
                         annotate("text", x = total.positions + 20, y = round(apply(posNEW[, 1:2], 1, mean)), label = unique(chains), parse = TRUE, color = "gray")+
                         scale_x_continuous(breaks = breaks, labels = labels) + scale_y_continuous(breaks = breaks, labels = labels)
  }
  
  if(Save){
    ggsave(paste(Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_map.png", sep = ""), plot = Graphic, width = 7, height = 6)
    cat(paste("Contact map is stored in ", Pdb$JobDir, "Images/", Pdb$PdbBase, "_", Pdb$Mode, "_map.png\n", sep = ""))
  }
  return(Graphic)
}
#view_frustration_pymol----
#' @title View Frustration Pymol.
#'
#' @description Generates a pymol session to observe the frustration patterns on top of the pdb protein structure.
#'
#' @param Pdb Pdb Frustration object.
#' 
#' @export
view_frustration_pymol <- function(Pdb){
  
  if(!dir.exists(paste(Pdb$JobDir, "/VisualizationScrips", sep = ""))) 
    stop("The representations in Pymol do not exist, re-exercise calculate_frustration() with Visualization = T")
  system(paste("cd ", Pdb$JobDir, "VisualizationScrips/ ; pymol ", Pdb$PdbBase, ".pdb_", Pdb$Mode, ".pml", sep = ""))
}
#plot_dynamic_res_5Adens_proportion----
#' @title Plot 5Adens proportion for specific residue.
#'
#' @description Generates plot to analyze the density of contacts around a sphere of 5 Armstrongs, centered in the C-alfa atom from the residue in Dynamic.
#'
#' @param Dynamic Dynamic frustration object.
#' @param Resno Specific residue. Type: numeric.
#' @param Chain Specific chain. Type: character.
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#' 
#' @import ggplot2
#'   
#' @export
plot_dynamic_res_5Adens_proportion <- function(Dynamic, Resno, Chain, Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  if(!is.null(Dynamic$ResiduesDynamic[[Chain]])){
    if(is.null(Dynamic$ResiduesDynamic[[Chain]][[paste("Res_", Resno, sep = "")]]))
      stop(paste("No analysis to ", Resno, " residue from ", Chain, " chain", sep = ""))
  }
  else  stop(paste("No analysis ", Chain, " chain", sep = ""))
  
  FrustrationResults <- read.table(Dynamic$ResiduesDynamic[[Chain]][[paste("Res_", Resno, sep = "")]], header = T)
  FrustrationResults <- as.data.frame(FrustrationResults)
  
  if(Dynamic$Mode == "configurational" | Dynamic$Mode == "mutational"){
    FrustrationResults <- cbind(seq(1, dim(FrustrationResults)[1]), FrustrationResults)
    colnames(FrustrationResults) <- c("Frame","Res", "ChainRes", "Total", "nHighlyFrst", "nNeutrallyFrst",
                                      "nMinimallyFrst", "relHighlyFrustrated", "relNeutralFrustrated", "relMinimallyFrustrated")
      
    MinimallyFrst = FrustrationResults$relMinimallyFrustrated
    NeutrallyFrst = FrustrationResults$relNeutralFrustrated
    MaximallyFrst = FrustrationResults$relHighlyFrustrated
    
    FrustrationData <- as.data.frame(rbind(MaximallyFrst, NeutrallyFrst, MinimallyFrst))
    FrustrationData$row <- seq_len(nrow(FrustrationData))
    FrustrationData <- melt(FrustrationData, id.vars = "row")
    FrustrationData$variable <- apply(as.matrix(FrustrationData$variable), c(1, 2),
                                      function(x){return(as.numeric(substr(start = 2, stop = nchar(x), x = x)))})
    FrustrationData$row <- as.factor(FrustrationData$row)
    
    Graphic <- ggplot(FrustrationData, aes(x = variable, y = value, fill = row), color = row) + geom_bar(stat = "identity", width = 1) 
    if(length(dim(FrustrationResults)[1]) * 0.05 > 1)
      Graphic <- Graphic + scale_x_continuous(breaks = seq(1 , length(dim(FrustrationResults)[1]), by = ceiling(length(dim(FrustrationResults)[1]) * 0.05)))
    Graphic <- Graphic + xlab("Frame") + ylab("Density arround 5A sphere (%)")
    Graphic <- Graphic + ggtitle(paste("Density arround 5A sphere (%) for residue ", Resno, " from chain ", Chain, " in dynamic", sep = ""))
    Graphic <- Graphic + scale_fill_manual(name = "", labels = c("Highly frustrated", "Neutral", "Minimally frustrated"), values = c("red", "gray", "green"))
    Graphic <- Graphic + scale_color_manual(name = "", values = c("red", "gray", "green"))
    Graphic <- Graphic + theme(plot.title = element_text(size = 11, hjust = 0.5), panel.background = element_blank())
    
    if(Save){
      ggsave(filename = paste(Dynamic$ResultsDir, "Dynamic_plots_res_", Resno, "_", Chain, "/dynamicHist5adens_", Dynamic$Mode, "_Res_", Resno, "_", Chain, ".png", sep = ""),
                 plot = Graphic, width = 10, height = 6)
      cat(paste("Dynamic_res 5Adens proportion plot is stored in ", Dynamic$ResultsDir, "Dynamic_plots_res_", Resno, "_", Chain, "/dynamicHist5adens_", Dynamic$Mode, "_Res_", Resno, "_", Chain, ".png\n", sep = ""))
    }
  }
  
  return(Graphic)
}
#gif_contact_map----
#' @title Contact map gif.
#'
#' @description Generate Gif of contact map of the structure under study along the dynamics.
#'
#' @param Dynamic Dynamic frustration object.
#' @param Show Print Gif on screen(TRUE), but(FALSE) store in corresponding directory. Type: bool. Default: FALSE.
#' 
#' @importFrom magick image_read image_join image_animate image_write
#' @importFrom bio3d basename.pdb
#' 
#' @export
gif_contact_map <- function(Dynamic, Show = FALSE){
  
  if(Show != T & Show != F)
    stop("Show must be a boolean value!")
  if(!requireNamespace("magick", quietly = TRUE)){
    stop("Please install magick package to continue!")
  }
  else library(magick)
  
  cat("----------------------------Getting paths-----------------------------\n")
  paths <- c()
  for(i in seq(1, length(Dynamic$OrderList))){
    paths[i] <- paste(Dynamic$ResultsDir, basename.pdb(Dynamic$OrderList[i]), ".done/Images/", basename.pdb(Dynamic$OrderList[i]), "_", Dynamic$Mode, "_map.png", sep = "")
  }
  cat("----------------------------Loading images-----------------------------\n")
  images <- list()
  for(i in seq(1, length(Dynamic$OrderList))){
    images[[i]] <- image_read(paths[i])
  }
  cat("----------------------------Concatenating images-----------------------------\n")
  images <- image_join(images)
  animation <- image_animate(images, fps = 2, optimize = TRUE)
  
  if(Show) print(animation)
  else{
    image_write(animation, paste(Dynamic$ResultsDir, "contactMap_", Dynamic$Mode, ".gif", sep = ""))
    cat(paste("Contact map gif in ", Dynamic$ResultsDir, "contactMap_", Dynamic$Mode, ".gif\n", sep = ""))
  }
}
#gif_5adens_proportions----
#' @title 5Adens proportions gif.
#'
#' @description Generate Gif of 5Adens proportions plot of the structure under study along the dynamics.
#'
#' @param Dynamic Dynamic frustration object.
#' @param Show Print Gif on screen(TRUE), but(FALSE) store in corresponding directory. Type: bool. Default: FALSE.
#' 
#' @importFrom magick image_read image_join image_animate image_write
#' @importFrom bio3d basename.pdb
#'
#' @export
gif_5adens_proportions <- function(Dynamic, Show = FALSE){
  
  if(Show != T & Show != F)
    stop("Show must be a boolean value!")
  if(!requireNamespace("magick", quietly = TRUE)){
    stop("Please install magick package to continue!")
  }
  else library(magick)
  
  paths <- c()
  for(i in seq(1, length(Dynamic$OrderList))){
    paths[i] <- paste(Dynamic$ResultsDir, basename.pdb(Dynamic$OrderList[i]), ".done/Images/", basename.pdb(Dynamic$OrderList[i]), "_", Dynamic$Mode, "_5Adens_around.png", sep = "")
  }
  images <- list()
  for(i in seq(1, length(Dynamic$OrderList))){
    images[[i]] <- image_read(paths[i])
  }

  images = image_join(images)

  animation <- image_animate(images, fps = 2, optimize = TRUE)
  
  if(Show) print(animation)
  else{
    image_write(image = animation, path = paste(Dynamic$ResultsDir, "5Adens_proportions_", Dynamic$Mode, ".gif", sep = ""))
    cat(paste("5Adens proportion gif in ", Dynamic$ResultsDir, "5Adens_proportions_", Dynamic$Mode, ".gif\n", sep = ""))
  }
}
#frustra_movie----
#' @title Frustration movie in Pymol.
#'
#' @description It generates visualization in pymol of the frames of the indicated dynamics and script for its execution in pymol of generation of the film.
#'
#' @param Dynamic Dynamic frustration object.
#' 
#' @return Stores pymol scripts representations_(Mode).pml and GenerateMovie_(Mode).pml
#' 
#' @importFrom bio3d basename.pdb
#' 
#' @export
frustra_movie <- function(Dynamic){

  ResultDir <- paste(Dynamic$ResultsDir, "FrustraMovie", sep = "")
  if(!dir.exists(ResultDir)) dir.create(ResultDir)

  fileChargePml <- paste(ResultDir, "/representations_", Dynamic$Mode, ".pml", sep = "")
  i <- 1
  for (Pdb in Dynamic$OrderList){
    write(paste("load ", Dynamic$ResultsDir, basename.pdb(Pdb), ".done/VisualizationScrips/", basename.pdb(Pdb), ".pdb_", Dynamic$Mode, ".pml", sep = ""),
          file = fileChargePml, append = TRUE)
    write(paste("group frame", i, ",*", basename.pdb(Pdb), sep = ""), file = fileChargePml, append = TRUE)
    i <- i + 1
  }

  fileMovie <- paste(ResultDir, "/GenerateMovie_", Dynamic$Mode, ".pml", sep = "")

  write(paste("mset 1x", length(Dynamic$OrderList), sep = ""), file = fileMovie, append = TRUE)
  for(Nframe in 1:length(Dynamic$OrderList)){
    write("hide all", file = fileMovie, append = TRUE)
    write(paste("show cgo, frame", Nframe, sep = ""), file = fileMovie, append = TRUE)
    write(paste("show dashes, frame", Nframe, sep = ""), file = fileMovie, append = TRUE)
    write(paste("show cartoon, frame", Nframe, sep = ""), file = fileMovie, append = TRUE)
    write(paste("scene F", Nframe, ",store", sep = ""), file = fileMovie, append = TRUE)
    write(paste("mview store,", Nframe, ",scene=F", Nframe, sep = ""), file = fileMovie, append = TRUE)
  }
  write("show all", file = fileMovie, append = TRUE)
  write("show cgo,all", file = fileMovie, append = TRUE)
  write("show dashes,all", file = fileMovie, append = TRUE)
  write("show cartoon,all", file = fileMovie, append = TRUE)
  
  cat(paste("Script that loads visualizations is stored in ", ResultDir, "/representations_", Dynamic$Mode, ".pml\n", sep = ""))
  cat(paste("Script that generates movie is stored in ", ResultDir, "/GenerateMovie_", Dynamic$Mode, ".pml\n", sep = ""))
  
}
#plot_delta_frus----
#' @title Delta frustration of mutated residue
#'
#' @description Generate a plot of the single residue frustration difference for the mutation of the specific residue given by mutate_res
#'
#' @param Pdb Pdb frustration object
#' @param Resno Specific residue. Type: numeric.
#' @param Chain Specific chain. Type: character.
#' @param Method Method indicates the method to use to perform the mutation (Threading or Modeller). Default: Threading
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#' 
#' @import ggplot2
#' @importFrom bio3d aa321
#'
#' @export
plot_delta_frus <- function(Pdb, Resno, Chain, Method = "threading", Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  if(Pdb$Mode != "singleresidue")
    stop("This graph is available for singleresidue index, run calculate_frustration() with Mode = 'singleresidue' and mutate_res()")
  Method <- tolower(Method)
  if(!any(Method %in% c("threading", "modeller")))
    stop(paste(Method, " it isn't a method of mutation. The methods are: threading or modeller!", sep = ""))
  if(!is.null(Pdb$Mutations[[Method]])){
    if(is.null(Pdb$Mutations[[Method]][[paste("Res_", Resno, "_", Chain, sep = "")]]))
      stop(paste("Not mutated to ", Resno, " residue from ", Chain, " chain", sep = ""))
  }
  else  stop(paste("It was not mutated with the ", Method, " method", sep = ""))

  Mutation <- Pdb$Mutations[[Method]][[paste("Res_", Resno, "_", Chain, sep = "")]]

  colClasses <- c("integer", "character", "character", "numeric")
  DataFrus <- read.table(Mutation$File, header = T, colClasses = colClasses)
  DataFrus <- as.data.frame(DataFrus)
  DataFrus <- cbind(DataFrus, seq(1, length(DataFrus[, 1])), seq(1, length(DataFrus[, 1])))

  colnames(DataFrus) <- c("Res1", "Chain1", "AA1", "FrstIndex", "FrstState", "Color")
  DataFrus[DataFrus$FrstIndex >= 0.58, 5] <- "minimally"
  DataFrus[DataFrus$FrstIndex <= (-1.0), 5] <- "highly"
  DataFrus[DataFrus$FrstIndex > (-1.0) & DataFrus$FrstIndex < 0.58, 5] <- "neutral"
  DataFrus$FrstState <- as.factor(DataFrus$FrstState)

  #Native residues
  Native <- atom.select(Pdb, resno = Mutation$Res, chain = Mutation$Chain, value = TRUE)
  Native <- aa321(unique(Native$atom$resid))

  DataFrus$Color[DataFrus$FrstState == "neutral"] <- "gray"
  DataFrus$Color[DataFrus$FrstState == "highly"] <- "red"
  DataFrus$Color[DataFrus$FrstState == "minimally"] <- "green"
  DataFrus$Color[DataFrus$AA1 == Native] <- "blue"
  #DataFrus$Color[DataFrus$AA1=='G']<-"orange"
  DataFrus$Color <- as.factor(DataFrus$Color)

  #DeltaFrus
  FrstIndexNative <- DataFrus$FrstIndex[DataFrus$AA1 == Native]
 
  DataFrus$FrstIndex <- DataFrus$FrstIndex - rep(FrstIndexNative, nrow(DataFrus))  
  x1 <- DataFrus$Res1[1] - 3
  x2 <- DataFrus$Res1[1] + 3
  y1 <- -4
  y2 <- 4
  if(min(DataFrus$FrstIndex) < (-4))  y1 <- min(DataFrus$FrstIndex)
  if(max(DataFrus$FrstIndex) > (4))  y2 <- max(DataFrus$FrstIndex)

  DataFrus <- rbind(DataFrus[DataFrus$Color != "blue", ], DataFrus[DataFrus$Color == "blue", ])
  
  Graphic <- ggplot(data = DataFrus) + geom_point(aes(x = Res1, y = FrstIndex, color = Color), shape = DataFrus$AA1, size = 3)
  Graphic <- Graphic + xlab("Residue Position") + ylab("Delta Frustration")
  Graphic <- Graphic + scale_x_continuous(breaks = unique(DataFrus$Res1), labels = paste(unique(DataFrus$Res1)))
  Graphic <- Graphic + scale_y_continuous(breaks = seq(y1, y2, 0.5), labels = as.character(seq(y1, y2, 0.5)), limits = c(y1, y2))
  Graphic <- Graphic + scale_color_manual("Frustration", breaks = levels(DataFrus$Color),
                                       labels = c("Native state", "Neutral", "Minimally frustrated", "Highly frustrated"), values = c("blue", "gray", "green", "red"))
  Graphic <- Graphic + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))
  
  if(Save){
    if(!dir.exists(paste(Pdb$JobDir, "MutationsData/Images", sep = "")))  
      dir.create(paste(Pdb$JobDir, "MutationsData/Images", sep = ""))
    ggsave(plot = Graphic, paste(Pdb$JobDir, "MutationsData/Images/Delta_frus_", Resno, "_", Chain, ".png", sep = ""), width = 10, height = 6)
    cat(paste("Delta frus plot is stored in ", Pdb$JobDir, "MutationsData/Images/Delta_frus_", Resno, "_", Chain, ".png\n", sep = ""))
  }
    
  return(Graphic)
}
#plot_mutate_res----
#' @title Frustration of mutated residue
#'
#' @description Plot the frustration for each of the 20 residue variants at a given position in the structure.
#'
#' @param Pdb Pdb frustration object
#' @param Resno Specific residue. Type: numeric.
#' @param Chain Specific chain. Type: character.
#' @param Method Method indicates the method to use to perform the mutation (Threading or Modeller). Default: Threading
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#' 
#' @import ggplot2
#' @importFrom bio3d aa321 aa123 atom.select
#'
#' @export
plot_mutate_res <- function(Pdb, Resno, Chain, Method = "threading", Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  Method <- tolower(Method)
  if(!any(Method %in% c("threading", "modeller")))
    stop(paste(Method, " it isn't a method of mutation. The methods are: threading or modeller!", sep = ""))
  if(!is.null(Pdb$Mutations[[Method]])){
    if(is.null(Pdb$Mutations[[Method]][[paste("Res_", Resno, "_", Chain, sep = "")]]))
      stop(paste("Not mutated to ", Resno, " residue from ", Chain, " chain", sep = ""))
  }
  else  stop(paste("It was not mutated with the ", Method, " method", sep = ""))
  
  colClasses <- c()
  if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational")
    colClasses <- c("integer", "integer", "character", "character",
                    "character", "character", "numeric", "character")
  else if(Pdb$Mode == "singleresidue")
    colClasses <- c("integer", "character", "character", "numeric")
    
  Mutation <- Pdb$Mutations[[Method]][[paste("Res_", Resno, "_", Chain, sep = "")]]
  
  DataFrus <- read.table(Mutation$File, header = T, colClasses = colClasses)
  DataFrus <- as.data.frame(DataFrus)
  DataFrus <- cbind(DataFrus,seq(1, length(DataFrus[, 1])))
  
  if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational"){
    
    colnames(DataFrus) <- c("Res1", "Res2", "Chain1", "Chain2", "AA1", "AA2", "FrstIndex", "FrstState", "Color")
    DataFrus[DataFrus$FrstIndex >= 0.78, 8] <- "minimally"
    DataFrus[DataFrus$FrstIndex <= (-1.0), 8] <- "highly"
    DataFrus[DataFrus$FrstIndex > (-1.0) & DataFrus$FrstIndex < 0.78, 8] <- "neutral"
    DataFrus$FrstState <- as.factor(DataFrus$FrstState)
    #Reacomodo teniendo en cuenta que los datos se almacenan en una matriz triangular superior
    DataFrus[DataFrus$Res2 == Mutation$Res, c(2, 4, 5, 6)] <- c(DataFrus$Res1[DataFrus$Res2 == Mutation$Res], DataFrus$Chain1[DataFrus$Res2 == Mutation$Res],
                                                                DataFrus$AA2[DataFrus$Res2 == Mutation$Res], DataFrus$AA1[DataFrus$Res2 == Mutation$Res])  
    DataFrus$Chain1[DataFrus$Chain1 != Mutation$Chain] <- Mutation$Chain
    DataFrus$Res1[DataFrus$Res1 != Mutation$Res] <- Mutation$Res
  }
  else if(Pdb$Mode == "singleresidue"){
    DataFrus <- cbind(DataFrus, seq(1, length(DataFrus[, 1])))
    colnames(DataFrus) <- c("Res1", "Chain1", "AA1", "FrstIndex", "FrstState", "Color")
    DataFrus[DataFrus$FrstIndex >= 0.58, 5] <- "minimally"
    DataFrus[DataFrus$FrstIndex <= (-1.0), 5] <- "highly"
    DataFrus[DataFrus$FrstIndex > (-1.0) & DataFrus$FrstIndex < 0.58, 5] <- "neutral"
    DataFrus$FrstState <- as.factor(DataFrus$FrstState)
  }
  
  #Se establece residuo nativo
  Native <- atom.select(Pdb, resno = Mutation$Res, chain = Mutation$Chain, value = TRUE)
  Native <- aa321(unique(Native$atom$resid))
  
  
  DataFrus$Color[DataFrus$FrstState == "neutral"] <- "gray"
  DataFrus$Color[DataFrus$FrstState == "highly"] <- "red"
  DataFrus$Color[DataFrus$FrstState == "minimally"] <- "green"
  DataFrus$Color[DataFrus$AA1 == Native] <- "blue"
  #DataFrus$Color[DataFrus$AA1=='G']<-"orange"
  DataFrus$Color <- as.factor(DataFrus$Color)
  
  if(Pdb$Mode == "configurational" | Pdb$Mode == "mutational"){
    
    Contacts <- unique(DataFrus[, c("Res2", "Chain2")])
    Contacts <- Contacts[order(Contacts[, 1]), ]
    Contacts <- cbind(Contacts, 1:length(Contacts[, 1]))
    Contacts <- as.data.frame(Contacts)
    colnames(Contacts) <- c("Res", "Chain", "Index")
    
    Contacts$Res <- as.numeric(Contacts$Res)
    Resid <- atom.select(Pdb, resno = Contacts$Res, elety = "CA", chain = Mutation$Chain, value = T)$atom$resid
    
    DataFrus$Res2 <- Contacts[match(DataFrus$Res2, Contacts$Res), 3]
    
    #DataFrus<-rbind(DataFrus[DataFrus$Color!="orange",],DataFrus[DataFrus$Color=="orange",])
    DataFrus <- rbind(DataFrus[DataFrus$Color != "blue", ], DataFrus[DataFrus$Color == "blue", ])
    
    y1 <- (-4)
    y2 <- 4
    Graphic <- ggplot(data = DataFrus) + geom_point(aes(x = Res2, y = FrstIndex, color = Color), shape = DataFrus$AA1, size = 3)
    Graphic <- Graphic + xlab("Contact residue") + ylab("Frustration Index")
    Graphic <- Graphic + scale_x_continuous(breaks = Contacts$Index, labels = paste(Resid, Contacts$Res, Contacts$Chain))
    Graphic <- Graphic + scale_y_continuous(breaks = seq(y1, y2, 0.5), labels = as.character(seq(y1, y2, 0.5)), limits = c(y1, y2))
    Graphic <- Graphic + geom_hline(yintercept = c(0.78, -1), color = "gray", linetype = "longdash")
    Graphic <- Graphic + ggtitle(paste("Contact Frustration ", Pdb$Mode, " of residue ", aa123(Native), "_", DataFrus$Res1[1], sep = ""))
    Graphic <- Graphic + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90), panel.background = element_blank())
    Graphic <- Graphic + scale_color_manual("", breaks = c("green", "gray", "red", "blue"), labels = c("Minimally frustrated", "Neutral", "Highly frustrated", "Native"),
                                            values = c("green", "gray", "red", "blue"))
    
  }
  else if( Pdb$Mode == "singleresidue"){
    y1 <- (-4)
    y2 <- 4
    Graphic <- ggplot(DataFrus, aes(x = aa123(AA1), y = FrstIndex, color = Color)) + geom_point(size = 2)
    Graphic <- Graphic + xlab("Residue") + ylab("Frustration Index")
    Graphic <- Graphic + scale_y_continuous(breaks = seq(y1, y2, 0.5), labels = as.character(seq(y1, y2, 0.5)), limits = c(y1, y2))
    Graphic <- Graphic + geom_hline(yintercept = c(0.58, -1), color = "gray", linetype = "longdash")
    Graphic <- Graphic + ggtitle(paste("Frustration of the 20 variants in position ", Mutation$Res, " of the structure ", sep = ""))
    Graphic <- Graphic + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90), panel.background = element_blank())
    Graphic <- Graphic + scale_color_manual("", breaks = c("green", "gray", "red", "blue"), labels = c("Minimally frustrated", "Neutral", "Highly frustrated", "Native"),
                                            values = c("green", "gray", "red", "blue"))
    
  }
  
  if(Save){
    if(!dir.exists(paste(Pdb$JobDir, "MutationsData/Images", sep = "")))  
      dir.create(paste(Pdb$JobDir, "MutationsData/Images", sep = ""))
    ggsave(plot = Graphic, paste(Pdb$JobDir, "MutationsData/Images/", Pdb$Mode, "_", DataFrus$Res1[1], "_", Mutation$Method, "_", Mutation$Chain, ".png", sep = ""), width = 10, height = 6)
    cat(paste("Mutate res plot is stored in ", Pdb$JobDir, "MutationsData/Images/", Pdb$Mode, "_", DataFrus$Res1[1], "_", Mutation$Method, "_", Mutation$Chain, ".png\n", sep = ""))
  }
  return(Graphic)  
}
#plot_dynamic_clusters_graph----
#' @title Plot Graph Dynamic Clusters
#'
#' @description Plot the graph obtained through detect_dynamic_clusters(),
#'  where each color represents a cluster
#'
#' @param Dynamic Dynamic Frustration Object
#' 
#' @importFrom igraph plot.igraph set_vertex_attr V
#' 
#' @export
plot_dynamic_clusters_graph <- function(Dynamic){
  
  if(is.null(Dynamic$Clusters[["Graph"]]))
    stop("Cluster detection failed, run detect_dynamic_clusters()")
  
  colorPalette = c(rgb(0.5, 1, 1), rgb(0.8, 0.1, 0.3), rgb(1, 0.7, 0.2), rgb(0.2, 1.0, 0.2), rgb(1, 1, 0),
                   rgb(1.00, 0.42, 0.42), rgb(0.2, 0.6, 0.2), rgb(0.75, 0.75, 1.0), rgb(1.0, 0.3, 0.3), rgb(0.55, 0.25, 0.60),
                   rgb(1.0, 0.2, 0.2), rgb(1.0, 0.55, 0.15), rgb(0.3, 0.3, 1.0), rgb(0.2, 1.0, 0.2), rgb(0.00, 0.75, 0.75),
                   rgb(0.6, 0.2, 0.2), rgb(0.70, 0.30, 0.40), rgb(0.24, 1.0, 0.00), rgb(0.47, 0.36, 0.89), rgb(1.0, 0.82, 0.13))
  Dynamic$Clusters$Graph <- set_vertex_attr(Dynamic$Clusters$Graph, index = V(Dynamic$Clusters$Graph) ,
                                            name = "color", value = colorPalette[Dynamic$Clusters$LeidenClusters$cluster])
  Dynamic$Clusters$Graph <- set_vertex_attr(Dynamic$Clusters$Graph, index = V(Dynamic$Clusters$Graph), 
                         name = "color", value = colorPalette[Dynamic$Clusters$LeidenClusters$cluster])
  
  #library(svglite)
  par(fig=c(0, 0.85, 0, 1))
  #svglite("/home/hacha/Documentos/grafo.svg",width = 20,height = 20)
  plot.igraph(Dynamic$Clusters$Graph, layout = layout_with_kk(Dynamic$Clusters$Graph, dim = 2), vertex.label.dist = 2, edge.width = 1,
              vertex.color = V(Dynamic$Clusters$Graph)$color, vertex.size = 10, edge.color = "black", vertex.label.cex = 1, label.color = "black")
  par(fig=c(0, 1, 0, 1), new=TRUE)
  legend(x = 0.8, y = 0.5, bty = "n", legend = paste(sort(unique(Dynamic$Clusters$LeidenClusters$cluster)), sep = ""), 
         fill = colorPalette[sort(unique(Dynamic$Clusters$LeidenClusters$cluster))], title = "Clusters")
}
#plot_res_dynamics----
#' @title Plot dynamics residues
#'
#' @description Plots the singleresidue frustration index for a residual of a specific chain along the dynamics and fitted model in the detect_dynamic_clusters() function.
#'
#' @param Dynamic Dynamic Frustration Object
#' @param Resno Specific residue. Type: numeric.
#' @param Chain Specific chain. Type: character.
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#' 
#' @import ggplot2
#' @importFrom bio3d basename.pdb atom.select
#'
#' @export
plot_res_dynamics <- function(Dynamic = Dynamic, Resno, Chain, Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  if(!is.null(Dynamic$ResiduesDynamic[[Chain]])){
    if(is.null(Dynamic$ResiduesDynamic[[Chain]][[paste("Res_", Resno, sep = "")]]))
      stop(paste("No analysis to ", Resno, " residue from ", Chain, " chain. Run dynamic_res().", sep = ""))
  }
  else  stop(paste("No analysis ",Chain," chain. Run dynamic_res().", sep = ""))
  
  Pdb <- read.pdb(paste(Dynamic$PdbsDir, Dynamic$OrderList[1], sep = ""), ATOM.only = T, rm.alt = T, rm.insert = T)
  AA <- atom.select(Pdb, resno = Resno, chain = Chain, elety = "CA", value = T)$atom$resid
  rm(Pdb)
  
  FrustrationResults <- read.table(Dynamic$ResiduesDynamic[[Chain]][[paste("Res_", Resno, sep = "")]], header = T)
  FrustrationResults <- as.data.frame(FrustrationResults)
  
  if(Dynamic$Mode == "configurational" | Dynamic$Mode == "mutational"){
    FrustrationResults <- cbind(seq(1, dim(FrustrationResults)[1]), FrustrationResults)
    colnames(FrustrationResults) <- c("Frame", "Res", "ChainRes", "Total", "nHighlyFrst", "nNeutrallyFrst",
                                      "nMinimallyFrst", "relHighlyFrustrated", "relNeutralFrustrated", "relMinimallyFrustrated")
    Maximum = max(c(max(FrustrationResults$nHighlyFrst), max(FrustrationResults$nMinimallyFrst),
                    max(FrustrationResults$nNeutrallyFrst), max(FrustrationResults$Total)))
    
    Graphic <- ggplot() + geom_line(aes(x = FrustrationResults$Frame, y = FrustrationResults$nHighlyFrst, colour = "1"))
    Graphic <- Graphic + geom_line(aes(x = FrustrationResults$Frame, y = FrustrationResults$nNeutrallyFrst, colour = "2"))
    Graphic <- Graphic + geom_line(aes(x = FrustrationResults$Frame, y = FrustrationResults$nMinimallyFrst, colour = "3"))
    Graphic <- Graphic + geom_line(aes(x =  FrustrationResults$Frame, y = FrustrationResults$Total, colour = "4" ))
    Graphic <- Graphic + ylab("Local frustration density (5A sphere)") + xlab("Frame") + ggtitle(paste("Frustration of resid", AA, Resno, Dynamic$Mode))
    Graphic <- Graphic + scale_colour_manual(name = "", labels = c("Highly frustrated", "Neutral", "Minimally frustrated", "Total"), values = c("red", "gray", "green", "black"))
    Graphic <- Graphic + scale_y_continuous(limits = c(0, Maximum), breaks = seq(0.0, Maximum, 5), labels = as.character(seq(0.0, Maximum, 5)))
    Graphic <- Graphic + scale_x_continuous(breaks = seq(1, length(FrustrationResults$Frame), ceiling(length(FrustrationResults$Frame) * 0.05)))
    Graphic <- Graphic + theme_classic() + theme(plot.title = element_text(size = 11, hjust = 0.5), axis.text.x = element_text(angle = 90))
    
    if(Save){
      ggsave(plot = Graphic, paste(Dynamic$ResultsDir, "Dynamic_plots_res_", Resno, "_", Chain, "/dynamic5adens_", Dynamic$Mode, "_Res", Resno, ".png", sep = ""), width = 10, height = 6)
      cat(paste("Dynamic_res 5Adens plot is stored in ", Dynamic$ResultsDir, "Dynamic_plots_res_", Resno, "_", Chain, "/dynamic5adens_", Dynamic$Mode, "_Res", Resno, ".png\n", sep = ""))
    }
  }
  else{
    FrustrationResults <- cbind(seq(1, dim(FrustrationResults)[1]), FrustrationResults[, 8], seq(1, dim(FrustrationResults)[1]))
    FrustrationResults <- as.data.frame(FrustrationResults)
    FrustrationResults[,2] <- as.numeric(FrustrationResults[,2])
    FrustrationResults[,1] <- as.numeric(FrustrationResults[,1])
    FrustrationResults[which(FrustrationResults[, 2] >= 0.58), 3] <- "Minimally frustrated"
    FrustrationResults[which((FrustrationResults[, 2] < 0.58) & (FrustrationResults[, 2] > (-1.0)) ), 3] <- "Neutral"
    FrustrationResults[which(FrustrationResults[, 2] <= (-1.0)), 3] <- "Highly frustrated"
    
    colnames(FrustrationResults) <- c("Frame", "IndexFrst", "Type")
    FrustrationResults$Type <- factor(x = FrustrationResults$Type, levels = c("Minimally frustrated", "Neutral", "Highly frustrated"))
    
    
    Graphic <- ggplot() + geom_point(data = FrustrationResults,aes(x = Frame, y = IndexFrst, colour = Type))
    Graphic <- Graphic + ylab("Index Frustration") + xlab("Frame") + ggtitle(paste("Index Frustration of resid", AA, Resno, "Single residue"))
    cols <- c("Minimally frustrated" = "green", "Neutral" = "gray", "Highly frustrated" = "red")
    Graphic <- Graphic + scale_colour_manual(name = "", values = scales::alpha(cols, 0.6)) + scale_y_reverse(limits = c(4, -4), breaks = seq(4.0, -4.0, -0.5), labels = as.character(seq(4.0, -4.0, -0.5)))
    Graphic <- Graphic + scale_x_continuous(breaks = seq(1, length(FrustrationResults$Frame), ceiling(length(FrustrationResults$Frame) * 0.05)))
    Graphic <- Graphic + theme_classic() + theme(plot.title = element_text(size = 11, hjust = 0.5), axis.text.x = element_text(angle = 90))
    
    if(is.null(Dynamic$Clusters[["Fitted"]])){
      modelo <- loess(FrustrationResults$IndexFrst ~ FrustrationResults$Frame, data = FrustrationResults, span = 0.05, degree = 1, family = "gaussian")
      DataFit <- cbind(seq(1, length(modelo$fitted)), modelo$fitted)
      cat("To make the graph, the Local Polynomial Regression Model was adjusted with Loesspan = 0.05. To use custom parameters, run previously detect_dynamic_clusters().\n\n")
    }else{
      DataFit <- cbind(seq(1, nrow(Dynamic$Clusters$Fitted)), Dynamic$Clusters$Fitted[,Resno])
    }
    Graphic <- Graphic + geom_line(aes(x = DataFit[,1], y = DataFit[,2]), color = "blue")
    
    if(Save){
      ggsave(paste(Dynamic$ResultsDir, "Dynamic_plots_res_", Resno, "_", Chain, "/dynamic_IndexFrustration_", Dynamic$Mode, "_Res", Resno, ".png", sep = ""), plot = Graphic, width = 10, height = 6)
      cat(paste("Dynamic_res plot is stored in ", Dynamic$ResultsDir, "Dynamic_plots_res_", Resno, "_", Chain, "/dynamic_IndexFrustration_", Dynamic$Mode, "_Res", Resno, ".png\n", sep = ""))
    }
  }
  
  return(Graphic)
}
#plot_variable_res_filter----
#' @title Plot variables residues filter
#'
#' @description Plot the filter in frustration dynamic range and the mean of the singleresidue index for each residue along the dynamics.
#'
#' @param Dynamic Dynamic Frustration Object
#' @param Save If it is TRUE it saves the graph, otherwise it does not. Type: bool. Default: FALSE.
#' 
#' @return ggplot2 object.
#' 
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom bio3d basename.pdb
#'
#' @export
plot_variable_res_filter <- function(Dynamic = Dynamic, Save = FALSE){
  
  if(Save != T & Save != F)
    stop("Save must be a boolean value!")
  if(is.null(Dynamic$Clusters[["Graph"]]))
    stop("Cluster detection failed, run detect_dynamic_clusters()")
  
  if(!requireNamespace("ggrepel", quietly = TRUE)){
    stop("Please install ggrepel package to continue!")
  }
  else library(ggrepel)
  
  ini <- read.table(paste(Dynamic$ResultsDir, basename.pdb(Dynamic$OrderList[1]),".done/FrustrationData/",
                          basename.pdb(Dynamic$OrderList[1]), ".pdb_singleresidue", sep = ""), header = T)
  residues <- ini$AA
  rm(ini)
  
  estadistics <- cbind(Dynamic$Clusters$FrstRange, Dynamic$Clusters$Means,
                       1:length(Dynamic$Clusters$Means), 1:length(Dynamic$Clusters$Means))
  estadistics <- as.data.frame(estadistics)
  colnames(estadistics)<- c("FrstRange", "Means", "Res", "Color")
  for(i in 1:length(Dynamic$Clusters$Means)){
    estadistics$Res[i] <- paste(residues[i], "_", i, sep = "")
    if((estadistics$FrstRange[i] > quantile(estadistics$FrstRange, probs = Dynamic$Clusters$MinFrstRange)) &
       (estadistics$Means[i] < -Dynamic$Clusters$FiltMean | estadistics$Means[i] > Dynamic$Clusters$FiltMean))
      estadistics$Color[i] <- "red"
    else estadistics$Color[i] <- "blue"
  }
  
  Graphic <- ggplot() + geom_point(aes(x = estadistics[estadistics$Color == "red",2], y = estadistics[estadistics$Color == "red",1], color = "red" ), size = 3)
  Graphic <- Graphic + geom_point(aes(x = estadistics[estadistics$Color == "blue",2], y = estadistics[estadistics$Color == "blue",1], color = "blue"), size = 3)
  Graphic <- Graphic + geom_hline(yintercept = quantile(estadistics$FrstRange, probs = Dynamic$Clusters$MinFrstRange), color = "red", linetype = "longdash")
  Graphic <- Graphic + geom_vline(xintercept = c(-Dynamic$Clusters$FiltMean,Dynamic$Clusters$FiltMean), color = "red", linetype = "longdash")
  Graphic <- Graphic + scale_color_manual("State", breaks = c("red", "blue"), labels = c("Variable", "Non variable"),
                                values = c("red", "gray35"))
  Graphic <- Graphic + xlab("Mean frustration") + ylab("Frustration dynamic range")
  Graphic <- Graphic + geom_label_repel(aes(x = estadistics[estadistics$Color == "red",2], 
                                            y = estadistics[estadistics$Color == "red",1],
                                            label = estadistics[estadistics$Color == "red",3]),
                                            box.padding   = 0.35,
                                            point.padding = 0.5,
                                            segment.color = 'grey50', size = 3)
  Graphic <- Graphic + theme(axis.line.x = element_line(color = "black"),
                             axis.line.y = element_line(color = "black"), 
                             text = element_text(size = 11), 
                             panel.background = element_blank())
  
  if(Save){
    if(!dir.exists(paste(Dynamic$ResultsDir, "DynamicClusters/", sep = "")))  
      dir.create(paste(Dynamic$ResultsDir, "DynamicClusters/", sep = ""))
    ggsave(plot = Graphic, paste(Dynamic$ResultsDir, "DynamicClusters/Filters_plot.png", sep = ""), width = 10, height = 6)
    cat(paste("Plot cluster detection filtering is stored in ", Dynamic$ResultsDir, "DynamicClusters/Filters_plot.png\n", sep = ""))
  }
  
  return(Graphic)
}
#plot_clusters_pymol----
#' @title Plot clusters in pymol
#'
#' @description Generates a pymol session to observe the residues belonging to the clusters indicate in Clusters.
#' 
#' @param Dynamic Dynamic Frustration Object.
#' @param Clusters Indicates the clusters, for example, c(1, 2, 3), clusters 1, 2 and 3. Default: "all".
#' 
#' @export
plot_clusters_pymol <- function(Dynamic, Clusters = "all"){
  
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
  
  colorPalette = c("aquamarine", "br7", "brightorange", "carbon", "dash",
                   "deepsalmon", "forest", "lightblue", "oxygen", "violetpurple",
                   "tv_red", "tv_orange", "tv_blue", "tv_green", "teal", 
                   "ruby", "raspberry", "calcium", "cerium", "gold")
  file <- paste(tempdir(), "/clustersPymol.pml", sep = "")
  
  write(paste("load ", Dynamic$PdbsDir, Dynamic$OrderList[1], ", structure", sep = ""), file = file)
  write("color grey, structure", file = file, append = T)
  
  clusters <- sort(unique(clusterData$Cluster))
  for (i in 1:length(clusters)) {
    residues <- paste(clusterData$Res[clusterData$Cluster == clusters[i]], collapse="+")
    write(paste("sele cluster", clusters[i], ",resi ", residues, sep = ""),
          file = file, append = T)
    write(paste("color ", colorPalette[i], ", cluster", clusters[i], sep = ""), file = file, append = T)
    write(paste("show sticks, cluster", clusters[i], sep = ""), file = file, append = T)
  }
  write("deselect", file = file, append = T)
  
  system(paste("pymol ", file, sep = ""))
  
  file.remove(file)
}
#save_res_dynamic_clusters----
#' @title Save plots of the dynamics of the residues belonging to a cluster 
#'
#' @description It graphs the frustration dynamics of each residue belonging to the indicated clusters (Clusters) and stores the results in files in PDF format.
#' 
#' @param Dynamic Dynamic Frustration Object
#' @param Clusters Indicates the clusters, for example, c(1, 2, 3), clusters 1, 2 and 3. Default: "all".
#' @param Ncol Number of columns on each page of the file output. Default: 2.
#' @param Nrow Number of rows on each page of the file output. Default: 2.
#' 
#' @importFrom gridExtra marrangeGrob
#' @import ggplot2
#' 
#' @export
save_res_dynamic_clusters <- function(Dynamic, Clusters = "all", Ncol = 2, Nrow = 2){
  
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("Please install gridExtra package to continue!")
  else library("gridExtra")
  
  ResultsDir <- paste(Dynamic$ResultsDir, "Dynamic_clusters_plots/", sep = "")
  if(!dir.exists(ResultsDir))
    dir.create(ResultsDir)
  if(Clusters[1] == "all"){
    Clusters <- unique(get_clusters(Dynamic = Dynamic)$Cluster)
  }
  for(cluster in Clusters)
  {
    plots<-list()
    for(i in 1:nrow(get_clusters(Dynamic, Clusters = cluster)))
    {
      Dynamic <- dynamic_res(Dynamic = Dynamic, Resno = get_clusters(Dynamic, Clusters = cluster)[i,"Res"], Chain = "A", Graphics = F)
      plot <- plot_res_dynamics(Dynamic = Dynamic, Resno = get_clusters(Dynamic, Clusters = cluster)[i,"Res"], Chain = "A")
      plots[[i]] <- plot
    }
    ml <- marrangeGrob(plots, nrow = Nrow, ncol = Ncol, top = "",width = 20, height = 10)
    ggsave(paste(ResultsDir, "cluster_", cluster, ".pdf", sep = ""), ml, width = 20, height = 10)
    cat(paste("The output file cluster_", cluster, ".pdf was stored in ", ResultsDir, "\n\n", sep = "" ))
    #Graphic <- gridExtra::grid.arrange(grobs = plots, ncol = Ncol)
  }
  cat(paste("All output files was stored in ", ResultsDir, "\n\n", sep = "" ))
}