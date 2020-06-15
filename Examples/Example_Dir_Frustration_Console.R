library(usethis)
library(devtools)
library(bio3d)
library(ggplot2)
library(reshape2)
library(magick)

args = commandArgs(trailingOnly=TRUE)

install_github("proteinphysiologylab/frustratometeR")
library(frustratometeR)

PdbDir<-args[1]
Modes<-args[2]
ResultsDir<-args[3]

dir_frustration(PdbDir=PdbDir, Modes=Modes, ResultsDir = ResultsDir)
