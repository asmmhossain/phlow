#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 2)
{
  stop('Usage: skyridePhylodyn.R <tree> <outname>')
}

suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(adephylo))
suppressPackageStartupMessages(require(pander))
suppressPackageStartupMessages(require(phylodyn))

trName <- args[1]
outName <- args[2]

tr <- read.tree(trName)

dates <- numeric(length(tr$tip.label))

for(i in 1:length(dates)){dates[i]<-as.numeric(lapply(strsplit(tr$tip.label[i],'_'),tail,n=1))}

present <- floor(max(dates))

tr_cond <- BNPR(tr)

gmax <- floor(max(tr_cond$grid))
gmin <- floor(min(tr_cond$grid))

tpos <- seq(gmax,gmin,by=-1)
tlabels <- present - tpos

axlabs <- list(x=tpos,labs=tlabels)

pdf(outName)

oName <- lapply(strsplit(trName,'/'),tail,n=1)

mst = paste('Skyride:',oName,sep=' ')

#plot_BNPR(tr_cond,axlabs=axlabs,main=mst)
plot_BNPR(tr_cond,main=mst)
mtext(paste('present=',present,sep=''),3,line=-1)

invisible(dev.off())