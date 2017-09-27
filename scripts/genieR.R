#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 2)
{
  stop('Usage: rtt2.R <tree> <outname>')
}

suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(pander))
suppressPackageStartupMessages(require(genieR))
#suppressPackageStartupMessages(require(adephylo))


trName <- args[1]
outName <- args[2]

tre <- read.tree(trName)

const <- Geniefit(tre,Model='const',start=c(100),upper=Inf,lower=0)

expo <- Geniefit(tre,Model='expo',start=c(100,.1),upper=Inf,lower=0)

log <- Geniefit(tre,Model='log',start=c(100,.1,.1),upper=Inf,lower=0)

popSize <- c(const$parr,expo$parr[1],log$parr[1])
growthRate <- c(NA,expo$parr[2],log$parr[2])
aic <- c(const$AIC,expo$AIC,log$AIC)
models <- c('Constant','Exponential','Logistic')

genie <- data.frame(models,popSize,growthRate,aic)
colnames(genie) <- c('Models','PopulationSize','GrowthRate','AIC')
pander(genie)

write.table(genie,outName,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)