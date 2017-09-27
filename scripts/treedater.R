#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 5)
{
  stop('Usage: treedater.R <tree> <alnName> <outTree> <plotName> <resultName>')
}


library('ape')
library('treedater')
library('seqinr')
library('pander')

# read the phylogenetic tree
tre <- read.tree(args[1])


# read the alignment and get length
seqs <- read.FASTA(args[2])
seqLen <- length(seqs[[1]])



# extract the dates from tip.labels
sts <- c()

nTips <- length(tre$tip.label)

for(i in 1:nTips){
  temp <- strsplit(tre$tip.label[i],'_')
  temp <- as.numeric(tail(temp[[1]],n=1))
  sts <- c(sts,temp)
}

# assign names to the dates
names(sts) <- tre$tip.label

# run basic treeDater function to get time-stamped tree
dtr <- dater(tre,sts,seqLen)

# write time-stamped tree
write.tree(dtr,args[3])

# run parametric bootstrap to get substitution rate and tMRCA
pb <- parboot.treedater(dtr)

# plot the estimated number of lineages through time
png(args[4])
plot.parboot.ltt( pb ) 
dev.off()

results <- data.frame(TMRCA=dtr$timeOfMRCA,SubstRate=dtr$meanRate,Clock=dtr$clock)
write.table(results,args[5],col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
pander(results)
