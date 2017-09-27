#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 2)
{
  stop('Usage: rtt2.R <tree> <outname>')
}

suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(adephylo))
suppressPackageStartupMessages(require(pander))

trName <- args[1]
outName <- args[2]

tr <- read.tree(trName)

dates <- numeric(length(tr$tip.label))

#for(i in 1:length(dates)){dates[i]<-as.numeric(strsplit(tr$tip.label[i],'_')[[1]][2])}
for(i in 1:length(dates)){dates[i]<-as.numeric(lapply(strsplit(tr$tip.label[i],'_'),tail,n=1))}

tr.rtt <- rtt(tr,dates)
write.tree(tr.rtt,paste(outName,".tre",sep=""))

# Now summarise results

rd <- distRoot(tr.rtt)
td <- dates[match(tr.rtt$tip.label,tr$tip.label)]
rtt.lm <- lm(rd~td)
root.time <- unname(-as.double(coef(rtt.lm)[1])/coef(rtt.lm)[2])
results <- data.frame(TMRCA=root.time,SubstRate=as.double(coef(rtt.lm)[2]))
write.table(results,paste(outName,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
pander(results)

pName <- paste(outName,'.rttRegression.pdf',sep='')
pdf(pName)
plot(rd~td,xlab="Time",ylab="Root to tip distance",ylim=c(0,max(rd)),xlim=c(root.time,max(td)),pch=16,col="red")
abline(rtt.lm)
abline(h=0,lty=2)
invisible(dev.off())