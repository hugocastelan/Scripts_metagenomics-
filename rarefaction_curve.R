

library("vegan")
args <- commandArgs(TRUE)
file_in <- as.character(args[1])
otu<-read.table(file_in, header=TRUE, row.names=1, sep=",")
otumatrix<-as.matrix(otu)
otumatrix_t<-t(otumatrix)
mycolor=c(rep("blue", 4), rep("yellow", 4), rep("red", 4))
Spec <- specnumber(otumatrix_t)
raremax <- min(rowSums(otumatrix_t)
out<-rarecurve(otumatrix_t, step = 200, sample = raremax, col=mycolor, cex = 0.6)
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)
pdf("curve.png",  width = 5*400, height = 5*300, res = 400, pointsize = 8)
dev.off()

