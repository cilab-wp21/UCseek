library("multtest")
library("DNAcopy")
library("aCGH")


lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

cbs.segment01 <- function(indir, outdir, bad.bins, varbin.gc, varbin.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width) {
	gc <- read.table(varbin.gc, header=T)
	bad <- read.table(bad.bins, header=F)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F) 
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	thisRatio$ratio <- a / mean(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)

	a <- quantile(gc$bin.length, 0.985)
	#thisRatioNobad <- thisRatio[which(bad[, 1] == 0), ] ##ûȥ
	thisRatioNobad <- thisRatio
	set.seed(25) 
	CNA.object <- CNA(log(thisRatio$lowratio, base=2), thisRatio$chrom, thisRatio$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatio$seg.mean.LOWESS <- m[, 1]

	n<- length(thisRatio$lowratio)
	k <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)
	 vecPred=log(thisRatio$seg.mean.LOWESS, base=2)
         vecObs=log(thisRatio$lowratio, base=2)
         k[,1]<- mergeLevels(vecObs,vecPred)$vecMerged
thisRatio$merge.level <- k[, 1]


set.seed(25) 
	CNA.object <- CNA(log(thisRatioNobad$lowratio, base=2), thisRatioNobad$chrom, thisRatioNobad$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatioNobad), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatioNobad$seg.mean.LOWESS <- m[, 1]
n<- length(thisRatioNobad$lowratio)
	k <- matrix(data=0, nrow=nrow(thisRatioNobad), ncol=1)
	 vecPred=log(thisRatioNobad$seg.mean.LOWESS, base=2)
         vecObs=log(thisRatioNobad$lowratio, base=2)
         k[,1]<- mergeLevels(vecObs,vecPred)$vecMerged
thisRatioNobad$merge.level <- k[, 1]


	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.short.txt", sep=""), quote=F, row.names=F) 
	write.table(thisRatioNobad, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.data.txt", sep=""), quote=F, row.names=F) 

}


cbs.segment01(indir="/indir",outdir="outdir",bad.bins="hg19.50k.k50.bad.bins.txt", varbin.gc="hg19.varbin.gc.content.50k.bowtie.k50.txt",varbin.data="sample_name", sample.name="smple_name_cbs", alt.sample.name="", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5)
