
ss <- function(TP,TN,FP,FN,bet=1){
	if (is.na(TP)) TP=0
	if (is.na(TN)) TN=0
	if (is.na(FP)) FP=0
	if (is.na(FN)) FN=0
	sen=TP/(TP+FN)
	spe=TN/(TN+FP)
	pre=TP/(TP+FP)
	f = (1+bet^2)*pre*sen/(bet^2*pre+sen)
	return(c(sen,spe,pre,f))
}

log_up <- function(x){
	ids = which(x<=0)
	x[ids] = 0.0000001
	return(log10(x))
}
library(ggplot2)

args = (commandArgs(trailingOnly = T))
print (args)
for (i in 1:length(args)){
	eval(parse(text=args[[i]]))
}

#dat <- read.table(file=finput,header=F,row.names=1,comment.char="")
dat <- read.table(file=finput,header=T,comment.char="")
rownames(dat) = paste(dat[,"chr"],dat[,"coor"],dat[,"strand"],sep="|")
dat = cbind(dat,0)
colnames(dat) <- c("chr","coordinate","strand","ifSNP","gene","reference_base","upstream_1base","downstream_1base","major_base","major_count","tot_count","major_ratio","if_MI","MI","pvalue_MI","estimated_allelic_ratio","ifNEG","if_GLM","pvalue_GLM","RNAE_t","A","C","G","T","ifRNAE")
ids = which(dat[,"major_base"]=="N")
if (length(ids)>0) dat = dat[-ids,] 
dat[rownames(subset(dat,pvalue_MI>0 & pvalue_MI<=0.05)),"ifRNAE"]= 1
dat.pos <- subset(dat, pvalue_MI>0 & pvalue_MI <= 0.01 & RNAE_t == "AG")
dat.neg <- subset(dat,ifNEG == 1)
#dat.neg <- subset(dat,ifNEG == 1 & RNAE_t=="AG")
dat.unknown = subset(dat,ifRNAE==0 & ifSNP ==0 & RNAE_t=="AG")

dat.pos <- cbind(dat.pos[,c("upstream_1base","downstream_1base")],abs(dat.pos[,"major_ratio"]-dat.pos[,"estimated_allelic_ratio"]),1)
dat.neg <- cbind(dat.neg[,c("upstream_1base","downstream_1base")],abs(dat.neg[,"major_ratio"]-dat.neg[,"estimated_allelic_ratio"]),0)
dat.unknown <- cbind(dat.unknown[,c("upstream_1base","downstream_1base")],abs(dat.unknown[,"major_ratio"]-dat.unknown[,"estimated_allelic_ratio"]))

colnames(dat.pos) = colnames(dat.neg) = c("upstream_1base","downstream_1base","dif","type")
colnames(dat.unknown) = c("upstream_1base","downstream_1base","dif")


dat.train = rbind(dat.pos[1:round(nrow(dat.pos)/2),],dat.neg[1:round(nrow(dat.neg)/2),])
dat.test = rbind(dat.pos[(round(nrow(dat.pos)/2)+1):nrow(dat.pos),],dat.neg[(round(nrow(dat.neg)/2)+1):nrow(dat.neg),])

glm.model <- glm(type ~ .,family=binomial,data = dat.train)
pre.test <- predict(glm.model, dat.test)

cuts=seq(min(pre.test)-0.05,max(pre.test),0.1)
cuts.f = matrix(0,nrow=0,ncol=5,dimnames=list(c(),c("cutoff","sen","spe","pre","f")))
for (i in cuts){
	sub.pos = which(pre.test >= i)
	sub.neg = which(pre.test < i)
	tb.pos = table(dat.test[names(sub.pos),"type"])
	tb.neg = table(dat.test[names(sub.neg),"type"])
	TP = as.numeric(tb.pos["1"])
	FP = as.numeric(tb.pos["0"])
	TN = as.numeric(tb.neg["0"])
	FN = as.numeric(tb.neg["1"])
	f=ss(TP,TN,FP,FN,bet=0.3)
	cuts.f=rbind(cuts.f,c(i,f))
}
cut.off = cuts.f[which.max(cuts.f[,"f"]),"cutoff"]

pre.unknown <- predict(glm.model, dat.unknown)
#dat[names(which(pre.unknown>=cut.off)),"ifRNAE_glm"]= 1
dat[names(which(pre.unknown>=cut.off)),"ifRNAE"]= 2 
dat = dat[,-c(13,18,19)]
write.table(file=fout,dat,quote=F,sep="\t")

