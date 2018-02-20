rm(list=ls())
require(reshape2)
require(lsa)
require(combinat)
require(cluster)

# ----------------------------------------------------------------------
# utility functions

# run limma to prune top hits
source("runLM.R")
source("silh.R") # silhouette
source("pcByClass.R"); # PCA by class
source("simFuns.R")

tis <- "KIRC"

# GBM
if (tis=="GBM") {
	xprFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/input/GBM_mRNA_core.txt"
	phenoFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/input/GBM_binary_survival.txt"
} else if (tis == "KIRC"){
	xprFile <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/input/KIRC_mRNA_core.txt"
	phenoFile <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/input/KIRC_binary_survival.txt"
} else if (tis == "OV") {
	xprFile <- "/home/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_OV/input/OV_mRNA_core.txt"
	phenoFile <- "/home/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_OV/input/OV_binary_survival.txt"
} else if (tis == "LUSC") {
	phenoFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/input/LUSC_binary_survival.txt"
	xprFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/input/LUSC_RPPA_core.txt"
	#xprFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/input/LUSC_mRNA_core.txt"
}

logFile <- sprintf("~/Desktop/%s.log",tis)
sink(logFile,split=TRUE)
tryCatch({
require(combinat)
xpr <- read.delim(xprFile,sep="\t",h=T,as.is=T)
sname <- xpr[,1]; xpr<- xpr[,-1]
xpr <- t(xpr)
xpr <- xpr[-nrow(xpr),]

#if (tis %in% c("KIRC","LUSC"))  xpr <- log(xpr+1) 
#before_num <- xpr
class(xpr) <- "numeric"
colnames(xpr) <- sname
pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
c1 <- pheno$feature[pheno$is_alive==1]
c2 <- pheno$feature[pheno$is_alive==0]

common <- intersect(colnames(xpr),pheno$feature)
xpr <- xpr[,which(colnames(xpr) %in% common)]
pheno <- pheno[which(pheno$feature %in% common),]

midx <- match(colnames(xpr),pheno$feature)
if (all.equal(pheno$feature[midx],colnames(xpr))!=TRUE) {
	cat("pheno/xpr don't match");browser()
}
pheno <- pheno[midx,]

require(RColorBrewer)
pal <- brewer.pal(n=2,name="Dark2")

# -----------------------------------------------------
# silhouette plot before/after LM-pruning
outFile <- sprintf("~/Desktop/%s_silhouette_before.pdf",tis)
pdf(outFile,width=8,height=8)
xbef <- silh(pheno$is_alive,xpr,title=sprintf("%s: Before lm-pruning",tis))
dev.off()
xpr_before <- xpr


# pca
outFile <- sprintf("~/Desktop/%s_PCAbeforePrune.pdf",tis)
pdf(outFile,width=9,height=4)
PCbyClass_simple(xpr,pheno$is_alive)
dev.off()

outFile <- sprintf("~/Desktop/%s_afterPrune.pdf",tis)
pdf(outFile,width=9,height=4)
tryCatch({
	cat("running lm\n")
	res <- runLM(xpr,pheno$is_alive,topVar=50) 
	res_pre <- res
	ct <- nrow(res)
	thresh_vec <- c(0.01,0.03,0.05,0.07,seq(0.1,0.7,0.1))

	sil_width <- matrix(NA,nrow=length(thresh_vec),ncol=3)
	colnames(sil_width) <- c("thresh","num_vars","avg_sil_width")
	sil_width[,1] <- thresh_vec
	ctr <- 1
	for (thresh in thresh_vec) {
		cat(sprintf("cutoff %1.2f\n", thresh))
		if (sum(res$adj.P.Val < thresh) < 5) {
			sil_width[ctr,2] <- 0
			cat("\t < 5 values left - ignore\n")
		} else {
		res_cur <- subset(res, adj.P.Val < thresh)
		cat(sprintf("\t%i of %i measures left\n",nrow(res_cur), ct,thresh))
		xpr_cur <- xpr[which(rownames(xpr) %in% rownames(res_cur)),]
		par(mfrow=c(1,1))
		x <- silh(pheno$is_alive, xpr_cur,
				  title=sprintf("LM prune:%1.2f",thresh))
		y <- summary(x)
		PCbyClass_simple(xpr_cur,pheno$is_alive)
		cat(sprintf("\tsilh = %1.2f,  %1.2f; avg = %1.2f\n",
				thresh, y$clus.avg.widths[1],y$clus.avg.widths[2],
				y$avg.width))
		sil_width[ctr,2:3] <- c(nrow(xpr_cur),y$avg.width)
		}
		cat("\n----------------\n")
		ctr <- ctr+1
	}
ybef <- summary(xbef)
sil_width <- rbind(c(1,nrow(xpr),ybef$avg.width),sil_width)
par(mfrow=c(1,1))
plot(sil_width[,2],sil_width[,3],xlab="num_vars",ylab="avg silh width",
	 main=sprintf("%s: Silhouette with LM pruning",tis),las=1,
	 ylim=c(-0.02,0.3), #max(sil_width[,3]*1.2)),
	 bty='n',cex.axis=1.3,pch=16,col="red")
abline(h=seq(0.1,1,0.1),lty=3,col='grey50')
abline(h=0,lwd=2,col='black')
text(sil_width[,2], 0.3-(0:(nrow(sil_width)-1) * 0.02),
	sil_width[,1],cex=1,font=3,col='grey20')

},error=function(ex) { print(ex)
},finally={
	dev.off()
})


# pick cutoff with best separation
write.table(sil_width,file=sprintf("%s_silhouette.pdf",tis),
			sep="\t",col=T,row=F,quote=F)
sil_width <- na.omit(sil_width)
idx <- which.max(sil_width[,3])
cat(sprintf("Best cutoff = %1.2f; %i measures; sil width = %1.2f\n",
		sil_width[idx,1], sil_width[idx,2],sil_width[idx,3]))
bestThresh <- sil_width[idx,1]

res <- subset(res, adj.P.Val < bestThresh)

# ----
# plot pairwise sim before/after filt
###pdf(sprintf("~/Desktop/%s_preFilt_sim.pdf",tis))
###cat("----------------\n")
###cat("Before\n")
###cat("----------------\n")
###plotSim(cor(xpr),name="Pearson")
###plotSim(sim.cos(xpr),name="cosine")
c1 <- pheno$feature[pheno$is_alive==1]
c2 <- pheno$feature[pheno$is_alive==0]
plotSim(sim.mi(xpr),name="mutinfo",c1,c2)
plotSim(sim.kern(xpr,"rbf"),name="rbf",c1,c2)
###dev.off()

xpr <- xpr[which(rownames(xpr) %in% rownames(res)),]
pdf(sprintf("~/Desktop/%s_postFilt_sim.pdf",tis))
cat("----------------\n")
cat("After\n")
cat("----------------\n")
plotSim(cor(xpr),name="Pearson",c1,c2)
#plotSim(cos.sim(xpr),name="cosine")
#plotSim(sim.dist(xpr,"euclidean"),name="euclidean")
#plotSim(sim.dist(xpr,"manhattan"),name="manhattan")
#plotSim(sim.dist(xpr,"minkowski"),name="minkowski")
plotSim(sim.dist(xpr,"mi"),name="mi",c1,c2)
for (sig in 0.05) { #seq(0.05,0.7,0.1)) {
	plotSim(sim.kern(xpr,"rbf",sig),name=sprintf("rbf,%1.2f",sig),c1,c2)
}
dev.off()

},error=function(ex){print(ex)
},finally={sink(NULL)
})
