#' ependymoma data - prepare and merge
#' merge German/Toronto cohort, correct for site with combat and save.
rm(list=ls())
require(GEOquery)

rootDir <- "/home/shraddhapai/BaderLab/2017_Ependymoma"
inDir <- sprintf("%s/input",rootDir)
outDir <- sprintf("%s/input/netDx_prepared",rootDir)

if (!file.exists(outDir)) dir.create(outDir)

# -------------------
# german samples
xpr <- read.delim(sprintf("%s/original_data/Germany-comparison-with-spinals/GER-ST-PFPURE-SPMIX-SEP16.gct",inDir),skip=2,h=T,as.is=T)
gplFile <- sprintf("%s/GPL6480-9577.txt",inDir)
gpl <- read.delim(gplFile,sep="\t",h=T,as.is=T,comment="#")
midx <- match(xpr$Description,gpl$ID)
if (all.equal(gpl$ID[midx],xpr$Description)!=TRUE) {
	cat("don't match"); browser()
}
xpr_genes <- gpl$GENE_SYMBOL[midx]
xpr <- xpr[,-(1:2)]
xpr_agg <- aggregate(xpr, by=list(GENE_SYMBOL=xpr_genes),FUN=mean,na.rm=T)
rownames(xpr_agg) <- xpr_agg[,1]
xpr_agg <- xpr_agg[,-1]
xpr <- xpr_agg


sampType <- scan(sprintf("%s/original_data/Germany-comparison-with-spinals/GER-ST-PFPURE-SPMIX-SEP16.cls",inDir),skip=2)
sampType <- as.integer(sampType)

# from Ruth
#  st = 0, PFPURE = 1 and PFMIX = 2
pheno <- data.frame(ID=colnames(xpr),INT_STATUS=sampType)
pheno$ID <- as.character(pheno$ID)
st <- c("ST","PFPURE","PFMIX")
pheno$STATUS <- st[sampType+1]

g_pheno <- pheno
g_xpr <- xpr
rm(pheno,xpr)

# ----------------------------------
# toronto samples
xpr <- read.delim(sprintf("%s/original_data/Toronto-comparison2-without-spinals/TOR-ST-PFPURE-PFMIX-SEP16.gct",inDir),skip=2,h=T,as.is=T)
rownames(xpr) <- xpr[,1]
xpr <- xpr[,-(1:2)]
sampType <- scan(sprintf("%s/original_data/Toronto-comparison2-without-spinals/TOR-ST-PFPURE-PFMIX-SEP16.cls",inDir),skip=2)
sampType <- as.integer(sampType)

# from Ruth
#  st = 0, PFPURE = 1 and PFMIX = 2
pheno <- data.frame(ID=colnames(xpr),INT_STATUS=sampType)
pheno$ID <- as.character(pheno$ID)
st <- c("ST","PFPURE","PFMIX")
pheno$STATUS <- st[sampType+1]

t_pheno <- pheno; 
t_xpr <- xpr
t_pheno$COHORT <- "Toronto"
t_xpr <- log(t_xpr);
rm(pheno,xpr)

g_pheno$COHORT <- "Germany"

common <- intersect(rownames(t_xpr),rownames(g_xpr))
t_xpr <- t_xpr[which(rownames(t_xpr)%in% common),]
g_xpr <- g_xpr[which(rownames(g_xpr)%in% common),]
midx <- match(rownames(t_xpr),rownames(g_xpr))
if (all.equal(rownames(g_xpr)[midx],rownames(t_xpr))!=TRUE) {
	cat("t and g don't match"); browser()
}
g_xpr <- g_xpr[midx,]

xpr <- cbind(t_xpr,g_xpr)
pheno <- rbind(t_pheno,g_pheno)

pheno$COHORT <- factor(pheno$COHORT)
pheno$STATUS <- factor(pheno$STATUS)

# apply combat to correct for cohort effect
require(sva)
mod0 <- model.matrix(~1,data=pheno)
combat_edata <- sva::ComBat(dat=xpr,batch=pheno$COHORT,mod=mod0,
	par.prior=TRUE,prior.plots=FALSE)

pdf(sprintf("%s/hclust_batchCorrect.pdf",outDir))
plotDendro_clr(xpr,pheno,groupPal=list(STATUS="Dark2",COHORT="Spectral"),topVar=round(0.5*nrow(xpr)),ttl="Ependymoma (before)")
plotDendro_clr(combat_edata,pheno,groupPal=list(STATUS="Dark2",COHORT="Spectral"),topVar=round(0.5*nrow(xpr)),ttl="Ependymoma (after)")
dev.off()

dt <- format(Sys.Date(),"%y%m%d")
save(pheno,xpr,file=sprintf("%s/Ependymoma_cohortMerged_%s.Rdata",outDir,dt))



