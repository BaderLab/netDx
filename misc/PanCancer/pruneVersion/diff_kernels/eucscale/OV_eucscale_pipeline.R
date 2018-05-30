#' PanCancer binarized survival: KIRC: Feature selection with one net per
# datatype
#' 10-fold CV predictor design
rm(list=ls())


inDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/OV/input"
outRoot <- "/home/shraddhapai/BaderLab/2017_PanCancer/OV/output"

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/eucscale_%s",outRoot,dt)

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/OV_clinical_core.txt",inDir),
	survival=sprintf("%s/OV_binary_survival.txt",inDir)
	)
datFiles <- list(
	rna=sprintf("%s/OV_mRNA_core.txt",inDir),
	prot=sprintf("%s/OV_RPPA_core.txt",inDir),
	mir=sprintf("%s/OV_miRNA_core.txt",inDir),
	dnam=sprintf("%s/OV_methylation_core.txt",inDir),
	cnv=sprintf("%s/OV_CNV_core.txt",inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"

surv <- read.delim(inFiles$survival,sep="\t",h=T,as.is=T)
colnames(surv)[1:2] <- c("ID","STATUS_INT")
survStr <- rep(NA,nrow(surv))
survStr[surv$STATUS_INT<1] <- "SURVIVENO"
survStr[surv$STATUS_INT>0] <- "SURVIVEYES"
surv$STATUS <- survStr
pheno <- merge(x=pheno,y=surv,by="ID")
pheno$X <- NULL
# pheno$gender <- ifelse(pheno$gender=="FEMALE",1, 0)
pheno_nosurv <- pheno[1:4]

cat("Collecting patient data:\n")
dats <- list() #input data in different slots
cat("\t* Clinical\n")
clin <- pheno
rownames(clin) <- clin[,1];
clin <- t(clin[,2,drop=FALSE])
dats$clinical <- clin; rm(clin)

# create master input net
for (nm in names(datFiles)) {
	cat(sprintf("\t* %s\n",nm))
	tmp <- read.delim(datFiles[[nm]],sep="\t",h=T,as.is=T)
	if (colnames(tmp)[ncol(tmp)]=="X") tmp <- tmp[,-ncol(tmp)]
	rownames(tmp) <- tmp[,1]
	tmp <- t(tmp[,-1])
	class(tmp) <- "numeric"
	dats[[nm]] <- tmp
}

cat("\t Ordering column names\n")
# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID), drop = FALSE]})
dats <- lapply(dats, function(x) {
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx, drop = FALSE]
	x
})


# confirm patient order the same for all input nets
pname <- colnames(dats[[1]])
for (k in 2:length(dats)) {
	if (all.equal(colnames(dats[[k]]),pname)!=TRUE) {
		cat(sprintf("Patient order doesn't match for %s\n",
			names(dats)[k]))
		browser()
	}
}

# input nets for each category
netSets <- lapply(dats, function(x) rownames(x))

# compile data
alldat <- do.call("rbind",dats)
pheno_all <- pheno

combList <- list(
    clinical="clinical_cont",
	mir="mir.profile",
	rna="rna.profile",
	prot="prot.profile",
	cnv="cnv.profile",
	dnam="dnam.profile",
    clinicalArna=c("clinical_cont","rna.profile"),
    clinicalAmir=c("clinical_cont","mir.profile"),
    clinicalAprot=c("clinical_cont","prot.profile"),
    clinicalAdnam=c("clinical_cont","dnam.profile"),
    clinicalAcnv=c("clinical_cont","cnv.profile"),
    all="all")

rm(pheno,pheno_nosurv)
rm(survStr,surv,tmp,nm,outRoot,inDir,dt,k,inFiles,datFiles,pname)

# -----------------------------------------------------------
# run predictor
source("PanCancer_eucscale.R")
runPredictor(mega_combList=combList,rngVals=1:20,netSets=netSets,
	dats=dats,pheno_all=pheno_all,megaDir=megaDir,
	cutoffSet=9)




