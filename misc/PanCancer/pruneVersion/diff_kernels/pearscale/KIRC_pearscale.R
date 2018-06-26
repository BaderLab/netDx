#' PanCancer binarized survival: KIRC: Feature selection with one net per
# datatype
#' 10-fold CV predictor design 
rm(list=ls())

inDir <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/input"
outRoot <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/output"

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/pearscale_%s",outRoot,dt)

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/KIRC_clinical_core.txt",inDir),
	survival=sprintf("%s/KIRC_binary_survival.txt",inDir)
	)
datFiles <- list(
	rna=sprintf("%s/KIRC_mRNA_core.txt",inDir),
	prot=sprintf("%s/KIRC_RPPA_core.txt",inDir),
	mir=sprintf("%s/KIRC_miRNA_core.txt",inDir),
	dnam=sprintf("%s/KIRC_methylation_core.txt",inDir),
	cnv=sprintf("%s/KIRC_CNV_core.txt",inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"

#======transform clinical data=========
pheno$grade <- as.vector(pheno$grade)
pheno$grade[pheno$grade=="G1"] <- "G2"
pheno$grade[pheno$grade=="GX"] <- "G2"
pheno$grade <- as.factor(pheno$grade)
pheno <- pheno[, -which(colnames(pheno)=="gender")]
#======================================

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
clinical <- pheno_nosurv
rownames(clinical) <- clinical[,1];
clinical$grade <- as.numeric(factor(clinical$grade))
clinical$stage <- as.numeric(factor(clinical$stage))
clinical$ID <- NULL
clinical <- t(clinical)
dats$clinical <- clinical; rm(clinical)

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
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID)]})
dats <- lapply(dats, function(x) { 
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx]
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
source("PanCancer_pearscale.R")
runPredictor(mega_combList=combList,rngVals=1:20,netSets=netSets,
	dats=dats,pheno_all=pheno_all,megaDir=megaDir,
	cutoffSet=9)


