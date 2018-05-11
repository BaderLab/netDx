#' PanCancer binarized survival: GBM: Feature selection with one net per
#' datatype
#' 10-fold CV predictor design 
#' multi cutoff evaluation
#' also pruning RNA before running

rm(list=ls())

rootDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/GBM"
inDir <- sprintf("%s/input",rootDir)
outRoot <- sprintf("%s/output",rootDir)

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/noprune_impute_sp0.3_%s",outRoot,dt)
cat(megaDir, file="test.txt",append=TRUE)

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/GBM_clinical_core.txt",inDir),
	survival=sprintf("%s/GBM_binary_survival.txt",inDir)
	)
datFiles <- list(
	rna=sprintf("%s/GBM_mRNA_core.txt",inDir),
	mir=sprintf("%s/GBM_miRNA_core.txt",inDir),
	dnam=sprintf("%s/GBM_methylation_core.txt",inDir),
	cnv=sprintf("%s/GBM_CNV_core.txt",inDir)
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
pheno_nosurv <- pheno[1:4]

cat("Collecting patient data:\n")
dats <- list() #input data in different slots
cat("\t* Clinical\n")
clinical <- pheno_nosurv
rownames(clinical) <- clinical[,1];
# =======================
# GBM-specific variables
clinical$performance_score[which(clinical$performance_score == "[Not Available]")] <- NA
clinical$performance_score <- strtoi(clinical$performance_score)
clinical$gender <- ifelse(pheno$gender=="FEMALE",1, 0)
# =======================
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
rm(pname)

# input nets for each category
netSets <- lapply(dats, function(x) rownames(x)) 

combList <- list(    
    clinicalAcnv=c("clinical_cont","cnv.profile"),    
    clinical="clinical_cont",    
	mir="mir.profile",
	rna="rna.profile",
	cnv="cnv.profile",
	dnam="dnam.profile",
    clinicalArna=c("clinical_cont","rna.profile"),
    clinicalAmir=c("clinical_cont","mir.profile"),    
    clinicalAdnam=c("clinical_cont","dnam.profile"),    
    all="all"  
)

pheno_all <- pheno

# cleanup
rm(pheno,pheno_nosurv)
rm(rootDir,survStr,surv,tmp,nm,outRoot,inDir,dt,k,inFiles,datFiles)

# -----------------------------------------------------------
# run predictor
source("PanCancer_noPrune_impute.R")
runPredictor(mega_combList=combList,rngVals=1:20,netSets=netSets,
	dats=dats,pheno_all=pheno_all,megaDir=megaDir,
	cutoffSet=9,maxEdge=6000,spCutoff=0.3)


