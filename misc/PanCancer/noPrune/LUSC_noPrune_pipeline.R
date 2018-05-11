#' PanCancer binarized survival: LUSC: Feature selection with one net per
#' datatype
#' 10-fold CV predictor design 

rm(list=ls())

inDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/input/"
outRoot <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/"

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/noPrune_sp0.3_%s",outRoot,dt)

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/LUSC_clinical_core.txt",inDir),
	survival=sprintf("%s/LUSC_binary_survival.txt",inDir)
	)
datFiles <- list(
	rna=sprintf("%s/LUSC_mRNA_core.txt",inDir),
	prot=sprintf("%s/LUSC_RPPA_core.txt",inDir),
 	mir=sprintf("%s/LUSC_miRNA_core.txt",inDir),
	cnv=sprintf("%s/LUSC_CNV_core.txt",inDir)
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

cat("Collecting patient data:\n")
dats <- list() #input data in different slots
cat("\t* Clinical\n")
clinical <- pheno
rownames(clinical) <- clinical[,1];
# =======================
# LUSC-specific variables
clinical$stage <- as.vector(clinical$stage)
clinical$stage[clinical$stage=="Stage IA"| clinical$stage=="Stage IB"] <- "I"
clinical$stage[clinical$stage=="Stage IIA"| clinical$stage=="Stage IIB"| clinical$stage=="Stage II"] <- "II"
clinical$stage[clinical$stage=="Stage IIIA"| clinical$stage=="Stage IIIB"] <- "III"
clinical$stage <- as.factor(clinical$stage)
clinical <- clinical[, -which(colnames(clinical)=="gender")]
clinical <- t(clinical[,c("age","stage")])
clinical[1,] <- as.integer(clinical[1,])
clinical[2,] <- as.integer(as.factor(clinical[2,]))
class(clinical) <- "numeric"
# =======================
dats$clinical <- clinical; rm(clinical)

# create master input net
for (nm in names(datFiles)) {
	cat(sprintf("\t* %s\n",nm))
	tmp <- read.delim(datFiles[[nm]],sep="\t",h=T,as.is=T)
	if (colnames(tmp)[ncol(tmp)]=="X") tmp <- tmp[,-ncol(tmp)]
	rownames(tmp) <- tmp[,1]
	tmp <- t(tmp[,-1])
	class(tmp) <- "numeric"
	if (nm == "rna") tmp <- log(tmp+1)
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
    clinicalArna=c("clinical_cont","rna.profile"),    
    clinicalAprot=c("clinical_cont","prot.profile"),
    clinical="clinical_cont",
	mir="mir.profile",
	rna="rna.profile",
	prot="prot.profile",
	cnv="cnv.profile",
    clinicalAmir=c("clinical_cont","mir.profile"),    
    clinicalAcnv=c("clinical_cont","cnv.profile"),    
    all="all"  
)

cat(sprintf("Clinical variables are: { %s }\n", 
	paste(rownames(dats$clinical),sep=",",collapse=",")))
rm(pheno)
rm(survStr,surv,tmp,nm,outRoot,inDir,dt,k,inFiles,datFiles,pname)

# -----------------------------------------------------------
# run predictor
source("PanCancer_noPrune.R")
runPredictor(mega_combList=combList,rngVals=1:20,netSets=netSets,
	dats=dats,pheno_all=pheno_all,megaDir=megaDir,
	cutoffSet=8,maxEdge=6000,spCutoff=0.3)
