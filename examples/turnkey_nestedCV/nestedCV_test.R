#' abstracted way of running nested CV with netDx
#' This function sets up the data for the particular application and 
#' calls runPredicted_nested() which is a generic function to run the 
#' predictor for binary classification

rm(list=ls())
require(netDx)
require(netDx.examples)

# VM1
inDir <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/input"
outRoot <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/output"

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/nestCV_%s",outRoot,dt)

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/KIRC_clinical_core.txt",inDir),
	survival=sprintf("%s/KIRC_binary_survival.txt",inDir),
	rna=sprintf("%s/KIRC_mRNA_core.txt",inDir),
	mut=sprintf("%s/from_firehose/KIRC_core_somatic_mutations.txt",
		inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"

#======transform clinical data=========
# KIRC-specific
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

dats <- list() #input data in different slots

# clinical
cat("\t* Clinical\n")
clinical <- pheno_nosurv
rownames(clinical) <- clinical[,1];
clinical$grade <- as.numeric(factor(clinical$grade))
clinical$stage <- as.numeric(factor(clinical$stage))
clinical$ID <- NULL
clinical <- t(clinical)
dats$clinical <- clinical; rm(clinical)

# RNA
cat("\t* RNA\n")
rna <- read.delim(inFiles$rna,sep="\t",h=T,as.is=T)
rna <- t(rna)
colnames(rna) <- rna[1,]; rna <- rna[-1,]; 
rna <- rna[-nrow(rna),]
class(rna) <- "numeric"
rownames(rna) <- sub("mRNA_","",rownames(rna))
rownames(rna) <- sub("\\..*","",rownames(rna))
dats$rna <- rna; rm(rna)

# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID)]})
dats <- lapply(dats, function(x) { 
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx]
	x
})



# Define net groupings
clinList <- list(age="age",grade="grade",stage="stage")
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
   path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
groupList <- list(rna=pathwayList,clinical=clinList)
rm(pheno_nosurv)

save(pheno,dats, groupList, file="nestedCV_input.rda")


source("runPredictor_nested.R")
source("makeNets.R")
source("KIRC_makeNets.R")
source("normDiff.R")
runPredictor_nested(pheno,dataList=dats, groupList=groupList,
		makeNetFunc=KIRC_makeNets,outDir=megaDir,numCores=8L,
		nFoldCV=10L, CVcutoff=9L,numSplits=3L)
