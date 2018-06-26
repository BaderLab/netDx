#' plot similarity for intra- and inter-class for netDx input net
rm(list=ls())

inDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/input"
inFiles <- list(
	clinical=sprintf("%s/LUSC_clinical_core.txt",inDir),
	survival=sprintf("%s/LUSC_binary_survival.txt",inDir)
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

netDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/pruneRBFsigma25_180201/networks"
netList <- dir(netDir,pattern="_cont.txt")
source("simFuns.R")
cur <- read.delim(sprintf("%s/%s",netDir,netList[1]),sep="\t",h=F,as.is=T)


