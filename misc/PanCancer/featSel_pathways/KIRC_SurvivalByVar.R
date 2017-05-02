#' KIRC see how age, stage, grade track with survival
rm(list=ls())

rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
inDir <- sprintf("%s/2017_TCGA_KIRC/input",rootDir)

inFiles <- list(
	clinical=sprintf("%s/KIRC_clinical_core.txt",inDir),
	survival=sprintf("%s/KIRC_binary_survival.txt",inDir)
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

yes_idx <- which(pheno$STATUS == "SURVIVEYES")
no_idx <- which(pheno$STATUS == "SURVIVENO")

pheno$orig_grade <- pheno$grade
pheno$orig_stage <- pheno$stage
pheno$grade <- as.numeric(factor(pheno$grade))
pheno$stage <- as.numeric(factor(pheno$stage))

require(ggplot2)
cat("Age\n")
p1 <- ggplot(pheno,aes(STATUS, age)) + geom_boxplot() 
wmw <- wilcox.test(pheno$age[yes_idx],pheno$age[no_idx])
p1 <- p1  + ggtitle(sprintf("age, p < %1.2e", wmw$p.value))
cat(sprintf("WMW: p < %1.2e\n", wmw$p.value))

cat("stage\n")
p2 <- ggplot(pheno,aes(STATUS, stage)) + geom_boxplot() 
wmw <- wilcox.test(pheno[yes_idx,"stage"],pheno$stage[no_idx])
p2 <- p2  + ggtitle(sprintf("stage, p < %1.2e", wmw$p.value))
cat(sprintf("WMW: p < %1.2e\n", wmw$p.value))

cat("grade\n")
p3 <- ggplot(pheno,aes(STATUS, grade)) + geom_boxplot() 
wmw <- wilcox.test(pheno[yes_idx,"grade"],pheno$grade[no_idx])
p3 <- p3  + ggtitle(sprintf("grade, p < %1.2e", wmw$p.value))
cat(sprintf("WMW: p < %1.2e\n", wmw$p.value))

source("multiplot.R")
plotList <- list(p1,p2,p3)
multiplot(plotlist=plotList,cols=3)

