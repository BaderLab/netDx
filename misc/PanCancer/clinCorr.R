#' troubleshoot clinical measure
rm(list=ls())

source("simFuns.R")
# -------------------------
datFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/input/LUSC_binary_survival.txt"
clinFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/input/LUSC_clinical_core.txt"

surv <- read.delim(datFile,sep="\t",h=T,as.is=T)
pheno <- read.delim(clinFile,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"

colnames(surv)[1:2] <- c("ID","STATUS_INT")
survStr <- rep(NA,nrow(surv))
survStr[surv$STATUS_INT<1] <- "SURVIVENO"
survStr[surv$STATUS_INT>0] <- "SURVIVEYES"
surv$STATUS <- survStr
pheno <- merge(x=pheno,y=surv,by="ID")

# =======================
# LUSC-specific variables
clinical <- pheno
clinical$stage <- as.vector(clinical$stage)
clinical$stage[clinical$stage=="Stage IA"| clinical$stage=="Stage IB"] <- "I"
clinical$stage[clinical$stage=="Stage IIA"| clinical$stage=="Stage IIB"| clinical$stage=="Stage II"] <- "II"
clinical$stage[clinical$stage=="Stage IIIA"| clinical$stage=="Stage IIIB"] <- "III"
#clinical$stage <- factor(clinical$stage,
#	levels=c("Stage IA","Stage IB",
#			 "Stage IIA","Stage IIB",
#			 "Stage IIIA","Stage IIIB"))
clinical <- clinical[, -which(colnames(clinical)=="gender")]
rownames(clinical) <- clinical$ID
clinical <- t(clinical[,c("age","stage")])
clinical[1,] <- as.integer(clinical[1,])
clinical[2,] <- as.integer(as.factor(clinical[2,]))
class(clinical) <- "numeric"
# =======================

ztrans <- function(x) (x-mean(x,na.rm=TRUE))/(sd(x,na.rm=TRUE))
clinical[1,] <- ztrans(clinical[1,])

yes <- pheno$ID[which(pheno$STATUS=="SURVIVEYES")]
no <- pheno$ID[which(pheno$STATUS=="SURVIVENO")]

pdf("~/Desktop/clinical_noDiff.pdf")
plotSim(sim.normDiff2(clinical),name="normDiff2",c1=yes,c2=no)
plotSim(sim.dist(clinical,"euclidean"),name="euclidean",c1=yes,c2=no)
plotSim(sim.dist(clinical,"manhattan"),name="manhattan",c1=yes,c2=no)
plotSim(sim.dist(clinical,"minkowski"),name="minkowski",c1=yes,c2=no)
plotSim(sim.cos(clinical),name="cosine",c1=yes,c2=no)
plotSim(sim.mi(clinical),name="mutinfo",c1=yes,c2=no)
for (sig in seq(0.05,0.8,0.05)) {
plotSim(sim.kern(clinical,"rbf",sigma=sig),name="rbf",c1=yes,c2=no)
}
dev.off()
