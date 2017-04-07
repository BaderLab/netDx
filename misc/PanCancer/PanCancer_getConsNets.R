# --------------------------------------------------------------
# data dirs for input
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/featSel_170222",rootDir),
	GBM=sprintf("%s/2017_TCGA_GBM/output/featSel_incMut_round2_170223",
		rootDir),
	OV=sprintf("%s/2017_TCGA_OV/output/OV_170227",rootDir),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/featSel_incMutRPPA_round2170223",
		rootDir)
)

source("PanCancer_featSel_writeConsensusNets.R")
for (curSet in names(dirList)){
	writeConsensusNets(dirList[[curSet]],outPfx=sprintf("./%s",curSet))
}
