
rootDir <- "/mnt/data2/BaderLab"
consNetDir <- sprintf("%s/PanCancer_common/pathwaysOnly_170502",rootDir)
netScoreFile <- list(
	SURVIVENO=sprintf("%s/pathways_thresh10_pctPass0.70_SURVIVENO_netScores.txt",
		consNetDir),
	SURVIVEYES=sprintf("%s/pathways_thresh10_pctPass0.70_SURVIVEYES_netScores.txt",
		consNetDir)
)

fsNets <- c() # nets feature selected in the original analysis and that
			  # we should exclude now.
for (nm in names(netScoreFile)) {
	netScores	<- read.delim(netScoreFile[[nm]],sep="\t",h=T,as.is=T)
	netNames 	<- netScores[,1]
	netScores <- netScores[,-1]

	wasHighScoring <- netScores >=7 
	wasHighScoring <- rowSums(wasHighScoring,na.rm=TRUE)
	idx <- which(wasHighScoring>=70)
	cat(sprintf("%i high-scoring nets removed\n", length(idx)))
	cat("-----\n"); print(netNames[idx]); cat("-----\n")
	fsNets <- c(fsNets, netNames[idx])
}
fsNets <- sub(".profile","",fsNets)

# now make sure these didn't get included in the list of networks.
# command to generate these.
# spai@debian:/mnt/data2/BaderLab/PanCancer_KIRC/output$ ls randomPathNotFS_randomNets_*170711/SURVIVE*/networks/*.profile | xargs -L1 sh -c 'basename $1 .profile' dummy | sort -k1,1 | uniq >  netNames.txt
netNames <- "/mnt/data2/BaderLab/PanCancer_KIRC/output/netNames.txt"
nn <- read.delim(netNames,sep="\t",h=F,as.is=T)[,1]

cat("number of fsNets overlapping those chosen in noFS test\n")
print(length(intersect(nn, fsNets)))
