suppressMessages(require(netDx))

# -----------------------------------------------------------
# prepare data
# -----------------------------------------------------------
cat("Preparing data\n")
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
curatedTCGAData(diseaseCode="BRCA", assays="*",dru.run=TRUE)

brca <- curatedTCGAData("BRCA",c("mRNAArray","Mutation"),FALSE)

# get subtype info
pID <- colData(brca)$patientID
pam50 <- colData(brca)$PAM50.mRNA
pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
pam50[which(pam50 %in% "Luminal A")] <- "LumA"
staget <- colData(brca)$pathology_T_stage
st2 <- rep(NA,length(staget))
st2[which(staget %in% c("t1","t1a","t1b","t1c"))] <- 1
st2[which(staget %in% c("t2","t2a","t2b"))] <- 2
st2[which(staget %in% c("t3","t3a"))] <- 3
st2[which(staget %in% c("t4","t4b","t4d"))] <- 4

colData(brca)$ID <- pID
colData(brca)$STATUS <- pam50
idx <- union(which(pam50 == "Normal-like"), which(is.na(st2)))
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]
brca <- brca[,,1] # keep only clinical and mRNA data

tmpDir <- tempdir()
outDir <- paste(tmpDir,"results",sep="/")
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

inDir <- sprintf("%s/extdata/example_output",
	path.package("netDx"))
all_rngs <- list.dirs(inDir, recursive = FALSE)
#print(head(basename(all_rngs)))

dir(all_rngs[1])

predClasses <- c("LumA","notLumA")
predFiles <- unlist(lapply(all_rngs, function(x) 
		paste(x, "predictionResults.txt", sep = "/")))
pdf("tmp.pdf")
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()

featScores <- getFeatureScores(inDir,predClasses=c("LumA","notLumA"))
dim(featScores[[1]])
head(featScores[[1]][,1:3])
featSelNet <- lapply(featScores, function(x) {
	callFeatSel(x, fsCutoff=2, fsPctPass=0.7)
})
pathwayList <- readPathways(getExamplePathways())

xpr_genes <- rownames(assays(brca)[[1]])
pathwayList <- lapply(pathwayList, function(x) x[which(x %in% xpr_genes)])

netInfoFile <- sprintf("%s/inputNets.txt",inDir)
netInfo <- read.delim(netInfoFile,sep="\t",h=TRUE,as.is=TRUE)

Emap_res <- writeEMapInput_many(featScores,pathwayList,
	maxScore=10,pctPass=0.7,netInfo,verbose=TRUE)

message("writing EMap")
# write emap results to file)
for (g in names(Emap_res)) {
	outFile <- sprintf("%s/%s_nodeAttrs.txt",outDir,g)
	write.table(Emap_res[[g]][["nodeAttrs"]],file=outFile,
		sep="\t",col=TRUE,row=FALSE,quote=FALSE)

	outFile <- sprintf("%s/%s.gmt",outDir,g)
	conn <- base::file(outFile,"w")
	tmp <- Emap_res[[g]][["featureSets"]]

	for (cur in names(tmp)) {
		curr <- sprintf("%s\t%s\t%s", cur,cur,
			paste(tmp[[cur]],collapse="\t"))
		writeLines(curr,con=conn)
	}
print(outFile)
close(conn)
}

sessionInfo()

