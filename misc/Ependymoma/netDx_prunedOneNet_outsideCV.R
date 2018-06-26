# Ependymoma
rm(list=ls())

require(netDx)
require(netDx.examples)
require(glmnet)

rootDir <- "/home/shraddhapai/BaderLab/2017_Ependymoma"
inDir <- sprintf("%s/input/netDx_prepared",rootDir)
outDir <- sprintf("%s/output",rootDir)
load(sprintf("%s/Ependymoma_cohortMerged_180125.Rdata",inDir))

# exclude ST
idx <- which(pheno$STATUS=="ST") 
pheno <- pheno[-idx,]
xpr <- xpr[,-idx]
    
makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c()
	# make RNA nets: group by pathway
	if (!is.null(groupList[["rna"]])) { 
	netList <- makePSN_NamedMatrix(dataList$rna, 
					rownames(dataList$rna),
			   	groupList[["rna"]],netDir,verbose=FALSE, 
			  	writeProfiles=TRUE,...) 
	netList <- unlist(netList)
	cat(sprintf("Made %i RNA pathway nets\n", length(netList)))
	}
	cat(sprintf("Total of %i nets\n", length(netList)))
	return(netList)
}

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/Epen_lassoOutsideCV_%s",outDir,dt)
if (!file.exists(megaDir)) dir.create(megaDir)

gps <- list(rna=list(rna=rownames(xpr)))
dats <- list(rna=xpr)
pheno$STATUS <- as.character(droplevels(pheno$STATUS))

	#### Test of doing something WRONG - lasso outside cv loop
	cat("Prefiltering enabled\n")
	for (nm in names(dats)) {
		if (nrow(dats[[nm]])<2)  # only has one var, take it.
			vars <- rownames(dats[[nm]])
		else { 
			set.seed(123)
			fit <- cv.glmnet(x=t(na.omit(dats[[nm]])),
					y=factor(pheno$STATUS), family="binomial", alpha=1) # lasso
			wt <- abs(coef(fit,s="lambda.min")[,1])
			vars <- setdiff(names(wt)[which(wt>.Machine$double.eps)],"(Intercept)")
			}
		if (length(vars)>0) {
			tmp <- dats[[nm]]
			tmp <- tmp[which(rownames(tmp) %in% vars),,drop=FALSE]
			dats[[nm]] <- tmp
			gps[[nm]] <- list(rna=rownames(dats[[nm]]))
		} 
		cat(sprintf("%s: %s pruned\n",nm,length(vars)))
		}

browser()
runPredictor_nestedCV(pheno,
   dataList=dats,groupList=gps,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred",megaDir),
   numCores=8L,nFoldCV=10L, CVcutoff=9L,numSplits=10L,startAt=1L,
	preFilter=FALSE)
