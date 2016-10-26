#' wrapper to build PSN + GM database + run query + get rankings for N-way
#'
#' @details This set of steps could be used for final model performance
#' once feature selection is complete. The steps are { build networks, 
#' build GM database, run query } once per class. Finally rank patients
#' based on best class. Currently supports only full matrix data and only
#' works for two-way classification (predClass/"other"). 
#' Extension to N-way classification
#' probably doesn't require much additional work
#' @param pheno (data.frame) patient metadata. ID, class labels (STATUS),
#' and whether train/test status (TT_STATUS): a patient is query if "TRAIN",
#' else is "TEST"
#' @param pdat (data.frame) patient data. Rows are unit variables, columns 
#' are patients. Rows and columns must be named.
#' @param predClass (char) class of interest to predict.
#' @param unitSets (list) groupings of unit variables; each grouping will
#' be converted into its own PSN. 
#' @param p_GR (GRanges) GRanges of patient CNVs. Has ID column in
#' metadata, containing patient IDs
#' @param unitSet_GR (list) sets of GRanges to group CNVs (e.g.
#' could have one GRanges per pathway, corresponding to regions in that 
#' pathway
#' @param patNets (list of chars) for each class (key), subset of "unitSets"
#' for which patient nets must be created. e.g. feature-selected networks
#' for the corresponding patient subtype.
#' @param outDir (char) directory in which results must be written
#' @param numCores (integer) number of cores for parallel processing
#' @return (list) 1. predRes: patient rankings for each class; 
#' 2. predClass: predicted patient labels 
#' 3. perfStats (list): model performance stats: tp,fp,tn,fn,accuracy,ppv
#' @export
GM_predClass_once <- function(pheno,pdat,predClass,unitSets,p_GR=NULL,
		unitSet_GR=NULL,patNets,outDir,numCores=1L) {
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)	

outRes <- list()
subtypes <- names(patNets)
for (g in subtypes) {
	pDir <- sprintf("%s/%s",outDir,g); dir.create(pDir)
	pTally <- patNets[[g]]

	cat(sprintf("%s: Final # nets = %i\n", g,length(pTally)))
	# create db
	profDir <- sprintf("%s/profiles",pDir)
	tmp <- makePSN_NamedMatrix(pdat,rownames(pdat),
		unitSets[which(names(unitSets)%in% pTally)],
		profDir,verbose=F,numCores=numCores,writeProfiles=TRUE) 
	if (!is.null(p_GR)) {
		netList2 <- makePSN_RangeSets(p_GR, unitSet_GR,
				profDir,verbose=FALSE)

	}
	dbDir <- GM_createDB(profDir,pheno$ID,pDir,numCores=numCores)

		# query of all training samples for this class
	qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
			 pheno$TT_STATUS%in%"TRAIN")]
	
	# run query , get rankings
	qFile <- sprintf("%s/%s_testQuery",pDir,g)
	pTally <- paste(pTally,".profile",sep="")
	GM_writeQueryFile(qSamps,pTally,nrow(pheno),qFile)

	resFile <- netDx::runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
	system(sprintf("unlink %s", resFile))
	outRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
}

outClass <- GM_OneVAll_getClass(outRes)
both <- merge(x=pheno,y=outClass,by="ID")
print(table(both[,c("STATUS","PRED_CLASS")]))

pos=(both$STATUS %in% predClass)
tp=sum(both$PRED_CLASS[pos]==predClass)
fp=sum(both$PRED_CLASS[!pos]==predClass)
tn=sum(both$PRED_CLASS[!pos]=="other")
fn=sum(both$PRED_CLASS[pos]=="other")

perfStats <- list(tp=tp,fp=fp,tn=tn,fn=fn,acc=(tp+tn)/nrow(both),
	ppv=tp/(tp+fp))
out <- list(predRes=outRes, predClass=outClass,perfStats=perfStats)

out
}
