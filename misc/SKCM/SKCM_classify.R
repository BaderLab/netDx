#' classify by survival for TCGA SKCM

inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_SKCM/input"
outRoot <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_SKCM/output"

dt <- format(Sys.Date(),"%y%m%d")
# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/cutaneous_melanoma_clinical.txt",inDir),
	rna=sprintf("%s/TCGA_SKCM_rnaseq.txt",inDir),
	mirna=sprintf("%s/TCGA_SKCM_mirna.txt",inDir),
	dnam=sprintf("%s/TCGA_SKCM_DNAm.txt",inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
pheno[,1] <- gsub("-",".",pheno[,1])
colnames(pheno)[1] <- "ID"

# rnaseq 
cat("\t* RNA\n")
rna <- read.delim(inFiles$rna,sep="\t",h=T,as.is=T)
xpr_name <- rna[,1]
rownames(rna) <- xpr_name
rna <- rna[,-1]
dpos <- gregexpr("\\.",colnames(rna))
sampName <- sapply(1:length(dpos), function(k) 
	substr(colnames(rna)[k],1,dpos[[k]][4]-2))
colnames(rna) <- sampName

# mirna
cat("\t* miR\n")
mir <- read.delim(inFiles$mirna,sep="\t",h=T,as.is=T)
rownames(mir) <- mir[,1]
mir_name <- mir[,1]
mir <- mir[,-1]
dpos <- gregexpr("\\.",colnames(mir))
sampName <- sapply(1:length(dpos), function(k) 
	substr(colnames(mir)[k],1,dpos[[k]][4]-2))
colnames(mir) <- sampName

# DNA methylation
cat("\t* DNA methylation\n")
dnam <- read.delim(inFiles$dnam,sep="\t",h=T,as.is=T)
dnam_name <- dnam[,1]
dnam <- dnam[,-1]
sum(colnames(dnam) %in% pheno$ID)

# prepare predictor variable
tmp <- pheno$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU
tmp[grep("Not ", tmp)] <- NA
tmp[grep("-", tmp)] <- NA
tmp <- as.integer(tmp)
survstat <- rep(NA,length(tmp))
survstat[which(tmp <= quantile(tmp,0.33,na.rm=T))] <- "LOW"
survstat[which(tmp >= quantile(tmp,0.66,na.rm=T))] <- "HIGH"
pheno$SURVIVAL <- survstat

# remove samples with NA values
pheno$STATUS <- pheno$SURVIVAL
pheno <- pheno[-which(is.na(pheno$STATUS)),]

# include only data for patients in classifier
rna <- rna[,which(colnames(rna)%in% pheno$ID)]
mir <- mir[,which(colnames(mir)%in% pheno$ID)]
dnam <- dnam[,which(colnames(dnam)%in% pheno$ID)]

common <- intersect(colnames(rna),colnames(mir))
common <- intersect(common, colnames(dnam))
pheno <- subset(pheno,ID %in% common)

# make column order the same to rbind
midx <- match(pheno$ID,colnames(rna))
rna <- rna[,midx]
midx <- match(pheno$ID,colnames(dnam))
dnam <- dnam[,midx]
midx <- match(pheno$ID,colnames(mir))
mir <- mir[,midx]

alldat <- rbind(rna,mir,dnam)
alldat_names <- c(xpr_name, mir_name, dnam_name)
rownames(alldat) <- alldat_names
rm(midx)

# ----------------------------------------------------------
# build classifier
outDir <- sprintf("%s/integrate_%s",outRoot,dt)
numCores <- 2L
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

netSets <- list(
	xpr=xpr_name,
	dnam=dnam_name,
	mir=mir_name)

# create similarity networks for all datatypes.
require(netDx)
profDir <- sprintf("%s/profiles",outDir)
dir.create(profDir)
tmp <- makePSN_NamedMatrix(alldat, rownames(alldat),netSets,
	   profDir,verbose=FALSE, numCores=numCores,
	   writeProfiles=TRUE)

# create a single GM database with training and test samples
dbDir 	<- GM_createDB(profDir,pheno$ID,outDir,numCores=numCores)

# datatype combinations to try.
combList <- list(xpr="xpr",dnam="dnam",mir="mir",
	xprdna=c("xpr","dnam"),
	all="all")

subtypes <- unique(pheno$STATUS)

# ---------------------------------------------------------------------
# run test for different train/test splits
sink(sprintf("%s/log.txt",outDir),split=TRUE)
tryCatch({
for (setSeed in seq(10,100,10)) {
	curd <- sprintf("%s/RNG%i",outDir,setSeed)
	dir.create(curd)
	
	# split into train/test for all 4 classes.
	pheno$TT_STATUS <- splitTestTrain(pheno,setSeed=setSeed)
	print(table(pheno[,c("STATUS","TT_STATUS")]))
	
	# run prediction using each data combination in turn
	for (cur in names(combList)) {
		t0 <- Sys.time()
		cat(sprintf("%s\n",cur)) 
		pDir <- sprintf("%s/%s",curd, cur)
		dir.create(pDir)
	
		outRes <- list()
		for (g in subtypes) {
			qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
				 pheno$TT_STATUS%in% "TRAIN")]
			qFile <- sprintf("%s/%s_testQuery",pDir,g)
	
			# use only selected nets
			if (cur == "all") nets <- "all" 
			else nets <- paste(combList[[cur]],".profile",sep="")
			GM_writeQueryFile(qSamps,nets,nrow(pheno),qFile)
	
			resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
			outRes[[g]] <- GM_getQueryROC(
				sprintf("%s.PRANK",resFile),pheno,g)
		}
		outClass <- GM_OneVAll_getClass(outRes)
		both <- merge(x=pheno,y=outClass,by="ID")
	
		#both <- both[-which(both$STATUS %in% "Normal"),]
		acc <- sum(both$STATUS == both$PRED_CLASS)/nrow(both)
		print(table(both[,c("STATUS","PRED_CLASS")]))
		
		save(outRes,file=sprintf("%s/outRes.Rdata",pDir))
		write.table(both,file=sprintf("%s/predictionResults.txt",pDir),
					sep="\t",col=TRUE,row=FALSE,quote=FALSE)
	
		cat(sprintf("%s complete\n", cur))
		print(Sys.time()-t0)
	}
}
},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})


