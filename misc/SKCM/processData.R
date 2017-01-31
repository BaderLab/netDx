#' compiles tables from patient-specific flatfiles downloaded from the 
#' tcga website

rootDir<-"/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_SKCM/input"


# --------------------------------------------------
# mirna
###inDir <- sprintf("%s/bcgsc.ca_SKCM.IlluminaHiSeq_miRNASeq.Level_3.1.11.0",
###	rootDir)
###
###fList <- dir(path=inDir,pattern="mirna.quantification.txt")
###out <- NA
###for (f in fList) {
###	nm <- sub(".mirna.quantification.txt","",f)
###	tmp <- read.delim(sprintf("%s/%s",inDir,f),sep="\t",h=T,as=T)
###	tmp <- tmp[,c(1,3)]; colnames(tmp)[2] <- nm
###
###	if (is.na(out)) out <- tmp
###	else out <- merge(x=out,y=tmp,by="miRNA_ID")
###}
###outFile <- sprintf("%s/TCGA_SKCM_mirna.txt",inDir)
###write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)

# mRNAseq
###inDir 	<- sprintf("%s/unc.edu_SKCM.IlluminaHiSeq_RNASeqV2.Level_3.1.10.0",
###	rootDir)
###fList 	<- dir(inDir, "genes.normalized_results")
###
#### maps extract names to patient IDs
###idFile	<- sprintf("%s/unc.edu_SKCM.IlluminaHiSeq_RNASeqV2.mage-tab.1.11.0/unc.edu_SKCM.IlluminaHiSeq_RNASeqV2.1.11.0.sdrf.txt",rootDir)
###ids <- read.delim(idFile,sep="\t",h=T,as.is=T)
###
###out <- NA
###for (f in fList) {
###	nm2 <- ids[which(ids$Derived.Data.File==f),2]
###	cat(sprintf("%s -> %s\n",nm,nm2))
###	tmp <- read.delim(sprintf("%s/%s",inDir,f),sep="\t",h=T,as=T)
###	colnames(tmp)[2] <- nm2
###
###	if (length(out)==1) { if (is.na(out)) out <- tmp} 
###	else out <- merge(x=out,y=tmp,by="gene_id")
###}
###outFile <- sprintf("%s/TCGA_SKCM_rnaseq.txt",inDir)
###write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)


# --------------------------------------------------
#  DNA methylation

tarPfx <- "jhu-usc.edu_SKCM.HumanMethylation450.Level_3"

for (k in 1:7) {
	curAr <- sprintf("%s.%i.9.0",tarPfx,k)
	cat(sprintf("k = %i - extracting\n",k))
	# get mid-level archive
	system(sprintf("tar xvfz %s/%s.tar.gz %s.tar.gz",rootDir,tarPfx,
		curAr))
	# get patient-level data
	system(sprintf("tar xvfz %s.tar.gz",curAr))
	system(sprintf("mv %s %s/.", curAr,rootDir))

	inDir <- sprintf("%s/%s",rootDir,curAr)
	fList <- dir(path=inDir,pattern=".txt$")
	fList <- fList[grep("jhu-usc.edu_SKCM.HumanMethylation450",
		fList)]
	
	out <- NA
	for (f in fList) {
		cat(sprintf("%s\n",f))
		dat <- read.delim(sprintf("%s/%s",inDir,f),sep="\t",
			skip=1,h=T,as.is=T)
		
		dpos <- gregexpr("\\.",f)[[1]]
		nm <- substr(f,dpos[5]+1,dpos[6]-1)

		tmp <- aggregate(dat$Beta_value,
			by=list(gene=dat$Gene_Symbol),FUN=mean,na.rm=T)
		colnames(tmp)[2] <- nm
		tmp[,2] <- signif(tmp[,2],digits=3)
		if (length(out)==1) {if (is.na(out)) out <- tmp}
		else out <- merge(x=out,y=tmp,by="gene")
	}
	outFile <- sprintf("%s/%s.compiled.txt",rootDir,curAr)
	write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)
	system(sprintf("rm -r %s/%s",rootDir,curAr))
	system(sprintf("unlink %s.tar.gz",curAr))
}


