# code to prepare input data for BRCA predictor using
# DNA methylation, miRNA and expression data
rm(list=ls())

inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/TCGA_breastCancer/input"

# match sample by XXXX-XX-XX-XX (up to sample ID, not including vial)

inFile <- list(
	xpr=sprintf("%s/BRCA.exp.547.med.txt", inDir),
	DNAm=sprintf("%s/BRCA.methylation.574probes.802.txt",inDir),
	mirna=sprintf("%s/BRCA.780.mimat.txt",inDir),
	subtype=sprintf("%s/BRCA.547.PAM50.SigClust.Subtypes.txt",inDir)
)

cat("Reading files...\n")
cat("subtypes\n")
pheno <- read.delim(inFile$subtype,sep="\t",h=T,as.is=T)
pheno$Sample <- gsub("-",".",pheno$Sample)
dpos <- gregexpr("\\.",pheno$Sample)
sampShort <- sapply(1:length(dpos),function(k) 
		substr(pheno$Sample[[k]],1,dpos[[k]][4]-2))
pheno$ID <- sampShort

# remove metastatic tumours
idx <- which(pheno$Type %in% "metastatic")
cat(sprintf("removing %i metastatic tumours\n", length(idx)))
pheno <- pheno[-idx,]

#xpr 
cat("xpr\n")
xpr <- read.delim(inFile$xpr,sep="\t",h=T,as.is=T)
xpr_name <- xpr[,1]
xpr <- xpr[,-1]
dpos <- gregexpr("\\.",colnames(xpr))
xprShort <- sapply(1:length(dpos),function(k) 
		substr(colnames(xpr)[k],1,dpos[[k]][4]-2))
colnames(xpr) <- xprShort
if (any(duplicated(xprShort))) stop("xpr: duplicate found\n")

#dnam
cat("dnam\n")
dnam <- read.delim(inFile$DNA,sep="\t",h=T,as.is=T)
colnames(dnam) <- substr(colnames(dnam),1,nchar(colnames(dnam))-1)
if (any(duplicated(colnames(dnam)))) stop("dnam: duplicate found\n")

#mirna
cat("mirna\n")
mirna <- read.delim(inFile$mirna,sep=",",h=T,as.is=T)
mirna_name <- mirna[,1]
mirna <- mirna[,-1]
dpos <- gregexpr("\\.",colnames(mirna))
mirShort <- sapply(1:length(dpos),function(k) 
		substr(colnames(mirna)[k],1,dpos[[k]][4]-2))
if (any(duplicated(mirShort))) stop("mirna: duplicate found\n")
colnames(mirna) <- mirShort

cat("find intersection\n")
x1 <- intersect(pheno$ID,xprShort)
x2 <- intersect(x1,mirShort)
x3 <- intersect(x2,colnames(dnam))
cat(sprintf("%i samples\n", length(x3)))
common <- x3

dnam	<- dnam[,which(colnames(dnam)%in% common)]
xpr 	<- xpr[,which(colnames(xpr)%in% common)]
mirna 	<- mirna[,which(colnames(mirna)%in% common)]
pheno	<- pheno[which(pheno$ID %in% common),]

dnam	<- dnam[,order(colnames(dnam))]
xpr		<- xpr[,order(colnames(xpr))]
mirna	<- mirna[,order(colnames(mirna))]
pheno	<- pheno[order(pheno$ID),]
if (all.equal(pheno$ID,colnames(xpr))!=TRUE) { 
		cat("pheno:xpr no match"); browser()}
if (all.equal(pheno$ID,colnames(dnam))!=TRUE) { 
		cat("pheno:dnam no match"); browser()}
if (all.equal(pheno$ID,colnames(mirna))!=TRUE) { 
		cat("pheno:mirna no match"); browser()}

rownames(xpr) <- xpr_name
rownames(mirna) <- mirna_name

cat("saving\n")
save(pheno,xpr,dnam,mirna,file="BRCA_3data_input_170117.Rdata")




