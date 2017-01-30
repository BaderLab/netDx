# code to prepare input data for BRCA predictor using
# DNA methylation, miRNA and expression data
rm(list=ls())

inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_DGG/input"
# match sample by XXXX-XX-XX-XX (up to sample ID, not including vial)

inFile <- list(
	xpr=sprintf("%s/LGG-GBM.gene_expression.normalized.txt", inDir),
	DNAm=sprintf("%s/LGG.GBM.meth.txt",inDir),
	subtype=sprintf("%s/DGG_pheno_170118.txt",inDir)
)

cat("Reading files...\n")
cat("subtypes\n")
pheno <- read.delim(inFile$subtype,sep="\t",h=T,as.is=T)
colnames(pheno)[6] <- "STATUS"
colnames(pheno)[1] <- "ID"
pheno$ID <- gsub("-",".",pheno$ID)

#xpr 
cat("xpr\n")
xpr <- read.delim(inFile$xpr,sep="\t",h=T,as.is=T)
xpr_name <- colnames(xpr)
xpr <- t(xpr)
colnames(xpr) <- gsub("-",".",colnames(xpr))
rownames(xpr) <- xpr_name

#dnam
cat("dnam\n")
dnam <- read.delim(inFile$DNA,sep="\t",h=T,as.is=T)
dnam_name <- dnam$Composite.Element.REF
dnam_anno <- dnam[,1:4]
dnam <- dnam[,-(1:4)]
dpos <- gregexpr("\\.",colnames(dnam))
tmp <- sapply(1:length(dpos),function(k) 
		substr(colnames(dnam)[k],1,dpos[[k]][3]-1))
if (any(duplicated(colnames(dnam)))) stop("dnam: duplicate found\n")
colnames(dnam) <- tmp
rownames(dnam) <- dnam_name

cat("find intersection\n")
x1 <- intersect(pheno$ID,colnames(xpr))
x3 <- intersect(x1,colnames(dnam))
cat(sprintf("%i samples\n", length(x3)))
common <- x3

dnam	<- dnam[,which(colnames(dnam)%in% common)]
xpr 	<- xpr[,which(colnames(xpr)%in% common)]
pheno	<- pheno[which(pheno$ID %in% common),]

dnam	<- dnam[,order(colnames(dnam))]
xpr		<- xpr[,order(colnames(xpr))]
pheno	<- pheno[order(pheno$ID),]
if (all.equal(pheno$ID,colnames(xpr))!=TRUE) { 
		cat("pheno:xpr no match"); browser()}
if (all.equal(pheno$ID,colnames(dnam))!=TRUE) { 
		cat("pheno:dnam no match"); browser()}

# remove missing
idx <- which(is.na(pheno$STATUS))
pheno <- pheno[-idx,]
xpr		<- xpr[,-idx]
dnam	<- dnam[,-idx]

cat("saving\n")
save(pheno,xpr,dnam,file=sprintf("%s/DGG_input_170118.Rdata",inDir))
