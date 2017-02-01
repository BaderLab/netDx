
fList <- dir(pattern="mirna.quantification.txt")
out <- NA
for (f in fList) {
		print(f)
	nm <- sub(".mirna.quantification.txt","",f)
	tmp <- read.delim(f,sep="\t",h=T,as=T)
	tmp <- tmp[,c(1,3)]; colnames(tmp)[2] <- nm

	if (is.na(out)) out <- tmp
	else out <- merge(x=out,y=tmp,by="miRNA_ID")
}




