#' check that train and test have no overlaps.

inRoot <- "/mnt/data2/BaderLab/PanCancer_LUSC/output/ownTrain_170205"

for (str in c("clnical", "clinicalArna","clinicalAmir","clinicalArppa")) {
ol_set <- numeric()
for (k in 1:100) {
	inDir <- sprintf("%s/run%i/all", inRoot,k)
	queryFile <- sprintf("%s/SURVIVEYES_testQuery",inDir)
	resFile <- sprintf("%s/predictionResults.txt",inDir)
	qry <- scan(what="character",file=queryFile,skip=1,nlines=1,quiet=TRUE)
	res <- read.delim(resFile,sep="\t",h=T,as.is=T)
	ol <- length(intersect(qry,res$ID))
	ol_set <- c(ol_set, ol)
}
print(str)
print(summary(ol_set))
}




