# save train/test split to check with guanlab.


require(netDx)
trainProp <- 0.67
inFile <- "/home/spai/BaderLab/DREAM_AD/input/GuanLab/C3/ADNI_Training_Q3_new.csv_matlab.csv"
outDir <- dirname(inFile)

dat <- read.delim(inFile,sep=",",h=F,as.is=T)
colnames(dat)[1] <- "ID"
colnames(dat)[7] <- "STATUS"
for (k in seq(5,100,5)) {
		cat(sprintf("Seed = %i\n",k))
TT_STATUS <- splitTestTrain(dat,pctT=trainProp,setSeed=k)
idx <- which(TT_STATUS %in% "TRAIN")
write.table(dat[idx,],sep=",",file=sprintf("%s.TRAIN_%i.txt",inFile,k),
			col=F,row=F)
idx <- which(TT_STATUS %in% "TEST")
write.table(dat[idx,],sep=",",file=sprintf("%s.TEST_%i.txt",inFile,k),
			col=F,row=F)
}


