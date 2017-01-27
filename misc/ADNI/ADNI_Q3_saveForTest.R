# save train/test split to check with guanlab.


require(netDx)
trainProp <- 0.67
inFile <- "/home/spai/BaderLab/DREAM_AD/input/GuanLab/C3/ADNI_Training_Q3_new.csv_matlab.csv"
outDir <- dirname(inFile)


dat <- read.delim(inFile,sep=",",h=F,as.is=T)
colnames(dat)[1] <- "ID"
colnames(dat)[7] <- "STATUS"
set.seed(102);
TT_STATUS <- splitTestTrain(dat,pctT=trainProp)
idx <- which(TT_STATUS %in% "TRAIN")
write.table(dat[idx,],sep=",",file=sprintf("%s.TRAIN.txt",inFile),col=F,row=F)
idx <- which(TT_STATUS %in% "TEST")
write.table(dat[idx,],sep=",",file=sprintf("%s.TEST.txt",inFile),col=F,row=F)


