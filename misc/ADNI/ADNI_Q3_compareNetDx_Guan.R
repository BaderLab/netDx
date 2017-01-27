# compare Q3 prediction results between guanlab and netDx.

guanFile <- "/home/spai/BaderLab/DREAM_AD/output/guan/Q3_guanlab_ADChallenge_netDxcompare.csv"
netFile <- "/home/spai/BaderLab/DREAM_AD/output/featSel_170125_test2/predClass.txt"

guan <- read.delim(guanFile,sep=",",h=T,as.is=T)
colnames(guan)[2:3] <- paste("guan",colnames(guan)[2:3],sep=".")
netDx <- read.delim(netFile,sep="\t",h=T,as.is=T)

both <- merge(x=guan,y=netDx,by="ID")
cat("guan accuracy\n")
print(sum(both$STATUS==both$guan.Diagnosis)/nrow(both))
cat("netdx accuracy\n")
print(sum(both$STATUS==both$PRED_CLASS)/nrow(both))


