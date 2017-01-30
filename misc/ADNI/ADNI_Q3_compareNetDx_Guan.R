# compare Q3 prediction results between guanlab and netDx.
rm(list=ls())

guanDir <- "/home/spai/BaderLab/DREAM_AD/output/guan"
netDir <- "/home/spai/BaderLab/DREAM_AD/output/featSel20RNG_170127"

netF <- 
guanD <- list()
setSeeds <- seq(5,100,5)
acc <- matrix(nrow=length(setSeeds),ncol=12)
colnames(acc) <- c("seed","n_netDx","n_guan","n_merge",
		"acc_netDx","acc_guan",
		"f1_netDx","f1_guan",
		"ppv_netDx","ppv_guan",
		"recall_netDx","recall_guan")
ctr <- 1
source("perfCalc_tmp.R")
for (setSeed in setSeeds) {
	guan <- read.delim(sprintf("%s/Q3_GuanLab_out_%i.csv",guanDir,setSeed),
		sep=",",h=T,as.is=T)
	guan[which(guan$Diagnosis %in% "MCI"),3] <- "LMCI"
	netDx <- read.delim(sprintf("%s/RNG%i/predClass.txt",netDir,setSeed),
		sep="\t",h=T,as.is=T)
	colnames(guan)[2:3] <- paste("guan",colnames(guan)[2:3],sep=".")

	both <- merge(x=netDx,y=guan,by="ID")
	
	# collect perf measures
	nperf	<- as.data.frame(perfCalc_tmp(both$STATUS, both$PRED_CLASS))
	gperf	<- as.data.frame(perfCalc_tmp(both$STATUS, both$guan.Diagnosis))
 
	# accuracy
	a_net <- sum(both$STATUS==both$PRED_CLASS)/nrow(both)
	a_guan <- sum(both$STATUS==both$guan.Diagnosis)/nrow(both)

	acc[ctr,] <- c(setSeed, nrow(netDx),nrow(guan),nrow(both),
			a_net, a_guan, 
			nperf$f1[4],gperf$f1[4], #f1
		    nperf$ppv[4],gperf$ppv[4],  #ppv
		    nperf$recall[4],gperf$recall[4]) #recall
	ctr <- ctr+1
}

# plot mean and sd
par(las=1,cex.axis=1.3,cex.lab=1.3,mfrow=c(2,2))
for (k in c(5,7,9,11)) {
	a <- k; b <- k+1
	mu <- colMeans(acc[,a:b])
	stdev <- unlist(sapply(a:b, function(x) sd(acc[,x])))

	ttl <- substr(colnames(acc)[a],1,regexpr("_",colnames(acc)[a])-1)
	x <- barplot(mu,ylim=c(0,1),
		ylab=sprintf("diagnosis %s", ttl))
	segments(x0=x,x1=x,y0=mu-stdev,y1=mu+stdev)
}
