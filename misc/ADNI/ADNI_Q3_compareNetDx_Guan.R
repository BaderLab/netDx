# compare Q3 prediction results between guanlab and netDx.
rm(list=ls())

guanDir <- "/home/spai/BaderLab/DREAM_AD/output/guan"
netDir <- "/home/spai/BaderLab/DREAM_AD/output/featSel20RNG_170127"

netF <- 
guanD <- list()
setSeeds <- seq(5,100,5)
acc <- matrix(nrow=length(setSeeds),ncol=6)
colnames(acc) <- c("seed","n_netDx","n_guan","n_merge",
		"acc_netDx","acc_guan")
ctr <- 1
for (setSeed in setSeeds) {
	guan <- read.delim(sprintf("%s/Q3_GuanLab_out_%i.csv",guanDir,setSeed),
		sep=",",h=T,as.is=T)
	netDx <- read.delim(sprintf("%s/RNG%i/predClass.txt",netDir,setSeed),
		sep="\t",h=T,as.is=T)
	colnames(guan)[2:3] <- paste("guan",colnames(guan)[2:3],sep=".")

	both <- merge(x=netDx,y=guan,by="ID")
	# accuracy
	a_net <- sum(both$STATUS==both$PRED_CLASS)/nrow(both)
	a_guan <- sum(both$STATUS==both$guan.Diagnosis)/nrow(both)

	acc[ctr,] <- c(setSeed, nrow(netDx),nrow(guan),nrow(both),
			a_net, a_guan)
	ctr <- ctr+1
}

# plot mean and sd
mu <- colMeans(acc[,5:6])
stdev <- unlist(sapply(5:6, function(x) sd(acc[,x])))

par(las=1,cex.axis=1.3,cex.lab=1.3)
x <- barplot(mu,ylim=c(0,1),ylab="diagnosis accuracy %")
segments(x0=x,x1=x,y0=mu-stdev,y1=mu+stdev)

