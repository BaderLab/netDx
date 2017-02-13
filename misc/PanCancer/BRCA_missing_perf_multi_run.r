# assess performance of BRCA classifier for different levels of missing
# data
rm(list=ls())
require(netDx)
require(RColorBrewer)

rngRange <- 1:10 # number of train/test split iterations run

inDir <- "/mnt/data2/BaderLab/TCGA_BRCA/output/msng_170213"

combSet <- paste("miss",c(10,50,70,85,90,95,99),sep="")
cols <- c(brewer.pal(n=length(combSet),name="Blues"), "darkblue")

cat(sprintf("Got %i combs\n", length(combSet)))

all_overall_acc <- list()
all_tot <- list()
all_f1 <- list()
all_acc <- list()
all_ppv <- list()
for (rngSeed in rngRange) {
	cat(sprintf("RNG %i\n",rngSeed))
    out <- list()
    overall_acc <- numeric()
    pctMiss <- numeric()
    numInt <- numeric()
    curRoc	<- list()

    for (cur in combSet) {
        inf <- sprintf("%s/rng%i/%s/predictionResults.txt",
                       inDir,rngSeed,cur)
		cat(sprintf("\t%s ", cur))
        dat <- read.delim(inf,sep="\t",h=T,as.is=T)
        dat <- dat[-which(dat$STATUS %in% "Normal"),]
        out[[cur]] <- perfCalc_multiClass(dat$STATUS,dat$PRED_CLASS)*100
        overall_acc <- c(overall_acc, 
                         sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)
    }
	cat("\n")
    names(overall_acc) <- combSet

    tot <- unlist(lapply(out,function(x) sum(x[1,1:4])/100))
    f1 <- unlist(lapply(out, function(x) x[nrow(x),7]))
    acc <- unlist(lapply(out, function(x) x[nrow(x),8]))
    ppv <- unlist(lapply(out, function(x) x[nrow(x),5]))
     
    all_tot[[rngSeed]] <- tot   
    all_f1[[rngSeed]] <- f1   
    all_acc[[rngSeed]] <- acc   
    all_ppv[[rngSeed]] <- ppv   
    all_overall_acc[[rngSeed]] <- overall_acc   

}

overall_avg_f1 <- c()
overall_avg_tot <- c()
overall_avg_acc <- c()
overall_avg_ppv <- c()
overall_avg_overall_acc <- c()


overall_sd_f1 <- c()
overall_sd_tot <- c()
overall_sd_acc <- c()
overall_sd_ppv <- c()
overall_sd_overall_acc <- c()

for (index in seq(1:length(combSet))) {
    avg_f1 <- c()
    avg_tot <- c()
    avg_acc <- c()
    avg_ppv <- c()
    avg_overall_acc <- c()
    for (rngSeed in rngRange) {
        avg_f1 <- c(avg_f1,all_f1[[rngSeed]][[index]])
        avg_tot <- c(avg_tot,all_tot[[rngSeed]][[index]])
        avg_acc <- c(avg_acc,all_acc[[rngSeed]][[index]])
        avg_ppv <- c(avg_ppv,all_ppv[[rngSeed]][[index]])
        avg_overall_acc <- c(avg_overall_acc,all_overall_acc[[rngSeed]][[index]])
    }
    current <- combSet[index]
    
    current_f1 <- paste(current, '.f1', sep = '')
    current_tot <- current
    current_acc <- paste(current, '.acc', sep = '')
    current_ppv <- paste(current, '.ppv', sep = '')
    current_overall_acc <- current
    
    overall_avg_f1[current_f1] <- mean(avg_f1)
    overall_avg_tot[current_tot] <- mean(avg_tot)
    overall_avg_acc[current_acc] <- mean(avg_acc)
    overall_avg_ppv[current_ppv] <- mean(avg_ppv)
    overall_avg_overall_acc[current_overall_acc] <- mean(avg_overall_acc)    
    
    overall_sd_f1[current_f1] <- sd(avg_f1)
    overall_sd_tot[current_tot] <- sd(avg_tot)
    overall_sd_acc[current_acc] <- sd(avg_acc)
    overall_sd_ppv[current_ppv] <- sd(avg_ppv)
    overall_sd_overall_acc[current_overall_acc] <- sd(avg_overall_acc)
}


pdf("BRCA_missing.pdf",width=8,height=4)
tryCatch({

# mean pairwise F1
f1_plot <- barplot(overall_avg_f1,main="mean F1",col=cols,
	ylab="F1",ylim=c(0,100),las =2, cex.names=0.75)

bottom_f1 <- overall_avg_f1 - overall_sd_f1
top_f1 <- overall_avg_f1 + overall_sd_f1
names(bottom_f1) <- NULL
names(top_f1) <- NULL       
arrows(f1_plot, bottom_f1, f1_plot,
       top_f1, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

# mean pairwise accuracy
acc_plot <- barplot(overall_avg_acc,main="mean accuracy",col=cols,
	ylab="accuracy",ylim=c(0,100),las =2, cex.names=0.75)
abline(h=c(50,85),col='red',lwd=2,lty=2)
bottom_acc <- overall_avg_acc - overall_sd_acc
top_acc <- overall_avg_acc + overall_sd_acc
names(bottom_acc) <- NULL
names(top_acc) <- NULL       
arrows(acc_plot, bottom_acc, acc_plot,
       top_acc, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

#ppv 
ppv_plot <- barplot(overall_avg_ppv,main="mean PPV",col=cols,ylab="PPV",
	ylim=c(0,100),las =2, cex.names=0.75)
abline(h=50,col='red',lwd=2,lty=2)
bottom_ppv<- overall_avg_ppv - overall_sd_ppv
top_ppv <- overall_avg_ppv + overall_sd_ppv
names(bottom_ppv) <- NULL
names(top_ppv) <- NULL       
arrows(ppv_plot, bottom_ppv, ppv_plot,
       top_ppv, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

# plot overall accuracy
overall_acc_plot <- barplot(overall_avg_overall_acc,main="Accuracy",col=cols,
	ylab="accuracy",ylim=c(0,100),las =2, cex.names=0.75)
abline(h=c(50,85),col='red',lwd=2,lty=2)
bottom_oall<- overall_avg_overall_acc - overall_sd_overall_acc
top_oall <- overall_avg_overall_acc + overall_sd_overall_acc
names(bottom_oall) <- NULL
names(top_oall) <- NULL       
arrows(overall_acc_plot, bottom_oall, overall_acc_plot,
       top_oall, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)


barplot(overall_avg_tot,main="total classified", las = 2, cex.names=0.75)
# names(numInt) <- names(tot)
# barplot(numInt/1000,main="# interactions (x1000)")
# names(pctMiss) <- names(tot)
# barplot(pctMiss*100,main="% missing in profile",ylim=c(0,100))

df <- do.call("rbind",all_overall_acc)
boxplot(df,ylab="% accuracy",main=sprintf("%% Overall accuracy (N=%i)",
	nrow(df)),col=cols,pars=list(boxwex=0.4),bty='n',ylim=c(0,100),las=1)
abline(h=25,col='red',lty=3,lwd=2)
},error=function(ex){
	print(ex)
},finally={
	dev.off()
})
