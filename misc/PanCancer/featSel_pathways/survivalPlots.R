rm(list=ls())
library(rms)
library(survival)

rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC"
outDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival"

dirList <- list(
	KIRC_oneClinNet="output/KIRC_pathway_clinOneNet_170428",
	KIRC_pathways="output/KIRC_featSel_pathway_170426"
	)

p_val_list = list(
	"KIRC_oneClinNet"=numeric(),
	"KIRC_pathways"=numeric()
)

numTest <- c()
for (cur in names(dirList)) {
    print(cur)
    for (cur_rng in c(1:100)){

        OS_dir <- sprintf("%s/input/KIRC_OS_core.txt", rootDir)
        clincore_dir <- sprintf("%s/input/KIRC_clinical_core.txt",rootDir)

        OS_dat <- read.delim(OS_dir,sep="\t",h=T,as.is=T)
        clincore_dat <- read.delim(clincore_dir,sep="\t",h=T,as.is=T)

        netdx_res_dir <- sprintf("%s/%s/rng%i/predictionResults.txt",rootDir,dirList[[cur]],cur_rng)
        netdx_res_dat <- read.delim(netdx_res_dir,sep="\t",h=T,as.is=T)

				if (cur_rng == 1) numTest <- c(numTest, nrow(netdx_res_dat))
        names(netdx_res_dat)[1] <- "feature"

        merged_dat <- merge(OS_dat, clincore_dat, by= "feature")
        netdx_merged <- merge(netdx_res_dat, merged_dat, by = "feature")

        #Need age, status in clinical and predicted status
        netdx_merged$SurvObj <- with(netdx_merged, Surv(OS_OS, STATUS_INT == 0))

        tester <- coxph(SurvObj ~ PRED_CLASS, data = netdx_merged)

        netdx_merged.status <- with(netdx_merged, data.frame(STATUS_INT=c(0, 1), 
					PRED_CLASS=c("SURVIVEYES", "SURVIVENO")))

        # plot(survfit(tester, newdata =  netdx_merged.status), lty=c(1,2))
        # legend("bottomleft", legend=c("survyes", "survno"), lty=c(1 ,2), inset=0.02)


        km.by.pred <- npsurv(SurvObj ~ PRED_CLASS, data = netdx_merged, 
						conf.type = "log-log")
        sdf <- survdiff(SurvObj ~ PRED_CLASS ,data=netdx_merged)
        p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        p.val <- round(p.val, 5)

        log_pval <- -log(p.val, 10)
        # p_val_list[[cur]] <- c(p_val_list[[cur]],p.val)
        p_val_list[[cur]] <- c(p_val_list[[cur]],log_pval)

        # survplot(km.by.pred, label.curves = FALSE, col = c('red','blue'), conf = "none", lwd = 1.2)

        # title(sprintf("Survival curves for %s\n p-value = %s", cur,toString(p.val)))
}
}


pdf(file=sprintf("%s/KIRC_pathway_survPlot.pdf", outDir),
	width = 6,height = 6)
ctr <- 1
par(bty='n')
boxplot(p_val_list,las=1,bty='n', cex.axis=1.3,
	cex.lab=1.5,
	ylab="-log(p),log-rank test for survival",
	main="KIRC: Survival prediction")
abline(h=-log10(c(0.05,0.05/100)),col='red')
###for(cur in names(p_val_list)){
###    log_5 <- -log(0.05, 10)
###    log_bon <- -log(0.0005, 10)
###    boxplot(p_val_list[[cur]], ylim = c(0,10), ylab = "-log10(p-value)")
###    abline(h= log_5, col = 'red')
###    abline(h=log_bon, col = 'red')
###    lowest_pval <- min(p_val_list[[cur]])
###
###    title(sprintf("Survival curves for %s\n num patients = %i", 
###			cur,numTest[ctr]))
###	ctr <- ctr+1
###}
dev.off()

