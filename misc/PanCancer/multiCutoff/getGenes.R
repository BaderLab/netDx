#' see gene scores. 
rm(list=ls())

dt <- format(Sys.Date(),"%y%m%d")
#outDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/output/eucscale_sp2max6000_180503"
#dNames <- sprintf("%s/rng%i", outDir,1:20)
#pdf(sprintf("GBM_topFeatures_%s_%s.pdf",setName,dt),width=6,height=13)
# ---
# lusc
###outDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/eucscale_sp26000_180503"
###dNames <- sprintf("%s/rng%i", outDir,1:11)
###setName <- "all"
###pdfFile <-sprintf("LUSC_topFeatures_%s_%s.pdf",setName,dt)
# ---
# kirc
outDir <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/output/eucscale_sp26000_180503"
dNames <- sprintf("%s/rng%i", outDir,1:5)
setName <- "all"
pdfFile <-sprintf("KIRC_topFeatures_%s_%s.pdf",setName,dt)

fSet <- list(
	SURVIVEYES=sprintf("%s/%s/SURVIVEYES/GM_results/SURVIVEYES_pathway_CV_score.txt",
		dNames,setName),
	SURVIVENO=sprintf("%s/%s/SURVIVENO/GM_results/SURVIVENO_pathway_CV_score.txt",
		dNames,setName)
)

require(netDx)

out <- list()
for (gp in names(fSet)) {
	fList <- fSet[[gp]]
	netColl <- list()
	for (scoreFile in fList) {
            tmp  <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
            colnames(tmp)[1] <- "PATHWAY_NAME"
			tmp[,1] <- sub("_cont","",tmp[,1])
                netColl[[scoreFile]] <- tmp

        }
	rpos <- regexpr("rng",fList)
	sname <- substr(fList,rpos,nchar(fList))
	spos <- regexpr("\\/",sname)
	s2name <- substr(sname,1,spos-1)
	names(netColl) <- s2name
			# filter for nets meeting cutoff criteria
			cat("* Computing consensus\n")
			cons <- getNetConsensus(netColl); x1 <- nrow(cons)
			na_sum <- rowSums(is.na(cons))
			cons <- cons[order(na_sum),]
			out[[gp]] <- cons
}

thresh <- 9
pctPassT <- lapply(out, function(x) {
	nm <- x[,1]; x <- x[,-1]
	y <- rowSums(x>=thresh,na.rm=T)
	y <- (y/ncol(x))*100
	names(y) <- nm
	y <- y[order(-y)]
	y <- y[1:20]
	y
})



pdf(sprintf("GBM_topFeatures_%s_%s.pdf",setName,dt),width=6,height=13)
pdf(pdfFile,width=8,height=13)
tryCatch({
	mar <- par("mar")
	par(mar=c(5,1,3,1))
	for (k in names(pctPassT)) {
	cur <- rev(pctPassT[[k]])
	y <- barplot(cur, horiz=TRUE,col="lightblue",border="white",
		main=sprintf("%s:%% splits with score>%i\n(N=%i splits)",
			k,thresh,length(dNames)),xlim=c(0,150))
	text(cur,y,names(cur),font=3,cex=1,pos=4)
}
},error=function(ex){print(ex)
},finally={
	dev.off()
})



