#' correlate feature with outcome
#'
#' @details Shows patient-level data for features of interest. For complex
#' features such as pathways, shows the patient-level PC projections, and
#' correlates these PC projects with outcome. It is a compact representation
#' of patient-level pathway activity. 
#' @param pheno (data.frame) must have ID and STATUS
#' @param datList (list) keys are datatypes, and values are data.frames or 
#' matrix with patients in columns and values in rows
#' @param groupList (list) unit groupings. keys are datatypes, values are lists
#' of unit groups. e.g. for pathway grouping in rna, groupList[["rna"]] would
#' be a list with pathway names as keys and genes as values.
#' @param inputNets (data.frame) contents of inputNets.txt. Two-column table
#' with column 1 having datatype and column 2 having net name. 
#' @param selFeatures (char) Features to correlate. Typically these would
#' be high-scoring features from the predictor. Names should match column 2
#' of inputNets
#' @param numPCs (integer) how many principal components to show data for
#' @param filePfx (char) prefix for output pdfs.
#' @import plotrix 
#' @return No value. Side effect of creating plots corr
#' @export
corrFeatWithOutcome <- function(pheno, datList, groupList,inputNets,
	selFeatures,numPCs=3,filePfx="corrFeat",plotLevels) {
		
	dt <- format(Sys.Date(),"%y%m%d")
	if (!class(pheno$STATUS)=="factor") pheno$STATUS <- factor(pheno$STATUS)
	
 	resMat <- matrix(0, nrow=length(selFeatures), ncol=numPCs*2)
    isDone  <- rep(FALSE, length(selFeatures))
	rownames(resMat) <- selFeatures
    #rownames(resMat) <- gsub("_"," ",sub(".profile$|_cont$","",selFeatures))
    #rownames(resMat) <- toTitleCase(rownames(resMat))
    plotList <- list();plotCtr <- 1

	for (idx in 1:length(selFeatures)) {
		curF <- selFeatures[idx]
		print(curF)
		netType <- inputNets[which(inputNets[,2]==curF),1]
		myDat <- datList[[netType]]
		mySet <- groupList[[netType]][[curF]]
		dat <- myDat[mySet,]
		dat <- na.omit(dat)

   		dat <- dat[,colSums(is.na(dat))==0]
        if (ncol(dat)>=3) {
      	  pr  <- prcomp(na.omit(t(dat)))
          pr  <- pr$x[,1:numPCs]
        } else {
        	cat(sprintf("\t\t%s: Has < 3 values!!\n",curF))
        	pr <- dat
       }

       maxDim <- min(ncol(pr),numPCs)
	 numPCs <- maxDim
       tmp <- data.frame(pr[,1:maxDim],STATUS=pheno$STATUS)
	   colnames(tmp)[1:maxDim] <- paste("PC",1:maxDim,sep="")

	combs <- combn(maxDim,2)

	for (k in 1:numPCs) {
		y <- cor.test(pr[,k],as.integer(pheno$STATUS),method="spearman")
		resMat[idx,k] <- y$estimate
		resMat[idx,3+k] <- -log10(y$p.value) # y$p.value)
	}

	for (k in 1:ncol(combs)) {
		i <- combs[1,k]; j <- combs[2,k]
		cat(sprintf("[%i %i]\n",i,j))

		# draw decision boundary in automated manner
		tmp2 <- tmp[,c("STATUS",sprintf("PC%i",i),sprintf("PC%i",j))]
		colnames(tmp2) <- c("y","x1","x2")
		mdl <- glm(y ~ ., data=tmp2,family=binomial)
		slope <- coef(mdl)[2]/(-coef(mdl)[3])
		intercept <- coef(mdl)[1]/(-coef(mdl)[3])
		
		showLeg_Flag <- (k == ncol(combs))
		if (ncol(combs) > 4) {
			pt <- 0.8; lwd <- 1 ;cex <-5 
		} else {
			pt <- 2; lwd <- 2; cex <- 12
		}

		p <- ggplot(tmp,aes_string(x=sprintf("PC%i",i),
			y=sprintf("PC%i",j)))
		p <- p + geom_point(aes(colour=factor(STATUS,levels=plotLevels)),alpha=0.6,
			size=pt,show.legend=showLeg_Flag)
		p <- p + geom_abline(intercept=intercept,slope=slope,
			colour="gray50",lwd=lwd)
		p <- p + ggtitle(sprintf("%s\ncor=%1.2f (p<%1.2e)",
			rownames(resMat)[idx],resMat[idx,i],10^-resMat[idx,3+i]))
		p <- p + theme(# legend.position="none",
			axis.ticks=element_blank(),
			axis.text=element_blank(),
			plot.title=element_text(size=cex))

		plotList[[plotCtr]] <- p
		plotCtr <- plotCtr+1
	}
}

# now plot the PC projections with colour-coded status


nr <- 3; nc <- choose(numPCs,2) 
pdf(sprintf("%s_PCview_%s.pdf", filePfx,dt),height=11,width=11)
tryCatch({
	for (sidx in seq(1,length(plotList),nr*nc)) {
	eidx <- sidx+((nr*nc)-1);
	cat(sprintf("%i-%i\n",sidx,eidx))
	if (eidx>length(plotList)) eidx <- length(plotList)
	multiplot(plotlist=plotList[sidx:eidx],
	layout=matrix(1:(nr*nc),ncol=nc,byrow=TRUE))
}
},error=function(ex) {print(ex)},finally={dev.off()})

resMat <- resMat[order(abs(resMat[,1]),decreasing=TRUE),]
pdf(sprintf("%s_PCtable_%s.pdf",filePfx,dt),height=11,width=11)
tryCatch({
#colnames(resMat)[1:3] <- paste("PC ",1:3,sep="")
par(mar = c(0.5, 35, 6.5, 0.5))
	plotrix::color2D.matplot(resMat[,1:3],show.values=TRUE,axes=F,
		xlab="",ylab="",vcex=2,vcol='black',
		cs1=c(1,1,0),cs2=c(0,1,0),cs3=c(0,1,1))
	axis(3,at=seq_len(3)-0.5,labels=colnames(resMat)[4:6],
			tick=F,cex.axis=1,line=-1)
	axis(2,at=seq_len(nrow(resMat))-0.5,
		labels=sub(".profile","",
			sub("_cont","",rev(rownames(resMat)))),tick=F,
		las=1,cex.axis=1)
},error=function(ex){print(ex)},finally={dev.off()})

outFile <- sprintf("corrFeat_table_%s.txt",dt)
write.table(resMat,file=outFile,sep="\t",col=T,row=F,quote=F)

return(resMat)
}

toTitleCase <- function(str) {
    str <- tolower(str)
    sp <- gregexpr(" ",str)
    str2 <- sapply(1:length(str), function(i) {
     z <- str[i]
    z <- paste(toupper(substr(z,1,1)),substr(z,2,nchar(z)),sep="")
     if (!sp[[i]][1]==-1) {
     for (idx in sp[[i]]) {
         z <- gsub(paste("^(.{",idx,"}).",sep=""),
            paste("\\1",toupper(substr(z,idx+1,idx+1)),sep=""),z);
        }
    }
    z
    })
    str2 <- unlist(str2)
    str2
}
