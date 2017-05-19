# for each feature selected variable check its concordance
# with the variable we want to predict (survival)
rm(list=ls())

require(plotrix) # to plot coloured correlation table
require(ggplot2) # geom_dotplot

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

### input directory for data
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
outDir <- sprintf("%s/2017_PanCancer_Survival",rootDir)

# must define universe of clinical variables
clinicalVars <- c("age","stage","grade","Karnofsky","gender")
clinicalVars <- paste(clinicalVars,"_cont",sep="")
consCutoff <- 10 
consPctPass <- 0.7
datSets <- c("KIRC") 

# to create pheno table
clinList <- list(
KIRC=sprintf("%s/2017_TCGA_KIRC/input/KIRC_clinical_core.txt",rootDir)
)
# binary survival
survList <- list(
KIRC=sprintf("%s/2017_TCGA_KIRC/input/KIRC_binary_survival.txt",rootDir)
)
# directory with consensus nets
consNetDir <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/pathOnly_consNets_170509", rootDir)
)


consListingDir <- outDir
dt <- format(Sys.Date(),"%y%m%d")

# ########## Loop over cancer datasets
for (curSet in c("KIRC")) { #LUSC","KIRC","OV","GBM")) {
	cat(sprintf("--------------------------\n"))
	cat(sprintf("%s\n--------------------------\n\n",curSet))
	
	dataDir <- consNetDir[[curSet]]
	profileDir <- list(
		SURVIVEYES=sprintf("%s/SURVIVEYES/networks",dataDir),
		SURVIVENO=sprintf("%s/SURVIVENO/networks",dataDir)
	)
	
	# the pathway scores for each class
	netScoreFile <- list(
		SURVIVEYES=sprintf("%s/pathwaysOnly_170502/%s_thresh10_pctPass0.70_SURVIVEYES_consNets.txt", outDir,curSet),
		SURVIVENO=sprintf("%s/pathwaysOnly_170502/%s_thresh10_pctPass0.70_SURVIVENO_consNets.txt", outDir,curSet)
	)
	
	# load survival data
	curSurv <- read.delim(survList[[curSet]],sep="\t",h=T,as.is=T)
	colnames(curSurv)[1] <- "ID"

	clrs <- plotrix::color.scale(x=curSurv$is_alive,alpha=0.5,
		extremes=c("orange","blue"))
	curSurv$patColor <- clrs; rm(clrs)
	
	# ----------------------------------------
	# compare with survival on a per-feature basis
	for (gps in names(netScoreFile)) {
		cat(sprintf("Group %s\n", gps))
		pTally <- read.delim(netScoreFile[[gps]],sep="\t",h=T,as.is=T)[,1]

	
		resMat <- matrix(0, nrow=length(pTally), ncol=6)
			isDone	<- rep(FALSE, length(pTally))
		rownames(resMat) <- gsub("_"," ",sub(".profile$|_cont$","",pTally))
		rownames(resMat) <- toTitleCase(rownames(resMat))
		plotList <- list();plotCtr <- 1
	
		#### Handle profiles
		idxSet <- grep("profile$",pTally)
		cat(sprintf("\t%i profiles\n", length(idxSet)))
		for (idx in idxSet) {
			#cat(sprintf("%s\n", pTally[idx]))
			pFile <- sprintf("%s/%s", profileDir[[gps]], pTally[idx])
			#if (file.exists(pFile)) {
			dat <- read.delim(pFile,sep="\t",h=F,as.is=T)
			rownames(dat) <-dat[,1] 
			pat_ID <-  dat[,1];
			dat <- dat[,-1]
		
			midx <- match(pat_ID, curSurv$ID)
			mysurv <- curSurv[midx,]

			dat <- dat[,colSums(is.na(dat))==0]
			if (ncol(dat)>=3) {
				pr 	<- prcomp(na.omit(dat))
				pr	<- pr$x[,1:3]
			} else {
				cat(sprintf("\t\t%s: Has < 3 values!!\n",pTally[idx]))
				pr <- dat
			}
			maxDim <- min(ncol(pr),3)
			tmp <- data.frame(pr[,1:maxDim],SURVIVE=mysurv$is_alive)
			colnames(tmp)[1:maxDim] <- paste("PC",1:maxDim,sep="")
			
			for (k in 1:maxDim) {
				y <- cor.test(pr[,k],mysurv$is_alive,method="spearman")
				resMat[idx,k] <- y$estimate
				resMat[idx,3+k] <- -log10(y$p.value) # y$p.value)
				if (k==1) { i <-1; j <- 2}
				else if (k==2) {
					if (maxDim < 3) {i <- 1; j <- 2}
					else { i <- 1; j <- 3}
				} else { i <- 2; j <- 3}
				
				p <- ggplot(tmp,aes_string(x=sprintf("PC%i",i),
					y=sprintf("PC%i",j)))
				p <- p + geom_point(aes(colour=factor(SURVIVE)),alpha=0.6,
						size=3)
				p <- p + ggtitle(sprintf("%s\ncor=%1.2f (p<%1.2e)",
					rownames(resMat)[idx],
					resMat[idx,k],10^-resMat[idx,3+k]))
				p <- p + theme(# legend.position="none",
						axis.ticks=element_blank(),
						axis.text=element_blank())
				plotList[[plotCtr]] <- p
				plotCtr <- plotCtr+1
				#plot(pr$x[,i],pr$x[,j],
				#	xlab=sprintf("PC%i",i),ylab=sprintf("PC%i",j),
				#	pch=16,col=mysurv$patColor,
			#		main=ttl,cex.axis=1.3,bty='n',cex=1.3,
			#		cex.main=1.3)
			}
			#}
			isDone[idx] <- TRUE
		}

		plotList2 <- list()
		plotCtr <- 1

		### Now handle mutation based nets
		idxSet <- grep("cont$",pTally)
		idx2 <- which(!(pTally[idxSet] %in%  clinicalVars))
		idxSet <- idxSet[idx2]
		cat(sprintf("\t%i binary non-clinical nets\n", length(idxSet)))
		for (idx in idxSet) {
			pFile 	<- sprintf("%s/%s.txt", profileDir[[gps]],pTally[idx])
			#if (file.exists(pFile)) {
			dat 	<- read.delim(pFile,sep="\t",h=F,as.is=T)
			in_net	<- unique(c(dat[,1],dat[,2]))
			mysurv	<- curSurv
			mysurv$IN_NET	 <- 0
			mysurv$IN_NET[mysurv$ID %in% in_net] <- 1 
			
			tmp <- matrix(NA,nrow=2,ncol=2)
			tmp[1,1] <- sum(mysurv$IN_NET>0 & mysurv$is_alive>0)
			tmp[1,2] <- sum(mysurv$IN_NET>0 & mysurv$is_alive<1)
			tmp[2,1] <- sum(mysurv$IN_NET<1 & mysurv$is_alive>0)
			tmp[2,2] <- sum(mysurv$IN_NET<1 & mysurv$is_alive<1)
			
			resMat[idx,] <- c(tmp[1,1]/sum(mysurv$is_alive>0),0,0,
				-log10(fisher.test(tmp)$p.value),1,1)
			
			tmp <- data.frame(SURV=mysurv$is_alive,VAR=mysurv$IN_NET)
			p <- ggplot(tmp, aes(x=factor(SURV),fill=factor(VAR),
				y=(..count..)/sum(..count..))) + geom_bar(position='dodge')
			p <- p + ggtitle(sprintf("%s\nFisher p < %1.2e",
				pTally[idx],10^-resMat[idx,4]))
			plotList2[[plotCtr]] <- p
			plotCtr <- plotCtr+1

			#}
			isDone[idx] <- TRUE
		}

		### Now handle clinical vars
		cat("\t done feature-wise analysis\n")
		cat(sprintf("\t%i clinical features\n", sum(!isDone)))
		curClin <- read.delim(clinList[[curSet]],sep="\t",h=T,as.is=T)
		colnames(curClin)[1] <- "ID"
		if (any(grep("performance_score",colnames(curClin)))) {
			colnames(curClin)[grep("performance_score",colnames(curClin))] <- "Karnofsky"
		}
		pheno <- merge(x=curClin,y=curSurv,by="ID")
			for (idx in which(!isDone)) {
				var <- pheno[,sub("_cont","",pTally[idx])]
				if (class(var)!="numeric") {
					cat(sprintf("\t\t%s -> converting to factor\n",
						pTally[idx]))
					var <- as.numeric(factor(var))
				}
				var_clean <- var[!is.na(var)]
				surv_clean <- pheno$is_alive[!is.na(var)]	
			
				x <- cor.test(var_clean,surv_clean,method="spearman")
				resMat[idx,] <- c(x$estimate,0,0,
					-log10(x$p.value),0,0)
				isDone[idx] <- TRUE

				binwidth <- max(diff(range(var_clean))/10,1)
				if (binwidth==1) dotwidth<-0.1 else dotwidth=0.2
 
				tmp <- data.frame(SURV=surv_clean,VAR=var_clean)
				p <- ggplot(tmp,aes(x=VAR,fill=factor(SURV)))+
					geom_dotplot(method="histodot", 
						binwidth=binwidth,
						dotsize=dotwidth,colour=NA,alpha=0.8,
						position="dodge"
						)
				p <- p + scale_y_continuous(NULL,breaks=NULL)
				p <- p + ggtitle(sprintf("%s: cor=%1.2f (p=%1.2e)",
						pTally[idx],x$estimate,x$p.value))
				p <- p + theme(legend.position="none",
						axis.ticks.y=element_blank(),
						axis.text.y=element_blank())
				plotList2[[plotCtr]] <- p
				plotCtr <- plotCtr+1
				#print(p)
		}

		pdf(sprintf("%s/%s_%s_PCA_%s.pdf",outDir,curSet,gps,dt),
			width=11,height=11)
		source("../multiplot.R")
		tryCatch({
			for (sidx in seq(1,length(plotList),9)) {
				eidx <- sidx+8; 
				cat(sprintf("%i-%i\n",sidx,eidx))
				if (eidx>length(plotList)) eidx <- length(plotList)	
				multiplot(plotlist=plotList[sidx:eidx],
					layout=matrix(1:9,ncol=3,byrow=TRUE))
			}
			if (length(plotList2)>0) {
			for (sidx in seq(1,length(plotList2),4)) {
				eidx <- sidx+3; 
				cat(sprintf("%i-%i\n",sidx,eidx))
				if (eidx>length(plotList2)) eidx <- length(plotList2)	
				multiplot(plotlist=plotList2[sidx:eidx],
					layout=matrix(1:4,ncol=2,byrow=TRUE))
			}
	}
		if (nrow(resMat)>1) {
		resMat <- resMat[order(apply(resMat[,4:6],1,max),decreasing=TRUE),,
			drop=F]
	}

		colnames(resMat) <-c("PC1","PC2","PC3",
			"-log(PC1p)","-log(PC2p)","-log(PC3p)")
		
		tmp <- resMat
		rownames(tmp) <- pTally
	write.table(tmp,
			file=sprintf("%s/%s_%s_correlations_%s.txt",
				outDir,curSet,gps,dt),
				sep="\t",col=T,row=T,quote=F)
		rm(tmp)
		tmp <- round(resMat[,4:6],1)
		if (nrow(tmp)>1) {
			isTop <- which(apply(tmp,1,max)>=2);
		} else {
			isTop <- which(max(tmp[1,])>=2)
		}
		cat(sprintf("%i nets with p < 0.01\n",length(isTop)))
	
		# write table twice - first, all results. then those with p < 0.01
		writeTables <- "all"
		if (length(isTop)>=1) writeTables <- c(writeTables,"top")
		for (writeVersion in writeTables) {
			if (writeVersion == "top") resMat <- resMat[isTop,,drop=FALSE]
		# plot correlation table
		if (nrow(resMat)>30) vcex <- 1 
		else if (nrow(resMat)>20) vcex <- 1.5
		else vcex <-2.5

		par(mar = c(0.5, 35, 5.5, 0.5))
		plotrix::color2D.matplot(resMat[,1:3],show.values=TRUE,axes=F,
			xlab="",ylab="",vcex=vcex,vcol='black',
			cs1=c(1,1,0),cs2=c(0,1,0), cs3=c(0,1,1),
			xrange=c(-1,0,1),main=sprintf("%s:%s features",curSet,gps)) 
		# column names
		axis(3,at=seq_len(3)-0.5,labels=colnames(resMat)[1:3],
				tick=F,cex.axis=2,line=-1)
		# rownames
		axis(2,at=seq_len(nrow(resMat))-0.5,
			labels=sub(".profile","",rev(rownames(resMat))),tick=F,
			las=1,cex.axis=1.7)
		# plot pvalues
		par(mar = c(0.5, 35, 6.5, 0.5))
		plotrix::color2D.matplot(resMat[,4:6],show.values=TRUE,axes=F,
			xlab="",ylab="",vcex=vcex,vcol='black',
			cs1=c(1,0),cs2=c(1,0),cs3=c(1,1),
			main=sprintf("%s:%s features",curSet,gps))
		axis(3,at=seq_len(3)-0.5,labels=colnames(resMat)[4:6],
				tick=F,cex.axis=1,line=-1)
		axis(2,at=seq_len(nrow(resMat))-0.5,
			labels=sub(".profile","",
				sub("_cont","",rev(rownames(resMat)))),tick=F,
			las=1,cex.axis=1)
		}

		},error=function(ex) {
			print(ex)
		},finally={
			dev.off()
		})
	
		}
}



