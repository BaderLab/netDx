#' Similarity functions.
#' these will be eventually moved into netDx
require(reshape2)

#' mutual information
sim.mi <- function(m) {
	tmp <- infotheo::discretize(m)
	infotheo::mutinformation(tmp)
}

#' radial basis function
#' @param m (matrix) data, columns are patients
#' @param nm (char) kernel to use, prefix to kernlab::*dot() functions.
#' e.g. rbf,tanh,laplace
sim.kern <- function(m,nm="rbf",sigma=0.05) {
	if (nm=="rbf") {
		func <- kernlab::rbfdot(sigma)
		cat(sprintf("Sigma = %1.2f\n", sigma))
	} else if (nm == "tanh") {
		cat("using tanh\n")
		func <- kernlab::tanhdot()
	}
	idx <- combinat::combn(1:ncol(m),2)
	out <- matrix(NA,nrow=ncol(m),ncol=ncol(m))
	for (comb in 1:ncol(idx)) {
		i <- idx[1,comb]; j <- idx[2,comb]
		x <- func(m[,i],m[,j])
		out[i,j] <- x; out[j,i] <- x
	}
	diag(out) <- 1
	colnames(out)<- colnames(m);
	rownames(out) <- colnames(m)
	out
}

#' cosine similarity
sim.cos <- function(m) {
	out <- matrix(NA,nrow=ncol(m),ncol=ncol(m))
	rownames(out) <- colnames(m); 
	colnames(out) <- colnames(m)
	idx <- combn(ncol(m),2)
	for (k in 1:ncol(idx)) {
		i <-idx[1,k]; j <- idx[2,k]
		out[i,j] <- lsa::cosine(m[,i],m[,j])
	}
	return(out)
}


# normalized difference 
# x is vector of values, one per patient (e.g. ages)
sim.normDiff <- function(x) {
    #if (nrow(x)>=1) x <- x[1,]
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX  <- max(x,na.rm=T)-min(x,na.rm=T)
    
    out <- matrix(NA,nrow=n,ncol=n);
    # weight between i and j is
    # wt(i,j) = 1 - (abs(x[i]-x[j])/(max(x)-min(x)))
    for (j in 1:n) out[,j] <- 1-(abs((x-x[j])/rngX))
    rownames(out) <- nm; colnames(out)<- nm
    out
}

# takes average of normdiff of each row in x
sim.normDiff2 <- function(x) {
	sim <- matrix(0,nrow=ncol(x),ncol=ncol(x))
	for (k in 1:nrow(x)) {
		tmp <- sim.normDiff(x[k,,drop=FALSE])
		sim <- sim + tmp
		rownames(sim) <- rownames(tmp)
		colnames(sim) <- colnames(tmp)
	}
	sim <- sim/nrow(x)
	sim
}

#' plot similarity matrix
#' 
#' @param x (matrix) similarity matrix
#' @return No value. Side effect of plotting the matrix with cell-colours
#' indicating strength of similarity. 
#' @examples 
#' rawDat <- matrix(rnorm(100),ncol=10)
#' sim <- abs(cor(rawDat))
#' plotSimMat(sim)
#' plotSimMat(NetPreProc::Sparsify.matrix(sim,k=5)) # sparsified
plotSimMat <- function(x) {
cellClrs <- matrix(NA, nrow=nrow(x),ncol=ncol(x))
plotrix::color2D.matplot(x,xrange=c(0,1),
	extremes=c("lightblue","darkblue"),
	na.color='grey50',show.values=(ncol(x)<15))
}

# given psn plot intra- and inter-class similarity
# matrix must have upper populated, lower can be empty
#' @param s1 (matrix) similarity matrix. If 3-column table provided, assumes
#' it's a SIF
#' @param c1,c2 (char) vector of patients in each of the two groups
plotSim <- function(s1,name="simfun",c1,c2,logT=FALSE) {
	if (logT) s1 <- log10(s1+.Machine$double.eps)
	if (ncol(s1) == 3) {
		cat("assuming SIF provided\n")
		colnames(s1) <- c("Var1","Var2","value")
		# do nothing
	} else {
		s1[lower.tri(s1,diag=TRUE)] <- NA
		s1 <- na.omit(melt(s1))
	}
	out <- list(
		pp=s1$value[which(s1$Var1 %in% c1 & s1$Var2 %in% c1)],
		mm=s1$value[which(s1$Var1 %in% c2 & s1$Var2 %in% c2)],
		pm=s1$value[union(which(s1$Var1 %in% c1 & s1$Var2 %in% c2),
				which(s1$Var1 %in% c2 & s1$Var2 %in% c1))]
	)
	cat(sprintf("Similarity by %s\n",name))
	cat("Median similarity\n")
	cat(sprintf("pp=%1.2f ; mm=%1.2f; pm = %1.2f\n",
		median(out$pp),median(out$mm),median(out$pm)))
	cat(sprintf("pp < pm: p < %1.2e\n", 
		wilcox.test(out$pp,out$pm,alternative="greater")$p.value))
	cat(sprintf("mm < pm: p < %1.2e\n", 
		wilcox.test(out$mm,out$pm,alternative="greater")$p.value))
	cat("------------\n")
	boxplot(out,main=name)
}





