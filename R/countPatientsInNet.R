#' Count number of patients in a network
#' 
#' @details This functionality is needed to count patient overlap when 
#' input data is in a form that results in highly missing data, rather than
#' when the same measures are available for almost all patients. An example
#' application is when patient networks are based on unique genomic events
#' in each patients (e.g. CNVs or indels), rather than 'full-matrix' data
#' (e.g. questionnaires or gene expression matrices). The former scenario
#' requires an update in the list of eligible networks each time some type
#' of patient subsetting is applied (e.g. label enrichment, or train/test
#' split). A matrix with patient/network membership serves as a lookup
#' table to prune networks as feature selection proceeds
#' @param netDir (char) dir with network set
#' @param fList (char) filenames of interaction networks to count in
#' @param ids (char) patient IDs to look for
#' @return (matrix) Size P by N, where P is num patients and N is 
#' number of networks networks; a[i,j] =1 if patient i in network j, else 0
#' @examples
#' d <- tempdir()
#' pids <- paste("P",1:5,sep="")
#' m1 <- matrix(c("P1","P1","P2","P2","P3","P4",1,1,1),
#' 	byrow=FALSE,ncol=3)
#' write.table(m1,file=sprintf("%s/net1.txt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#' m2 <- matrix(c("P3","P4",1),nrow=1)
#' write.table(m2,file=sprintf("%s/net2.txt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#' x <- countPatientsInNet(d,c("net1.txt","net2.txt"), pids)
#' @export
countPatientsInNet <- function(netDir, fList, ids) {
    
    outmat <- matrix(0, nrow = length(ids), ncol = length(fList))
    colnames(outmat) <- fList
    rownames(outmat) <- ids
    
    ctr <- 1
    for (f in fList) {
        dat <- read.delim(sprintf("%s/%s", netDir, f), sep = "\t", 
						header = FALSE, 
            as.is = TRUE)
        memb <- c(dat[, 1], dat[, 2])  # patients in this network
        outmat[which(ids %in% memb), ctr] <- 1
        
        ctr <- ctr + 1
    }
    
    return(outmat)
}
