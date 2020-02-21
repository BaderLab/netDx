#' moves interaction networks when compiling database for sparse genetic
#' workflow
#'
#' @param netDir (char) source directory
#' @param outDir (char) target directory
#' @param pheno (data.frame) contains patient ID and STATUS
#' @param fileSfx (char) suffix to strip from network file names before
#' registering in metadata tables
#' @return No value. Side effect of moving interaction nets to target
#' directory and creating network-related metadata files used to compile
#' feature database
#' @importFrom utils write.table
moveInteractionNets <- function(netDir,outDir,pheno,fileSfx="_cont.txt") {
netList <- dir(path=netDir,pattern=fileSfx)
	netID <- data.frame(ID = seq_len(length(netList)),
                name = netList, ID = seq_len(length(netList)),
        name2 = netList, 0, 1, stringsAsFactors = TRUE)	
		dir.create(sprintf("%s/INTERACTIONS",netDir))
        for (p in netList) {
			dat <- read.delim(sprintf("%s/%s",netDir,p),sep="\t",
				header=FALSE,as.is=TRUE)
			dat2 <- dat
			dat2[,1] <- pheno$INTERNAL_ID[match(dat[,1],pheno$ID)]
			dat2[,2] <- pheno$INTERNAL_ID[match(dat[,2],pheno$ID)]
			write.table(dat2,file=sprintf("%s/INTERACTIONS/1.%i.txt", 
                netDir, netID$ID[which(netID$name == p)]),
				sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
        }
    
    # write NETWORKS.txt
    write.table(netID, file = sprintf("%s/NETWORKS.txt", netDir), 
			sep = "\t", col.names = FALSE, 
        row.names = FALSE, quote = FALSE)
    
    # write NETWORK_GROUPS.txt
    con <- file(sprintf("%s/NETWORK_GROUPS.txt", netDir), "w")
    write(paste(1, "dummy_group", "geneset_1", "dummy_group", 1, sep = "\t"),
	file = con)
    close(con)
    
    con <- file(sprintf("%s/NETWORK_METADATA.txt", netDir), "w")
    tmp <- paste(netID$ID, "", "", "", "", "", "", "", 
			"", "", 0, "", "", 0, "", 
        "", "", "", "", sep = "\t")
    write.table(tmp, file = con, sep = "\t", col.names = FALSE, 	
			row.names = FALSE, 
        quote = FALSE)
    close(con)
}
