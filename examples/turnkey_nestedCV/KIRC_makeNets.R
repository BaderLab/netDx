# User-defined function to make nets
KIRC_makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c()
	# make RNA nets: group by pathway
	if (!is.null(groupList[["rna"]])) { 
	netList <- makePSN_NamedMatrix(dataList$rna, 
					rownames(dataList$rna),
			   	groupList[["rna"]],netDir,verbose=FALSE, 
			  	writeProfiles=TRUE,...) 
	netList <- unlist(netList)
	cat(sprintf("Made %i RNA pathway nets\n", length(netList)))
	}
	
	# make clinical nets
	netList2 <- c()
	if (!is.null(groupList[["clinical"]])) {
	netList2 <- makePSN_NamedMatrix(dataList$clinical, 
		rownames(dataList$clinical),
		groupList[["clinical"]],netDir, simMetric="custom",customFunc=normDiff,
		sparsify=TRUE,verbose=TRUE,append=TRUE,...)
	}
	netList2 <- unlist(netList2)
	cat(sprintf("Made %i clinical nets\n", length(netList2)))
	netList <- c(netList,netList2) 
	cat(sprintf("Total of %i nets\n", length(netList)))
	return(netList)
}
