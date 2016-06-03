#' Parse GMT file and return pathways as list
#' 
#' @details The GMT file format currently supported should match the ones
#' found at http://downloads.baderlab.org. The original GMT file format is:
#' <set name><set description><member 1><member 2>...<member N>, 
#' one row per set, with values tab-delimited.
#' The version at baderlab.org has additional unique formatting of the
#' <set name> column as follows:
#' <pathway_full_name>%<pathway_source>%<pathway_source_id>
#' 
#' Example:
#' UREA CYCLE%HUMANCYC%PWY-4984	urea cycle	ASS1	ARG1	CPS1	ASL	OTC
#' ABACAVIR METABOLISM%REACTOME DATABASE ID RELEASE 55%2161541	Abacavir metabolism	ADH1A	GUK1	ADAL	PCK1	NT5C2
#' 
#' This function requires the specific formatting of the first column
#' to assign the key name of the output list (see \code{useIDasName} 
#' argument).
#' @param fname (char) path to pathway file in gmt format
#'	pathway score to include pathway in the filter list
#' @param MIN_SIZE (integer) min num genes allowed in a pathway. Pathways
#' with fewer number of genes are excluded from the output list
#' @param MAX_SIZE (integer) max num genes allowed in a pathway. Pathways
#' with gene counts greater than this are excluded from the output list
#' @param EXCLUDE_KEGG: (boolean) If TRUE exclude KEGG pathways. Our
#' experience has been that some KEGG gene sets are to broad to be 
#' physiologically relevant
#' @param IDasName: (boolean) Value for key in output list. 
#' If TRUE, uses db name and ID as name (e.g.  KEGG:hsa04940)
#' If FALSE, pathway name. If TRUE, 
#' @param getOrigNames (logical) when TRUE also returns a mapping of the
#' cleaned pathway names to the original names
#' @return Depends on value of getOrigNames. If FALSE (Default), list with
#' pathway name as key, vector of genes as value. If TRUE, returns list of
#' length two, (1) geneSets: pathway-gene mappings as default, 
#' (2) pNames: data.frame with original and cleaned names.
#' @examples pathFile <- sprintf("%s/extdata/pathways.gmt",
#'		path.package("netDx"))
#'	pathwayList    <- readPathways(pathFile)
#' 
#' @export
readPathways <- function(fname,MIN_SIZE=10L, MAX_SIZE=500L, 
	EXCLUDE_KEGG=TRUE,IDasName=FALSE,verbose=TRUE,getOrigNames=FALSE) {

# change locale to accommodate nonstandard chars in pathway names
oldLocale   <- Sys.getlocale("LC_ALL")
Sys.setlocale("LC_ALL","C")
    out <- list()
    # read list of master pathways
	if (verbose) cat("---------------------------------------\n")
	if (verbose) cat(sprintf("File: %s\n\n", basename(fname)))
    f	<- file(fname,"r")
    # TODO: deal with duplicate pathway names
    
    #pName <- list()
	ctr <- 0
	options(warn=1)
    repeat {
        s	<- scan(f, what="character",nlines=1,quiet=TRUE,sep="\t")
        if (length(s)==0) break;

        pPos<- gregexpr("%",s[1])[[1]];
		src <- ""
		src_id	<- ""
        if (pPos[1]==-1) {
            #cat("\n\n% symbol not found in pathway name")
			s[1]	<- s[1]
        } else {

			src		<- substr(s[1],pPos[1]+1,pPos[2]-1)
			src_id	<- substr(s[1],pPos[2]+1,nchar(s[1]))
			if (IDasName) 
				s[1]	<- paste(src,src_id,sep=":")
			else 
        		s[1]	<- substr(s[1],1,pPos[1]-1)
		}
		if (!EXCLUDE_KEGG || (src!="KEGG")) {
			idx <- which(s=="") # remove trailing blank rows.
			if (any(idx)) s <- s[-idx]
        	out[[s[1]]]	<- s[3:length(s)]
        	#pName[[s[1]]] <- s[2] # stores pathway source - prob not needed
		} 
		ctr <- ctr+1
    }
    close(f)
	if (verbose) {
		cat(sprintf("Read %i pathways in total, internal list has %i entries\n",
				ctr, length(out)))
    	cat(sprintf("\tFILTER: sets with num genes in [%i, %i]\n",
				MIN_SIZE,MAX_SIZE))
	}
    ln <- unlist(lapply(out, length))
    idx	<- which(ln < MIN_SIZE | ln >= MAX_SIZE)
    out[idx] <- NULL
    #pName[idx] <- NULL
    if (verbose) cat(sprintf("\t  => %i pathways excluded\n\t  => %i left\n", 
				length(idx),length(out)))

	# clean pathway names
	nm	<- cleanPathwayName(names(out))
	if (getOrigNames) {
		pnames <- cbind(names(out), nm)
    	names(out) <- nm
		out <- list(geneSets=out,pNames=pnames)
	} else {
		names(out) <- nm
	}
    return(out)
}
