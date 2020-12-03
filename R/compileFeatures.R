#' Create GeneMANIA database
#'
#' @details Creates a generic_db for use with GeneMania QueryRunner.
#' The database is in tab-delimited format, and indexes are built using Apache 
#' lucene.
#' NOTE: This pipeline expects input in the form of interaction networks
#' and not profiles.
#' Profile tables have patient-by-datapoint format (e.g. patient-by-genotype)
#' Interaction networks have pairwise similarity measures:
#' <PatientA> <PatientB><similarity>
#' Documentation: https://github.com/GeneMANIA/pipeline/wiki/GenericDb
#' @param netDir (char) path to dir with input networks/profiles. All
#' networks in this directory will be added to the GM database. Note:
#' This needs to be an absolute path, not relative.
#' @param outDir (char) path to dir in which GeneMANIA database is created. 
#' The database will be under \code{outDir/dataset}.
#' @param simMetric (char) similarity measure to use in converting 
#' @param netSfx (char) pattern for finding network files in \code{netDir}.
#' profiles to interaction networks. 
#' @param verbose (logical) print messages
#' @param numCores (integer) num cores for parallel processing
#' @param P2N_threshType (char) Most users shouldn't have to change this.
#' ProfileToNetworkDriver's threshold option. One of 'off|auto'. 
#' unit testing
#' @param P2N_maxMissing (integer 5-100)
#' @param JavaMemory (integer) Memory for GeneMANIA (in Gb)
#' @param altBaseDir (char) Only use this if you're developing netDx. Used in
#' unit tests
#' @param debugMode (logical) when TRUE runs jobs in serial instead of parallel and 
#' prints verbose messages. Also prints system Java calls and prints all standard out
#' and error output associated with these calls.
#' @param ... params for \code{writeQueryBatchFile()}
#' @return (list). 'dbDir': path to GeneMANIA database 
#' 'netDir': path to directory with interaction networks. If profiles
#' are provided, this points to the INTERACTIONS/ subdirectory within 
#' the text-based GeneMANIA generic database
#' If the DB creation process results in an erorr, these values return 
#' NA
#' @examples
#' data(xpr,pheno)
#' pathwayList <- list(pathA=rownames(xpr)[1:10],
#'	pathB=rownames(xpr)[21:50])
#' 
#' dataList <- list(rna=xpr)  #only one layer type
#' groupList <- list(rna=pathwayList) # group genes by pathways
#' 
#' makeNets <- function(dataList, groupList, netDir,...) {
#'     netList <- makePSN_NamedMatrix(dataList[['rna']],
#'					rownames(dataList[['rna']]),
#'          groupList[['rna']],netDir,verbose=FALSE,
#'					writeProfiles=TRUE,...)
#'     unlist(netList)
#' }
#' tmpDir <- tempdir(); netDir <- paste(tmpDir,"nets",
#'	sep=getFileSep())
#' dir.create(netDir,recursive=TRUE)
#' 
#' pheno_id <- setupFeatureDB(pheno,netDir)
#' netList <- createPSN_MultiData(dataList=dataList, groupList=groupList,
#'     pheno=pheno_id,netDir=netDir,customFunc=makeNets,verbose=TRUE)
#' 
#' outDir <- paste(tmpDir,'dbdir',sep=getFileSep()); 
#'	dir.create(outDir)
#' dbDir <- compileFeatures(netDir,outDir)
#' @import doParallel
#' @export
compileFeatures <- function(netDir, outDir = tempdir(), 
		simMetric = "pearson", 
		netSfx = "txt$", verbose = TRUE, numCores = 1L, 
		P2N_threshType = "off", P2N_maxMissing = 100, 
    JavaMemory = 4L, altBaseDir = NULL, debugMode=FALSE,...) {
    
    dataDir <- paste(outDir,"dataset",sep=getFileSep())
    GM_jar <- getGMjar_path()
    
    if (P2N_maxMissing < 5) 
        PSN_maxMissing <- 5
    if (P2N_maxMissing > 100) 
        PSN_maxMissing <- 100
    if (!P2N_threshType %in% c("off", "auto")) 
        P2N_threshType <- "off"
    
    if (!file.exists(dataDir)) dir.create(dataDir)
    curwd <- getwd()
    setwd(netDir)
    
    netList1 <- dir(path = paste(netDir,"profiles",sep=getFileSep()),
				pattern = "profile$")
    netList2 <- dir(path = paste(netDir,"INTERACTIONS",sep=getFileSep()),
				pattern = netSfx)
    netList <- c(netList1, netList2)

   tryCatch({
        .jcheck(silent=FALSE)
    },error=function(ex){
        .jinit()
        .jaddClassPath(getGMjar_path())
    })
    x <- getGMjar_path()
    if (!x %in% .jclassPath()) .jaddClassPath(x)
    jObj <- .jnew("org.genemania.engine.core.evaluation.ProfileToNetworkDriver")    
    
    if (verbose) 
        message(sprintf("Got %i networks", length(netList)))
    idFile <- paste(outDir,"ids.txt", sep=getFileSep())
    writeQueryBatchFile(netDir, netList, netDir, idFile, ...)
    
    if (length(netList1) > 0) {
        if (verbose) 
            message("\t* Converting profiles to interaction networks")
        
        cl <- makeCluster(numCores, 
			outfile = paste(netDir,"P2N_log.txt",
			sep=getFileSep()))
        registerDoParallel(cl)
        
        if (simMetric == "pearson") {
            corType <- "PEARSON"
        } else if (simMetric == "MI") {
            corType <- "MUTUAL_INFORMATION"
        }
        
        args <- c("-proftype", "continuous", "-cor", corType)
        args <- c(args, c("-threshold", P2N_threshType, 
							"-maxmissing", 
							sprintf("%1.1f", P2N_maxMissing)))
        profDir <- paste(netDir,"profiles",sep=getFileSep())
        netOutDir <- paste(netDir,"INTERACTIONS",sep=getFileSep())
        tmpsfx <- sub("\\$", "", netSfx)
        
        curProf <- ""
		`%myinfix%` <- ifelse(debugMode, `%do%`, `%dopar%`)
        foreach(curProf = dir(path = profDir, pattern = "profile$")) %myinfix% {
            args2 <- c("-in", paste(profDir,curProf,sep=getFileSep()))
            args2 <- c(args2, "-out", 
		paste(netOutDir,sub(".profile", ".txt", curProf),
			sep=getFileSep()))
            args2 <- c(args2, "-syn", 
		paste(netDir,"1.synonyms",sep=getFileSep()),
			"-keepAllTies", "-limitTies")

    
      tryCatch({
        .jcheck(silent=FALSE)
    },error=function(ex){
        .jinit()
        .jaddClassPath(getGMjar_path())
    })
    jObj <- .jnew("org.genemania.engine.core.evaluation.ProfileToNetworkDriver")
	.jcall(jObj,"V",method="main",c(args,args2))
        }
        stopCluster(cl)
        netSfx = ".txt"
        netList2 <- dir(path = netOutDir, pattern = netSfx)
				msg2 <- paste("This problem usually occurs because of a failed",
							"Java call. Try upgrading to Java 11. If that doesn't",
							"work, contact shraddha.pai@utoronto.ca with a copy",
							"of the complete log file after running buildPredict()",
							"with debugMode=TRUE",sep="\n")

				if (length(netList2)<length(netList)) {
				  warnings(paste("",
						"---------------------------------",
						"One or more profiles did not successfully convert to PSNs!",
						"This usually happens because of the underlying call to a",
						"Java library failed. Upgrading to Java 11 usually fixes",
						"this problem. If not, please send a copy of the detailed",
						"error message from the call below to",
						"shraddha.pai@utoronto.ca","",sep="\n")
					)
					
					curProf <- dir(profDir,"profile$")[1]
	        args2 <- c("-in", paste(profDir, curProf,sep=getFileSep()))
	        args2 <- c(args2, "-out", paste(netOutDir, 
				sub(".profile", ".txt", curProf),
				sep=getFileSep()))
	        args2 <- c(args2, "-syn", 
			paste(netDir,"1.synonyms",sep=getFileSep()),
				"-keepAllTies", "-limitTies")
		tmp <- paste(c(args,args2),collapse=" ")
        .jcall(jObj,"V",method="main",c(args,args2))       
		stop("Stopping netDx now. See error message above.")
	 }  
        if (verbose) 
            message(sprintf("Got %i networks from %i profiles", 
		length(netList2), length(netList)))
        netList <- netList2
        rm(netOutDir, netList2)
    }
    
    #### Build GeneMANIA index
    if (verbose) 
        message("\t* Build GeneMANIA index")
    setwd(dataDir)
    args <- c("-Xmx10G", "-cp", GM_jar)
    args <- c(args, paste("org.genemania.mediator.lucene.",
			"exporter.Generic2LuceneExporter",sep=""))

     args <- c(args, 
        paste(netDir,"db.cfg",sep=getFileSep()), 
        netDir, # base path
	 	paste(netDir,"colours.txt",sep=getFileSep()),
        "none", # profile - set null value per code
        paste(dataDir,"lucene_index",sep=getFileSep())) # need this for .jcall because absolute path
      tryCatch({
        .jcheck(silent=FALSE)
    },error=function(ex){
        .jinit()
        .jaddClassPath(getGMjar_path())
    })
     jObj3 <- .jnew(args[4],check=TRUE)
     .jcall(jObj3,"V",method="main",args[-(1:4)])
    
    olddir <- paste(dataDir,"lucene_index", sep=getFileSep())
    flist <- list.files(olddir, recursive = TRUE)
    dirs <- list.dirs(olddir, recursive = TRUE, full.names = FALSE)
    dirs <- setdiff(dirs, "")
    for (d in dirs) dir.create(paste(dataDir, d, sep = getFileSep()))
    file.copy(from = paste(olddir, flist, sep = getFileSep()), 
	to = paste(dataDir, flist,sep =getFileSep()))
    unlink(olddir)

	# Check: need to replace commas used as decimal separators, into periods
	tmp <- dir(path=sprintf("%s/INTERACTIONS",netDir),pattern="txt$")[1]
	tmp <- sprintf("%s/INTERACTIONS/%s",netDir,tmp)
 	if (sum(grepl(pattern=",",readLines(tmp,n=1))>0)) { # detect comma
		replacePattern(path=sprintf("%s/INTERACTIONS",netDir))
	}
    
    # Build GeneMANIA cache
    if (verbose) 
        message("\t* Build GeneMANIA cache")

    args <- c('-cachedir', paste(getwd(),'cache',sep=getFileSep()),
                 '-indexDir', getwd(), 
				'-networkDir', 
			paste(netDir,"INTERACTIONS",sep=getFileSep()), 
				"-log",             
			paste(netDir,"test.log",sep=getFileSep()))

    tryCatch({
        .jcheck(silent=FALSE)
    },error=function(ex){
        .jinit()
        .jaddClassPath(getGMjar_path())
    })
    jObj2 <- .jnew("org.genemania.engine.apps.CacheBuilder")
    .jcall(jObj2,"V",method="main",args)
    
    # Cleanup
    if (verbose) 
        message("\t * Cleanup")
    GM_xml <- system.file("extdata","genemania.xml",package="netDx")
    file.copy(from = GM_xml, to = paste(dataDir,".",sep=getFileSep()))
    
    setwd(curwd)
    return(list(dbDir = dataDir, netDir = netDir))
}

#' Replace pattern in all files in dir
#' @description find/replace pattern in all files of specified file type
#' in specified directory. Needed to modify number format when intefacing
#' with GeneMANIA, on  French locale machines. Without this step,
#' CacheBuilder throws error with commas. 
#' @param pattern (char) pattern to find
#' @param target (char) pattern to replace
#' @param path (char) dir to replace pattern in
#' @param fileType (char) pattern for files to replace pattern in
#' @return No value. Files have patterns replaced in place.
#' @export
replacePattern <- function(pattern=",",target=".",path=getwd(),fileType="txt$") {
fList <- dir(path,fileType)
for (currF in fList) {
        fFull <- sprintf("%s/%s",path,currF)
        tx <- readLines(fFull)
        tx2 <- gsub(",",".",tx)
        writeLines(tx2,con=fFull)
}
}
