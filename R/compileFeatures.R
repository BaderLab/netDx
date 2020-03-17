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
#' profiles to interaction networks. 
#' @param netSfx (char) pattern for finding network files in \code{netDir}.
#' @param verbose (logical) print messages
#' @param numCores (integer) num cores for parallel processing
#' @param P2N_threshType (char) Most users shouldn't have to change this.
#' ProfileToNetworkDriver's threshold option. One of 'off|auto'. 
#' unit testing
#' @param P2N_maxMissing (integer 5-100)
#' @param JavaMemory (integer) Memory for GeneMANIA (in Gb)
#' @param altBaseDir (char) Only use this if you're developing netDx. Used in
#' unit tests
#' @param ... params for \code{writeQueryBatchFile()}
#' @return (list). 'dbDir': path to GeneMANIA database 
#' 'netDir': path to directory with interaction networks. If profiles
#' are provided, this points to the INTERACTIONS/ subdirectory within 
#' the text-based GeneMANIA generic database
#' If the DB creation process results in an erorr, these values return 
#' NA
#' @examples
#' data(xpr,pheno)
#' pathwayList <- list(pathA=rownames(xpr)[1:10],pathB=rownames(xpr)[21:50])
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
#' tmpDir <- tempdir(); netDir <- sprintf('%s/nets',tmpDir)
#' dir.create(netDir,recursive=TRUE)
#' 
#' pheno_id <- setupFeatureDB(pheno,netDir)
#' netList <- createPSN_MultiData(dataList=dataList, groupList=groupList,
#'     pheno=pheno_id,netDir=netDir,customFunc=makeNets,verbose=TRUE)
#' 
#' outDir <- sprintf('%s/dbdir',tmpDir); dir.create(outDir)
#' dbDir <- compileFeatures(netDir,outDir)
#' @import doParallel
#' @export
compileFeatures <- function(netDir, outDir = tempdir(), 
		simMetric = "pearson", 
		netSfx = "txt$", verbose = TRUE, numCores = 1L, 
		P2N_threshType = "off", P2N_maxMissing = 100, 
    JavaMemory = 4L, altBaseDir = NULL, ...) {
    
    dataDir <- sprintf("%s/dataset", outDir)
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
    
    netList1 <- dir(path = sprintf("%s/profiles", netDir), 
				pattern = "profile$")
    netList2 <- dir(path = sprintf("%s/INTERACTIONS", netDir), 
				pattern = netSfx)
    netList <- c(netList1, netList2)
    
    if (verbose) 
        message(sprintf("Got %i networks", length(netList)))
    idFile <- sprintf("%s/ids.txt", outDir)
    writeQueryBatchFile(netDir, netList, netDir, idFile, ...)
    
    if (length(netList1) > 0) {
        if (verbose) 
            message("\t* Converting profiles to interaction networks")
        
        cl <- makeCluster(numCores, 
			outfile = sprintf("%s/P2N_log.txt", netDir))
        registerDoParallel(cl)
        
        if (simMetric == "pearson") {
            corType <- "PEARSON"
        } else if (simMetric == "MI") {
            corType <- "MUTUAL_INFORMATION"
        }
        
        args <- c(sprintf("-Xmx%iG", JavaMemory), "-cp", GM_jar)
        args <- c(args, 
					paste("org.genemania.engine.core.",
					"evaluation.ProfileToNetworkDriver",sep=""))
        args <- c(args, c("-proftype", "continuous", "-cor", corType))
        args <- c(args, c("-threshold", P2N_threshType, 
							"-maxmissing", 
							sprintf("%1.1f", P2N_maxMissing)))
        profDir <- sprintf("%s/profiles", netDir)
        netOutDir <- sprintf("%s/INTERACTIONS", netDir)
        tmpsfx <- sub("\\$", "", netSfx)
        
        curProf <- ""
        foreach(curProf = dir(path = profDir, pattern = "profile$")) %dopar% {
            args2 <- c("-in", sprintf("%s/%s", profDir, curProf))
            args2 <- c(args2, "-out", sprintf("%s/%s", netOutDir, 
								sub(".profile", ".txt", curProf)))
            args2 <- c(args2, "-syn", sprintf("%s/1.synonyms", netDir), 
								"-keepAllTies", "-limitTies")
            system2("java", args = c(args, args2), wait = TRUE, 
				stdout = NULL)
        }
        stopCluster(cl)
        netSfx = ".txt"
        netList2 <- dir(path = netOutDir, pattern = netSfx)
        
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
    args <- c(args, sprintf("%s/db.cfg", netDir), netDir, 
				sprintf("%s/colours.txt", netDir))
    system2("java", args, wait = TRUE, stdout = NULL)
    
    olddir <- sprintf("%s/lucene_index", dataDir)
    flist <- list.files(olddir, recursive = TRUE)
    dirs <- list.dirs(olddir, recursive = TRUE, full.names = FALSE)
    dirs <- setdiff(dirs, "")
    for (d in dirs) dir.create(paste(dataDir, d, sep = "/"))
    file.copy(from = paste(olddir, flist, sep = "/"), 
							to = paste(dataDir, flist,sep = "/"))
    unlink(olddir)
    
    # Build GeneMANIA cache
    if (verbose) 
        message("\t* Build GeneMANIA cache")
    args <- c("-Xmx10G", "-cp", GM_jar, 
				"org.genemania.engine.apps.CacheBuilder")
    args <- c(args, "-cachedir", "cache", "-indexDir", ".", 
				"-networkDir", sprintf("%s/INTERACTIONS",netDir), 
				"-log", sprintf("%s/test.log", netDir))
    system2("java", args = args, stdout = NULL)
    
    # Cleanup
    if (verbose) 
        message("\t * Cleanup")
    GM_xml <- system.file("extdata","genemania.xml",package="netDx")
    file.copy(from = GM_xml, to = sprintf("%s/.", dataDir))
    
    setwd(curwd)
    return(list(dbDir = dataDir, netDir = netDir))
}
