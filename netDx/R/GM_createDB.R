#' Create GeneMANIA database
#'
#' @details Creates a generic_db for use with GeneMania QueryRunner.
#' The database is in tab-delimited format, and indexes are built using Apache lucene.
#' NOTE: This pipeline expects input in the form of interaction networks
#' and not profiles.
#' Profile tables have patient-by-datapoint format (e.g. patient-by-genotype)
#' Interaction networks have pairwise similarity measures:
#' <PatientA>	<PatientB>	<similarity>
#' Documentation: https://github.com/GeneMANIA/pipeline/wiki/GenericDb
#' @param netDir (char) path to dir with input networks/profiles. All
#' networks in this directory will be added to the GM database. Note:
#' This needs to be an absolute path, not relative.
#' @param patientID (char) vector of patient IDs.
#' @param outDir (char) path to dir in which GeneMANIA database is created. 
#' The database will be under \code{outDir/dataset}.
#' @param simMetric (char) similarity measure to use in converting 
#' profiles to interaction networks. 
#' @param netSfx (char) pattern for finding network files in \code{netDir}.
#' @param verbose (logical) print messages
#' @param numCores (integer) num cores for parallel processing
#' @param P2N_threshType (char) Most users shouldn't have to change this.
#' ProfileToNetworkDriver's threshold option. One of "off|auto". 
#' @param P2N_maxMissing (integer 5-100)
#' @param GMmemory (integer) Memory for GeneMANIA (in Gb)
#' @param ... params for \code{GM_writeBatchFile()}
#' @return (list). "dbDir": path to GeneMANIA database 
#' 	"netDir": path to directory with interaction networks. If profiles
#' are provided, this points to the INTERACTIONS/ subdirectory within 
#' the text-based GeneMANIA generic database
#' If the DB creation process results in an erorr, these values return 
#' NA
#' @examples
#' data(TCGA_mini,pathwayList);
#' # note: the paths in the calls below need to be absolute. If you 
#' # do not have write access to /tmp, change to a different directory.
#'	n <- makePSN_NamedMatrix(xpr,rownames(xpr),pathwayList,"/tmp/nets/",
#'		writeProfiles=TRUE); 
#'	db <- GM_createDB("/tmp/nets/",pheno$ID,"/tmp")
#' @export
GM_createDB <- function(netDir,patientID,outDir,simMetric="pearson",
	netSfx="_cont.txt$",verbose=TRUE,numCores=1L, P2N_threshType="off",
	P2N_maxMissing=100,GMmemory=4L, ...) {
	# tmpDir/ is where all the prepared files are stored.
	# GeneMANIA uses tmpDir as input to create the generic database. 
	# The database itself will be in outDir/
	tmpDir <- sprintf("%s/tmp",outDir)
	dataDir <- sprintf("%s/dataset",outDir)

	if (P2N_maxMissing < 5) PSN_maxMissing <- 5
	if (P2N_maxMissing >100) PSN_maxMissing <- 100
	if (!P2N_threshType %in% c("off","auto")) P2N_threshType <- "off"

	if (file.exists(tmpDir))  unlink(tmpDir,recursive=TRUE)
	if (file.exists(dataDir)) unlink(dataDir,recursive=TRUE)
	dir.create(dataDir)

	curwd <- getwd()
	tryCatch( {
	
	system(sprintf("cp -r %s %s", netDir,tmpDir))
	setwd(tmpDir)

	# write batch.txt and ids.txt, currently required 
	# for these scripts to work
	netList1 <- dir(path=netDir,pattern="profile$")
	# must copy networks dir to tmp instead of copying contents via
	# cp -r networks/* tmp
	# latter gives "argument list too long" error in OS/X.
	netList2	<- dir(path=netDir, pattern=netSfx)
	netList <- c(netList1,netList2)

	if (verbose) cat(sprintf("Got %i networks\n",length(netList)))
	idFile	<- sprintf("%s/ids.txt",outDir)
	GM_writeBatchFile(netDir,netList,netDir,idFile,...)
	system(sprintf("cp %s/batch.txt .", netDir))
	write.table(patientID,file=idFile,sep="\t",
				row=FALSE,col=FALSE,quote=FALSE)

	#### Step 1. placeholder files
	if (verbose) cat("\t* Creating placeholder files\n")

	# move files to tmpDir
	#file.copy(netDir,tmpDir,recursive=TRUE)
	file.copy(idFile, sprintf("%s/ids.txt",tmpDir))
	system("chmod u+w *.*")

	fBasic <- c("ATTRIBUTES.txt", "ATTRIBUTE_GROUPS.txt",
				"ONTOLOGY_CATEGORIES.txt","ONTOLOGIES.txt", "TAGS.txt", 
				"NETWORK_TAG_ASSOC.txt", "INTERACTIONS.txt")
	for (f in fBasic) {
		file.create(f)
	}

	#### Step 2. recoding networks in a format GeneMANIA prefers
	# TODO this step uses a Python script. The script is simple enough 
	# that it should be written in R. Avoid calls to other languages 
	# unless there is additional value in 
	# doing so.
	if (verbose) cat("\t* Populating database files, recoding identifiers\n")
	dir.create("profiles")
	procNet <- paste(path.package("netDx"),
					 "python/process_networks.py",sep="/")
	cmd <- sprintf("python2 %s batch.txt",procNet)
	system(cmd,wait=TRUE)
	# using new jar file
	GM_jar	<- sprintf("%s/java/genemania-cytoscape-plugin-3.5.0.jar",
						 path.package("netDx"))
	#### Step 3. (optional). convert profiles to interaction networks.
	### TODO. This step is currently inefficient. We are writing all the
	### profile files (consumes disk space) in makePSN_NamedMatrix.R
	### and then converting them 
	### enmasse to interaction networks here. 
	### Necessary because process_networks.py is prereq for ProfileToNetwork
	### Driver and that doesn't get called until the step above.
	if (length(netList1)>0) {
		cat("\t* Converting profiles to interaction networks\n")

		cl	<- makeCluster(numCores,outfile=sprintf("%s/P2N_log.txt",tmpDir))
		registerDoParallel(cl)
	
		if (simMetric=="pearson") {
			corType <- "PEARSON"
		} else if (simMetric == "MI") {
			corType <- "MUTUAL_INFORMATION"
		}

		cmd1 <- sprintf("java -Xmx%iG -cp %s org.genemania.engine.core.evaluation.ProfileToNetworkDriver", GMmemory,GM_jar)
		cmd3 <- sprintf("-proftype continuous -cor %s",corType)
		cmd5 <- sprintf("-threshold %s -maxmissing %1.1f", P2N_threshType,
			P2N_maxMissing)
		profDir <- sprintf("%s/profiles",tmpDir)
		netOutDir <- sprintf("%s/INTERACTIONS",tmpDir)
		tmpsfx <- sub("\\$","",netSfx)
		print(system.time(
		foreach (curProf=dir(path=profDir,pattern="profile$")) %dopar% {
			cmd2 <- sprintf("-in %s/%s -out %s/%s",
				profDir, curProf,netOutDir, sub(".profile",".txt",curProf))
			cmd4 <- sprintf("-syn %s/1.synonyms -keepAllTies -limitTies",
							tmpDir)
			cmd <- sprintf("%s %s %s %s %s", cmd1,cmd2,cmd3,cmd4,cmd5)
			print(cmd)
			system(cmd)
		}
		))
		stopCluster(cl)
		netSfx=".txt"
		netList2 <- dir(path=netOutDir,pattern=netSfx)

		cat(sprintf("Got %i networks from %i profiles\n", length(netList2),
			length(netList)))

		netDir <- netOutDir
		netList <- netList2; rm(netOutDir,netList2)
	}

	#### Step 4. Build GeneMANIA index
	if (verbose) cat("\t* Build GeneMANIA index\n")
	setwd(dataDir)
	cmd1 <- sprintf("java -Xmx10G -cp %s org.genemania.mediator.lucene.exporter.Generic2LuceneExporter",GM_jar)
	cmd2 <- sprintf("%s/db.cfg %s %s/colours.txt",tmpDir,tmpDir,tmpDir)
	cmd	<- sprintf("%s %s", cmd1,cmd2)
	system(cmd)

	cmd <- sprintf("mv %s/lucene_index/* %s/.",dataDir,dataDir)
	system(cmd)

	#### Step 5. Build GeneMANIA cache
	if (verbose) cat("\t* Build GeneMANIA cache\n")
	cmd1<- sprintf("java -Xmx10G -cp %s", GM_jar)
	cmd2<- "org.genemania.engine.apps.CacheBuilder"
	cmd3<- sprintf("-cachedir cache -indexDir . -networkDir %s/INTERACTIONS"
				   , tmpDir)
	cmd <- sprintf("%s %s %s",cmd1,cmd2,cmd3)
	system(cmd)

	#### Step 6. Cleanup.
	if (verbose) cat("\t * Cleanup")
	GM_xml	<- sprintf("%s/java/genemania.xml",path.package("netDx"))
	system(sprintf("cp %s %s/.", GM_xml, dataDir))

	}, error=function(ex) {
		print(ex)
		return(NA)
	}, finally={
		setwd(curwd)
	})

	return(list(dbDir=dataDir,netDir=netDir))
}
