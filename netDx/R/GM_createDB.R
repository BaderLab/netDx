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
#' networks in this directory will be added to the GM database
#' @param patientID (char) vector of patient IDs.
#' @param outDir (char) path to dir in which GeneMANIA database is created. 
#' The database will be under \code{outDir/dataset}.
#' @param simMetric (char) similarity measure to use in converting 
#' profiles to interaction networks. 
#' @param netSfx (char) pattern for finding network files in \code{netDir}.
#' @param verbose (logical) print messages
#' @param numCores (integer) num cores for parallel processing
#' @param attrib_DF (data.frame) attribute info. NULL for no attributes.
#' Column names: 1) attribute_group (char) attribute group name
#' 2) attribute_name (char) attribute name
#' 3) ID (integer) patient ID
#' 4) attribute_value (char) attribute value
#' @param ... params for \code{GM_writeBatchFile()}
#' @return (list). "dbDir": path to GeneMANIA database 
#' 	"netDir": path to directory with interaction networks. If profiles
#' are provided, this points to the INTERACTIONS/ subdirectory within 
#' the text-based GeneMANIA generic database
#' @export
GM_createDB <- function(netDir,patientID,outDir,simMetric="cor_pearson",
		netSfx="_cont.txt$",verbose=TRUE,numCores=1L,
		attrib_DF=NULL,...) {
	# tmpDir/ is where all the prepared files are stored.
	# GeneMANIA uses tmpDir as input to create the generic database. 
	# The database itself will be in outDir/
	tmpDir <- sprintf("%s/tmp",outDir)
	dataDir <- sprintf("%s/dataset",outDir)

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
	cmd <- sprintf("python %s batch.txt",procNet)
	system(cmd,wait=TRUE)
	
	GM_jar	<- sprintf("%s/java/GeneMANIA-3.2B7.jar",
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

		cl	<- makeCluster(numCores)
		registerDoParallel(cl)

		cmd1 <- sprintf("java -Xmx10G -cp %s org.genemania.engine.core.evaluation.ProfileToNetworkDriver", GM_jar)
		cmd3 <- "-proftype continuous -cor PEARSON"
		cmd5 <- "-threshold off -maxmissing 100.0"
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
			system(cmd)
		}
		))
		stopCluster(cl)
		netList2 <- dir(path=netOutDir,pattern=netSfx)
		netSfx=".txt"

		cat(sprintf("Got %i networks from %i profiles\n", length(netList2),
			length(netList)))

		netDir <- netOutDir
		netList <- netList2; rm(netOutDir,netList2)
	}

	#### Step X. Writing attributes
	# format for file obtained from:
	# https://github.com/GeneMANIA/pipeline/wiki/GenericDb

	# ATTRIBUTE_GROUPS.txt file - current default is to select attribute
	if (!is.null(attrib_DF)) {
	cat("Attributes provided. Processing\n")
	attrib_group <- unique(attrib_DF$attribute_group)
	attgroup_id <- 1:length(attrib_group)
	names(attgroup_id) <- attrib_group
	n <- length(attrib_group)
	tmp <- data.frame(attgroup_id, 1,
					  attrib_group,"",
					  attrib_group,"","",1,"","")
	write.table(tmp,file="ATTRIBUTE_GROUPS.txt",
				sep="\t",col=F,row=F,quote=F)
	tmp <- tmp[,c(1,3)]
	colnames(tmp)<- c("group_id","attribute_group")
	attrib_DF <- merge(x=tmp,y=attrib_DF,by="attribute_group")
	# ATTRIBUTES.txt file
	attrib_name <- attrib_DF[!duplicated(attrib_DF[,
				c("group_id","attribute_group","attribute_name")]),]
	n <- nrow(attrib_name)
	attset <- data.frame(id=1:nrow(attrib_name),org=1,
			group_id=attrib_name$group_id,aname=attrib_name$attribute_name,
			attrib_name$attribute_name,attrib_name$attribute_name)
	write.table(attset,file="ATTRIBUTES.txt",sep="\t",col=F,row=F,quote=F)

	# convert to internal GENES.txt ID before writing.
	gn <- read.delim("GENES.txt",sep="\t",h=F,as.is=T)[,c(1:2)]
	colnames(gn) <- c("GM_ID","ID")
	attrib_DF <- merge(x=attrib_DF,y=gn,by="ID")

	# ATTRIBUTES/ dir
	dir.create("ATTRIBUTES")
	for (nm in 1:nrow(attset)) {
		tmp <- subset(attrib_DF, group_id %in% attset$group_id[nm]&
					  			 attribute_name%in% attset$aname[nm])
		if (verbose) cat(sprintf("\t%s:%s: %i entries\n",attset[nm,3],
			attset[nm,4],nrow(tmp)))
		write.table(tmp[,c("GM_ID","attribute_value")],
					file=sprintf("ATTRIBUTES/%i.txt",nm),sep="\t",col=F,
					row=F,quote=F)
	}
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
	print(cmd)
	system(cmd)

	#### Step 6. Cleanup.
	if (verbose) cat("\t * Cleanup")
	GM_xml	<- sprintf("%s/java/genemania.xml",path.package("netDx"))
	system(sprintf("cp %s %s/.", GM_xml, dataDir))

	}, error=function(ex) {
		print(ex)
	}, finally={
		setwd(curwd)
	})

	return(list(dbDir=dataDir,netDir=netDir))
}
