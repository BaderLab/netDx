runPredictor_nestedCV <- function(pheno,dataList,groupList,outDir,makeNetFunc,
	nFoldCV=10L,trainProp=0.8,numSplits=10L,numCores,CVmemory=4L,CVcutoff=9L,
	keepAllData=FALSE) {
    ### tests# pheno$ID and $status must exist
    if (missing(dataList)) stop("dataList must be supplied.\n")
    if (missing(groupList)) stop("groupList must be supplied.\n")
    if (trainProp <= 0 | trainProp >= 1)
    		stop("trainProp must be greater than 0 and less than 1")

    megaDir <- outDir
    if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
    dir.create(megaDir)

    # set aside for testing within each split
    pheno_all <- pheno;

    # run featsel once per subtype
    subtypes <- unique(pheno$STATUS)

    cat(sprintf("-------------------------------\n"))
    cat(sprintf("# patients = %i\n", nrow(pheno)))
    cat(sprintf("# classes = %i { %s }\n", length(subtypes),
    	paste(subtypes,collapse=",")))
    cat("Sample breakdown by class\n")
    print(table(pheno$STATUS))
    cat(sprintf("Nested CV design = %i CV x %i splits\n", nFoldCV, numSplits))
    cat(sprintf("Datapoints:\n"))
    for (nm in names(dataList)) {
    	cat(sprintf("\t%s: %i units\n", nm, nrow(dataList[[nm]])))
    }

    # create master list of possible networks
    cat("# input nets provided:\n")
    netFile <- sprintf("%s/inputNets.txt", megaDir)
    cat("NetType\tNetName\n",file=netFile)
    for (nm in names(groupList)) {
    	curNames <- names(groupList[[nm]])
    	for (nm2 in curNames) {
    		cat(sprintf("%s\t%s\n",nm,nm2),file=netFile,append=TRUE)
    	}
    }

    cat("\n\nCustom function to generate input nets:\n")
    print(makeNetFunc)
    cat(sprintf("-------------------------------\n\n"))

    for (rngNum in 1:numSplits) {
    	cat(sprintf("-------------------------------\n"))
    	cat(sprintf("RNG seed = %i\n", rngNum))
    	cat(sprintf("-------------------------------\n"))
    	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
    	dir.create(outDir)

      pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
    											  setSeed=rngNum*5)
    	pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")

    	if(typeof(dataList[[1]]) == 'S4'){
    		temp_list <- list()
    		for(cur_dat in names(dataList)){
    			dats_train_temp <- dataList[[cur_dat]][mcols(dataList[[cur_dat]])$ID %in% pheno$ID]
    			temp_list[[cur_dat]] <- dats_train_temp
    		}
    		dats_train <- temp_list
    	}else{
    		dats_train <- lapply(dats,function(x) {
    							 x[,which(colnames(x) %in% pheno$ID)]})
    		}

    	netDir <- sprintf("%s/networks",outDir)
    	netList <- createPSN_MultiData(dataList=dats_train,groupList=groupList,
    			netDir=netDir,customFunc=makeNetFunc,numCores=numCores)
    	if(typeof(dataList[[1]]) == 'S4'){
    		p <- countPatientsInNet(netDir,netList, pheno$ID)
    		tmp	<- updateNets(p,pheno,writeNewNets=FALSE)
    		pheno	<- tmp[[2]]
        p_FULL <- p; pheno_FULL <- pheno;
    	}

      Nway_netSum_nestedCV(p_FULL,pheno_FULL,predictClass,outDir,netDir,
      			cliqueReps=2500L,numCores=8L)


      # run query for this class
  		qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS%in%"TRAIN")]
  		qFile <- sprintf("%s/%s_query",pDir,g)
  		GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
  		resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
  		predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)



    }






  }
