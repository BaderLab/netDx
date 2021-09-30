#' returns valid built-in similarity functions
allowedSims <- function(){
  c("pearson_corr","normDiff","avgNormDiff",
        "sim.pearscale","sim.eucscale")
}

#' checks if provided similarity functions are valid. Returns error if not
#'
#' @param sims (list) keys are layer names, values are functions or characters 
#' (names of built-in similarity functions)
#' @return TRUE if all pass check. Else throws error.
checkSimValid <- function(sims){
    allowed <- allowedSims
    for (k in names(sims)){
        if (class(sims[[k]])!="function"){
            if class(sims[[k]])!="character"{
                stop(paste("Invalid sims datatype. ",
                    "sims entries must be functions or keywords (characters) ",
                    "for built-in similarity functions.",sep=""))
            } else {
                if (!sims[[k]] %in% allowed){
                    stop(paste(
                            sprintf("sims[[%s]] has invalid similarity type:",k),
                            sims[[k]],". ",
                            "Allowed values are: {%s}",
                            paste(allowed,collapse=",")))
                }
            }
        }
    }
    return(TRUE)
}

makeNetFunc <- function(dataList, groupList, netDir, sims,...){    
    settings <- list(dataList=dataList,groupList=groupList,
                    netDir=netDir,sims=sims)
    netList <- c()    
    for (nm in names(sims)){
        csim <- sims[[nm]]
        netList_cur <- NULL

        cur_set <- settings; 
        cur_set[["name"]] <- nm; cur_set[["similarity"]] <- csim

        if (!is.null(groupList[[nm]])){
            if (class(csim)=="function") {# custom function
                netList_cur <- builtInPSN(cur_set,csim,...)
            } else if (csim == "pearson_corr") {
                netList_cur <- corrPSN(cur_set,...)
            } else {
                netList_cur <- builtInPSN(cur_set,...)
            }
            netList <- c(netList,netList_cur)
        }
    }
    unlist(netList)
}

#' make PSN for built-in similarity functions
#' @param settings (list) from makeNetFunc
builtInPSN <- function(settings,...){
funcs <- list(
    "normDiff"=normDiff,
    "avgNormDiff"=avgNormDiff,
    "sim.pearscale"=sim.pearscale,
    "sim.eucscale"=sim.eucscale
)

    message(sprintf("Layer %s: Function %s",settings$name,settings$similarity))

    nm <- settings$name
    netList <- makePSN_NamedMatrix(
        settings$dataList[[nm]],
		rownames(settings$dataList[[nm]]),
		settings$groupList[[nm]],
        settings$netDir,
		simMetric="custom",
        customFunc=funcs[[settings$similarity]], # custom function
		writeProfiles=FALSE,
		sparsify=TRUE,...
    )
    netList
}

#' make PSN for custom similarity functions
#' @param settings (list) from makeNetFunc
customPSN <- function(settings,fn, ...){
    nm <- settings$name
    message(sprintf("Layer %s: CUSTOM FUNCTION",settings$name))
    netList <- makePSN_NamedMatrix(
        settings$dataList[[nm]],
		rownames(settings$dataList[[nm]]),
		settings$groupList[[nm]],
        settings$netDir,
		simMetric="custom",customFunc=fn, # custom function
		writeProfiles=FALSE,
		sparsify=TRUE,...
    )
    netList
}

#' wrapper for PSNs using Pearson correlation
#' @param settings (list) from makeNetFunc
corrPSN <- function(settings,...){
    message(sprintf("Layer %s: PEARSON CORR",settings$name))
    nm <- settings$name
    netList <- makePSN_NamedMatrix(
				settings$dataList,
				rownames(settings$dataList[[nm]]),	## names of measures (e.g. genes, CpGs)
				settings$groupList[[nm]],			## how to group measures in that layer
				settings$netDir,						## leave this as-is, netDx will figure out where this is.
				verbose=FALSE, 			
				writeProfiles=TRUE,   		## use Pearson correlation-based similarity
				...
				)
    return(netList)
}