#' Write batch.txt file required to create GeneMANIA database
#'
#' @details Currently creates a single entry for network group. Future versions of this
#' package may want to 
#' @param netDir (char) path to dir with networks
#' @param netList (char) vector of network names
#' @param outDir (char) directory to write batch file
#' @param idFile (char) path to file with patient IDs
#' @param orgName (char) organism name. Don't change the default unless
#' you know what you are doing.
#' @param orgDesc (char) organism description. Similar to \code{orgName},
#' don't change the default
#' @param orgAlias (char) organism alias. Similar to \code{orgName}, don't
#' change the default.
#' @param taxID (integer) taxonomyID required for GeneMANIA . Similar to 
#' \code{orgName}, don't change the default.
#' @return No value. Side effect of writing batch file to 
#' \code{<outDir>/batch.txt}.
#' @export
#' @examples
#' data(npheno)
#' netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
#' netList <- dir(netDir,pattern="txt$")
#' GM_writeBatchFile(netDir,netList, "~/tmp", npheno$ID)
GM_writeBatchFile <- function(netDir,netList,outDir,idFile,
	orgName="predictor",orgDesc="my_predictor",orgAlias="my_predictor",
	taxID=1339) {
	
	outF <- sprintf("%s/batch.txt",outDir)
	if (file.exists(outF)) unlink(outF)

	# organism info
	tmp 	<- c("#organism","id","file","name","description","alias","taxonomyid")
	tmp2	<- c("organism",basename(idFile),orgName,orgDesc,orgAlias,
				 as.character(taxID))
	cat(sprintf("%s\n",paste(tmp,collapse="\t")),file=outF)
	cat(sprintf("%s\n\n",paste(tmp2,collapse="\t")),file=outF,append=TRUE)
	rm(tmp,tmp2)

	# group info
	groupName	<- "dummy_group"
	groupCode	<- "geneset_1"
	groupDesc	<- "dummy_group"
	tmp		<- c("#group","name","code","description","RRGGBB colour","organism")
	tmp2	<- c("group",groupName,groupCode,groupDesc,"ff00ff",orgName)
	cat(sprintf("%s\n",paste(tmp,collapse="\t")),file=outF,append=TRUE)
	cat(sprintf("%s\n\n",paste(tmp2,collapse="\t")),file=outF,append=TRUE)
	rm(tmp,tmp2)

	# network info - header
	tmp		<- c("#network","filename","name","description","group code")
	cat(sprintf("%s\n",paste(tmp,collapse="\t")),file=outF,append=TRUE)
	rm(tmp)

	# write networks
	net_DF	<- data.frame(type="network",
			filename=netList,name=sub(".txt","",netList),
			description=netList,groupCode=groupCode)
	write.table(net_DF,file=outF,sep="\t",
			col=FALSE,row=FALSE,quote=FALSE,append=TRUE)
}
