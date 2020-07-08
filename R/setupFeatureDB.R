#' setup database of features for feature selection
#'
#' Creates all the input files for the collection of features used in 
#' feature selection.
#' @param pheno (data.frame) patient metadata. Must contain ID column
#' @param prepDir (char) directory in which to setup database
#' @return (data.frame) internal numerical id for patients (INTERNAL_ID) and
#' user-provided ID (ID)
#' @examples
#' data(xpr,pheno)
#' pathwayList <- list(pathA=rownames(xpr)[1:10],pathB=rownames(xpr)[21:50])
#' 
#' dataList <- list(rna=xpr)  #only one layer type
#' groupList <- list(rna=pathwayList) # group genes by pathways
#' 
#' makeNets <- function(dataList, groupList, netDir,...) {
#'     netList <- makePSN_NamedMatrix(dataList[['rna']],
#'			rownames(dataList[['rna']]),
#'      groupList[['rna']],netDir,verbose=FALSE,
#' 			writeProfiles=TRUE,...)
#'     unlist(netList)
#' }
#' tmpDir <- tempdir(); netDir <- paste(tmpDir,"nets",sep=.Platform$file.sep)
#' dir.create(netDir,recursive=TRUE)
#' 
#' pheno_id <- setupFeatureDB(pheno,netDir)
#' @export
setupFeatureDB <- function(pheno, prepDir=tempdir()) {
    
    curd <- getwd()
    setwd(prepDir)
    
    pheno$INTERNAL_ID <- seq_len(nrow(pheno))
    
    # ORGANISMS.txt
    con <- file("ORGANISMS.txt", "w")
    write(paste("1", "predictor", "my_predictor", "my_predictor", -1, 1339, 
				sep = "\t"), 
        file = con, append = FALSE)
    close(con)
    
    # NODES.txt
    tmp <- pheno[, c("INTERNAL_ID", "ID", "INTERNAL_ID")]
    tmp$dummy <- 1
    write.table(tmp, file = "NODES.txt", sep = "\t", col.names = FALSE, 
				row.names = FALSE, 
        quote = FALSE)
    
    # GENES.txt
    tmp <- paste(pheno$INTERNAL_ID, pheno$ID, "N/A", 1, pheno$INTERNAL_ID, 1, 
					0, sep = "\t")
    write.table(tmp, file = "GENES.txt", sep = "\t", col.names = FALSE, 
				row.names = FALSE, 
        quote = FALSE)
    
    # GENE_DATA.txt
    write.table(pheno[, c("INTERNAL_ID", "ID")], file = "GENE_DATA.txt", 
				sep = "\t", 
        col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # synonym file
    write.table(cbind(pheno$INTERNAL_ID, pheno$INTERNAL_ID), 
				file = "1.synonyms", 
        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    # GENE_NAMING_SOURCES.txt
    con <- file("GENE_NAMING_SOURCES.txt", "w")
    write("1\tOther\t1\tOther", file = con, append = TRUE)
    write("2\tEntrez Gene ID\t0\tEntrez Gene ID", file = con, append = TRUE)
    close(con)
    
    # TAGS.txt
    file.create("TAGS.txt")
    
    # NETWORK_TAG_ASSOC.txt
    file.create("NETWORK_TAG_ASSOC.txt")
    
    # ONTOLOGIES.txt
    file.create("ONTOLOGIES.txt")
    
    # ONTOLOGY_CATEGORIES
    file.create("ONTOLOGY_CATEGORIES.txt")
    
    # colours.txt
    con <- file("colours.txt", "w")
    write("geneset_1\tff00ff", file = con, append = TRUE)
    close(con)
    
    # ATTRIBUTES.txt
    file.create("ATTRIBUTES.txt")
    
    # ATTRIBUTE_GROUPS.txt
    file.create("ATTRIBUTE_GROUPS.txt")
    
    # db.cfg
    con <- file("db.cfg", "w")
    write("[FileLocations]", file = con, append = TRUE)
    write("generic_db_dir = .", file = con, append = TRUE)
    write("[Organisms]", file = con, append = TRUE)
    write("organisms = org_1", file = con, append = TRUE)
    write("[org_1]", file = con, append = TRUE)
    write("gm_organism_id = 1", file = con, append = TRUE)
    write("short_name = predictor", file = con, append = TRUE)
    write("common_name = my_predictor", file = con, append = TRUE)
    write("", file = con, append = TRUE)
    close(con)
    
    setwd(curd)
    
    return(pheno)
}
