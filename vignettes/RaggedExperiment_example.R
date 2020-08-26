require(curatedTCGAData)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TCGAutils)
require(rtracklayer)
require(BiocFileCache)
require(rappdirs)


tcga <- curatedTCGAData("OV","Mutation",dry.run=FALSE)
clin=tcga@colData
cnames <- colnames(clin)
smp <- sampleMap(tcga);
samps <- smp[which(!duplicated(smp$primary)),]

stage <- clin[,"patient.stage_event.clinical_stage"]

#Prepare genomic map for encoding mutations not in the current version
lifturl <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz"

message("* fetch chain file")
bfc <- BiocFileCache(ask=FALSE)
cache <- rappdirs::user_cache_dir(appname = "netDx_vignette")
    BiocFileCache::BiocFileCache(cache,ask=FALSE)
qfile <- bfcquery(bfc, "18to19chain", exact = TRUE)[["rpath"]]
cfile <-
  if (length(qfile) && file.exists(qfile)) {
    bfcquery(bfc, "18to19chain", exact = TRUE)[["rpath"]]
  } else {
    bfcadd(bfc, "18to19chain", lifturl)
  }
chainfile <- file.path(cache, gsub("\\.gz", "", basename(cfile)))
if (!file.exists(chainfile)) {
	R.utils::gunzip(cfile, destname = chainfile, remove = FALSE)
}

chain <- suppressMessages(
  rtracklayer::import.chain(chainfile)
)

message("* run liftOver")
genome(tcga[[1]]) <- translateBuild(genome(tcga[[1]]))
seqlevelsStyle(tcga[[1]]) <- c("UCSC","Ensembl")
ranges19 <- liftOver(rowRanges(tcga[[1]]), chain)

message("* Convert to SummarizedExperiment")
re19 <- tcga[[1]][as.logical(lengths(ranges19))]
ranges19 <- unlist(ranges19)
genome(ranges19) <- "hg19"
rowRanges(re19) <- ranges19
# replace with hg19 ranges
tcga[[1]] <- re19
x=tcga[, , 1]
x1=qreduceTCGA(x)

#Harmonization of the mutation matrix and clinical information
assay_df<-as.data.frame(assay(x1[[1]][,samps$colname]))
colnames(assay_df)<-samps$primary

#Build of the clinical information dataframe
info_pos<-match(colnames(assay_df),clin$patientID)
clin_df<-clin[info_pos,]


