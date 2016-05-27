#' Create patient similarity interaction networks based on range sets
#'
#' @details Creates patient similarity networks when data consist of 
#' genomic events associated with patients. Examples include CNV or 
#' indel data for patients. To generate networks from full matrices such
#' gene expression data, use \code{makePSN_NamedMatrix} instead.
#' Genomic ranges corresponding to events in patients (gr) should be named.
#' One network is created per named range set (rangeSet). Each set
#' reflects a group of related loci ; for example, genomic ranges associated
#' with genes in the same cellular pathway. 
#' Currently, the only similarity measure supported is binary; two patients
#' are related in a network N if they both overlap elements of set N.
#' @param gr (GRanges) patient ranges. Metadata should contain:
#'	ID: (char) unique patient ID
#'	name: (comma-separated char) named ranges overlapped
#' @param rangeSet (list) list of GRanges, one entry per range set.
#' 	Key is the name of the range set, and value is a GRanges object with
#' corresponding ranges
#' @param netDir (char) path to directory where networks should be written
#' @param simMetric (char) Similarity metric. Currently only 'coincide' 
#' is supported; two patients share an edge if they overlap elements in the
#' the same gene set. E.g. Two patients with CNVs that overlap different
#' genes of the same pathway would be related, but patients overlapping
#' genes that don't share a pathway (or, more accurately, a named-set 
#' grouping) would not be related. The edge weight is therefore binary.
#' @param verbose (logical) print detailed messages
#' @param quorum (integer) minimum number of patients in a network for the 
#' network to be constructed
#' @return Vector of network filenames
#' @export
#' @import GenomicRanges
#' @import bigmemory
#' @import foreach
#' @import parallel
makePSN_RangeSets <- function(gr, rangeSet, netDir, simMetric="coincide",
  quorum=2L,verbose=TRUE) {
if (!file.exists(netDir)) dir.create(netDir)

TEST_MODE <- FALSE # for debugging

# num patients per network
netCountFile	<- sprintf("%s/patient_count.txt",netDir)
# IDs of patients with 1+ interaction in set of networks
incPatientFile	<- sprintf("%s/inc_patients.txt", netDir)

uq_loci    	<- unique(unlist(lapply(rangeSet, function(x) { x$name })))
uq_patients	<- unique(gr$ID)

if (!simMetric %in% "coincide")
		stop("Only value supported for simMetric is 'coincide'");

# LOCUS_NAMES not provided? Compute these
if (!"LOCUS_NAMES" %in% names(elementMetadata(gr))) {
	cat("\tLOCUS_NAMES column not provided; computing overlap of patients
		with regions\n")
	gr <- getRegionOL(gr, rangeSet)
}

# set up a matrix of patient by locus.
# a[i,j] = 1 if patient i has a CNV affecting locus j
# else 0
cat("* Preparing patient-locus matrix\n")
cat(sprintf("\t%i unique patients, %i unique locus symbols\n",
            length(uq_patients), length(uq_loci)))
pgMat   <- big.matrix(0,
		 nrow=length(uq_patients), ncol=length(uq_loci),
         type="integer")
pgDesc  <- describe(pgMat)
cat("\n")

x <- foreach (k=1:length(uq_patients)) %do% {
    inner_mat 	<- attach.big.matrix(pgDesc)
	idx			<- which(gr$ID %in% uq_patients[k])
    myloci		<- unlist(strsplit(gr$LOCUS_NAMES[idx],","))
	myloci		<- setdiff(myloci,"")

	if (length(myloci)>0) {
        inner_mat[k, which(uq_loci %in% myloci)] <- 1L
    }
    if (k %% 100==0) cat(".")
}

hit_p 		<- integer(length(rangeSet))

# binary vector indicating if this network was included
inc_status	<-  integer(length(rangeSet))
# patients that have at least one interaction in one network
inc_patients <- integer(length(uq_patients));
names(inc_patients) <- uq_patients

# now group set-by-set
outFiles	<- character()
for (idx in 1:length(rangeSet)) {
	curP 		<- names(rangeSet)[idx]
    if (verbose) cat(sprintf("\t%s: ", curP))

    locus_idx 	<- which(uq_loci %in% rangeSet[[idx]]$name)
	if (length(locus_idx)>=2) {
    	hit_pathway <- rowSums(pgMat[,locus_idx])
	} else { # one-locus pathway, probably risk locus
		hit_pathway	<- pgMat[,locus_idx]
	}

    hit_p[idx]	<- sum(hit_pathway>0)
    if (verbose) cat(sprintf("%i patients with interactions",
		hit_p[idx]))

	pScore	<- 1	# similarity score for default binary option
	# pathway included in analysis
    if (hit_p[idx]>=quorum) {
		if (verbose) cat(sprintf("\n\t\tlength=%i; score = %1.2f",
								 length(rangeSet[[idx]]), pScore))
		if (!TEST_MODE) {
		inc_patients[hit_pathway>0] <-  inc_patients[hit_pathway>0]+1;
		pat_pairs <- t(combn(uq_patients[hit_pathway>0],2))
		pat_pairs <- cbind(pat_pairs,pScore);

		# write network for pathway
		outFile		<- sprintf("%s/%s_cont.txt",netDir,curP)
		write.table(pat_pairs, file=outFile,sep="\t",
					col=FALSE,row=FALSE,quote=FALSE)
		outFiles 	<- c(outFiles, basename(outFile))
		}
		status	<- 1;
    } else {
		status	<- 0;
	}
    if (idx %% 100==0) cat(".")
    if (verbose) cat("\n")
	inc_status[idx] <- status
	
	if (!verbose) {
		if (idx %% 100 == 0) cat(".")
		if (idx %% 1000 == 0) cat("\n")
	}
}

outFiles

}
