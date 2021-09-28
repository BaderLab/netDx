# test convertToMAE.R

test_that("convertToMAE works", {
    # 20 patients, 10 case, 10 control
	pheno <- data.frame(ID=sprintf("PAT%i",1:20),
		STATUS=rep(c("case","control"),each=10))
	# 100 dummy genes
	rna <- matrix(rnorm(100*20),nrow=100); 
	colnames(rna) <- pheno$ID
	rownames(rna) <- sprintf("gene%i",1:100)	
	# 2 dummy clin variables
	clin <- t(data.frame(AGE=runif(20,10,50)))
	colnames(clin) <- pheno$ID
    clin <- t(clin)

	# netDx files
	dataList <- list(rna=rna,pheno=clin)
    
    x <- convertToMAE(dataList)
    expect_is(x, "MultiAssayExperiment")
})

test_that("convertToMAE works with more than one assay", {
	# 20 patients, 10 case, 10 control
	pheno <- data.frame(ID=sprintf("PAT%i",1:20),
		STATUS=rep(c("case","control"),each=10))
	# 100 dummy genes
	rna <- matrix(rnorm(100*20),nrow=100); 
	colnames(rna) <- pheno$ID
	rownames(rna) <- sprintf("gene%i",1:100)
	# 200 dummy proteins
	prot <- matrix(rnorm(200*20), nrow = 200);
	colnames(prot) <- pheno$ID
	rownames(prot) <- sprintf("protein%i",1:200) 	
	# 2 dummy clin variables
	clin <- t(data.frame(AGE=runif(20,10,50)))
	colnames(clin) <- pheno$ID
    clin <- t(clin)

	# netDx files
	dataList <- list(rna = rna, proteomics = prot, pheno = clin)

	x <- convertToMAE(dataList)
	expect_is(x, "MultiAssayExperiment")
})

test_that ("convertToMAE removes duplicated sample", {
	# 20 patients, 10 case, 10 control
	pheno <- data.frame(ID=sprintf("PAT%i",1:20),
                    STATUS=rep(c("case","control"),each=10))
	# 100 dummy genes, with first sample duplicated
	rna <- matrix(rnorm(100*20),nrow=100); 
	colnames(rna) <- pheno$ID
	rownames(rna) <- sprintf("gene%i",1:100)
	rna <- cbind(rna, rna[,1])
	colnames(rna)[21] <- colnames(rna)[1]
	# 200 dummy proteins
	prot <- matrix(rnorm(200*20), nrow = 200);
	colnames(prot) <- pheno$ID
	rownames(prot) <- sprintf("protein%i",1:200) 	
	# 2 dummy clin variables
	clin <- t(data.frame(AGE=runif(20,10,50)))
	colnames(clin) <- pheno$ID
	clin <- t(clin)
	
	# netDx files
	dataList <- list(rna = rna, proteomics = prot, pheno = clin)

	x <- convertToMAE(dataList)
	expect_is(x, "MultiAssayExperiment")
	# number of samples in rna assay vs colData should differ by 1
	expect_equal((dim(rna)[2] - dim(colData(x))[1]), 1)
	# number of samples in metadata should agree with colData
	expect_equal((dim(clin)[1] - dim(colData(x))[1]), 0)
})