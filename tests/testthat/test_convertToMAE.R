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
	dataList <- list(rna=rna,patientPheno=clin)
    
    x <- convertToMAE(dataList)
    expect_is(x, "MultiAssayExperiment")
})