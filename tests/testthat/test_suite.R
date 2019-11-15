# test utilities

test_that("readPathways works", {
	  x <- readPathways(getExamplePathways(),MIN_SIZE=10L, MAX_SIZE=200L)
		ln <- unlist(lapply(x,length))
    expect_that(x,is_a("list"))
	  expect_that(x[[1]],is_a("character"))
		expect_gt(min(ln),9)
		expect_lt(max(ln),201)
})


###test_that("lasso filtering works", {
###	# make own subroutine
###})
###
###test_that("imputation works", {
###	# make own subroutine
###})
###
#### ------------------------------------------
#### feature creation
#### ------------------------------------------
###
#### similarity methods
###test_that("similarity works: normDiff", {
###})
###test_that("similarity works: AvgNormDiff", {
###})
###
###test_that("similarity works: euc + exp scaling", {
###})
###
###test_that("similarity works: Pearson", {
###})
###
#### sparsification methods
###test_that("sparsification works: sparsify2", {
###})
###
###test_that("sparsification works: sparsify3", {
###})
###
