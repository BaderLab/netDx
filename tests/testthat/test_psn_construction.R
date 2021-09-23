

test_that("invalid sims is flagged",{

    expect_equal(TRUE, checkSimValid(list(a="pearson_corr")))
    expect_identical(TRUE, checkSimValid(list(a="pearson_corr",b=function(x) 2+4)))
    expect_identical(TRUE, checkSimValid(list(a="normDiff")))
    expect_error(checkSimValid(list(a="normDifff")))
    expect_error(checkSimValid(list(a=list(a=2))))
})