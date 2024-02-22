testthat::context("colnames rownames dimnames")
library(chiptsne2)
library(testthat)

ct2 = exampleChIPtsne2.with_meta()
exp_row_cn = c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF",  "peak_MCF10CA1_CTCF")
exp_col_cn = c("sample", "cell",  "mark")
exp_min_cn = c("id", "x", "y", "sample")
test_that("getTidyProfile basic", {
    prof_dt = getTidyProfile(ct2)
    expect_setequal(colnames(prof_dt), exp_min_cn)
    expect_equal(nrow(prof_dt), 4200)
})

test_that("getTidyProfile all", {
    prof_dt = getTidyProfile(ct2, meta_VARS = TRUE)
    expect_setequal(colnames(prof_dt), c(exp_row_cn, exp_col_cn, exp_min_cn))
    expect_equal(nrow(prof_dt), 4200)
})

test_that("getTidyProfile peak_MCF10A_CTCF", {
    prof_dt = getTidyProfile(ct2, meta_VARS = "peak_MCF10A_CTCF")
    expect_setequal(colnames(prof_dt), c("peak_MCF10A_CTCF", exp_min_cn))
    expect_equal(nrow(prof_dt), 4200)
})

test_that("getTidyProfile bad", {
    expect_error({
        prof_dt = getTidyProfile(ct2, meta_VARS = "not good")
    }, regexp = "Invalid meta_VARS specified.")

})
