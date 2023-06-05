testthat::context("converters")
# flipping viewGranges
library(chiptsne2)
library(testthat)

query_gr = exampleQueryGR()
prof_dt = exampleProfDT()

meta_dt = prof_dt %>%
    dplyr::select(sample) %>%
    unique %>%
    tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)

prof_dt1 = getTidyProfile(ct2)

prof_dt1.cell = getTidyProfile(ct2, sample_meta_VARS = "cell")
prof_dt1.true = getTidyProfile(ct2, sample_meta_VARS = TRUE)

test_that("Conversion", {
    expect_setequal(colnames(prof_dt1), c("id", "x", "y", "sample"))
    expect_setequal(colnames(prof_dt1.cell), c("id", "x", "y", "sample", "cell"))
    expect_setequal(colnames(prof_dt1.true), c("id", "x", "y", "sample", "cell", "mark"))
})
