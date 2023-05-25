testthat::context("converters")
# flipping viewGranges
library(chiptsne2)
library(testthat)

query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
prof_dt = seqsetvis::CTCF_in_10a_profiles_dt
meta_dt = prof_dt %>%
    dplyr::select(sample) %>%
    unique %>%
    tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)

colData(ct2)
rowRanges(ct2)
rowData(ct2)

.prof_dt_from_chiptsne2 = chiptsne2:::.prof_dt_from_chiptsne2
prof_dt1 = .prof_dt_from_chiptsne2(ct2)
prof_dt1.cell = .prof_dt_from_chiptsne2(ct2, sample_meta_VARS = "cell")
prof_dt1.true = .prof_dt_from_chiptsne2(ct2, sample_meta_VARS = TRUE)

test_that("Constructors - invalid", {
    expect_setequal(colnames(prof_dt1), c("id", "x", "y", "name"))
    expect_setequal(colnames(prof_dt1.cell), c("id", "x", "y", "name", "cell"))
    expect_setequal(colnames(prof_dt1.true), c("id", "x", "y", "name", "cell", "mark"))
})
