testthat::context("ChIPtsne_class")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(SummarizedExperiment)
library(data.table)

query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
prof_dt = seqsetvis::CTCF_in_10a_profiles_dt


meta_dt = prof_dt %>% dplyr::select(sample) %>% unique
meta_dt = meta_dt %>% tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
colData(ct2)

ct2.no_meta = ChIPtsne2.from_tidy(prof_dt, query_gr)

prof_dt2 = copy(prof_dt)
prof_dt2 = prof_dt2 %>% tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)
prof_dt2$extra = 1

ct2.auto = ChIPtsne2.from_tidy(prof_dt2, query_gr)
ct2.no_auto = ChIPtsne2.from_tidy(prof_dt2, query_gr, auto_sample_metadata = FALSE)

test_that("Meta", {
    expect_setequal(colnames(colData(ct2)), c("cell", "mark"))
    expect_setequal(colnames(colData(ct2.no_meta)), character())
    expect_setequal(colnames(colData(ct2.auto)), c("cell", "mark", "extra"))
    expect_setequal(colnames(colData(ct2.no_auto)), character())
})

