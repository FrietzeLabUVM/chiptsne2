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

ChIPtsne2.from_tidy(prof_dt, query_gr) %>% colData

prof_dt2 = copy(prof_dt)
prof_dt2 = prof_dt2 %>% tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

ChIPtsne2.from_tidy(prof_dt2, query_gr) %>% colData
ChIPtsne2.from_tidy(prof_dt2, query_gr, auto_sample_metadata = FALSE) %>% colData

test_that("Constructors - valid", {
})

