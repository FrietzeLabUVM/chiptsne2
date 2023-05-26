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
meta_dt$mapped_reads = seq(nrow(meta_dt))

ct2.no_mr = ChIPtsne2.from_tidy(prof_dt, query_gr)
ct2.mr = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)


test_that("normalizeSignalRPM", {
    expect_error(normalizeSignalRPM(ct2.no_mr), "mapped_reads_data was not supplied and mapped_reads_VAR")
    normalizeSignalRPM(ct2.mr)
})
