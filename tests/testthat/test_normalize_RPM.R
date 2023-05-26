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
meta_dt$test_mapped_reads = seq(nrow(meta_dt))*1000
map_read = meta_dt$mapped_reads
names(map_read) = meta_dt$sample

ct2.no_mr = ChIPtsne2.from_tidy(prof_dt, query_gr)
ct2.mr = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)

# debug(normalizeSignalRPM, signature = "ChIPtsne2")

test_that("normalizeSignalRPM", {
    expect_error(normalizeSignalRPM(ct2.no_mr), "mapped_reads_data was not supplied and mapped_reads_VAR")
    ct2.norm = normalizeSignalRPM(ct2.mr)
    ct2.norm2 = normalizeSignalRPM(ct2.mr, mapped_reads_VAR = "test_mapped_reads")
    m_norm2 = rowToRowMat(ct2.norm2) %>% mean
    m_norm1 = rowToRowMat(ct2.norm) %>% mean
    m_raw = rowToRowMat(ct2.mr) %>% mean
    expect_gt(m_norm1, m_raw * 1e5)
    expect_gt(m_norm1, m_norm2 * 900)

    normalizeSignalRPM(ct2.no_mr, mapped_reads_data = meta_dt)
    expect_error(normalizeSignalRPM(ct2.no_mr, mapped_reads_data = meta_dt, mapped_reads_VAR = "bad_VAR"), "mapped_reads_data does not contain mapped_reads_VAR: bad_VAR")
    normalizeSignalRPM(ct2.no_mr, mapped_reads_data = meta_dt, mapped_reads_VAR = "test_mapped_reads")


    normalizeSignalRPM(ct2.no_mr, mapped_reads_data = map_read)
    map_read.bad = map_read
    names(map_read.bad) = NULL
    expect_error(normalizeSignalRPM(ct2.no_mr, mapped_reads_data = map_read.bad), "When mapped_reads_data is supplied, names must be set.")
})
