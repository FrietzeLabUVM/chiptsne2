testthat::context("norm cap")
# flipping viewGranges
library(chiptsne2)
library(testthat)

query_gr = exampleQueryGR()
prof_dt = exampleProfDT()
meta_dt = prof_dt %>%
    dplyr::select(sample) %>%
    unique %>%
    tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)
meta_dt$cap_value = seq(nrow(meta_dt))*5
meta_dt$test_cap = seq(nrow(meta_dt))
cap_values = meta_dt$cap_value
names(cap_values) = meta_dt$sample

ct2.no_mr = ChIPtsne2.from_tidy(prof_dt, query_gr)
ct2.mr = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)

meta_dt2 = meta_dt
meta_dt2$cap_value = NULL
ct2.mr2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt2)

# debug(normalizeSignalCapValue, signature = "ChIPtsne2")

test_that("normalizeSignalCapValue", {
    expect_error(normalizeSignalCapValue(ct2.no_mr), "signal_cap_data was not supplied and signal_cap_VAR: cap_value is not in colData")
    ct2.norm = normalizeSignalCapValue(ct2.mr)
    ct2.norm2 = normalizeSignalCapValue(ct2.mr2, signal_cap_VAR = "test_cap")
    expect_error(normalizeSignalCapValue(ct2.mr2))
    m_norm2 = rowToRowMat(ct2.norm2) %>% mean
    m_norm1 = rowToRowMat(ct2.norm) %>% mean
    m_raw = rowToRowMat(ct2.mr) %>% mean
    expect_lt(m_norm1, m_raw * 10)
    expect_gt(m_norm2, m_norm1)


    normalizeSignalCapValue(ct2.no_mr, signal_cap_data = meta_dt, signal_cap_VAR = "cap_value")
    expect_error(
        normalizeSignalCapValue(
            ct2.no_mr,
            signal_cap_data = meta_dt,
            signal_cap_VAR = "bad_VAR"),
        "Supplied signal_cap_data does not contain bad_VAR")
    normalizeSignalCapValue(ct2.no_mr, signal_cap_data = meta_dt, signal_cap_VAR = "test_cap")


    normalizeSignalCapValue(ct2.no_mr, signal_cap_data = cap_values)
    cap_values.bad = cap_values
    names(cap_values.bad) = NULL
    expect_error(normalizeSignalCapValue(ct2.no_mr, signal_cap_data = cap_values.bad), "vector signal_cap_data must be named.")
})

test_that("calculateSignalCapValue", {
    ct2.capped = ct2.no_mr %>%
        calculateSignalCapValue %>%
        normalizeSignalCapValue
    expect_equal(rowToRowMat(ct2.capped) %>% max, 1)

    ct2.no_norm1 = ct2.no_mr %>%
        calculateSignalCapValue %>%
        normalizeSignalCapValue(norm_to_1 = FALSE)
    expect_gt(rowToRowMat(ct2.no_norm1) %>% max, 50)
    expect_lt(rowToRowMat(ct2.no_norm1) %>% max, 100)

    ct2.no_trim = ct2.no_mr %>%
        calculateSignalCapValue %>%
        normalizeSignalCapValue(trim_values_to_cap = FALSE)
    expect_gt(rowToRowMat(ct2.no_trim) %>% max, 1)
    expect_lt(rowToRowMat(ct2.no_trim) %>% max, 2)

    expect_warning(ct2.no_mr %>%
                       calculateSignalCapValue %>%
                       normalizeSignalCapValue(trim_values_to_cap = FALSE, norm_to_1 = FALSE),
                   regexp = "With do_not_cap and do_not_scaleTo1, only cap_value will be appended")
})

