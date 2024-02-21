testthat::context("history")
# flipping viewGranges
library(chiptsne2)
library(testthat)

query_gr = exampleQueryGR()
query_gr = seqsetvis::prepare_fetch_GRanges_width(query_gr, win_size = 50)
prof_dt = exampleProfDT()
meta_dt = prof_dt %>%
    dplyr::select(sample) %>%
    unique %>%
    tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
ct2.dupe = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
#birthday is to the second so we need to override
ct2.dupe@metadata$birthday = "test_value"

# undebug(centerProfilesAndTrim, signature = c("ChIPtsne2", "numeric"))
# debug(getTidyProfile, signature = c("ChIPtsne2"))
ct2.c = centerProfilesAndTrim(ct2, view_size = 200)
ct2.c = centerProfilesAndTrim(ct2.c, view_size = 100)
ct2.c = centerProfilesAndTrim(ct2.c, view_size = 50)

ChIPtsne2.history(ct2.c)
ct2.rerun = rerun_history(ct2.dupe, ct2.c)

hist.1 = ChIPtsne2.history(ct2)
hist.c = ChIPtsne2.history(ct2.c)
hist.rerun = ChIPtsne2.history(ct2.rerun)

test_that("History names", {
    expect_setequal(names(hist.1), c("birthday", "session_info", "chiptsne2_version"))
    expect_setequal(names(hist.c), c("birthday", "session_info", "chiptsne2_version", "centerProfilesAndTrim", "centerProfilesAndTrim", "centerProfilesAndTrim"))
    expect_setequal(names(hist.rerun), c("birthday", "session_info", "chiptsne2_version", "centerProfilesAndTrim", "centerProfilesAndTrim", "centerProfilesAndTrim"))
})

test_that("ct2 rowRanges widths", {
    expect_equal(width((ct2)) %>% unique, 700)
    expect_equal(width((ct2.c)) %>% unique, 300)
    expect_equal(width((ct2.rerun)) %>% unique, 300)
})


test_that("ct2 profile result", {
    original_dig = rowToRowMat(ct2) %>% digest::digest()
    expect_failure(expect_equal(
        rowToRowMat(ct2.c) %>% digest::digest(), original_dig
    ))
    expect_failure(expect_equal(
        rowToRowMat(ct2.rerun) %>% digest::digest(), original_dig
    ))
    #result from rerunning history should equal orginal result
    expect_equal(
        rowToRowMat(ct2.c) %>% digest::digest(), rowToRowMat(ct2.rerun) %>% digest::digest()
    )
})

test_that("ct2 birthday", {
    #birthday is the same for prerun and postrun
    expect_equal(
        hist.c$birthday, hist.1$birthday
    )
    #birthday is different for prerun and rerun
    expect_failure(expect_equal(
        hist.1$birthday, hist.rerun$birthday
    ))

})
