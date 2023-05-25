#' .centerSignal
#'
#' @param ct
#'
#' @return A chiptsne2 object updated to reflect centering procecure.
#'
#' @importFrom seqsetvis centerGRangesAtMax
#'
#' @examples
#' library(tidyverse)
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' query_gr = seqsetvis::prepare_fetch_GRanges_width(query_gr, win_size = 50)
#' prof_dt = seqsetvis::CTCF_in_10a_profiles_dt
#' meta_dt = prof_dt %>%
#'   select(sample) %>%
#'   unique %>%
#'   separate(sample, c("cell", "mark"), remove = FALSE)
#' ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
#' chiptsne2:::.centerSignalGRanges(ct2)
.centerSignalGRanges = function(ct2){
    w = rowRanges(ct2) %>% GenomicRanges::width() %>% unique
    prof_dt = .prof_dt_from_chiptsne2(ct2)
    win_size = prof_dt$x %>% unique %>% diff %>% unique
    center_gr = rowRanges(ct2) %>% GenomicRanges::resize(1, fix = "center")
    prof_dt$strand = as.character(GenomicRanges::strand(center_gr[prof_dt$id]))
    prof_dt$seqnames = as.character(GenomicRanges::seqnames(center_gr[prof_dt$id]))
    prof_dt$start = GenomicRanges::start(center_gr[prof_dt$id]) + prof_dt$x - win_size/2
    prof_dt$end = GenomicRanges::start(center_gr[prof_dt$id]) + prof_dt$x + win_size/2
    seqsetvis::centerGRangesAtMax(prof_dt, rowRanges(ct2), width = w)
}

#'
#' @param ct2
#' @param view_size
#'
#' @return
#'
#' @examples
#' library(tidyverse)
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' query_gr = seqsetvis::prepare_fetch_GRanges_width(query_gr, win_size = 50)
#' prof_dt = seqsetvis::CTCF_in_10a_profiles_dt
#' meta_dt = prof_dt %>%
#'   select(sample) %>%
#'   unique %>%
#'   separate(sample, c("cell", "mark"), remove = FALSE)
#' ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
#' chiptsne2:::.centerSignalProfile(ct2, view_size = 500)
.centerSignalProfile = function(ct2, view_size){
    w = rowRanges(ct2) %>% GenomicRanges::width() %>% unique
    prof_dt = .prof_dt_from_chiptsne2(ct2)
    # win_size = prof_dt$x %>% unique %>% diff %>% unique
    # center_gr = rowRanges(ct2) %>% GenomicRanges::resize(1, fix = "center")
    # prof_dt$strand = as.character(GenomicRanges::strand(center_gr[prof_dt$id]))
    # prof_dt$seqnames = as.character(GenomicRanges::seqnames(center_gr[prof_dt$id]))
    # prof_dt$start = GenomicRanges::start(center_gr[prof_dt$id]) + prof_dt$x - win_size/2
    # prof_dt$end = GenomicRanges::start(center_gr[prof_dt$id]) + prof_dt$x + win_size/2
    new_prof_dt = seqsetvis::centerAtMax(prof_dt, trim_to_valid = TRUE, view_size = view_size, check_by_dupes = FALSE)
    rng = new_prof_dt$x %>% range
    rng_min = min(abs(rng))
    new_prof_dt = filter(new_prof_dt, x <= rng_min & x >= -rng_min)
    new_w = new_prof_dt$x %>% range %>% diff
    new_query_gr = GenomicRanges::resize(query_gr, new_w, fix = "center")
    ChIPtsne2.from_tidy(new_prof_dt, new_query_gr, sample_metadata = colData(ct2), name_VAR = "name")
}
