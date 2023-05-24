#' .centerSignal
#'
#' @param ct
#'
#' @return A chiptsne2 object updated to reflect centering procecure.
#'
#' @examples
#' library(tidyverse)
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' query_gr = prepare_fetch_GRanges_width(query_gr, win_size = 50)
#' prof_dt = seqsetvis::CTCF_in_10a_profiles_dt
#' meta_dt = prof_dt %>%
#'   select(sample) %>%
#'   unique %>%
#'   separate(sample, c("cell", "mark"), remove = FALSE)
#' ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_meta_dt = meta_dt)

#'
.centerSignal = function(ct2){
    w = rowRanges(ct2) %>% width %>% unique
    centerGRangesAtMax()
}
