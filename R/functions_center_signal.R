#' .centerSignal
#'
#' @param ct2 A ChIPtsne2 object
#'
#' @return A chiptsne2 object updated to reflect centering procedure.
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
    win_size = prof_dt[[ct2@position_VAR]] %>% unique %>% diff %>% unique
    center_gr = rowRanges(ct2) %>% GenomicRanges::resize(1, fix = "center")
    prof_dt$strand = as.character(GenomicRanges::strand(center_gr[prof_dt[[ct2@region_VAR]]]))
    prof_dt$seqnames = as.character(GenomicRanges::seqnames(center_gr[prof_dt[[ct2@region_VAR]]]))
    prof_dt$start = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] - win_size/2
    prof_dt$end = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] + win_size/2
    seqsetvis::centerGRangesAtMax(prof_dt, rowRanges(ct2), width = w)
}

#' .centerSignalProfile
#'
#' @param ct2 A ChIPtsne2 object
#' @param view_size bp range to search for max
#'
#' @return A chiptsne2 object updated to reflect centering procedure.
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
#' debug(chiptsne2:::.centerSignalProfile)
#' ct2.c = chiptsne2:::.centerSignalProfile(ct2, view_size = 500)
.centerSignalProfile = function(ct2, view_size){
    args = get_args()

    w = rowRanges(ct2) %>% GenomicRanges::width() %>% unique
    prof_dt = .prof_dt_from_chiptsne2(ct2)
    # win_size = prof_dt[[ct2@position_VAR]] %>% unique %>% diff %>% unique
    # center_gr = rowRanges(ct2) %>% GenomicRanges::resize(1, fix = "center")
    # prof_dt$strand = as.character(GenomicRanges::strand(center_gr[prof_dt[[ct2@region_VAR]]]))
    # prof_dt$seqnames = as.character(GenomicRanges::seqnames(center_gr[prof_dt[[ct2@region_VAR]]]))
    # prof_dt$start = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] - win_size/2
    # prof_dt$end = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] + win_size/2
    new_prof_dt = seqsetvis::centerAtMax(prof_dt, trim_to_valid = TRUE, view_size = view_size, check_by_dupes = FALSE, x_ = ct2@position_VAR, y_ = ct2@value_VAR)
    rng = new_prof_dt[[ct2@position_VAR]] %>% range
    rng_min = min(abs(rng))
    new_prof_dt = filter(new_prof_dt, x <= rng_min & x >= -rng_min)
    new_w = new_prof_dt[[ct2@position_VAR]] %>% range %>% diff
    new_query_gr = GenomicRanges::resize(query_gr, new_w, fix = "center")

    history_item = list(centerSignalProfile = list(FUN = .centerSignalProfile, ARG = args))
    ChIPtsne2.from_tidy(new_prof_dt,
                        new_query_gr,
                        sample_metadata = colData(ct2),
                        obj_history = c(ChIPtsne2.history(ct2), history_item),
                        name_VAR = ct2@name_VAR,
                        position_VAR = ct2@position_VAR,
                        value_VAR = ct2@value_VAR,
                        region_VAR = ct2@region_VAR
    )
}

# hist = ChIPtsne2.history(ct2.c)
# names(hist)
history_list = ct2.c

.test_history = function(x){
    if(is.list(x)){
        if(!is.null(x$FUN) & !is.null(x$ARG)){
            return(TRUE)
        }
    }
    FALSE
}

#' Title
#'
#' @param ct2
#' @param history_list
#'
#' @return
#' @export
#'
#' @examples
rerun_history = function(ct2, history_list){
    if(!is.list(history_list)){
        history_list = ChIPtsne2.history(history_list)
    }
    keep = sapply(history_list, .test_history)
    history_list = history_list[keep]
    for(hist_i in history_list){
        ct2 = do.call(hist_i$FUN, c(list(ct2 = ct2), hist_i$ARG))
    }
    ct2
}
