#### centerProfilesAndRefetch ####


.centerProfilesAndRefetch = function(ct2){
    message("centerProfilesAndRefetch ...")
    args = get_args()
    if(isFetchConfigNull(ct2@fetch_config)){
        stop("FetchConfig must valid and not NULL. Use centerProfilesAndTrim or create ChIPtsne2 with ChIPtsne2.from_FetchConfig.")
    }
    w = ct2@fetch_config$view_size
    prof_dt = getTidyProfile(ct2)
    win_size = prof_dt[[ct2@position_VAR]] %>% unique %>% diff %>% unique
    center_gr = rowRanges(ct2) %>% GenomicRanges::resize(1, fix = "center")
    prof_dt$strand = as.character(GenomicRanges::strand(center_gr[prof_dt[[ct2@region_VAR]]]))
    prof_dt$seqnames = as.character(GenomicRanges::seqnames(center_gr[prof_dt[[ct2@region_VAR]]]))
    prof_dt$start = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] - win_size/2
    prof_dt$end = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] + win_size/2
    new_query_gr = seqsetvis::centerGRangesAtMax(prof_dt, rowRanges(ct2), width = w, x_ = ct2@position_VAR, y_ = ct2@value_VAR, by_ = ct2@region_VAR)

    history_item = list(centerProfilesAndRefetch = list(FUN = .centerProfilesAndRefetch, ARG = args))

    ChIPtsne2.from_FetchConfig(ct2@fetch_config,
                               new_query_gr,
                               obj_history = c(ChIPtsne2.history(ct2), history_item),
                               init = FALSE)
}

#' centerSignal
#'
#' @param ct2 A ChIPtsne2 object
#'
#' @return A chiptsne2 object updated to reflect centering procedure. Width will be the same as original but this requires a second fetch.
#'
#' @importFrom seqsetvis centerGRangesAtMax
#' @export
#'
#' @examples
#' bam_cfg_f = exampleBamConfigFile()
#' fetch_config = FetchConfig.load_config(bam_cfg_f)
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' ct2 = ChIPtsne2.from_FetchConfig(fetch_config, query_gr)
#' ct2.c = centerProfilesAndRefetch(ct2)
#' ct2.c
setGeneric("centerProfilesAndRefetch", function(ct2) standardGeneric("centerProfilesAndRefetch"))

#' @export
setMethod("centerProfilesAndRefetch", c("ChIPtsne2"), .centerProfilesAndRefetch)

#### centerProfilesAndTrim ####

.centerProfilesAndTrim = function(ct2, view_size){
    message("centerProfilesAndTrim ...")
    args = get_args()
    w = rowRanges(ct2) %>% GenomicRanges::width() %>% unique
    prof_dt = getTidyProfile(ct2)
    new_prof_dt = seqsetvis::centerAtMax(prof_dt, trim_to_valid = TRUE, view_size = view_size, check_by_dupes = FALSE, x_ = ct2@position_VAR, y_ = ct2@value_VAR, by_ = ct2@region_VAR)
    rng = new_prof_dt[[ct2@position_VAR]] %>% range
    rng_min = min(abs(rng))
    new_prof_dt = dplyr::filter(new_prof_dt, get(ct2@position_VAR) <= rng_min & get(ct2@position_VAR) >= -rng_min)
    new_w = new_prof_dt[[ct2@position_VAR]] %>% range %>% diff
    new_query_gr = GenomicRanges::resize(rowRanges(ct2), new_w, fix = "center")

    history_item = list(centerProfilesAndTrim = list(FUN = .centerProfilesAndTrim, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        prof_dt = new_prof_dt,
        query_gr = new_query_gr,
        obj_history = c(ChIPtsne2.history(ct2), history_item),
        init = FALSE
    )
}


#' centerProfilesAndTrim
#'
#' @param ct2 A ChIPtsne2 object
#' @param view_size bp range to search for max
#'
#' @return A chiptsne2 object updated to reflect centering procedure. Some x
#'   values will have been lost.
#' @importFrom GenomicRanges resize width
#' @export
#'
#' @examples
#' library(magrittr)
#' library(tidyr)
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' query_gr = seqsetvis::prepare_fetch_GRanges_width(query_gr, win_size = 50)
#' prof_dt = seqsetvis::CTCF_in_10a_profiles_dt
#' meta_dt = prof_dt %>%
#'   select(sample) %>%
#'   unique %>%
#'   separate(sample, c("cell", "mark"), remove = FALSE)
#' ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
#' ct2.c = centerProfilesAndTrim(ct2, view_size = 500)
setGeneric("centerProfilesAndTrim", function(ct2, view_size) standardGeneric("centerProfilesAndTrim"))

#' @export
setMethod("centerProfilesAndTrim", c("ChIPtsne2_no_rowRanges", "numeric"), .centerProfilesAndTrim)

