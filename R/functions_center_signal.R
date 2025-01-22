#### centerProfilesAndRefetch ####


.centerProfilesAndRefetch = function(ct2, view_size = NULL, use_cache = TRUE){
    #visible binding NOTE
    `.` = `:=` = original_center = start = end = strand_multiplier = strand = NULL
    message("centerProfilesAndRefetch ...")
    args = get_args()
    if(isFetchConfigNull(ct2@fetch_config)){
        stop("FetchConfig must valid and not NULL. Use centerProfilesAndTrim or create ChIPtsne2 with ChIPtsne2.from_FetchConfig.")
    }

    prof_dt = getTidyProfile(ct2)
    win_size = prof_dt[[ct2@position_VAR]] %>% unique %>% diff %>% round(., digits = 5) %>% unique

    center_gr = rowRanges(ct2) %>% GenomicRanges::resize(., 1, fix = "center")
    prof_dt$strand = as.character(GenomicRanges::strand(center_gr[prof_dt[[ct2@region_VAR]]]))
    prof_dt$seqnames = as.character(GenomicRanges::seqnames(center_gr[prof_dt[[ct2@region_VAR]]]))
    #win_size less than 1 is appropriate if data is fetched as summaries of regions instead of sample
    #a win_size of 1 would not be appropriate for sample
    #this handles ifelse block handles the summary type (position is fraction of region) and sample type (position is bp relative to region)
    if(win_size <= 1){
        w_dt = as.data.table(as.data.frame(rowRanges(ct2)))
        data.table::set(w_dt, j = ct2@region_VAR, value = rownames(ct2))
        w_dt[, original_center := round((start + end) / 2)]
        w_dt = w_dt[, c(ct2@region_VAR, "width", "original_center"), with = FALSE]
        prof_dt = merge(prof_dt, w_dt, by = ct2@region_VAR)
        if(any(GenomicRanges::strand(rowRanges(ct2)) == "-")){
            prof_dt[, strand_multiplier := ifelse(strand == "-", -1, 1)]
            prof_dt[, start := original_center + strand_multiplier*((width * get(ct2@position_VAR)) - (win_size / 2 * width))]
            prof_dt[, end := original_center + strand_multiplier*((width * get(ct2@position_VAR)) + (win_size / 2 * width))-1]
        }else{
            prof_dt[, start := original_center + ((width * get(ct2@position_VAR)) - (win_size / 2 * width))]
            prof_dt[, end := original_center + ((width * get(ct2@position_VAR)) + (win_size / 2 * width))-1]
        }
    }else{
        if(any(GenomicRanges::strand(rowRanges(ct2)) == "-")){
            prof_dt[, strand_multiplier := ifelse(strand == "-", -1, 1)]
            prof_dt[, start := GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + strand_multiplier*prof_dt[[ct2@position_VAR]] - win_size/2]
            # prof_dt$start = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + strand_multiplier*prof_dt[[ct2@position_VAR]] - win_size/2
            prof_dt[, end := GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + strand_multiplier*prof_dt[[ct2@position_VAR]] + win_size/2]
            # prof_dt$end = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + strand_multiplier*prof_dt[[ct2@position_VAR]] + win_size/2
        }else{
            prof_dt[, start := GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] - win_size/2]
            # prof_dt$start = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] - win_size/2
            prof_dt[, end := GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] + win_size/2]
            # prof_dt$end = GenomicRanges::start(center_gr[prof_dt[[ct2@region_VAR]]]) + prof_dt[[ct2@position_VAR]] + win_size/2
        }

    }
    new_query_gr = seqsetvis::centerGRangesAtMax(prof_dt,
                                                 rowRanges(ct2),
                                                 width = ct2@fetch_config$view_size,
                                                 x_ = ct2@position_VAR,
                                                 y_ = ct2@value_VAR,
                                                 by_ = ct2@region_VAR,
                                                 view_size = view_size)

    history_item = list(centerProfilesAndRefetch = list(FUN = .centerProfilesAndRefetch, ARG = args))

    ChIPtsne2.from_FetchConfig(ct2@fetch_config,
                               new_query_gr,
                               obj_history = c(ChIPtsne2.history(ct2), history_item),
                               init = FALSE,
                               use_cache = use_cache)
}

#' centerSignal
#'
#' @param ct2 A ChIPtsne2 object
#' @param view_size bp range to search for max
#' @param use_cache  If TRUE, default [BiocFileCache::BiocFileCache] will be used. If FALSE, no caching will be done. You may also supply a user created [BiocFileCache::BiocFileCache].
#'
#' @return A chiptsne2 object updated to reflect centering procedure. Width will be the same as original but this requires a second fetch.
#'
#' @importFrom seqsetvis centerGRangesAtMax
#' @export
#' @rdname ct2-center-refetch
#'
#' @examples
#' bam_cfg_f = exampleBamConfigFile()
#' fetch_config = FetchConfig.load_config(bam_cfg_f)
#' query_gr = exampleQueryGR()
#' ct2 = ChIPtsne2.from_FetchConfig(fetch_config, query_gr)
#' ct2.c = centerProfilesAndRefetch(ct2)
#' ct2.c
setGeneric("centerProfilesAndRefetch", function(ct2, view_size = NULL, use_cache = TRUE) standardGeneric("centerProfilesAndRefetch"))

#' @export
#' @rdname ct2-center-refetch
setMethod("centerProfilesAndRefetch", c("ChIPtsne2"), .centerProfilesAndRefetch)

#### centerProfilesAndTrim ####

.centerProfilesAndTrim = function(ct2, view_size){
    message("centerProfilesAndTrim ...")
    args = get_args()
    prof_dt = getTidyProfile(ct2)
    new_prof_dt = seqsetvis::centerAtMax(prof_dt, trim_to_valid = TRUE, view_size = view_size, check_by_dupes = FALSE, x_ = ct2@position_VAR, y_ = ct2@value_VAR, by_ = ct2@region_VAR)
    rng = new_prof_dt[[ct2@position_VAR]] %>% range
    rng_min = min(abs(rng))
    new_prof_dt = dplyr::filter(new_prof_dt, get(ct2@position_VAR) <= rng_min & get(ct2@position_VAR) >= -rng_min)
    new_w = new_prof_dt[[ct2@position_VAR]] %>% range %>% diff
    new_rowRanges = GenomicRanges::resize(rowRanges(ct2), new_w, fix = "center")

    history_item = list(centerProfilesAndTrim = list(FUN = .centerProfilesAndTrim, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        new_prof_dt = new_prof_dt,
        new_rowRanges = new_rowRanges,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
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
#' @rdname ct2-center-trim
#'
#' @examples
#' library(magrittr)
#' query_gr = exampleQueryGR()
#' query_gr = seqsetvis::prepare_fetch_GRanges_width(query_gr, win_size = 50)
#' prof_dt = exampleProfDT()
#' meta_dt = prof_dt %>%
#'   dplyr::select(name) %>%
#'   unique %>%
#'   tidyr::separate(name, c("cell", "mark"), remove = FALSE)
#' ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
#' ct2.c = centerProfilesAndTrim(ct2, view_size = 500)
setGeneric("centerProfilesAndTrim", function(ct2, view_size) standardGeneric("centerProfilesAndTrim"))

#' @export
#' @rdname ct2-center-trim
setMethod("centerProfilesAndTrim", c("ChIPtsne2_no_rowRanges", "numeric"), .centerProfilesAndTrim)

