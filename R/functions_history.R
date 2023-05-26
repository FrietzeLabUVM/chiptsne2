
.test_history = function(x){
    if(is.list(x)){
        if(!is.null(x$FUN) & !is.null(x$ARG)){
            return(TRUE)
        }
    }
    FALSE
}

#' rerun_history
#'
#' @param ct2  A ChIPtsne2 object to apply history too.
#' @param history_source list of history items or a ChIPtsne2 object to pull history from.
#'
#' @return Input ct2 with all functions from history_source applied.
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
#' ct2.c = chiptsne2:::.centerSignalProfile(ct2, view_size = 500)
#' rerun_history(ct2, ct2.c)
#' rerun_history(ct2, ChIPtsne2.history(ct2.c))
.rerun_history = function(ct2, history_source){
    if(!is.list(history_source)){
        history_source = ChIPtsne2.history(history_source)
    }
    keep = sapply(history_source, .test_history)
    history_source = history_source[keep]
    for(hist_i in history_source){
        ct2 = do.call(hist_i$FUN, c(list(ct2 = ct2), hist_i$ARG))
    }
    ct2
}

#' @export
setGeneric("rerun_history", function(ct2, history_source) standardGeneric("rerun_history"))

#' @export
setMethod("rerun_history", c("ChIPtsne2", "list"), .rerun_history)

#' @export
setMethod("rerun_history", c("ChIPtsne2", "ChIPtsne2"), .rerun_history)

