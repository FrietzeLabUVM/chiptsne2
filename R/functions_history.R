
.test_history = function(x){
    if(is.list(x)){
        if(!is.null(x$FUN) & !is.null(x$ARG)){
            return(TRUE)
        }
    }
    FALSE
}


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


#' rerun_history
#'
#' @param ct2  A ChIPtsne2 object to apply history too.
#' @param history_source list of history items or a ChIPtsne2 object to pull history from.
#'
#' @return Input ct2 with all functions from history_source applied.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2.c = centerProfilesAndTrim(ct2, view_size = 2e3)
#' ChIPtsne2.history(ct2.c)
#' #histories can be rerun just by supplying a source ct2 object
#' ct2.c_re1 = rerun_history(ct2, history_source = ct2.c)
#' ChIPtsne2.history(ct2.c_re1)
#' hist_src = ChIPtsne2.history(ct2.c)
#' #you can also supply the history items list, potentially modifying entries like so
#' hist_src$centerProfilesAndTrim$ARG$view_size = 1000
#' ct2.c_re2 = rerun_history(ct2, history_source = hist_src)
#' ChIPtsne2.history(ct2.c_re2)
setGeneric("rerun_history", function(ct2, history_source) standardGeneric("rerun_history"))

#' @export
setMethod("rerun_history", c("ChIPtsne2_no_rowRanges", "list"), .rerun_history)

#' @export
setMethod("rerun_history", c("ChIPtsne2_no_rowRanges", "ChIPtsne2_no_rowRanges"), .rerun_history)

.add_history_entry = function(ct2, entry_name, FUN, ARG){
    history_item = list(list(FUN = FUN, ARG = ARG))
    names(history_item) = entry_name
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)
    ct2
}
