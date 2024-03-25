.setSeed = function(ct2, seed = NULL){
    message("setSeed ...")
    if(is.null(seed)){
        seed = as.numeric(Sys.time())
    }
    args = get_args()
    set.seed(seed)
    history_item = list(setSeed = list(FUN = .setSeed, ARG = args))
    cloneChIPtsne2(ct2, new_metadata = c(ChIPtsne2.history(ct2), history_item))
}

#' setSeed
#'
#' Calls [set.seed].
#'
#' Recorded in history to allow reproducibilty with [rerurun_history]
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param seed New set.seed value. Default of NULL will use Sys.time. This will still be saved in history.
#'
#' @return `r doc_ct2_nrr()` with seed stored in histroy
#' @export
#' @rdname setSeed
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = setSeed(ct2)
#' ChIPtsne2.history(ct2)$setSeed$ARG$seed
setGeneric("setSeed", function(ct2, seed = NULL) standardGeneric("setSeed"), signature = "ct2")

#' @export
#' @rdname setSeed
setMethod("setSeed", c("ChIPtsne2_no_rowRanges"), .setSeed)
