.setSeed = function(ct2, seed = NULL){
    if(is.null(seed)){
        seed = as.numeric(Sys.time())
    }
    args = get_args()
    set.seed(seed)
    history_item = list(setSeed = list(FUN = .setSeed, ARG = args))
    cloneChIPtsne2(ct2, .metadata = c(ChIPtsne2.history(ct2), history_item))
}

#' @export
setGeneric("setSeed", function(ct2, seed = NULL) standardGeneric("setSeed"), signature = "ct2")

#' @export
setMethod("setSeed", c("ChIPtsne2"), .setSeed)
