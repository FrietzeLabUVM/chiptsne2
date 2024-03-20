#### Constructor ####

#' ChIPtsne2List
#'
#' @param ... A number of ChIPtsne2_no_rowRanges or ChIPtsne2 objects.
#'
#' @export
#'
#' @importFrom S4Vectors SimpleList
ChIPtsne2List <- function(
        ...)
{
    ctl <- S4Vectors::SimpleList(...)
    list_classes = unique(sapply(ctl, class))
    if(!length(list_classes) == 1){
        stop("ChIPtsne2List must contain only either ChIPtsne2 or ChIPtsne2_no_rowRanges objects.")
    }
    if(list_classes == "ChIPtsne2"){
        out = .ChIPtsne2List(listData = ctl@listData, elementType = "ChIPtsne2", elementMetadata = ctl@elementMetadata, metadata = ctl@metadata)
    }else if(list_classes == "ChIPtsne2_no_rowRanges"){
        out = .ChIPtsne2List(listData = ctl@listData, elementType = "ChIPtsne2_no_rowRanges", elementMetadata = ctl@elementMetadata, metadata = ctl@metadata)
    }else{
        stop("ChIPtsne2List must contain only either ChIPtsne2 or ChIPtsne2_no_rowRanges objects.")
    }
    out
}

#' cbind-ChIPtsne2List
#'
#' @param ChIPtsne2List `r doc_ct2list()`
#'
#' @return a ChIPtsne2 object of concatenated columns/samples of all items in input ChIPtsne2List
#' @rdname cbind
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2.by_cell = split(ct2, "cell")
#' ct2.cbind = cbind(ct2.by_cell)
setMethod("cbind", "ChIPtsne2List", function(..., deparse.level=1) {
    args <- list(...)
    if(length(args) > 1){
        stop("Only one ChIPtsne2List may be supplied.")
    }
    do.call(cbind, as.list(args[[1]]))
})


#### rbind ####
#' rbind-ChIPtsne2List
#'
#' @param ChIPtsne2List `r doc_ct2list()`
#'
#' @return a ChIPtsne2 object of concatenated rows/regions of all items in input ChIPtsne2List
#' @rdname rbind
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2.by_peak = split(ct2, "peak_MCF10AT1_CTCF")
#' ct2.rbind = rbind(ct2.by_peak)
setMethod("rbind", "ChIPtsne2List", function(..., deparse.level=1) {
    args <- list(...)
    if(length(args) > 1){
        stop("Only one ChIPtsne2List may be supplied.")
    }
    do.call(rbind, as.list(args[[1]]))
})


#### lapply shortcuts for ChIPtsne2 ####
#' lapply-ChIPtsne2List
#'
#' @export
#'
#' @rdname lapply
#' @return A ChIPtsne2List if FUN returns a ChIPtsne2 object. Otherwise a simple list.
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2.by_peak = split(ct2, "peak_MCF10AT1_CTCF")
#' getNameVariable(ct2.by_peak)
#' ct2.by_peak = lapply(ct2.by_peak, setNameVariable, "name")
#' getNameVariable(ct2.by_peak)
setMethod("lapply", c("ChIPtsne2List"), function(X, FUN, ...){
    # when the result of lapply on a ChIPtsne2List is a list of ChIPtsne2_no_rowRanges,
    # force it to be an actual ChIPtsne2List
    res = lapply(as.list(X), FUN, ...)
    is_ct2 = sapply(res, is, "ChIPtsne2_no_rowRanges")
    if(all(is_ct2)){
        res = ChIPtsne2List(res)
    }
    res
})

#' @export
#' @rdname ct2-getset
setMethod("setNameVariable", c("ChIPtsne2List"), function(ct2, new_name_VAR){
    lapply(ct2, setNameVariable, new_name_VAR)
})

#' @export
#' @rdname ct2-getset
setMethod("swapNameVariable", c("ChIPtsne2List"), function(ct2, new_name_VAR){
    lapply(ct2, swapNameVariable, new_name_VAR)
})

#' @export
#' @rdname ct2-getset
setMethod("getNameVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getNameVariable)
})

#' @export
#' @rdname ct2-getset
setMethod("setValueVariable", c("ChIPtsne2List"), function(ct2, new_value_VAR){
    lapply(ct2, setValueVariable, new_value_VAR)
})
#' @export
#' @rdname ct2-getset
setMethod("getValueVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getValueVariable)
})

#' @export
#' @rdname ct2-getset
setMethod("setRegionVariable", c("ChIPtsne2List"), function(ct2, new_region_VAR){
    lapply(ct2, setRegionVariable, new_region_VAR)
})
#' @export
#' @rdname ct2-getset
setMethod("getRegionVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getRegionVariable)
})

#' @export
#' @rdname ct2-getset
setMethod("setPositionVariable", c("ChIPtsne2List"), function(ct2, new_position_VAR){
    lapply(ct2, setPositionVariable, new_position_VAR)
})
#' @export
#' @rdname ct2-getset
setMethod("getPositionVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getPositionVariable)
})

