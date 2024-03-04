#### Constructor ####

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
        out = new("ChIPtsne2List", listData = ctl@listData, elementType = "ChIPtsne2", elementMetadata = ctl@elementMetadata, metadata = ctl@metadata)
    }else if(list_classes == "ChIPtsne2_no_rowRanges"){
        out = new("ChIPtsne2List", listData = ctl@listData, elementType = "ChIPtsne2_no_rowRanges", elementMetadata = ctl@elementMetadata, metadata = ctl@metadata)
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
setMethod("setNameVariable", c("ChIPtsne2List"), function(ct2, new_name_VAR){
    lapply(ct2, setNameVariable, new_name_VAR)
})

#' @export
setMethod("swapNameVariable", c("ChIPtsne2List"), function(ct2, new_name_VAR){
    lapply(ct2, swapNameVariable, new_name_VAR)
})

#' @export
setMethod("getNameVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getNameVariable)
})

#' @export
setMethod("setValueVariable", c("ChIPtsne2List"), function(ct2, new_value_VAR){
    lapply(ct2, setValueVariable, new_value_VAR)
})
#' @export
setMethod("getValueVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getValueVariable)
})

#' @export
setMethod("setRegionVariable", c("ChIPtsne2List"), function(ct2, new_region_VAR){
    lapply(ct2, setRegionVariable, new_region_VAR)
})
#' @export
setMethod("getRegionVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getRegionVariable)
})

#' @export
setMethod("setPositionVariable", c("ChIPtsne2List"), function(ct2, new_position_VAR){
    lapply(ct2, setPositionVariable, new_position_VAR)
})
#' @export
setMethod("getPositionVariable", c("ChIPtsne2List"), function(ct2){
    sapply(ct2, getPositionVariable)
})

#
# ct2.l = split(ct2, "cell")
# res = setNameVariable(ct2.l, "mark")
# is(class(res), "ChIPtsne2List")
# lapply(ct2.l, colnames)
# lapply(ct2.l, setNameVariable, "mark")
#
# mutateCol = function(.data, ...){
#     colData(.data) = S4Vectors::DataFrame(dplyr::mutate(as.data.frame(colData(.data)), ...))
#     .data
# }
# ct2.mut = mutateCol(ct2, group = paste(cell, mark, sep = "\n"))
# colData(ct2.mut)
#
# mutateRow = function(.data, ...){
#     rowData(.data) = S4Vectors::DataFrame(dplyr::mutate(as.data.frame(rowData(.data)), ...))
#     .data
# }
#
# ct2 = mutateRow(ct2, new_col = peak_MCF10A_CTCF | peak_MCF10AT1_CTCF)
# rowData(ct2)
# separateCol = function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
#                         convert = FALSE, extra = "warn", fill = "warn", ...){
#     colData(data) = S4Vectors::DataFrame(tidyr::separate(
#         as.data.frame(colData(data)),
#         col = col,
#         into = into,
#         sep = sep,
#         remove = remove,
#         convert = convert,
#         extra = extra,
#         fill = fill,
#         ...))
#     data
# }
#
# separateRow = function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
#                         convert = FALSE, extra = "warn", fill = "warn", ...){
#     rowData(data) = S4Vectors::DataFrame(tidyr::separate(
#         as.data.frame(rowData(data)),
#         col = col,
#         into = into,
#         sep = sep,
#         remove = remove,
#         convert = convert,
#         extra = extra,
#         fill = fill,
#         ...))
#     data
# }
