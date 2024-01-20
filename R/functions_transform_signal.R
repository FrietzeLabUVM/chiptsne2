.transformSignal = function(ct2, transform_FUN){
    args = get_args()


    for(cols in ct2@colToRowMatCols){
        ct2@rowToRowMat[, cols] = transform_FUN(ct2@rowToRowMat[, cols])
    }
    ct2 = .recalculateMax_ct2(ct2)

    history_item = list(transformSignal = list(FUN = .transformSignal, ARG = args))
    ct2@metadata = c(ct2@metadata, history_item)
    ct2
}

#' @export
setGeneric("transformSignal",
           function(ct2, transform_FUN)
               standardGeneric("transformSignal"),
           signature = "ct2")

#' @export
setMethod("transformSignal", c("ChIPtsne2_no_rowRanges"), .transformSignal)
