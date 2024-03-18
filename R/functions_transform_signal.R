.recalculateMax_ct2 = function(ct2){
    r2rm = ct2@rowToRowMat
    c2rmc = ct2@colToRowMatCols
    .recalculateMax(r2rm, c2rmc)
}

.recalculateMax = function(r2rm, c2rmc){
    abs_max = function(x){
        x[which.max(abs(x))]
    }
    resl = lapply(c2rmc, function(x){
        apply(r2rm[,x,drop = FALSE], 1, abs_max)
    })
    df = as.data.frame(resl)
    colnames(df) = names(c2rmc)
    rownames(df) = rownames(r2rm)
    as.matrix(df)
}

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

#' transformSignal method for objects of class \code{ChIPtsne2_no_rowRanges}.
#'
#' @docType methods
#' @name transformSignal-ChIPtsne2_no_rowRanges
#' @rdname transformSignal-ChIPtsne2_no_rowRanges
#' @aliases transformSignal-ChIPtsne2_no_rowRanges transformSignal,ChIPtsne2_no_rowRanges-method
#'
#' @param x A \code{ChIPtsne2_no_rowRanges} object.
#'
#' @export
setMethod("transformSignal", c("ChIPtsne2_no_rowRanges"), .transformSignal)
