.recalculateMax_ct2 = function(ct2){
    r2rm = ct2@rowToRowMat
    c2rmc = ct2@colToRowMatCols
    assay(ct2, "max") = .recalculateMax(r2rm, c2rmc)
    ct2
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
    for(cn in names(ct2@colToRowMatCols)){
        sel_cols = ct2@colToRowMatCols[[cn]]
        ct2@rowToRowMat[, sel_cols] = transform_FUN(ct2@rowToRowMat[, sel_cols], getSampleMetaData(ct2)[cn,])
    }
    ct2 = .recalculateMax_ct2(ct2)

    history_item = list(transformSignal = list(FUN = .transformSignal, ARG = args))
    ct2@metadata = c(ct2@metadata, history_item)
    ct2
}

#' transformSignal method for objects of class \code{ChIPtsne2_no_rowRanges}.
#'
#' @docType methods
#' @rdname ct2-trans
#' @aliases transformSignal-ChIPtsne2_no_rowRanges transformSignal,ChIPtsne2_no_rowRanges-method
#'
#' @param x A \code{ChIPtsne2_no_rowRanges} object.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ct2 = exampleChIPtsne2.with_meta()
#' plotSignalLinePlot(ct2) + expand_limits(y = 30)
#' # a transformation function takes a matrix and returns the same size matrix
#' trans_fun_1 = function(mat, ...){
#'     mat/2
#' }
#'
#' ct2.t1 = transformSignal(ct2, trans_fun_1)
#' plotSignalLinePlot(ct2.t1) + expand_limits(y = 30)
#'
#' # as a second argument, sample metadata is passed in per sample
#' trans_fun_2 = function(mat, sample_meta){
#'   if(sample_meta$cell == "MCF10CA1"){
#'     mat / 4
#'   }else{
#'     mat
#'   }
#' }
#' ct2.t2 = transformSignal(ct2, trans_fun_2)
#' plotSignalLinePlot(ct2.t2) + expand_limits(y = 30)
#'
setGeneric("transformSignal",
           function(ct2, transform_FUN)
               standardGeneric("transformSignal"),
           signature = "ct2")

#' @rdname ct2-trans
#' @export
setMethod("transformSignal", c("ChIPtsne2_no_rowRanges"), .transformSignal)
