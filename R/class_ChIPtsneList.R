#### Constructor ####

#' @export
#'
#' @importFrom S4Vectors SimpleList
ChIPtsne2List <- function(
        ...)
{
    ctl <- S4Vectors::SimpleList(...)
    # .ChIPtsne2List(ctl)
    new("ChIPtsne2List", listData = ctl@listData, elementType = "ChIPtsne2", elementMetadata = ctl@elementMetadata, metadata = ctl@metadata)
}

#' @param ChIPtsne2
#'
#' @return
#' @export
#'
#' @examples
setMethod("cbind", "ChIPtsne2List", function(..., deparse.level=1) {
    args <- list(...)
    if(length(args) > 1){
        stop("Only one ChIPtsne2List may be supplied.")
    }
    do.call(cbind, as.list(args[[1]]))
})


#### rbind ####
#' @param ChIPtsne2
#'
#' @return
#' @export
#'
#' @examples
setMethod("rbind", "ChIPtsne2List", function(..., deparse.level=1) {
    args <- list(...)
    if(length(args) > 1){
        stop("Only one ChIPtsne2List may be supplied.")
    }
    do.call(rbind, as.list(args[[1]]))
})

unsplit = BiocGenerics::unsplit
#### unsplit ####
#' @param ChIPtsne2
#'
#' @return
#' @export
#'
#' @examples
setMethod("unsplit", "ChIPtsne2List", function (value, f, drop = FALSE){
    mode = "by_column"
    x = value[[1]]
    sample_meta_data = getSampleMetaData(x)
    region_meta_data = getRegionMetaData(x)
    if(is.null(f)){
        f = colnames(x)
        names(f) = f
    }
    if(length(f) == 1){
        if(f %in% colnames(sample_meta_data)){
            f = split(rownames(sample_meta_data), sample_meta_data[[f]])
        }else if(f %in% colnames(region_meta_data)){
            f = split(region_meta_data[[x@region_VAR]], region_meta_data[[f]])
            mode = "by_row"
        }
    }else{
        if(ncol(x) == nrow(x)){
            stop("Cannot unambiguously split ChIPtsne2 using a vector when ncol == nrow. Try using a row or column attribute name.")
        }
        if(length(f) == ncol(x)){
            f = split(colnames(x), f)
        }else if(length(f) == nrow(x)){
            f = split(rownames(x), f)
            mode = "by_row"
        }
    }

    if(mode == "by_column"){
        x.unsplit = cbind(value)
    }else{
        x.unsplit = rbind(value)
    }
    x.unsplit
})
