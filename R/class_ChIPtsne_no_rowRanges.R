


#### Constructor ####

#' @export
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges mcols
ChIPtsne2_no_rowRanges <- function(
        rowToRowMat=matrix(0,0,0),
        colToRowMatCols=list(),
        name_VAR = "sample",
        position_VAR = "x",
        value_VAR = "y",
        region_VAR = "id",
        fetch_config = FetchConfig.null(),
        ...)
{
    se <- SummarizedExperiment(...)
    .ChIPtsne2_no_rowRanges(se,
               rowToRowMat = rowToRowMat,
               colToRowMatCols = colToRowMatCols,
               name_VAR = name_VAR,
               position_VAR = position_VAR,
               value_VAR = value_VAR,
               region_VAR = region_VAR,
               fetch_config = fetch_config
    )
}

#### Getters ####
#' @export
setGeneric("rowToRowMat", function(x, ...) standardGeneric("rowToRowMat"))

ct2_nrr_get_rowToRowMat = function(x, withDimnames=TRUE) {
    out <- x@rowToRowMat
    if (withDimnames)
        rownames(out) <- rownames(x)
    out
}

#' @export
setMethod("rowToRowMat", "ChIPtsne2_no_rowRanges", ct2_nrr_get_rowToRowMat)

#' @export
setGeneric("colToRowMatCols", function(x, ...) standardGeneric("colToRowMatCols"))

ct2_nrr_get_colToRowMatCols =  function(x, withDimnames=TRUE) {
    out <- x@colToRowMatCols
    out
}
#' @export
setMethod("colToRowMatCols", "ChIPtsne2_no_rowRanges", ct2_nrr_get_colToRowMatCols)

#### Validity ####

ct2_nrr_validity = function(object) {
    NR <- NROW(object)
    NC <- NCOL(object)
    msg <- NULL

    # 2D
    if (NROW(rowToRowMat(object, withDimnames=FALSE)) != NR) {
        msg <- c(
            msg, "'nrow(rowToRowMat)' should be equal to the number of rows"
        )
    }
    # list
    if (length(colToRowMatCols(object)) != NC) {
        msg <- c(
            msg, "'length(colToRowMatCols)' should be equal to the number of columns"
        )
    }
    if (length(msg)) {
        msg
    } else TRUE
}

#' @importFrom S4Vectors setValidity2
#' @importFrom BiocGenerics NCOL NROW
S4Vectors::setValidity2("ChIPtsne2_no_rowRanges", ct2_nrr_validity)

#### Show ####

ct2_nrr_show = function(object) {
    callNextMethod()
    cat(
        "rowToRowMat has ", ncol(rowToRowMat(object)), " columns\n",
        "colToRowMatCols has ", length(colToRowMatCols(object)), " items\n",
        sep=""
    )
}

#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "ChIPtsne2_no_rowRanges", ct2_nrr_show)

#### Setter ####

#' @export
setGeneric("rowToRowMat<-", function(x, ..., value)
    standardGeneric("rowToRowMat<-")
)

ct2_nrr_set_rowToRowMat = function(x, value) {
    x@rowToRowMat <- value
    validObject(x)
    x
}

#' @export
setReplaceMethod("rowToRowMat", "ChIPtsne2_no_rowRanges", ct2_nrr_set_rowToRowMat)

#' @export
setGeneric("colToRowMatCols<-", function(x, ..., value)
    standardGeneric("colToRowMatCols<-")
)

ct2_nrr_set_colToRowMatCols = function(x, value) {
    x@colToRowMatCols <- value
    validObject(x)
    x
}

#' @export
setReplaceMethod("colToRowMatCols", "ChIPtsne2_no_rowRanges", ct2_nrr_set_colToRowMatCols)

#### Subsetting by index ####

ct2_nrr_index_accessor = function(x, i, j, drop=TRUE) {
    rrm <- rowToRowMat(x, withDimnames=FALSE)
    c2rrm = colToRowMatCols(x)

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                j, colnames(x), fmt
            )
        }
        j <- as.vector(j)
        c2rrm = c2rrm[j]
        rrm <- rrm[, unlist(c2rrm),drop=FALSE]
    }

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(x), fmt
            )
        }
        i <- as.vector(i)
        rrm <- rrm[i,,drop=FALSE]
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out,
                                rowToRowMat = rrm,
                                colToRowMatCols = c2rrm,
                                check=FALSE)
}
#' @export
setMethod("[", "ChIPtsne2_no_rowRanges", ct2_nrr_index_accessor)

#### split ####

ct2_nrr_split = function(x, f = NULL, drop=FALSE, ...){
    mode = "by_column"
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
            stop("Cannot unambiguously split ChIPtsne2_no_rowRanges using a vector when ncol == nrow. Try using a row or column attribute name.")
        }
        if(length(f) == ncol(x)){
            f = split(colnames(x), f)
        }else if(length(f) == nrow(x)){
            f = split(rownames(x), f)
            mode = "by_row"
        }
    }

    if(mode == "by_column"){
        x.split = lapply(f, function(split_val)x[, split_val])
    }else{
        x.split = lapply(f, function(split_val)x[split_val, ])
    }
    ChIPtsne2List(x.split)
}

#' @param ChIPtsne2_no_rowRanges
#'
#' @export
#' @examples
#' split(ct2, "sample")
#' split(ct2, colnames(ct2))
#' split(ct2, "cell")
#' split(ct2, "peak_MCF10CA1_CTCF")
#' split(ct2, ct2$cell)
#'
#' sample_meta_data = getSampleMetaData(ct2)
#' region_meta_data = getRegionMetaData(ct2)
#'
#' split(ct2, sample_meta_data$mark)
#' split(ct2, region_meta_data$peak_MCF10A_CTCF)
#'
setMethod("split", "ChIPtsne2_no_rowRanges", ct2_nrr_split)

.validate_names_match = function(args, dim_FUN, str){
    ref = args[[1]]
    for(test in args[-1]){
        is_match = dim_FUN(ref) == dim_FUN(test)
        if(!all(is_match)){
            a = dim_FUN(ref)[!is_match]
            b = dim_FUN(test)[!is_match]
            stop(paste(c(paste0(str, " names must be identical for all ChIPtsne2_no_rowRanges objects. Example mismatches: "), head(paste(a, b, sep = " != "))), collapse = "\n"))
        }
    }
}

.validate_names_unique = function(args, dim_FUN, str){
    cns = unname(unlist(lapply(args, dim_FUN)))
    cn_dupes = duplicated(cns)
    if(any(cn_dupes)){
        stop(paste0("Duplicated ", str, " names are not allowed when combining ChIPtsne2_no_rowRanges objects. You may need to use setNameVariable to differentiate names between ChIPtsne2_no_rowRanges objects. Duplicated examples:\n"),
             paste(head(unique(cns[cn_dupes])), collapse = "\n"))
    }
}

#### cbind ####

ct2_nrr_cbind = function(..., deparse.level=1) {
    args <- list(...)
    .validate_names_unique(args, colnames, "Column")
    .validate_names_match(args, rownames, "Row")


    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)
    all.c2rrm <- lapply(args, colToRowMatCols, withDimnames=FALSE)

    all.rrm <- do.call(cbind, all.rrm)
    names(all.c2rrm) = NULL
    all.c2rrm <- do.call(c, all.c2rrm)

    # Checks for identical column state.
    ref <- args[[1]]
    ref.rrm <- rowToRowMat(ref, withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(rownames(ref.rrm), rownames(rowToRowMat(x, withDimnames=FALSE))))
        {
            stop("per-row values are not compatible")
        }
    }

    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(
        out,
        rowToRowMat=all.rrm,
        colToRowMatCols=all.c2rrm,
        check=FALSE)
}

#' @param ChIPtsne2_no_rowRanges
#'
#' @return a ChIPtsne2 object of concatenated columns/samples of all items in input ChIPtsne2List
#' @export
#' @rdname cbind
#'
#' @examples
#' ct2.left = exampleChIPtsne2()
#' ct2.right = exampleChIPtsne2()
#' colnames(ct2.right) = paste0("right_", colnames(ct2.right))
#' ct2 = cbind(ct2.left, ct2.right)
#' dim(ct2)
#' colnames(ct2)
setMethod("cbind", "ChIPtsne2_no_rowRanges", ct2_nrr_cbind)

ct2_nrr_rbind = function(..., deparse.level=1) {
    args <- list(...)
    .validate_names_unique(args, rownames, "Row")
    .validate_names_match(args, colnames, "Column")

    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)

    all.rrm <- do.call(rbind, all.rrm)

    # Checks for identical column state.
    ref <- args[[1]]
    ref.rrm <- rowToRowMat(ref, withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(colnames(ref.rrm), colnames(rowToRowMat(x, withDimnames=FALSE))))
        {
            stop("per-row values are not compatible")
        }
    }

    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(
        out,
        rowToRowMat=all.rrm,
        check=FALSE)
}


#### rbind ####
#' @param ChIPtsne2_no_rowRanges
#'
#' @return a ChIPtsne2 object of concatenated rows/regions of all items in input ChIPtsne2List
#' @export
#' @aliases rbind cbind
#'
#' @examples
#' ct2_a = exampleChIPtsne2()
#' ct2_b = exampleChIPtsne2()
#' # duplicated rownames are not allowed so we need to modify before rbind
#' rownames(ct2_b) = paste0("b_", rownames(ct2_b))
#' ct2_rbind = rbind(ct2_a, ct2_b)
#' rownames(ct2_rbind)
setMethod("rbind", "ChIPtsne2_no_rowRanges", ct2_nrr_rbind)

ct2_nrr_set_dimnames =  function(x, value){
    x = .update_ct2_rownames(x, new_names = value[[1]])
    x = .update_ct2_colnames(x, new_names = value[[2]])
    x
}

setReplaceMethod("dimnames", c("ChIPtsne2_no_rowRanges", "list"), ct2_nrr_set_dimnames)

