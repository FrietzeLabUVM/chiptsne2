


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

#' @export
setMethod("rowToRowMat", "ChIPtsne2_no_rowRanges", function(x, withDimnames=TRUE) {
    out <- x@rowToRowMat
    if (withDimnames)
        rownames(out) <- rownames(x)
    out
})

#' @export
setGeneric("colToRowMatCols", function(x, ...) standardGeneric("colToRowMatCols"))

#' @export
setMethod("colToRowMatCols", "ChIPtsne2_no_rowRanges", function(x, withDimnames=TRUE) {
    out <- x@colToRowMatCols
    out
})

#### Validity ####
#' @importFrom S4Vectors setValidity2
#' @importFrom BiocGenerics NCOL NROW
S4Vectors::setValidity2("ChIPtsne2_no_rowRanges", function(object) {
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
})

#### Show ####

#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "ChIPtsne2_no_rowRanges", function(object) {
    callNextMethod()
    cat(
        "rowToRowMat has ", ncol(rowToRowMat(object)), " columns\n",
        "colToRowMatCols has ", length(colToRowMatCols(object)), " items\n",
        sep=""
    )
})

#### Setter ####

#' @export
setGeneric("rowToRowMat<-", function(x, ..., value)
    standardGeneric("rowToRowMat<-")
)

#' @export
setReplaceMethod("rowToRowMat", "ChIPtsne2_no_rowRanges", function(x, value) {
    x@rowToRowMat <- value
    validObject(x)
    x
})

#' @export
setGeneric("colToRowMatCols<-", function(x, ..., value)
    standardGeneric("colToRowMatCols<-")
)

#' @export
setReplaceMethod("colToRowMatCols", "ChIPtsne2_no_rowRanges", function(x, value) {
    x@colToRowMatCols <- value
    validObject(x)
    x
})

#' @export
rowData = SummarizedExperiment::rowData
#' @export
colData = SummarizedExperiment::colData
#' @export
assay = SummarizedExperiment::assay

#### Subsetting by index ####

#' @export
setMethod("[", "ChIPtsne2_no_rowRanges", function(x, i, j, drop=TRUE) {
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
})

hasDimReduce = function(ct2){
    meta_dt = getRegionMetaData(ct2)
    all(c("tx", "ty") %in% colnames(meta_dt))
}

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

#### split ####

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
setMethod("split", "ChIPtsne2_no_rowRanges", function(x, f = NULL, drop=FALSE, ...){
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
})

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
cbind = SummarizedExperiment::cbind
# rowToRowMat = rowToRowMat,
# colToRowMatCols = colToRowMatCols,
# name_VAR = name_VAR,
# position_VAR = position_VAR,
# value_VAR = value_VAR,
# region_VAR = region_VAR,
# fetch_config = fetch_config

#' @param ChIPtsne2_no_rowRanges
#'
#' @return
#' @export
#'
#' @examples
setMethod("cbind", "ChIPtsne2_no_rowRanges", function(..., deparse.level=1) {
    args <- list(...)
    # cns = unname(unlist(lapply(args, colnames)))
    # cn_dupes = duplicated(cns)
    # if(any(cn_dupes)){
    #     stop("Duplicated colnames are not allowed when combining ChIPtsne2_no_rowRanges objects. You may need to use setNameVariable to differentiate names between ChIPtsne2_no_rowRanges objects. Duplicated examples:\n",
    #          paste(head(unique(cns[cn_dupes])), collapse = "\n"))
    # }
    .validate_names_unique(args, colnames, "Column")
    .validate_names_match(args, rownames, "Row")


    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)
    all.c2rrm <- lapply(args, colToRowMatCols, withDimnames=FALSE)

    # cns = unlist(lapply(all.rrm, colnames))
    # cn_dupes = duplicated(cns)
    # if(any(cn_dupes)){
    #     stop("Duplicated colnames are not allowed when combining rowToRowMat. You may need to use setNameVariable to differentiate names between ChIPtsne2_no_rowRanges objects. Duplicated examples:\n",
    #          paste(head(unique(cns[cn_dupes])), collapse = "\n"))
    # }

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
})


#### rbind ####
#' @export
rbind = SummarizedExperiment::rbind
#' @param ChIPtsne2_no_rowRanges
#'
#' @return
#' @export
#'
#' @examples
setMethod("rbind", "ChIPtsne2_no_rowRanges", function(..., deparse.level=1) {
    args <- list(...)
    .validate_names_unique(args, rownames, "Row")
    .validate_names_match(args, colnames, "Column")
    # ref = args[[1]]
    # for(test in args[-1]){
    #     is_match = colnames(ref) == colnames(test)
    #     if(!all(is_match)){
    #         a = colnames(ref)[!is_match]
    #         b = colnames(test)[!is_match]
    #         stop(paste(c("Column names must be identical for all ChIPtsne2_no_rowRanges objects. Example mismatches: ", head(paste(a, b, sep = " != "))), collapse = "\n"))
    #     }
    # }

    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)

    # cns = unlist(lapply(all.rrm, colnames))
    # cn_dupes = duplicated(cns)
    # if(any(cn_dupes)){
    #     stop("Duplicated colnames are not allowed when combining rowToRowMat. You may need to use setNameVariable to differentiate names between ChIPtsne2_no_rowRanges objects. Duplicated examples:\n",
    #          paste(head(unique(cns[cn_dupes])), collapse = "\n"))
    # }

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
})

# getMethod("colnames", "SummarizedExperiment")
# getMethod("rownames", "SummarizedExperiment")
# getMethod(`rownames<-`, "SummarizedExperiment")
# showMethods(`rownames<-`)
# getMethod(`rownames<-`, "DFrame")
setMethod("colnames", "ChIPtsne2_no_rowRanges", function(x, do.NULL = TRUE, prefix = "col") {
    callNextMethod()
})
setMethod("rownames", "ChIPtsne2_no_rowRanges", function(x, do.NULL = TRUE, prefix = "row") {
    callNextMethod()
})
setMethod("colnames<-", "ChIPtsne2_no_rowRanges", function(x, value) {
    out = callNextMethod()
    .update_ct2_colnames(x, new_names = colnames(out))
})
setMethod("rownames<-", "ChIPtsne2_no_rowRanges", function(x, value) {
    stop("NYI")
    out = callNextMethod()
    # head(x@rowToRowMat)
    browser()
})
