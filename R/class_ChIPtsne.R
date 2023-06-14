


#### Constructor ####

#' @export
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges mcols
ChIPtsne2 <- function(
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
    if(!is.null(rowRanges(se))){
        k = colnames(GenomicRanges::mcols(rowRanges(se))) %in% colnames(se)
        if(any(k)){
            old_names = colnames(GenomicRanges::mcols(se@rowRanges))[k]
            new_names = paste0("region_", colnames(GenomicRanges::mcols(rowRanges(se)))[k])
            warning("Modifying region metadata column names to prevent collision with ChIPtsne2 colnames.\n",
                    paste(paste(old_names, "->", new_names), collapse = "\n"))
            colnames(GenomicRanges::mcols(se@rowRanges))[k] =
                new_names
        }
    }
    .ChIPtsne2(se,
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
setMethod("rowToRowMat", "ChIPtsne2", function(x, withDimnames=TRUE) {
    out <- x@rowToRowMat
    if (withDimnames)
        rownames(out) <- rownames(x)
    out
})

#' @export
setGeneric("colToRowMatCols", function(x, ...) standardGeneric("colToRowMatCols"))

#' @export
setMethod("colToRowMatCols", "ChIPtsne2", function(x, withDimnames=TRUE) {
    out <- x@colToRowMatCols
    out
})

#### Validity ####
#' @importFrom S4Vectors setValidity2
#' @importFrom BiocGenerics NCOL NROW
S4Vectors::setValidity2("ChIPtsne2", function(object) {
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
    #rowRanges mcols cannot include colnames
    if(any(colnames(GenomicRanges::mcols(rowRanges(object))) %in% colnames(object))){
        msg <- c(
            msg, "'mcols(rowRanges(object))' may not have entries equal to colnames of ChIPtsne2 object."
        )
    }

    if (length(msg)) {
        msg
    } else TRUE
})

#### Show ####

#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "ChIPtsne2", function(object) {
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
setReplaceMethod("rowToRowMat", "ChIPtsne2", function(x, value) {
    x@rowToRowMat <- value
    validObject(x)
    x
})

#' @export
setGeneric("colToRowMatCols<-", function(x, ..., value)
    standardGeneric("colToRowMatCols<-")
)

#' @export
setReplaceMethod("colToRowMatCols", "ChIPtsne2", function(x, value) {
    x@colToRowMatCols <- value
    validObject(x)
    x
})

#' @export
rowData = SummarizedExperiment::rowData
#' @export
rowRanges = SummarizedExperiment::rowRanges
#' @export
colData = SummarizedExperiment::colData
#' @export
assay = SummarizedExperiment::assay

#' getSampleMetaData
#'
#' @param ct2 A ChIPtsne object
#'
#' @return data.frame with sample meta data, similar to colData but suitable for tidyverse operations.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' getSampleMetaData(ct2)
getSampleMetaData = function(ct2){
    cd = colData(ct2)
    df = as.data.frame(cd)
    df[[ct2@name_VAR]] = rownames(cd)
    rownames(df)
    df
}

#' getRegionMetaData
#'
#' @param ct2 A ChIPtsne object
#'
#' @return data.frame with region meta data, similar to rowRanges but suitable for tidyverse operations.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' getRegionMetaData(ct2)
getRegionMetaData = function(ct2){
    gr = rowRanges(ct2)
    df = GenomicRanges::mcols(gr) %>% as.data.frame
    df[[ct2@region_VAR]] = names(gr)
    df
}

#' exampleQueryGR
#'
#' @return GRanges example
#' @export
#'
#' @importFrom GenomicRanges mcols
#' @examples
#' exampleQueryGR()
exampleQueryGR = function(){
    query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
    colnames(GenomicRanges::mcols(query_gr)) = paste0(
        "peak_",
        colnames(GenomicRanges::mcols(query_gr))
    )
    query_gr
}
#' exampleProfDT
#'
#' @return data.table example
#' @export
#'
#' @examples
#' exampleProfDT()
exampleProfDT = function(){
    seqsetvis::CTCF_in_10a_profiles_dt
}

#' query_gr
#'
#' @return ChIPtsne object for testing
#' @export
#'
#' @examples
#' exampleChIPtsne2()
exampleChIPtsne2 = function(){
    query_gr = exampleQueryGR()
    prof_dt = exampleProfDT()

    ChIPtsne2.from_tidy(prof_dt, query_gr)
}

#' exampleChIPtsne2.with_meta
#'
#' @return ChIPtsne object for testing, includes meta data
#' @export
#'
#' @examples
#' exampleChIPtsne2.with_meta()
exampleChIPtsne2.with_meta = function(){
    query_gr = exampleQueryGR()
    prof_dt = exampleProfDT()
    meta_dt = prof_dt %>%
        dplyr::select(sample) %>%
        unique %>%
        tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

    ChIPtsne2.from_tidy(prof_dt, query_gr, sample_metadata = meta_dt)
}

# For SummarizedExperiment slots
# Again, we can use the setter methods defined in SummarizedExperiment to modify slots in the base class. These should generally not require any re-defining. However, if it is necessary, the methods should use callNextMethod internally:
#
#     #' @export
#     #' @importMethodsFrom SummarizedExperiment "rowData<-"
#     setReplaceMethod("rowData", "ChIPtsne2", function(x, ..., value) {
#         y <- callNextMethod() # returns a modified ChIPtsne2
#
#         # Do something extra here.
#         message("hi!\n")
#
#         y
#     })

#### Subsetting by index ####

#' @export
setMethod("[", "ChIPtsne2", function(x, i, j, drop=TRUE) {
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

