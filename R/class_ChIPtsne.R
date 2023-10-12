


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

#' addRegionAnnotation
#'
#' An alternative to setRegionMetaData that adds region metadata information based on overlaps with supplied anno_gr.
#'
#' @param ct2 A ChIPtsne object
#' @param anno_gr A GenomicRanges object to annotate ct2 based on overlap with rowRanges of ct2.
#' @param anno_VAR Attribute in mcols of anno_gr to pull values from.
#' @param anno_VAR_renames Matched vector to anno_VAR specificying final names in rowRanges of ct2. Essentially renames anno_VAR.
#' @param no_overlap_value Value for when there is no overlap with anno_gr. Default is "no_hit".
#' @param overlap_value Value for when there is an overlap, only relevant if anno_VAR is not in mcols of anno_gr. I.e. adding a single "hit" "no hit" annotation.
#'
#' @return A modified ChIPtsne2 object with added/overwritten sample metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' anno_gr = rowRanges(ct2)
#' anno_gr = subset(anno_gr, seqnames == "chr1")
#' anno_gr$start_pos = start(anno_gr)
#' anno_gr$end_pos = end(anno_gr)
#' anno_gr$chr_name = as.character(seqnames(anno_gr))
#' ct2.anno = addRegionAnnotation(ct2, anno_gr, anno_VAR = c("start_pos", "end_pos", "chr_name"))
#' rowRanges(ct2.anno)
addRegionAnnotation = function(ct2,
                               anno_gr,
                               anno_VAR = "annotation",
                               anno_VAR_renames = anno_VAR,
                               no_overlap_value = "no_hit",
                               overlap_value = "hit"){
    message("addRegionAnnotation ...")
    args = get_args()
    valid_anno_names = colnames(GenomicRanges::mcols(anno_gr))
    if(length(anno_VAR) != length(anno_VAR_renames)){
        stop("Length of anno_VAR and anno_VAR_renames must be identical")
    }
    if(length(anno_VAR) > 1){
        if(!all(anno_VAR %in% valid_anno_names)){
            stop("Supplied anno_VARS must present in mcols of anno_gr. Missing: ", paste(setdiff(anno_VAR, valid_anno_names), collapse = ", "))
        }
    }else{
        if(!anno_VAR %in% valid_anno_names){
            GenomicRanges::mcols(anno_gr)[[anno_VAR]] = overlap_value
        }
    }
    if(length(no_overlap_value) != 1 | length(no_overlap_value) != length(anno_VAR)){
        stop("length(no_overlap_value) must be 1 or match anno_VAR length. Was: ", length(no_overlap_value), ".")
    }
    if(length(no_overlap_value) == 1){
        no_overlap_value = rep(no_overlap_value, length(anno_VAR))
    }

    olaps = GenomicRanges::findOverlaps(query = ct2, subject = anno_gr)
    new_gr = rowRanges(ct2)
    for(i in seq_along(anno_VAR)){
        av = anno_VAR[i]
        av_new = anno_VAR_renames[i]
        anno_vals = GenomicRanges::mcols(anno_gr)[[av]][S4Vectors::subjectHits(olaps)]
        GenomicRanges::mcols(new_gr)[[av_new]] = no_overlap_value[i]
        GenomicRanges::mcols(new_gr)[[av_new]][S4Vectors::queryHits(olaps)] = anno_vals
    }

    ct2@rowRanges = new_gr
    history_item = list(addRegionAnnotation = list(FUN = addRegionAnnotation, ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)
    ct2
}

#### split ####

#' @param ChIPtsne2
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
setMethod("split", "ChIPtsne2", function(x, f = NULL, drop=FALSE, ...){
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
            stop(paste(c(paste0(str, " names must be identical for all ChIPtsne2 objects. Example mismatches: "), head(paste(a, b, sep = " != "))), collapse = "\n"))
        }
    }
}

.validate_names_unique = function(args, dim_FUN, str){
    cns = unname(unlist(lapply(args, dim_FUN)))
    cn_dupes = duplicated(cns)
    if(any(cn_dupes)){
        stop(paste0("Duplicated ", str, " names are not allowed when combining ChIPtsne2 objects. You may need to use setNameVariable to differentiate names between ChIPtsne2 objects. Duplicated examples:\n"),
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

#' @param ChIPtsne2
#'
#' @return
#' @export
#'
#' @examples
setMethod("cbind", "ChIPtsne2", function(..., deparse.level=1) {
    args <- list(...)
    # cns = unname(unlist(lapply(args, colnames)))
    # cn_dupes = duplicated(cns)
    # if(any(cn_dupes)){
    #     stop("Duplicated colnames are not allowed when combining ChIPtsne2 objects. You may need to use setNameVariable to differentiate names between ChIPtsne2 objects. Duplicated examples:\n",
    #          paste(head(unique(cns[cn_dupes])), collapse = "\n"))
    # }
    .validate_names_unique(args, colnames, "Column")
    .validate_names_match(args, rownames, "Row")


    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)
    all.c2rrm <- lapply(args, colToRowMatCols, withDimnames=FALSE)

    # cns = unlist(lapply(all.rrm, colnames))
    # cn_dupes = duplicated(cns)
    # if(any(cn_dupes)){
    #     stop("Duplicated colnames are not allowed when combining rowToRowMat. You may need to use setNameVariable to differentiate names between ChIPtsne2 objects. Duplicated examples:\n",
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
#' @param ChIPtsne2
#'
#' @return
#' @export
#'
#' @examples
setMethod("rbind", "ChIPtsne2", function(..., deparse.level=1) {
    args <- list(...)
    .validate_names_unique(args, rownames, "Row")
    .validate_names_match(args, colnames, "Column")
    # ref = args[[1]]
    # for(test in args[-1]){
    #     is_match = colnames(ref) == colnames(test)
    #     if(!all(is_match)){
    #         a = colnames(ref)[!is_match]
    #         b = colnames(test)[!is_match]
    #         stop(paste(c("Column names must be identical for all ChIPtsne2 objects. Example mismatches: ", head(paste(a, b, sep = " != "))), collapse = "\n"))
    #     }
    # }

    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)

    # cns = unlist(lapply(all.rrm, colnames))
    # cn_dupes = duplicated(cns)
    # if(any(cn_dupes)){
    #     stop("Duplicated colnames are not allowed when combining rowToRowMat. You may need to use setNameVariable to differentiate names between ChIPtsne2 objects. Duplicated examples:\n",
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
