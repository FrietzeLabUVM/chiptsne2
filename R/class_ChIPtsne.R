


#### Constructor ####
#https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/Extensions.html
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

#### Validity ####
ct2_validity = function(object) {
    NR <- NROW(object)
    NC <- NCOL(object)
    msg <- NULL

    # 2D
    if (NROW(rowToRowMat(object)) != NR) {
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
}
#' @importFrom S4Vectors setValidity2
#' @importFrom BiocGenerics NCOL NROW
S4Vectors::setValidity2("ChIPtsne2", ct2_validity)

#### addRegionAnnotation ####


#' addRegionAnnotation
#'
#' An alternative to setRegionMetaData that adds region metadata information based on overlaps with supplied anno_gr.
#'
#' @param ct2 A ChIPtsne object
#' @param anno_gr A GenomicRanges object to annotate ct2 based on overlap with rowRanges of ct2.
#' @param anno_VAR Attribute in mcols of anno_gr to pull values from.
#' @param anno_VAR_renames Matched vector to anno_VAR specifying final names in rowRanges of ct2. Essentially renames anno_VAR.
#' @param no_overlap_value Value for when there is no overlap with anno_gr. Default is "no_hit".
#' @param overlap_value Value for when there is an overlap, only relevant if anno_VAR is not in mcols of anno_gr. I.e. adding a single "hit" "no hit" annotation.
#'
#' @return A modified ChIPtsne2 object with added/overwritten sample metadata.
#' @export
#'
#' @examples
#' library(GenomicRanges)
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
    if(length(no_overlap_value) != 1 & length(no_overlap_value) != length(anno_VAR)){
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


#### replace rowRanges, names, dimnames ####

ct2_replace_rowRanges = function(x, ..., value){
    if(is.null(value)){
        ChIPtsne2_no_rowRanges(
            assays = assays(x),
            rowData = rowData(x),
            colData = colData(x),
            rowToRowMat = x@rowToRowMat,
            colToRowMatCols = x@colToRowMatCols,
            name_VAR = x@name_VAR,
            position_VAR = x@position_VAR,
            value_VAR = x@value_VAR,
            region_VAR = x@region_VAR)
    }else{
        stop("Manipulating rowRanges is not supported except for removal by NULL assignment.")
    }
}

setReplaceMethod("rowRanges", c("ChIPtsne2", "NULL"), ct2_replace_rowRanges)


