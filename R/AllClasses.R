#' .FetchConfig
#'
#' @slot view_size Consistent size to use when viewing assessment regions. Uses
#'   CT_VIEW_SIZE option or 3kb as default.
#' @slot read_mode Read mode of signal data, one of bam_SE, bam_PE, or bigwig.
#'   Use CT_READ_MODES$.
#' @slot fetch_options Named list of additional arguments to pass to signal
#'   fetch function.
#'
#' @rdname FetchConfig
#' @export
#'
.FetchConfig = setClass("FetchConfig",
                        representation = list(
                            meta_data = "data.frame",
                            is_null = "logical",
                            view_size = "numeric",
                            window_size = "numeric",
                            read_mode = "character",
                            fetch_options = "list",
                            name_VAR = "character"
                        )
)

#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.ChIPtsne2_no_rowRanges <- setClass("ChIPtsne2_no_rowRanges",
                       slots= representation(
                           rowToRowMat="matrix",
                           colToRowMatCols="list",
                           name_VAR = "character",
                           position_VAR = "character",
                           value_VAR = "character",
                           region_VAR = "character",
                           fetch_config = "FetchConfig"
                       ),
                       contains="SummarizedExperiment"
)

#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.ChIPtsne2 <- setClass("ChIPtsne2",
                       contains= c("RangedSummarizedExperiment", "ChIPtsne2_no_rowRanges")
)

#' @export
#' @importClassesFrom S4Vectors List
.ChIPtsne2List = setClass("ChIPtsne2List",
                          contains="SimpleList")

