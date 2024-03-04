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
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
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
                       contains= c("ChIPtsne2_no_rowRanges", "RangedSummarizedExperiment")
)

#' @export
#' @importClassesFrom S4Vectors List
.ChIPtsne2List = setClass("ChIPtsne2List",
                          contains="SimpleList")

#### Generics ####

#' getters and setters for ChIPtsne2 and ChIPtsne2_no_rowRanges objects
#'
#' @export
#' @rdname ct2-getset
setGeneric("setNameVariable",
           function(ct2, new_name_VAR){
               standardGeneric("setNameVariable")
           })
#' @export
#' @rdname ct2-getset
setGeneric("swapNameVariable",
           function(ct2, new_name_VAR){
               standardGeneric("swapNameVariable")
           })
#' @export
#' @rdname ct2-getset
setGeneric("getNameVariable",
           function(ct2){
               standardGeneric("getNameVariable")
           })

#' @export
#' @rdname ct2-getset
setGeneric("setValueVariable",
           function(ct2, new_value_VAR){
               standardGeneric("setValueVariable")
           })
#' @export
#' @rdname ct2-getset
setGeneric("getValueVariable",
           function(ct2){
               standardGeneric("getValueVariable")
           })
#' @export
#' @rdname ct2-getset
setGeneric("setRegionVariable",
           function(ct2, new_region_VAR){
               standardGeneric("setRegionVariable")
           })
#' @export
#' @rdname ct2-getset
setGeneric("getRegionVariable",
           function(ct2){
               standardGeneric("getRegionVariable")
           })
#' @export
#' @rdname ct2-getset
setGeneric("setPositionVariable",
           function(ct2, new_position_VAR){
               standardGeneric("setPositionVariable")
           })
#' @export
#' @rdname ct2-getset
setGeneric("getPositionVariable",
           function(ct2){
               standardGeneric("getPositionVariable")
           })

#' @export
#' @rdname ct2-getset
setGeneric("getRegionMetaData", function(ct2, select_VARS = NULL) standardGeneric("getRegionMetaData"))
