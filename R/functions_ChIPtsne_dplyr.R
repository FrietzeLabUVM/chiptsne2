#longterm should use:
#https://dplyr.tidyverse.org/reference/dplyr_extending.html

#' subsetRegions
#'
#' @param ct2 `doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from region/row metadata.
#'
#' @return A subsetted `doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' subsetRegions(ct2, peak_MCF10AT1_CTCF == TRUE)
subsetRegions = function(ct2, expr){
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, eval(ps))
}
#' subsetSamples
#'
#' @param ct2 `doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from sample/column metadata.
#'
#' @return A subsetted `doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' subsetSamples(ct2, cell %in% c("MCF10A", "MCF10AT1"))
subsetSamples = function(ct2, expr){
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, TRUE, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, TRUE, eval(ps))
}

#' subsetRow
#'
#' @param ct2 `doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from region/row metadata.
#'
#' @return A subsetted `doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' #deprecated, use subsetRegions instead
#' subsetRegions(ct2, peak_MCF10AT1_CTCF == TRUE)
subsetRow = function(ct2, expr){
    .Deprecated("subsetRegions")
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, eval(ps))
}

#' subsetCol
#'
#' @param ct2 `doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from sample/column metadata.
#'
#' @return A subsetted `doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' #deprecated, use subsetRegions instead
#' subsetSamples(ct2, cell %in% c("MCF10A", "MCF10AT1"))
subsetCol = function(ct2, expr){
    .Deprecated("subsetSamples")
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, TRUE, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, TRUE, eval(ps))
}

# #' @export
# filter.ChIPtsne2 = function(.data, ...){
#     subset(.data, ...)
# }

