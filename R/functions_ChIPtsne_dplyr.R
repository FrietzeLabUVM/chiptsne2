#longterm should use:
#https://dplyr.tidyverse.org/reference/dplyr_extending.html

#' subsetRegions
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from region/row metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
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
#' @param ct2 `r doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from sample/column metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
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
#' @param ct2 `r doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from region/row metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
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
#' @param ct2 `r doc_ct2_nrr()`
#' @param expr expression, indicating columns to select from sample/column metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
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

#

#' mutateSamples
#'
#' Should work just like [dplyr::mutate] but applied to colData/sample metadata of ChIPtsne2 objects.
#'
#' @param .data See [dplyr::mutate]
#' @param ... Passed to [dplyr::mutate]
#'
#' @return A `r doc_ct2_nrr()` with modified colData/sample metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' colData(ct2)
#' getSampleMetaData(ct2)
#' ct2 = mutateSamples(ct2, cell_mark = paste(cell, mark))
#' colData(ct2)
#' getSampleMetaData(ct2)
mutateSamples = function(.data, ...){
    colData(.data) = S4Vectors::DataFrame(dplyr::mutate(as.data.frame(colData(.data)), ...))
    .data
}

#' mutateRegions
#'
#' Should work just like [dplyr::mutate] but applied to rowData/region metadata of ChIPtsne2 objects.
#'
#' @param .data See [dplyr::mutate]
#' @param ... Passed to [dplyr::mutate]
#'
#' @return A `r doc_ct2_nrr()` with modified rowData/region metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' rowData(ct2)
#' getRegionMetaData(ct2)
#' ct2 = mutateRegions(ct2, new_col = peak_MCF10A_CTCF | peak_MCF10AT1_CTCF)
#' rowData(ct2)
#' getRegionMetaData(ct2)
mutateRegions = function(.data, ...){
    rowData(.data) = S4Vectors::DataFrame(dplyr::mutate(as.data.frame(rowData(.data)), ...))
    .data
}

#' separateSamples
#'
#' Should work just like [tidyr::separate] but applied to colData/sample metadata of ChIPtsne2 objects.
#'
#' @param data See [tidyr::separate]
#' @param col See [tidyr::separate]
#' @param into See [tidyr::separate]
#' @param sep See [tidyr::separate]
#' @param remove See [tidyr::separate]
#' @param convert See [tidyr::separate]
#' @param extra See [tidyr::separate]
#' @param fill See [tidyr::separate]
#' @param ... Passed to [tidyr::separate]
#'
#' @return A `r doc_ct2_nrr()` with modified colData/sample metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' colData(ct2)
#' getSampleMetaData(ct2)
#' ct2 = separateSamples(ct2, "cell_mark", c("v1", "v2"), sep = " ", remove = FALSE)
#' colData(ct2)
#' getSampleMetaData(ct2)
separateSamples = function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
                        convert = FALSE, extra = "warn", fill = "warn", ...){
    colData(data) = S4Vectors::DataFrame(tidyr::separate(
        as.data.frame(colData(data)),
        col = col,
        into = into,
        sep = sep,
        remove = remove,
        convert = convert,
        extra = extra,
        fill = fill,
        ...))
    data
}

#' separateRegions
#'
#' Should work just like [tidyr::separate] but applied to rowData/region metadata of ChIPtsne2 objects.
#'
#' @param data See [tidyr::separate]
#' @param col See [tidyr::separate]
#' @param into See [tidyr::separate]
#' @param sep See [tidyr::separate]
#' @param remove See [tidyr::separate]
#' @param convert See [tidyr::separate]
#' @param extra See [tidyr::separate]
#' @param fill See [tidyr::separate]
#' @param ... Passed to [tidyr::separate]
#'
#' @return A `r doc_ct2_nrr()` with modified rowData/region metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' rowData(ct2)
#' getRegionMetaData(ct2)
#' ct2 = mutateRegions(ct2, col_ab = paste("a", "b"))
#' ct2 = separateRegions(ct2, "col_ab", c("a", "b"), sep = " ", remove = FALSE)
#' rowData(ct2)
#' getRegionMetaData(ct2)
separateRegions = function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
                        convert = FALSE, extra = "warn", fill = "warn", ...){
    rowData(data) = S4Vectors::DataFrame(tidyr::separate(
        as.data.frame(rowData(data)),
        col = col,
        into = into,
        sep = sep,
        remove = remove,
        convert = convert,
        extra = extra,
        fill = fill,
        ...))
    data
}

