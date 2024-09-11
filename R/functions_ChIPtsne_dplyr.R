#longterm should use:
#https://dplyr.tidyverse.org/reference/dplyr_extending.html

#' subsetRegions
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param subset_expression expression, indicating columns to select from region/row metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' subsetRegions(ct2, peak_MCF10AT1_CTCF == TRUE)
subsetRegions = function(ct2, subset_expression){
    #because we can't store an expression, we need to convert to character for history
    test_expr = substitute(subset_expression)
    if(is.call(test_expr)){
        subset_expression = deparse(test_expr)
    }
    remove("test_expr")

    message("subsetRegions ...")
    args = get_args()
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    meta_data = getRegionMetaData(ct2, include_value_max = TRUE)
    meta_data = eval(substitute(subset(meta_data, eval(parse(text = subset_expression)))))
    ct2 = ct2[meta_data[[ct2@region_VAR]],]

    history_item = list(subsetRegions  = list(FUN = subsetRegions , ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)

    ct2
}
#' subsetSamples
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param subset_expression expression, indicating columns to select from sample/column metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' subsetSamples(ct2, cell %in% c("MCF10A", "MCF10AT1"))
subsetSamples = function(ct2, subset_expression){
    #because we can't store an expression, we need to convert to character for history
    test_expr = substitute(subset_expression)
    if(is.call(test_expr)){
        subset_expression = deparse(test_expr)
    }
    remove("test_expr")

    message("subsetSamples ...")
    args = get_args()
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    meta_data = getSampleMetaData(ct2)
    meta_data = eval(substitute(subset(meta_data, eval(parse(text = subset_expression)))))
    ct2 = ct2[, rownames(meta_data)]

    history_item = list(subsetSamples  = list(FUN = subsetSamples , ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)

    ct2
}

# .ids_from_value_selection(ct2, MCF10A_CTCF > 20 & MCF10AT1_CTCF > 30)
.ids_from_value_selection = function(ct2, subset_expression){
    values = as.data.frame(assays(ct2)$max)
    ps = substitute(subset_expression)
    values_filtered = subset(values, eval(ps))
    sel_ids = rownames(values_filtered)
    sel_ids
}

#' subsetValues
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param value_test expression, indicating column values to filter from assay slot.
#'
#' @return A subsetted `r doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' subsetValues(ct2, MCF10A_CTCF > 20 & MCF10AT1_CTCF > 20)
#' subsetValues(ct2, MCF10A_CTCF > Inf)
#'
#' min_signal = 20
#' ct2 = subsetValues(ct2, MCF10A_CTCF > min_signal & MCF10AT1_CTCF > min_signal)
#' ct2
subsetValues = function(ct2, value_test){
    #because we can't store an expression, we need to convert to character for history
    test_expr = substitute(value_test)
    if(is.call(test_expr)){
        value_test = deparse(test_expr)
    }
    remove("test_expr")

    message("subsetValues ...")
    args = get_args()

    sel_ids = eval(substitute(.ids_from_value_selection(ct2, eval(parse(text = value_test)))))
    ct2 = ct2[sel_ids,]

    history_item = list(subsetValues  = list(FUN = subsetValues , ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)

    ct2
}

#' subsetRow
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param value_test expression, indicating columns to select from region/row metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' #deprecated, use subsetRegions instead
#' subsetRegions(ct2, peak_MCF10AT1_CTCF == TRUE)
subsetRow = function(ct2, value_test){
    .Deprecated("subsetRegions")
    stop()
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(value_test))
    # subset(ct2, eval(parse(text = ssubset)))
    ps <- substitute(value_test)
    subset(ct2, eval(ps))
}

#' subsetCol
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param subset_expression expression, indicating columns to select from sample/column metadata.
#'
#' @return A subsetted `r doc_ct2_nrr()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' #deprecated, use subsetRegions instead
#' subsetSamples(ct2, cell %in% c("MCF10A", "MCF10AT1"))
subsetCol = function(ct2, subset_expression){
    .Deprecated("subsetSamples")
    stop()
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(subset_expression))
    # subset(ct2, TRUE, eval(parse(text = ssubset)))
    ps <- substitute(subset_expression)
    subset(ct2, TRUE, eval(ps))
}

#' mutateSamples
#'
#' Should work similarly to [dplyr::mutate] but applied to colData/sample metadata of ChIPtsne2 objects.
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param mutate_name  Name of new variable created by `mutate_expression`.
#' @param mutate_expression Expression to derive new variable values.
#' @param .by See [dplyr::mutate]
#' @param .keep  See [dplyr::mutate]
#' @param .before  See [dplyr::mutate]
#' @param .after  See [dplyr::mutate]
#'
#' @return A `r doc_ct2_nrr()` with modified colData/sample metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' colData(ct2)
#' getSampleMetaData(ct2)
#' ct2 = mutateSamples(ct2, "cell_mark", paste(cell, mark), )
#' colData(ct2)
#' getSampleMetaData(ct2)
mutateSamples = function(ct2,
                         mutate_name,
                         mutate_expression,
                         .by = NULL,
                         .keep = c("all", "used", "unused", "none")[1],
                         .before = NULL,
                         .after = NULL){
    #because we can't store an expression, we need to convert to character for history
    test_expr = substitute(mutate_expression)
    if(is.call(test_expr)){
        mutate_expression = deparse(test_expr)
    }
    remove("test_expr")

    message("mutateSamples ...")
    args = get_args()

    meta_data = getSampleMetaData(ct2)
    meta_data = eval(substitute(
        dplyr::mutate(meta_data,
                      eval(parse(text = mutate_expression)),
                      .by = .by,
                      .keep = .keep,
                      .before = .before,
                      .after = .after)
    ))
    k = grepl("eval.parse.text", colnames(meta_data))
    if(sum(k) != 1){
        stop("Something has gone wrong evaluating the supplied expression. There may be something screwy with supplied sample metadata. If not, please report this issue.")
    }
    colnames(meta_data)[k] = mutate_name

    ct2 = setSampleMetaData(ct2, new_meta = meta_data, silent = TRUE)

    history_item = list(mutateSamples  = list(FUN = mutateSamples , ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)

    ct2
}

#' mutateRegions
#'
#' Should work similarly to [dplyr::mutate] but applied to rowData/region metadata of ChIPtsne2 objects.
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param mutate_name  Name of new variable created by `mutate_expression`.
#' @param mutate_expression Expression to derive new variable values.
#' @param .by See [dplyr::mutate]
#' @param .keep  See [dplyr::mutate]
#' @param .before  See [dplyr::mutate]
#' @param .after  See [dplyr::mutate]
#'
#' @return A `r doc_ct2_nrr()` with modified rowData/region metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' rowData(ct2)
#' getRegionMetaData(ct2)
#' ct2 = mutateRegions(ct2, "either_10a_or_at1", peak_MCF10A_CTCF | peak_MCF10AT1_CTCF)
#' rowData(ct2)
#' getRegionMetaData(ct2)
mutateRegions = function(
        ct2,
        mutate_name,
        mutate_expression,
        .by = NULL,
        .keep = c("all", "used", "unused", "none")[1],
        .before = NULL,
        .after = NULL){
    #because we can't store an expression, we need to convert to character for history
    test_expr = substitute(mutate_expression)
    if(is.call(test_expr)){
        mutate_expression = deparse(test_expr)
    }
    remove("test_expr")

    message("mutateRegions ...")
    args = get_args()

    meta_data = getRegionMetaData(ct2, include_value_max = TRUE)
    meta_data = eval(substitute(
        dplyr::mutate(
            meta_data,
            eval(parse(text = mutate_expression)),
            .by = .by,
            .keep = .keep,
            .before = .before,
            .after = .after)
    ))
    k = grepl("eval.parse.text", colnames(meta_data))
    if(sum(k) != 1){
        stop("Something has gone wrong evaluating the supplied expression. There may be something screwy with supplied region metadata. If not, please report this issue.")
    }
    colnames(meta_data)[k] = mutate_name

    ct2 = setRegionMetaData(ct2, new_meta = meta_data[, c(colnames(getRegionMetaData(ct2)), mutate_name)], silent = TRUE)

    history_item = list(mutateRegions  = list(FUN = mutateRegions , ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)

    ct2
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
#' ct2 = mutateSamples(ct2, "cell_mark", paste(cell, mark))
#' ct2 = separateSamples(ct2, "cell_mark", c("v1", "v2"), sep = " ", remove = FALSE)
#' colData(ct2)
#' getSampleMetaData(ct2)
separateSamples = function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
                           convert = FALSE, extra = "warn", fill = "warn", ...){
    message("separateSamples ...")
    args = get_args()
    new_meta_data = tidyr::separate(
        as.data.frame(getSampleMetaData(data)),
        col = col,
        into = into,
        sep = sep,
        remove = remove,
        convert = convert,
        extra = extra,
        fill = fill,
        ...)
    data = setSampleMetaData(data, new_meta_data, silent = TRUE)

    history_item = list(separateSamples  = list(FUN = separateSamples , ARG = args))
    data@metadata = c(ChIPtsne2.history(data), history_item)

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
#' ct2 = mutateRegions(ct2,
#'   mutate_name = "either_10a_or_at1",
#'   mutate_expression = paste(peak_MCF10A_CTCF, id))
#' ct2 = separateRegions(ct2, "either_10a_or_at1", c("a", "b"), sep = " ")
#' rowData(ct2)
#' getRegionMetaData(ct2)
separateRegions = function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
                           convert = FALSE, extra = "warn", fill = "warn", ...){
    message("separateRegions ...")
    args = get_args()
    new_meta_data = tidyr::separate(
        as.data.frame(getRegionMetaData(data)),
        col = col,
        into = into,
        sep = sep,
        remove = remove,
        convert = convert,
        extra = extra,
        fill = fill,
        ...)
    data = setRegionMetaData(data, new_meta = new_meta_data, silent = TRUE)

    history_item = list(separateRegions  = list(FUN = separateRegions , ARG = args))
    data@metadata = c(ChIPtsne2.history(data), history_item)

    data
}

