#### name_VAR ####


#' @export
#' @rdname ct2-getset
setMethod("setNameVariable", c("ChIPtsne2_no_rowRanges"), function(ct2, new_name_VAR){
    args = get_args()
    history_item = list(setNameVariable = list(FUN = setNameVariable, ARG = args))
    ct2@metadata = c(ct2@metadata, history_item)
    if(new_name_VAR %in% colnames(colData(ct2))){
        .Deprecated("swapNameVariable", msg = "Specified sample metadata variable already exists. In the future, swapping to this variable must be explicit using 'swapNameVariable'.")
        cn = colData(ct2)[[new_name_VAR]]
        if(any(duplicated(cn))){
            bad_cn = unique(cn[duplicated(cn)])
            stop("New name variable is already present in colData(ct2). This is only allowed when all values are unique. Offending values: ", paste(bad_cn, collapse = ", "))
        }else{
            #swap old and new
            old_name_VAR = ct2@name_VAR
            ct2 = .update_ct2_colnames(ct2, old_name_VAR = old_name_VAR, new_name_VAR = new_name_VAR)
        }
    }else{
        #no special considerations if name_VAR isn't otherwise present.
    }
    ct2@name_VAR = new_name_VAR
    ct2
})

#' @export
#' @rdname ct2-getset
setMethod("swapNameVariable", c("ChIPtsne2_no_rowRanges"), function(ct2, new_name_VAR){
    args = get_args()
    history_item = list(swapNameVariable = list(FUN = swapNameVariable, ARG = args))
    ct2@metadata = c(ct2@metadata, history_item)
    if(new_name_VAR %in% colnames(colData(ct2))){
        cn = colData(ct2)[[new_name_VAR]]
        if(any(duplicated(cn))){
            bad_cn = unique(cn[duplicated(cn)])
            stop("All values of new name variable are not unique. Offending values: ", paste(bad_cn, collapse = ", "))
        }else{
            #swap old and new
            old_name_VAR = ct2@name_VAR
            ct2 = .update_ct2_colnames(ct2, old_name_VAR = old_name_VAR, new_name_VAR = new_name_VAR)
        }
    }else{
        stop("New name variable is not present in colData(ct2). Specify an existing variable or did you mean to use setNameVariable?")
    }
    ct2@name_VAR = new_name_VAR
    ct2
})

.update_ct2_rownames = function(ct2, new_names = NULL, old_name_VAR = NULL, new_name_VAR = NULL){
    old_names = rownames(ct2)
    if(is.null(new_names) & is.null(new_name_VAR)){
        stop("One of new_names or new_name_VAR required.")
    }
    if(!is.null(new_names) & !is.null(new_name_VAR)){
        stop("Only one of new_names or new_name_VAR is allowed.")
    }
    if(!is.null(old_name_VAR)){
        rowData(ct2)[[old_name_VAR]] = rownames(ct2)
    }
    if(!is.null(new_name_VAR)){
        new_names = rowData(ct2)[[new_name_VAR]]
        rowData(ct2)[[new_name_VAR]] = NULL
    }
    names(new_names) = old_names
    #updating internal data
    r2rm = rowToRowMat(ct2)
    new_rn = new_names[rownames(r2rm)]
    names(new_rn) = NULL

    rownames(r2rm) = new_rn
    rowToRowMat(ct2) = r2rm

    old_assays = ct2@assays@data
    for(i in seq_along(old_assays)){
        rownames(old_assays[[i]]) = new_rn
    }
    ct2@assays = SummarizedExperiment::Assays(old_assays)


    if(is(ct2, "ChIPtsne2")){#unclear if equivalent operation of ChIPtsne2_no_rowRanges required
        names(ct2@rowRanges) = new_rn
    }else{
        rownames(ct2@elementMetadata) = new_rn
        ct2@NAMES = new_rn
    }
    if(!is.null(new_name_VAR)){
        ct2@region_VAR = new_name_VAR
    }
    ct2
}

.update_ct2_colnames = function(ct2, new_names = NULL, old_name_VAR = NULL, new_name_VAR = NULL){
    old_names = rownames(colData(ct2))
    if(is.null(new_names) & is.null(new_name_VAR)){
        stop("One of new_names or new_name_VAR required.")
    }
    if(!is.null(new_names) & !is.null(new_name_VAR)){
        stop("Only one of new_names or new_name_VAR is allowed.")
    }
    if(!is.null(old_name_VAR)){
        ct2@colData[[old_name_VAR]] = rownames(colData(ct2))
    }
    if(!is.null(new_name_VAR)){
        rownames(ct2@colData) = colData(ct2)[[new_name_VAR]]
        new_names = rownames(ct2@colData)
        ct2@colData[[new_name_VAR]] = NULL
    }
    names(new_names) = old_names
    #updating internal data
    r2rm = rowToRowMat(ct2)
    c2rmc = colToRowMatCols(ct2)
    if(is.null(new_name_VAR)){
        rownames(colData(ct2)) = new_names[rownames(colData(ct2))]
    }
    all_old_cn = unlist(c2rmc)
    for(nam in names(c2rmc)){
        match_names = c2rmc[[nam]]
        i = which(old_names == nam)
        new_nam = new_names[i]
        c2rmc[[nam]] = sub(nam, new_nam, match_names, fixed = TRUE)
    }
    all_new_cn = unlist(c2rmc)
    names(all_new_cn) = all_old_cn
    #assigning named vector seems problematic so go through tmp and remove names
    tmp = all_new_cn[colnames(r2rm)]
    names(tmp) = NULL
    colnames(r2rm) = tmp
    tmp = new_names[names(c2rmc)]
    names(tmp) = NULL
    names(c2rmc) = tmp

    old_assays = ct2@assays@data
    for(i in seq_along(old_assays)){
        new_cn = new_names[colnames(old_assays[[i]])]
        names(new_cn) = NULL
        colnames(old_assays[[i]]) = new_cn
    }
    ct2@assays = SummarizedExperiment::Assays(old_assays)

    rowToRowMat(ct2) = r2rm
    colToRowMatCols(ct2) = c2rmc
    if(!is.null(new_name_VAR)){
        ct2@name_VAR = new_name_VAR
    }
    ct2
}



#' @export
setMethod("getNameVariable", c("ChIPtsne2_no_rowRanges"), function(ct2){
    ct2@name_VAR
})



#### value_VAR ####

#' @export
#' @rdname ct2-getset
setMethod("setValueVariable", c("ChIPtsne2_no_rowRanges"), function(ct2, new_value_VAR){
    args = get_args()
    history_item = list(setValueVariable = list(FUN = setValueVariable, ARG = args))
    ct2@metadata = c(ct2@metadata, history_item)
    ct2@value_VAR = new_value_VAR
    ct2
})



#' @export
#' @rdname ct2-getset
setMethod("getValueVariable", c("ChIPtsne2_no_rowRanges"), function(ct2){
    ct2@value_VAR
})

#### region_VAR ####


#' @export
#' @rdname ct2-getset
setMethod("setRegionVariable", c("ChIPtsne2_no_rowRanges"), function(ct2, new_region_VAR){
    args = get_args()
    history_item = list(setRegionVariable = list(FUN = setRegionVariable, ARG = args))
    ct2@metadata = c(ct2@metadata, history_item)
    ct2@region_VAR = new_region_VAR
    ct2
})



#' @export
#' @rdname ct2-getset
setMethod("getRegionVariable", c("ChIPtsne2_no_rowRanges"), function(ct2){ct2@region_VAR})

#### position_VAR ####


#' @export
#' @rdname ct2-getset
setMethod("setPositionVariable", c("ChIPtsne2_no_rowRanges"), function(ct2, new_position_VAR){
    args = get_args()
    history_item = list(setPositionVariable = list(FUN = setPositionVariable, ARG = args))
    ct2@metadata = c(ct2@metadata, history_item)
    ct2@position_VAR = new_position_VAR
    ct2
})


#' @export
#' @rdname ct2-getset
setMethod("getPositionVariable", c("ChIPtsne2_no_rowRanges"), function(ct2){ct2@position_VAR})

#### SampleMetaData ####

#' getSampleMetaData
#'
#' @param ct2 A ChIPtsne object
#' @param select_VARS character vector of variables to select from sample
#'   metadata. Default of NULL will select all available sample metadata
#'   variables.
#'
#' @return data.frame with sample meta data, similar to colData but suitable for
#'   tidyverse operations.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' getSampleMetaData(ct2)
#' getSampleMetaData(ct2, "sample")
getSampleMetaData = function(ct2, select_VARS = NULL){
    cd = colData(ct2)
    df = as.data.frame(cd)
    df[[ct2@name_VAR]] = factor(rownames(cd), levels = colnames(ct2))
    if(!is.null(select_VARS)){
        if(!all(select_VARS %in% colnames(df))){
            stop(
                paste(collapse = "\n",
                      c(
                          "select_VARS not found in sample metadata:",
                          setdiff(select_VARS, colnames(df))
                      )
                )
            )
        }
        df = df[, union(ct2@name_VAR, select_VARS), drop = FALSE]
    }
    df
}

#' setSampleMetaData
#'
#' @param ct2 A ChIPtsne object
#' @param new_meta A data.frame with new metadata information. Must include same name_VAR as ct2 or have equivalent rownames. Variables already present in ct2 will result in overiting those variables.
#'
#' @return A modified ChIPtsne2 object with added/overwritten sample metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' getSampleMetaData(ct2)
#' new_meta = data.frame(sample = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF"), id = seq(3))
#' ct2 = setSampleMetaData(ct2, new_meta)
#' getSampleMetaData(ct2)
#'
#' #metadata may be overriden
#' new_meta2 = data.frame(id = LETTERS[seq(3)], id2 = LETTERS[seq(3)])
#' rownames(new_meta2) = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF")
#' ct2 = setSampleMetaData(ct2, new_meta2)
#' getSampleMetaData(ct2)
setSampleMetaData = function(ct2, new_meta){
    message("setSampleMetaData ...")
    args = get_args()
    cd = getSampleMetaData(ct2)
    if(!ct2@name_VAR %in% colnames(new_meta)){
        if(is.null(rownames(new_meta))){
            stop("new_meta must contain name_VAR or have rownames.")
        }else{
            new_meta[[ct2@name_VAR]] = rownames(new_meta)
        }
    }

    if(!setequal(new_meta[[ct2@name_VAR]], rownames(cd))){
        stop(paste(sep = "\n",
                   "name_VAR is not equivalent in new metadata.",
                   paste(c("Extra entries in new_meta:", setdiff(new_meta[[ct2@name_VAR]], rownames(cd))), collapse = "\n"),
                   paste(c("Missing entries from new_meta:", setdiff(rownames(cd), new_meta[[ct2@name_VAR]])), collapse = "\n")
        ))
    }
    retained_cn = setdiff(colnames(cd),
                          setdiff(colnames(new_meta), ct2@name_VAR)
    )
    new_cd = merge(cd[, retained_cn, drop = FALSE], new_meta, by = ct2@name_VAR)
    rownames(new_cd) = new_cd[[ct2@name_VAR]]
    new_cd[[ct2@name_VAR]] = NULL
    new_cd = S4Vectors::DataFrame(new_cd)

    ct2@colData = new_cd
    history_item = list(setSampleMetaData = list(FUN = setSampleMetaData, ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)
    ct2
}

#### RegionMetaData ####

#### Getters ####

.getRegionMetaData = function(ct2, select_VARS = NULL){
    if(is(ct2, "ChIPtsne2")){
        gr = rowRanges(ct2)
        df = GenomicRanges::mcols(gr) %>% data.frame(check.names = FALSE)
        df[[ct2@region_VAR]] = names(gr)
    }else if (is(ct2, "ChIPtsne2_no_rowRanges")){
        df = rowData(ct2) %>% data.frame(check.names = FALSE)
        df[[getRegionVariable(ct2)]] = rownames(df)

    }
    df = df[, c(getRegionVariable(ct2), setdiff(colnames(df), getRegionVariable(ct2))), drop = FALSE]
    if(!is.null(select_VARS)){
        if(!all(select_VARS %in% colnames(df))){
            stop(
                paste(collapse = "\n",
                      c(
                          "select_VARS not found in region metadata:",
                          setdiff(select_VARS, colnames(df))
                      )
                )
            )
        }
        df = df[, union(ct2@region_VAR, select_VARS), drop = FALSE]
    }
    df
}

#' getRegionMetaData
#'
#' @param ct2 A ChIPtsne object
#' @param select_VARS character vector of variables to select from region
#'   metadata. Default of NULL will select all available region metadata
#'   variables.
#'
#' @return data.frame with region meta data, similar to rowRanges but suitable
#'   for tidyverse operations.
#'
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' getRegionMetaData(ct2)
#' getRegionMetaData(ct2, c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
setMethod("getRegionMetaData", "ChIPtsne2_no_rowRanges", .getRegionMetaData)

# internals of setRegionMetaData
# does not impact history
setRegionMetaData.no_history = function(ct2, new_meta){
    cd = getRegionMetaData(ct2)
    if(!ct2@region_VAR %in% colnames(new_meta)){
        if(is.null(rownames(new_meta))){
            stop("new_meta must contain region_VAR or have rownames.")
        }else{
            new_meta[[ct2@region_VAR]] = rownames(new_meta)
        }
    }

    if(!setequal(new_meta[[ct2@region_VAR]], rownames(cd))){
        stop(paste(sep = "\n",
                   "region_VAR is not equivalent in new metadata.",
                   paste(c("Extra entries in new_meta:", setdiff(new_meta[[ct2@region_VAR]], rownames(cd))), collapse = "\n"),
                   paste(c("Missing entries from new_meta:", setdiff(rownames(cd), new_meta[[ct2@region_VAR]])), collapse = "\n")
        ))
    }
    retained_cn = setdiff(colnames(cd),
                          setdiff(colnames(new_meta), ct2@region_VAR)
    )
    new_cd = merge(cd[, retained_cn, drop = FALSE], new_meta, by = ct2@region_VAR)
    if(!is(ct2, "ChIPtsne2")){# this should be a generic method dispatch
        overwrite = TRUE
        old_rd = rowData(ct2)
        new_cd
        conflicting_cn = intersect(colnames(old_rd), colnames(new_cd))
        if(length(conflicting_cn) > 0){
            if(overwrite){
                old_rd = old_rd[, setdiff(colnames(old_rd), conflicting_cn), drop = FALSE]
            }else{
                stop(paste(c("Conflicting colnames in region_metadata already present in query_gr:", conflicting_cn), collapse = "\n"))
            }
        }
        new_rd = cbind(
            old_rd,
            as.data.frame(new_cd %>% dplyr::select(!dplyr::all_of(c(ct2@region_VAR))))
        )
        rowData(ct2) = new_rd
    }else{
        new_gr = .add_region_metadata(rowRanges(ct2), region_metadata = new_cd, region_VAR = ct2@region_VAR, overwrite = TRUE)
        ct2@rowRanges = new_gr
    }
    ct2
}

#' setRegionMetaData
#'
#' @param ct2 A ChIPtsne object
#' @param new_meta A data.frame with new metadata information. Must include same
#'   name_VAR as ct2 or have equivalent rownames. Variables already present in
#'   ct2 will result in overiting those variables.
#' @param silent If TRUE, no messages will be displayed. Default is FALSE.
#'
#' @return A modified ChIPtsne2 object with added/overwritten sample metadata.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' new_meta = getRegionMetaData(ct2)
#' new_meta[["10A_and_AT1"]] = ifelse(
#'   new_meta$peak_MCF10A_CTCF & new_meta$peak_MCF10AT1_CTCF,
#'   "yes",
#'   "no"
#' )
#' ct2 = setRegionMetaData(ct2, new_meta)
#' getRegionMetaData(ct2)
#'
#' ct2.nrr = exampleChIPtsne2.with_meta()
#' rowRanges(ct2.nrr) = NULL
#' getRegionMetaData(ct2.nrr)
#' ct2.nrr = setRegionMetaData(ct2.nrr, new_meta)
#' getRegionMetaData(ct2.nrr)
setRegionMetaData = function(ct2, new_meta, silent = FALSE){
    if(!silent){
        message("setRegionMetaData ...")
    }
    if(!is(ct2, "ChIPtsne2_no_rowRanges")){
        stop("ct2 must be a ChIPtsne2 object")
    }
    args = get_args()
    ct2 = setRegionMetaData.no_history(ct2, new_meta)
    history_item = list(setRegionMetaData =
                            list(FUN = setRegionMetaData, ARG = args))
    ct2@metadata = c(ChIPtsne2.history(ct2), history_item)
    ct2
}


