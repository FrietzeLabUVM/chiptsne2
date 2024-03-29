

#' ChipTSNE2 arithmetic
#' Arithmetic operators are supports for ChIPtsne2 and ChIPtsne2_no_rowRanges objects.
#'
#' +, -, *, and /
#' @rdname ct2-op
#' @importFrom utils sessionInfo
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = ct2[, c("MCF10A_CTCF", "MCF10AT1_CTCF")]
#' ct2.by_cell = split(ct2, "cell")
#' plotSignalLinePlot(ct2)
#'
#' #operators can be used on 2 ChIPtsne2_no_rowRanges objects or 1 object and a number
#' ct2.diff = ct2.by_cell$MCF10A - ct2.by_cell$MCF10AT1
#' plotSignalLinePlot(ct2.diff)
#'
#' # numeric operators work in either order
#' ct2.2div10a = cbind(2 / (ct2.by_cell$MCF10A+1), ct2.by_cell$MCF10AT1)
#' plotSignalLinePlot(ct2.2div10a)
#' ct2.10adiv2 = cbind(ct2.by_cell$MCF10A / 2, ct2.by_cell$MCF10AT1)
#' plotSignalLinePlot(ct2.10adiv2)
#'
apply_ChIPtsne2_operator = function(e1, e2, operator = "-"){
    if(!all(dim(e1) == dim(e2))){
        stop("ChIPtsne2 objects are not dimensionally compatible. Check dim() of each.")
    }
    if(!all(rowRanges(e1) == rowRanges(e2))){
        stop('ChIPtsne2 objects must have identical rowRanges. Check rowRanges() or each.')
    }
    new_history = list(
        list(
            e1 = e1@metadata,
            e2 = e2@metadata
        )
    )
    switch (operator,
            `-` = {names(new_history) = "subtraction"},
            `+` = {names(new_history) = "addition"},
            `/` = {names(new_history) = "division"},
            `*` = {names(new_history) = "multiplication"},
            stop("unrecognized operation for ChIPtsne2: ", operator)
    )
    #locate comparative variables, validate, generate new values
    cd1 = colData(e1)
    cd2 = colData(e2)
    common_cn = intersect(colnames(cd1), colnames(cd2))
    cd1 = cd1[, common_cn, drop = FALSE]
    cd2 = cd2[, common_cn, drop = FALSE]
    is_numeric = sapply(cd1, is.numeric) & sapply(cd2, is.numeric)
    cd1 = cd1[, !is_numeric, drop = FALSE]
    cd2 = cd2[, !is_numeric, drop = FALSE]

    is_diff = sapply(seq(ncol(cd1)), function(i){
        !(all(cd1[,i] == cd2[,i]))
    })
    new_colData = cd1
    # arithmetic operators make sense to apply in 2 modes
    if(all(colnames(e1) == colnames(e2))){
        # 1) when all name_VAR values are equal
        if(!all(c(
            all(colnames(e1@rowToRowMat) == colnames(e2@rowToRowMat)),
            all(names(e1@colToRowMatCols) == names(e2@colToRowMatCols)),
            all(unlist(e1@colToRowMatCols) == unlist(e2@colToRowMatCols)),
            all(rownames(cd1) == rownames(cd2))
        ))){
            stop("There is an internal issue with column names. Unless you've manipulated slots directly, this is a bug and should be reported.")
        }
        for(i in which(is_diff)){
            new_colData[, i] = paste(cd1[, i], operator, cd2[,i])
        }
    }else{
        # 2) when only unique values are present in colData for both
        for(i in which(is_diff)){
            if(length(unique(cd1[, i])) != 1){
                stop("Error in e1, too many values. Use subset to select single value.")
            }
            if(length(unique(cd2[, i])) != 1){
                stop("Error in e2, too many values. Use subset to select single value.")
            }
            new_colData[, i] = paste(cd1[, i], operator, cd2[,i])
        }
        if(!all(rownames(cd1) == rownames(cd2))){
            rownames(new_colData) = paste(rownames(cd1), operator, rownames(cd2))
        }
    }

    #generate new ChIPtsne2 internal data structures
    new_rowToRowMat = get(operator)(rowToRowMat(e1), rowToRowMat(e2))
    new_colToRowMatCols = colToRowMatCols(e1)

    #update colnames of new_rowToRowMat and names/values in colToRowMatCols
    for(i in seq(nrow(cd1))){
        old_rn = rownames(cd1)[i]
        new_rn = rownames(new_colData)[i]
        new_colToRowMatCols[[old_rn]] = sub(old_rn, new_rn, new_colToRowMatCols[[old_rn]], fixed = TRUE)
        names(new_colToRowMatCols)[names(new_colToRowMatCols) == old_rn] = new_rn
    }
    new_cn = unlist(new_colToRowMatCols)
    names(new_cn) = NULL
    colnames(new_rowToRowMat) = new_cn

    new_prof_max_mat = .recalculateMax(new_rowToRowMat, new_colToRowMatCols)
    init_history = list(birthday = date(), session_info = utils::sessionInfo(), chiptsne2_version = utils::packageDescription("chiptsne2")$Version)
    if("ChIPtsne2_no_rowRanges" %in% c(class(e1), class(e2))){
        ct2_result = ChIPtsne2_no_rowRanges(
            assay = list(max = new_prof_max_mat),
            rowData = rowData(e1),
            rowToRowMat = new_rowToRowMat,
            colToRowMatCols = new_colToRowMatCols,
            colData = new_colData,
            name_VAR = e1@name_VAR,
            position_VAR = e1@position_VAR,
            value_VAR = e1@value_VAR,
            region_VAR = e1@region_VAR,
            fetch_config = e1@fetch_config,
            metadata = init_history
        )
    }else{
        ct2_result = ChIPtsne2(
            assay = list(max = new_prof_max_mat),
            rowRanges = rowRanges(e1),
            rowToRowMat = new_rowToRowMat,
            colToRowMatCols = new_colToRowMatCols,
            colData = new_colData,
            name_VAR = e1@name_VAR,
            position_VAR = e1@position_VAR,
            value_VAR = e1@value_VAR,
            region_VAR = e1@region_VAR,
            fetch_config = e1@fetch_config,
            metadata = init_history
        )
    }
    ct2_result@metadata = c(ct2_result@metadata, new_history)
    ct2_result
}

#' @importFrom utils sessionInfo
#' @rdname ct2-op
apply_ChIPtsne2_operator.num = function(e1, e2, operator = "-"){
    reverse_mode = FALSE
    if(is.numeric(e1)){
        tmp = e1
        e1 = e2
        e2 = tmp
        # subtraction and division are not transitive
        if(operator %in% c("-", "/")){
            reverse_mode = TRUE
        }

    }
    new_history = list(
        list(
            e1 = e1@metadata,
            e2 = e2
        )
    )
    switch (operator,
            `-` = {names(new_history) = "subtraction"},
            `+` = {names(new_history) = "addition"},
            `/` = {names(new_history) = "division"},
            `*` = {names(new_history) = "multiplication"},
            stop("unrecognized operation for ChIPtsne2: ", operator)
    )
    #locate comparative variables, validate, generate new values
    cd1 = colData(e1)
    new_colData = cd1
    rownames(new_colData) = paste0("(", rownames(cd1), ") ", operator, " ", e2)

    #generate new ChIPtsne2 internal data structures
    new_colToRowMatCols = colToRowMatCols(e1)
    if(!reverse_mode){
        new_rowToRowMat = get(operator)(rowToRowMat(e1), e2)
        # new_prof_max_mat = get(operator)(assay(e1, "max"), e2)
    }else{
        if(operator == "-"){
            new_rowToRowMat = (-1 * rowToRowMat(e1)) + e2
            # new_prof_max_mat = (-1 * assay(e1, "max")) + e2
        }else if(operator == "/"){
            new_rowToRowMat = (1 / rowToRowMat(e1)) * e2
            # new_prof_max_mat = (1 / assay(e1, "max")) * e2
        }else{
            stop("Invalid operator/reverse_mode, please report this error.")
        }
    }


    #update colnames of new_rowToRowMat and names/values in colToRowMatCols
    for(i in seq(nrow(cd1))){
        old_rn = rownames(cd1)[i]
        new_rn = rownames(new_colData)[i]
        new_colToRowMatCols[[old_rn]] = sub(old_rn, new_rn, new_colToRowMatCols[[old_rn]], fixed = TRUE)
        names(new_colToRowMatCols)[names(new_colToRowMatCols) == old_rn] = new_rn
    }
    new_cn = unlist(new_colToRowMatCols)
    names(new_cn) = NULL
    colnames(new_rowToRowMat) = new_cn


    new_prof_max_mat = .recalculateMax(new_rowToRowMat, new_colToRowMatCols)



    init_history = list(birthday = date(), session_info = utils::sessionInfo(), chiptsne2_version = utils::packageDescription("chiptsne2")$Version)
    if("ChIPtsne2_no_rowRanges" %in% c(class(e1))){
        ct2_result = ChIPtsne2_no_rowRanges(
            assay = list(max = new_prof_max_mat),
            rowData = rowData(e1),
            rowToRowMat = new_rowToRowMat,
            colToRowMatCols = new_colToRowMatCols,
            colData = new_colData,
            name_VAR = e1@name_VAR,
            position_VAR = e1@position_VAR,
            value_VAR = e1@value_VAR,
            region_VAR = e1@region_VAR,
            fetch_config = e1@fetch_config,
            metadata = init_history
        )
    }else{
        ct2_result = ChIPtsne2(
            assay = list(max = new_prof_max_mat),
            rowRanges = rowRanges(e1),
            rowToRowMat = new_rowToRowMat,
            colToRowMatCols = new_colToRowMatCols,
            colData = new_colData,
            name_VAR = e1@name_VAR,
            position_VAR = e1@position_VAR,
            value_VAR = e1@value_VAR,
            region_VAR = e1@region_VAR,
            fetch_config = e1@fetch_config,
            metadata = init_history
        )
    }
    ct2_result@metadata = c(ct2_result@metadata, new_history)
    ct2_result
}



#' @rdname ct2-op
setMethod("-", signature = c("ChIPtsne2_no_rowRanges", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator(e1, e2, "-")
})
#' @rdname ct2-op
setMethod("+", signature = c("ChIPtsne2_no_rowRanges", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator(e1, e2, "+")
})
#' @rdname ct2-op
setMethod("/", signature = c("ChIPtsne2_no_rowRanges", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator(e1, e2, "/")
})
#' @rdname ct2-op
setMethod("*", signature = c("ChIPtsne2_no_rowRanges", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator(e1, e2, "*")
})
#' @rdname ct2-op
setMethod("-", signature = c("ChIPtsne2_no_rowRanges", "numeric"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "-")
})
#' @rdname ct2-op
setMethod("+", signature = c("ChIPtsne2_no_rowRanges", "numeric"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "+")
})
#' @rdname ct2-op
setMethod("/", signature = c("ChIPtsne2_no_rowRanges", "numeric"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "/")
})
#' @rdname ct2-op
setMethod("*", signature = c("ChIPtsne2_no_rowRanges", "numeric"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "*")
})
#' @rdname ct2-op
setMethod("-", signature = c("numeric", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "-")
})
#' @rdname ct2-op
setMethod("+", signature = c("numeric", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "+")
})
#' @rdname ct2-op
setMethod("/", signature = c("numeric", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "/")
})
#' @rdname ct2-op
setMethod("*", signature = c("numeric", "ChIPtsne2_no_rowRanges"), definition = function(e1, e2){
    apply_ChIPtsne2_operator.num(e1, e2, "*")
})

