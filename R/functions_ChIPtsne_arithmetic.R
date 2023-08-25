apply_chiptsne2_operator = function(e1, e2, operator = "-"){
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
    is_diff = sapply(seq(ncol(cd1)), function(i){
        cd1[,i] != cd2[,i]
    })
    new_colData = cd1
    subtract_names = colnames(cd1)[is_diff]
    for(i in which(is_diff)){
        if(length(unique(cd1[, i])) != 1){
            stop("Error in e1, too many values. Use subset to select single value.")
        }
        if(length(unique(cd2[, i])) != 1){
            stop("Error in e2, too many values. Use subset to select single value.")
        }
        new_colData[, i] = paste(cd1[, i], operator, cd2[,i])
    }
    rownames(new_colData) = paste(rownames(cd1), operator, rownames(cd2))

    #generate new ChIPtsne2 internal data structures
    new_rowToRowMat = get(operator)(rowToRowMat(e1), rowToRowMat(e2))
    new_colToRowMatCols = colToRowMatCols(e1)

    #update colnames of new_rowToRowMat and names/values in colToRowMatCols
    for(i in seq(nrow(cd1))){
        old_rn = rownames(cd1)[i]
        new_rn = rownames(new_colData)[i]
        new_colToRowMatCols[[old_rn]] = sub(old_rn, new_rn, new_colToRowMatCols[[old_rn]])
        names(new_colToRowMatCols)[names(new_colToRowMatCols) == old_rn] = new_rn
    }
    new_cn = unlist(new_colToRowMatCols)
    names(new_cn) = NULL
    colnames(new_rowToRowMat) = new_cn

    new_prof_max_mat = .recalculateMax(new_rowToRowMat, new_colToRowMatCols)
    init_history = list(birthday = date(), session_info = sessionInfo(), chiptsne2_version = utils::packageDescription("chiptsne2")$Version)
    ct2_result = ChIPtsne2(assay = list(max = new_prof_max_mat),
                           rowRanges = rowRanges(e1),
                           rowToRowMat = new_rowToRowMat,
                           colToRowMatCols = new_colToRowMatCols,
                           colData = new_colData,
                           name_VAR = e1@name_VAR,
                           position_VAR = e1@position_VAR,
                           value_VAR = e1@value_VAR,
                           region_VAR = e1@region_VAR,
                           fetch_config = e1@fetch_config,
                           metadata = init_history)
    ct2_result@metadata = c(ct2_result@metadata, new_history)
    ct2_result
}
apply_chiptsne2_operator.num = function(e1, e2, operator = "-"){
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

    new_prof_max_mat = .recalculateMax(new_rowToRowMat, new_colToRowMatCols)

    init_history = list(birthday = date(), session_info = sessionInfo(), chiptsne2_version = utils::packageDescription("chiptsne2")$Version)
    ct2_result = ChIPtsne2(assay = list(max = new_prof_max_mat),
                           rowRanges = rowRanges(e1),
                           rowToRowMat = new_rowToRowMat,
                           colToRowMatCols = new_colToRowMatCols,
                           colData = new_colData,
                           name_VAR = e1@name_VAR,
                           position_VAR = e1@position_VAR,
                           value_VAR = e1@value_VAR,
                           region_VAR = e1@region_VAR,
                           fetch_config = e1@fetch_config,
                           metadata = init_history)
    ct2_result@metadata = c(ct2_result@metadata, new_history)
    ct2_result
}

setMethod("-", signature = c("ChIPtsne2", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator(e1, e2, "-")
})
setMethod("+", signature = c("ChIPtsne2", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator(e1, e2, "+")
})
setMethod("/", signature = c("ChIPtsne2", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator(e1, e2, "/")
})
setMethod("*", signature = c("ChIPtsne2", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator(e1, e2, "*")
})

setMethod("-", signature = c("ChIPtsne2", "numeric"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "-")
})
setMethod("+", signature = c("ChIPtsne2", "numeric"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "+")
})
setMethod("/", signature = c("ChIPtsne2", "numeric"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "/")
})
setMethod("*", signature = c("ChIPtsne2", "numeric"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "*")
})

setMethod("-", signature = c("numeric", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "-")
})
setMethod("+", signature = c("numeric", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "+")
})
setMethod("/", signature = c("numeric", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "/")
})
setMethod("*", signature = c("numeric", "ChIPtsne2"), definition = function(e1, e2){
    apply_chiptsne2_operator.num(e1, e2, "*")
})

