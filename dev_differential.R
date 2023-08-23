#method to subtract
#method to filter
library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()
subset(ct2, peak_MCF10CA1_CTCF, cell == "MCF10A")

filter.ChIPtsne2 = function(.data, ...){
    subset(.data, ...)
}

library(dplyr)
filter(ct2, peak_MCF10CA1_CTCF, cell == "MCF10A")

setMethod("-", signature = c("ChIPtsne2", "ChIPtsne2"), definition = function(e1, e2){
    message(colData(e1)$cell, " - ", colData(e2)$cell)
    message("yes")
})

ct2_a = filter(ct2, peak_MCF10CA1_CTCF, cell == "MCF10A")

ct2_b = filter(ct2, peak_MCF10CA1_CTCF, cell == "MCF10AT1")

ct2_a - ct2_b
e1 = ct2_a
e2 = ct2_b

operator = "+"
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

    new_prof_max_mat = get(operator)(assay(e1, "max"), assay(e2, "max"))
    colnames(new_prof_max_mat) = paste(colnames(assay(e1, "max")), operator, colnames(assay(e2, "max")))


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
ct2_sub = apply_chiptsne2_operator(ct2_a, ct2_b, "-")
ct2_plus = apply_chiptsne2_operator(ct2_a, ct2_b, "+")
ct2_div = apply_chiptsne2_operator(ct2_a, ct2_b, "/")
ct2_mult = apply_chiptsne2_operator(ct2_a, ct2_b, "*")
apply_chiptsne2_operator(ct2_a, ct2_b, "asdf")

plotSignalHeatmap(ct2_sub, heatmap_colors = scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-40, 40)))

# debug(getTidyProfile, "ChIPtsne2")
getTidyProfile(ct2_sub)
getTidyProfile(ct2_plus)
getTidyProfile(ct2_div)
getTidyProfile(ct2_mult)
# debug(plotSignalHeatmap, "ChIPtsne2")
plotSignalHeatmap(ct2_plus)
plotSignalHeatmap(ct2_div)
plotSignalHeatmap(ct2_mult)
