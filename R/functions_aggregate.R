
#' aggregateRegionsByGroup
#'
#' @param ct2 `r doc_ct2()`
#' @param group_VAR Attribute name to aggregate regions to.  There will be 1 meta-region per unique entry in `group_VAR`. `group_VAR` may specify multiple attributes, in which case there will 1 meta-region per combination of entries in all `group_VAR`.
#' @param new_region_VAR The new region variable of the resulting ChIPtsne2 object. Defaults to `group_VAR` for single `group_VAR` but then defaults to "grouped_regions" if there are multiple.
#' @param sep Works just like `sep` parameter for base::paste. Used to combine multiple group_VAR values.
#'
#' @return A ChIPtsne2_no_rowRanges object with meta-regions for combinations of `group_VAR` values.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#'
#' ct2.agg1 = aggregateRegionsByGroup(ct2, "peak_MCF10A_CTCF")
#' rowData(ct2.agg1)
#' getRegionVariable(ct2.agg1)
#'
#' ct2.agg2 = aggregateRegionsByGroup(
#'   ct2,
#'   c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"),
#'   new_region_VAR = "peak_overlap"
#' )
#' rowData(ct2.agg2)
#' getRegionVariable(ct2.agg2)
aggregateRegionsByGroup = function(ct2, group_VAR, new_region_VAR = ifelse(length(group_VAR) == 1, group_VAR, "grouped_regions"), sep = " "){
    centroid = calculateGroupCentroid(ct2, group_VAR, sep = sep)

    df = do.call(rbind,
                 lapply(names(ct2@colToRowMatCols), function(nam){
                     i = ct2@colToRowMatCols[[nam]]
                     df = reshape2::melt(centroid[, i, drop = FALSE])
                     df$Var2 = factor(df$Var2)
                     levels(df$Var2) = as.numeric(sub(paste0(nam, "_"), "", levels(df$Var2), fixed = TRUE))
                     df[[ct2@name_VAR]] = nam
                     df
                 })
    )
    colnames(df) = c(new_region_VAR, ct2@position_VAR, ct2@value_VAR, ct2@name_VAR)
    df[[new_region_VAR]] = factor(df[[new_region_VAR]], levels = rownames(centroid))

    rd = unique(rowData(ct2)[, group_VAR, drop = FALSE])
    rownames(rd) = apply(as.data.frame(rd), 1, paste, collapse = sep)
    rd[[new_region_VAR]] = rownames(rd)

    ct2.meta = ChIPtsne2.from_tidy(
        df,
        query_gr = NULL,
        sample_metadata = colData(ct2),
        region_metadata = rd,
        position_VAR = ct2@position_VAR,
        name_VAR = ct2@name_VAR,
        value_VAR = ct2@value_VAR,
        region_VAR = new_region_VAR)
    ct2.meta
}

#' aggregateSamplesByGroup
#'
#' @param ct2 `r doc_ct2()`
#' @param group_VAR Attribute name to aggregate samples to.  There will be 1 meta-sample per unique entry in `group_VAR`. `group_VAR` may specify multiple attributes, in which case there will 1 meta-region per combination of entries in all `group_VAR`.
#' @param new_name_VAR The new name variable of the resulting ChIPtsne2 object. Defaults to `group_VAR` for single `group_VAR` but then defaults to "grouped_regions" if there are multiple.
#' @param sep Works just like `sep` parameter for base::paste. Used to combine multiple group_VAR values.
#'
#' @return A ChIPtsne2_no_rowRanges object with meta-regions for combinations of `group_VAR` values.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' colData(ct2)
#'
#' ct2.agg1 = aggregateSamplesByGroup(ct2, "mark")
#' colData(ct2.agg1)
#' getNameVariable(ct2.agg1)
#'
#' #pretend that we have reps to aggregate
#' ct2.r1 = exampleChIPtsne2.with_meta()
#' colData(ct2.r1)$rep = "rep1"
#' colnames(ct2.r1) = paste0(colnames(ct2.r1), "_rep1")
#' ct2.r2 = exampleChIPtsne2.with_meta()
#' colData(ct2.r2)$rep = "rep2"
#' colnames(ct2.r2) = paste0(colnames(ct2.r2), "_rep2")
#' ct2.reps = cbind(ct2.r1, ct2.r2)
#' ct2.agg2 = aggregateSamplesByGroup(ct2.reps,  c("cell", "mark"), "group")
#' colData(ct2.agg2)
#' getNameVariable(ct2.agg2)
#'
#' ct2.agg3 = aggregateSamplesByGroup(ct2.reps, c("cell", "mark"))
#' colData(ct2.agg3)
#' getNameVariable(ct2.agg3)
aggregateSamplesByGroup = function(ct2, group_VAR, new_name_VAR = ifelse(length(group_VAR) == 1, group_VAR, "group_name"), sep = " "){
    .validate_allowed_input(group_VAR, colnames(colData(ct2)), "Some values of group_VAR are not present in colData:")
    # carried_VARS = unique(c(group_VAR, new_name_VAR))
    carried_VARS = setdiff(unique(c(group_VAR)), new_name_VAR)
    ct2@colData[[new_name_VAR]] =  apply(
        as.data.frame(colData(ct2))[, group_VAR, drop = FALSE],
        1,
        function(x){
            paste(as.character(x), collapse = sep)
        }
    )
    ct2.sp = split(ct2, new_name_VAR)
    ct2.sp = lapply(ct2.sp, function(x){rowData(x) = NULL; x})
    ct2.parts = list()
    for(name in names(ct2.sp)){
        x = ct2.sp[[name]]
        if(ncol(x) == 1){
            ct2.new = x
        }else{
            ct2.new = x[,1]
            for(i in seq(2, ncol(x))){
                ct2.new = ct2.new + x[,i]
            }
            ct2.new = ct2.new / ncol(x)
        }
        colnames(ct2.new) = name

        ct2.new@colData = ct2.new@colData[, carried_VARS, drop = FALSE]
        ct2.parts[[name]] = ct2.new
    }
    ct2.meta = do.call(cbind, ct2.parts)
    rowData(ct2.meta) = rowData(ct2)
    ct2.meta = setNameVariable(ct2.meta, new_name_VAR)
    ct2.meta
}

#' aggregateByGroup
#'
#' @param ct2 `r doc_ct2()`
#' @param group_VARS Attribute name to aggregate samples and/or regions to.  There will be 1 meta-sample per unique sampled entry in `group_VARS` and 1 meta-region per unique region entry. `group_VARS` may specify multiple attributes, in which case there will 1 meta-region and/or sample per combination of entries in all `group_VARS`.
#' @param sample_sep Separator character to use when pasting multiple sample grouping variables.
#' @param region_sep Separator character to use when pasting multiple region grouping variables.
#'
#' @return Either a ChIPtsne2 or ChIPtsne2_no_rowRanges object, with meta-regions for combinations of `group_VARS` values if region aggregation occurred.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#'
#' aggregateRegionsByGroup(ct2, c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
#' aggregateSamplesByGroup(ct2, c("mark"))
#' aggregateByGroup(ct2, c("mark", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
#'
#' ct2.r1 = exampleChIPtsne2.with_meta()
#' colData(ct2.r1)$rep = "rep1"
#' colnames(ct2.r1) = paste0(colnames(ct2.r1), "_rep1")
#' ct2.r2 = exampleChIPtsne2.with_meta()
#' colData(ct2.r2)$rep = "rep2"
#' colnames(ct2.r2) = paste0(colnames(ct2.r2), "_rep2")
#' ct2.reps = cbind(ct2.r1, ct2.r2)
#' ct2.agg3 = aggregateByGroup(ct2.reps, c("cell", "mark", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
#' colData(ct2.agg3)
#' rowData(ct2.agg3)
aggregateByGroup = function(ct2, group_VARS, sample_sep = " ", region_sep = " "){
    group_VARS.col = group_VARS[group_VARS %in% colnames(colData(ct2))]
    group_VARS.row = group_VARS[group_VARS %in% colnames(rowData(ct2))]
    if(length(group_VARS.col) > 0){
        ct2.meta = aggregateSamplesByGroup(ct2, group_VARS.col, sep = sample_sep)
    }else{
        ct2.meta = ct2
    }
    if(length(group_VARS.row) > 0){
        ct2.meta = aggregateRegionsByGroup(ct2.meta, group_VARS.row, sep = region_sep)
    }

    ct2.meta
}
