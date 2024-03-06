
#'
#' @importFrom dplyr select
#' @importFrom data.table as.data.table
#' @importFrom reshape2 melt
.getTidyProfile.with_meta = function(ct2,
                                   sample_meta_VARS = NULL,
                                   region_meta_VARS = NULL){
    df = do.call(rbind,
                 lapply(names(ct2@colToRowMatCols), function(nam){
                     i = ct2@colToRowMatCols[[nam]]
                     df = reshape2::melt(ct2@rowToRowMat[, i, drop = FALSE])
                     df$Var2 = factor(df$Var2)
                     levels(df$Var2) = as.numeric(sub(paste0(nam, "_"), "", levels(df$Var2), fixed = TRUE))
                     df[[ct2@name_VAR]] = nam
                     df
                 })
    )
    colnames(df) = c(ct2@region_VAR, ct2@position_VAR, ct2@value_VAR, ct2@name_VAR)
    df[[ct2@region_VAR]] = factor(df[[ct2@region_VAR]], levels = rownames(ct2))
    df[[ct2@position_VAR]] = as.numeric(as.character(df[[ct2@position_VAR]]))
    df[[ct2@name_VAR]] = factor(df[[ct2@name_VAR]], levels = colnames(ct2))
    if(!is.null(sample_meta_VARS) & length(sample_meta_VARS) > 0){
        if(length(sample_meta_VARS) == 1 && sample_meta_VARS == TRUE){
            sample_meta_VARS = colnames(colData(ct2))
        }
        meta_df = colData(ct2) %>%
            as.data.frame %>%
            dplyr::select(dplyr::all_of(sample_meta_VARS))
        meta_df[[ct2@name_VAR]] = rownames(meta_df)
        df = merge(meta_df, df, by = ct2@name_VAR)
    }
    if(!is.null(region_meta_VARS) & length(region_meta_VARS) > 0){
        if(length(region_meta_VARS) == 1 && region_meta_VARS == TRUE){
            region_meta_VARS = colnames(rowData(ct2))
        }
        reg_df = getRegionMetaData(ct2) %>%
            dplyr::select(dplyr::all_of(c(ct2@region_VAR, region_meta_VARS)))
        df = merge(df, reg_df, by = ct2@region_VAR)
    }
    df[[ct2@name_VAR]] = factor(df[[ct2@name_VAR]], levels = colnames(ct2))
    data.table::as.data.table(df)
}

setGeneric("getTidyProfile.with_meta", function(ct2, sample_meta_VARS = NULL, region_meta_VARS = NULL) standardGeneric("getTidyProfile.with_meta"), signature = "ct2")

setMethod("getTidyProfile.with_meta", c("ChIPtsne2_no_rowRanges"), .getTidyProfile.with_meta)

.getTidyProfile = function(ct2, meta_VARS = NULL){
    if(length(meta_VARS) == 1 && meta_VARS == TRUE){
        sample_meta_VARS = TRUE
        region_meta_VARS = TRUE
    }else{
        missed = setdiff(meta_VARS, c(colnames(getSampleMetaData(ct2)), colnames(getRegionMetaData(ct2))))
        if(length(missed) > 0){
            stop(paste0("Invalid meta_VARS specified. Not found in region or sample metadata:\n",
                        paste(missed, collapse = "\n")))
        }
        sample_meta_VARS = intersect(meta_VARS, colnames(getSampleMetaData(ct2)))
        sample_meta_VARS = setdiff(sample_meta_VARS, ct2@name_VAR)
        region_meta_VARS = intersect(meta_VARS, colnames(getRegionMetaData(ct2)))
        if(length(intersect(sample_meta_VARS, region_meta_VARS)) > 0){
            stop("Ambiguous meta_VARS between sample and region.")
        }
    }

    .getTidyProfile.with_meta(ct2 = ct2,
                          sample_meta_VARS = sample_meta_VARS,
                          region_meta_VARS = region_meta_VARS)
}

#' @export
setGeneric("getTidyProfile", function(ct2, meta_VARS = NULL) standardGeneric("getTidyProfile"), signature = "ct2")

#' @export
setMethod("getTidyProfile", c("ChIPtsne2_no_rowRanges"), .getTidyProfile)

