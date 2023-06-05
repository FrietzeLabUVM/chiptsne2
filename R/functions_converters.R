#'
#' @importFrom dplyr select
#' @importFrom data.table as.data.table
.prof_dt_from_chiptsne2 = function(ct2,
                                   sample_meta_VARS = NULL,
                                   region_meta_VARS = NULL){
    df = do.call(rbind,
                 lapply(names(ct2@colToRowMatCols), function(nam){
                     i = ct2@colToRowMatCols[[nam]]
                     df = reshape2::melt(ct2@rowToRowMat[, i])
                     df$Var2 = factor(df$Var2)
                     levels(df$Var2) = as.numeric(sub(paste0(nam, "_"), "", levels(df$Var2)))
                     df[[ct2@name_VAR]] = nam
                     df
                 })
    )
    colnames(df) = c(ct2@region_VAR, ct2@position_VAR, ct2@value_VAR, ct2@name_VAR)
    df[[ct2@region_VAR]] = factor(df[[ct2@region_VAR]], levels = names(rowRanges(ct2)))
    df[[ct2@position_VAR]] = as.numeric(as.character(df[[ct2@position_VAR]]))
    df[[ct2@name_VAR]] = factor(df[[ct2@name_VAR]], levels = colnames(ct2))
    if(!is.null(sample_meta_VARS)){
        if(length(sample_meta_VARS) == 1 && sample_meta_VARS == TRUE){
            sample_meta_VARS = colnames(colData(ct2))
        }
        meta_df = colData(ct2) %>%
            as.data.frame %>%
            dplyr::select(all_of(sample_meta_VARS))
        meta_df[[ct2@name_VAR]] = rownames(meta_df)
        df = merge(meta_df, df, by = ct2@name_VAR)
    }
    if(!is.null(region_meta_VARS)){
        browser()
    }
    data.table::as.data.table(df)
}



#' @export
setGeneric("getTidyProfile", function(ct2, sample_meta_VARS = NULL, region_meta_VARS = NULL) standardGeneric("getTidyProfile"), signature = "ct2")

#' @export
setMethod("getTidyProfile", c("ChIPtsne2"), .prof_dt_from_chiptsne2)

