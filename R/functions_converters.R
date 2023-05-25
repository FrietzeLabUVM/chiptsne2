#'
#' @importFrom dplyr select
#' @importFrom data.table as.data.table
.prof_dt_from_chiptsne2 = function(ct, sample_meta_VARS = NULL){
    df = do.call(rbind,
                 lapply(names(ct@colToRowMatCols), function(nam){
                     i = ct@colToRowMatCols[[nam]]
                     df = reshape2::melt(ct@rowToRowMat[, i])
                     # df$Var1 = NULL
                     df$Var2 = factor(df$Var2)
                     # df$Var2 = as.numeric(sapply(strsplit(levels(df$Var2), "_"), function(x)x[length(x)]))
                     levels(df$Var2) = as.numeric(sub(paste0(nam, "_"), "", levels(df$Var2)))
                     df$name = nam
                     df
                 })
    )
    colnames(df) = c("id", "x", "y", "name")
    df$id = factor(df$id, levels = names(rowRanges(ct)))
    df$x = as.numeric(as.character(df$x))
    df$name = factor(df$name, levels = colnames(ct))
    if(!is.null(sample_meta_VARS)){
        if(length(sample_meta_VARS) == 1 && sample_meta_VARS == TRUE){
            sample_meta_VARS = colnames(colData(ct))
        }
        meta_df = colData(ct) %>%
            as.data.frame %>%
            dplyr::select(all_of(sample_meta_VARS))
        meta_df$name = rownames(meta_df)
        df = merge(meta_df, df, by = "name")
    }
    data.table::as.data.table(df)
}
