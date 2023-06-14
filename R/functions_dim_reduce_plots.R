

#' plotDimReducePoints
#'
#' @param ct2 valid ChIPtsne2 after dimReduce has been run
#' @param color_VAR Control color assignment in plot. Can match entries in *either* sample metadata (colnames) or region metadata (rowRanges). Default of NULL will plot max signal for all sample profiles. NA will perform no color mapping.
#' @param point_size Size of points in plot.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2() %>%
#'    dimReduceUMAP() %>%
#'    groupRegionsByDimReduceCluster(group_VAR = "umap_cluster") %>%
#'    groupRegionsBySignalCluster(group_VAR = "signal_cluster")
#'
#' plotDimReducePoints(ct2, NA)
#' plotDimReducePoints(ct2)
#' plotDimReducePoints(ct2, "umap_cluster")
#' plotDimReducePoints(ct2, c("umap_cluster", "signal_cluster"))
#' plotDimReducePoints(ct2, c("MCF10A_CTCF", "MCF10AT1_CTCF"))
.plotDimReducePoints = function(ct2, color_VAR = NULL, point_size = NULL){
    xy_df = getRegionMetaData(ct2) %>%
        select(all_of(c("tx", "ty", ct2@region_VAR)))
    if(is.null(point_size)){
        nr = nrow(xy_df)
        point_size = 1/nr*100
        if(point_size < .05) point_size = .05
        if(point_size > 1) point_size = 1
    }
    if(is.null(color_VAR)){
        color_VAR = colnames(ct2)
    }
    if(all(is.na(color_VAR))){
        p = ggplot(xy_df, aes(x = tx, y = ty)) +
            geom_point(size = point_size)
    }else if(all(color_VAR %in% colnames(getRegionMetaData(ct2)))){
        xy_df = getRegionMetaData(ct2) %>%
            select(all_of(c("tx", "ty", ct2@region_VAR, color_VAR)))
        xy_df = pivot_longer(xy_df, setdiff(colnames(xy_df), c("id", "tx", "ty")), names_to = "group")
        p = ggplot(xy_df, aes(x = tx, y = ty, color = value)) +
            geom_point(size = point_size) +
            facet_wrap(paste0("~", "group"))
    }else if(all(color_VAR %in% colnames(ct2))){
        signal_df = SummarizedExperiment::assay(ct2, "max") %>%
            as.data.frame
        signal_df = signal_df[, color_VAR]
        signal_df[[ct2@region_VAR]] = rownames(signal_df)
        xy_df = merge(xy_df, signal_df, by = ct2@region_VAR)
        xy_df = pivot_longer(xy_df, setdiff(colnames(xy_df), c("id", "tx", "ty")), names_to = ct2@name_VAR, values_to = "max")
        p = ggplot(xy_df, aes(x = tx, y = ty, color = max)) +
            geom_point(size = point_size) +
            facet_wrap(paste0("~", ct2@name_VAR))
    }

    p
}

#' @export
setGeneric("plotDimReducePoints", function(ct2, color_VAR = NULL, point_size = NULL) standardGeneric("plotDimReducePoints"), signature = "ct2")

#' @export
setMethod("plotDimReducePoints", c("ChIPtsne2"), .plotDimReducePoints)