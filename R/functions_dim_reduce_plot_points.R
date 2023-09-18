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
#' ct2 = exampleChIPtsne2.with_meta() %>%
#'    dimReduceUMAP() %>%
#'    groupRegionsByDimReduceCluster(group_VAR = "umap_cluster") %>%
#'    groupRegionsBySignalCluster(group_VAR = "signal_cluster")
#'
#' plotDimReducePoints(ct2, NA)
#' plotDimReducePoints(ct2)
#' ct2_diff = subsetCol(ct2, cell == "MCF10A") - subsetCol(ct2, cell == "MCF10AT1")
#' plotDimReducePoints(ct2_diff)
#' plotDimReducePoints(ct2, "umap_cluster")
#' plotDimReducePoints(ct2, c("umap_cluster", "signal_cluster"))
#' plotDimReducePoints(ct2, c("MCF10A_CTCF", "MCF10AT1_CTCF"))
.plotDimReducePoints = function(ct2,
                                color_VAR = NULL,
                                point_size = NULL,
                                point_color_limits = c(NA, NA),
                                has_symmetrical_limits = NULL,
                                point_colors = NULL){
    if(!hasDimReduce(ct2)){
        stop("No dimensional reduction data present in this ChIPtsne2 object. Run dimReduceTSNE/PCA/UMAP first then try again.")
    }
    xy_df = getRegionMetaData(ct2) %>%
        dplyr::select(all_of(c("tx", "ty", ct2@region_VAR)))
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
            dplyr::select(all_of(c("tx", "ty", ct2@region_VAR, color_VAR)))
        xy_df = tidyr::pivot_longer(xy_df, setdiff(colnames(xy_df), c("id", "tx", "ty")), names_to = "group")
        p = ggplot(xy_df, aes(x = tx, y = ty, color = value)) +
            geom_point(size = point_size) +
            facet_wrap(paste0("~", "group"))
    }else if(all(color_VAR %in% colnames(ct2))){
        signal_df = SummarizedExperiment::assay(ct2, "max") %>%
            as.data.frame
        signal_df = signal_df[, color_VAR, drop = FALSE]
        signal_df[[ct2@region_VAR]] = rownames(signal_df)
        xy_df = merge(xy_df, signal_df, by = ct2@region_VAR)
        xy_df = tidyr::pivot_longer(xy_df, setdiff(colnames(xy_df), c("id", "tx", "ty")), names_to = ct2@name_VAR, values_to = "max")

        point_colors = .prep_color_scale(xy_df$max, point_colors)
        point_color_limits = .prep_symmetrical(xy_df$max, has_symmetrical_limits, point_color_limits)
        xy_df$max = .apply_limits(xy_df$max, point_color_limits)

        p = ggplot(xy_df, aes(x = tx, y = ty, color = max)) +
            geom_point(size = point_size) +
            facet_wrap(paste0("~", ct2@name_VAR)) +
            labs(color = paste("max", ct2@value_VAR, "\nper", ct2@region_VAR))
        p = .apply_scale(p, point_colors, point_color_limits, fill = FALSE)
    }

    p
}

generic_plotDimReducePoints = function(ct2,
                                       color_VAR = NULL,
                                       point_size = NULL,
                                       point_color_limits = c(NA, NA),
                                       has_symmetrical_limits = NULL,
                                       point_colors = NULL){
    standardGeneric("plotDimReducePoints")
}

#' @export
setGeneric("plotDimReducePoints",
           generic_plotDimReducePoints,
           signature = "ct2")

#' @export
setMethod("plotDimReducePoints", c("ChIPtsne2"), .plotDimReducePoints)
