
.plotDimReducePoints = function(ct2,
                                color_VAR = NULL,
                                point_size = NULL,
                                point_color_limits = c(NA, NA),
                                has_symmetrical_limits = NULL,
                                point_colors = NULL,
                                extra_VARS = NULL,
                                background_annotation_color = NULL,
                                underlayer_FUN = function(p)p){
    #visible binding NOTE
    tx = ty = value = NULL
    if(!hasDimReduce(ct2)){
        stop("No dimensional reduction data present in this ChIPtsne2 object. Run dimReduceTSNE/PCA/UMAP first then try again.")
    }
    xy_df = getRegionMetaData(ct2, select_VARS = c("tx", "ty"))
    if(is.null(point_size)){
        nr = nrow(xy_df)
        point_size = 1/nr*100
        if(point_size < .05) point_size = .05
        if(point_size > 1) point_size = 1
    }
    if(is.null(color_VAR)){
        color_VAR = colnames(ct2)
    }

    enforce_extra_VARS = function(ct2, df, ev){
        ev_missed = setdiff(ev, colnames(df))
        if(length(ev_missed) > 0){
            if(ct2@region_VAR %in% colnames(df) & ct2@name_VAR %in% colnames(df)){
                full_region_df = getRegionMetaData(ct2)
                region_df = getRegionMetaData(ct2, select_VARS = intersect(extra_VARS, colnames(full_region_df)))
                full_sample_df = getSampleMetaData(ct2)
                sample_df = getSampleMetaData(ct2, select_VARS = intersect(extra_VARS, colnames(full_sample_df)))
                df = merge(df, region_df, by = ct2@region_VAR)
                df = merge(df, sample_df, by = ct2@name_VAR)
            }else if(ct2@region_VAR %in% colnames(df)){
                full_region_df = getRegionMetaData(ct2)
                region_df = getRegionMetaData(ct2, select_VARS = intersect(extra_VARS, colnames(full_region_df)))
                df = merge(df, region_df, by = ct2@region_VAR)
            }else if(ct2@name_VAR %in% colnames(df)){
                full_sample_df = getSampleMetaData(ct2)
                sample_df = getSampleMetaData(ct2, select_VARS = intersect(extra_VARS, colnames(full_sample_df)))
                df = merge(df, sample_df, by = ct2@name_VAR)
            }else{
                stop("confusing")
            }
        }
        df
    }

    background_FUN = function(p){
        if(!is.null(background_annotation_color)){
            bg_df = unique(xy_df[, c("tx", "ty")])
            p = p + annotate("point", x = bg_df$tx, y = bg_df$ty, color = background_annotation_color, size = .7*point_size)
        }
        p
    }

    if(all(is.na(color_VAR))){
        xy_df = enforce_extra_VARS(ct2, xy_df, extra_VARS)
        p = ggplot(xy_df, aes(x = tx, y = ty))
        p = underlayer_FUN(p)
        p = background_FUN(p)
        p = p +
            geom_point(size = point_size)
    }else if(all(color_VAR %in% colnames(getRegionMetaData(ct2)))){
        xy_df = getRegionMetaData(ct2) %>%
            dplyr::select(dplyr::all_of(c("tx", "ty", ct2@region_VAR, color_VAR)))
        xy_df = tidyr::pivot_longer(xy_df, setdiff(colnames(xy_df), c(ct2@region_VAR, "tx", "ty")), names_to = "group")
        xy_df = enforce_extra_VARS(ct2, xy_df, extra_VARS)
        point_colors = .prep_color_scale(xy_df$value, color_scale = point_colors)
        p = ggplot(xy_df, aes(x = tx, y = ty, color = value)) +
            scale_color_manual(values = point_colors)
        p = underlayer_FUN(p)
        p = background_FUN(p)
        p = p +
            geom_point(size = point_size) +
            facet_wrap(paste0("~", "group"))
    }else if(all(color_VAR %in% colnames(ct2))){
        signal_df = SummarizedExperiment::assay(ct2, "max") %>%
            as.data.frame
        signal_df = signal_df[, color_VAR, drop = FALSE]
        signal_df[[ct2@region_VAR]] = rownames(signal_df)
        xy_df = merge(xy_df, signal_df, by = ct2@region_VAR)
        xy_df = tidyr::pivot_longer(xy_df, setdiff(colnames(xy_df), c(ct2@region_VAR, "tx", "ty")), names_to = ct2@name_VAR, values_to = "max")
        point_colors = .prep_color_scale(xy_df$max, has_symmetrical_limits, point_colors)
        point_color_limits = .prep_symmetrical(xy_df$max, has_symmetrical_limits, point_color_limits)
        xy_df$max = .apply_limits(xy_df$max, point_color_limits)
        xy_df = enforce_extra_VARS(ct2, xy_df, extra_VARS)
        p = ggplot(xy_df, aes(x = tx, y = ty, color = max))
        p = underlayer_FUN(p)
        p = background_FUN(p)
        p = p +
            geom_point(size = point_size) +
            facet_wrap(paste0("~", ct2@name_VAR)) +
            labs(color = paste("max", ct2@value_VAR, "\nper", ct2@region_VAR))
        p = .apply_scale(p, point_colors, point_color_limits, fill = FALSE)
    }else{
        stop("color_VAR: \"", color_VAR, "\" was not recognized. Check vs colnames of ct2 object or colnames of rowData(ct2).")
    }
    p
}

generic_plotDimReducePoints = function(ct2,
                                       color_VAR = NULL,
                                       point_size = NULL,
                                       point_color_limits = c(NA, NA),
                                       has_symmetrical_limits = NULL,
                                       point_colors = NULL,
                                       extra_VARS = NULL,
                                       background_annotation_color = NULL,
                                       underlayer_FUN = function(p)p){
    standardGeneric("plotDimReducePoints")
}


#' plotDimReducePoints
#'
#' @param ct2 valid ChIPtsne2 after dimReduce has been run
#' @param color_VAR Control color assignment in plot. Can match entries in
#'   *either* sample metadata (colnames) or region metadata (rowRanges). Default
#'   of NULL will plot max signal for all sample profiles. NA will perform no
#'   color mapping.
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
#' ct2_diff = subsetSamples(ct2, cell == "MCF10A") -
#'   subsetSamples(ct2, cell == "MCF10AT1")
#' plotDimReducePoints(ct2_diff)
#' plotDimReducePoints(ct2, "umap_cluster")
#' plotDimReducePoints(ct2, c("umap_cluster", "signal_cluster"))
#' plotDimReducePoints(ct2, c("MCF10A_CTCF", "MCF10AT1_CTCF"))
#'
#' # layer plot elements beneath the final plot with a function like this:
#' base_plot = function(p){
#'   p +
#'     annotate("rect",
#'              xmin = 0, xmax = .3,
#'              ymin = -.05, ymax = .13,
#'              fill = "gray80", color = "red")
#' }
#'
#' plotDimReducePoints(ct2, extra_VARS = "peak_MCF10CA1_CTCF",
#'   background_annotation_color = "gray50",
#'   underlayer_FUN = base_plot) +
#'     facet_grid(peak_MCF10CA1_CTCF~sample) +
#'     labs(title = "signal per sample facetted by peak_MCF10CA1_CTCF")
setGeneric("plotDimReducePoints",
           generic_plotDimReducePoints,
           signature = "ct2")

#' @export
setMethod("plotDimReducePoints", c("ChIPtsne2_no_rowRanges"), .plotDimReducePoints)
