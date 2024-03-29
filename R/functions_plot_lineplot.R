
.plotSignalLinePlot = function(ct2,
                               group_VAR = NULL,
                               color_VAR = NULL,
                               facet_VAR = ct2@name_VAR,
                               extra_VARS = character(),
                               linewidth = 1.5,
                               moving_average_window = 1,
                               n_splines = 1,
                               return_data = FALSE){
    all_VARS = unique(c(group_VAR, color_VAR, facet_VAR, extra_VARS))
    meta_VARS = setdiff(all_VARS, ct2@name_VAR)
    row_VARS = meta_VARS[meta_VARS %in% colnames(rowData(ct2))]
    if(length(row_VARS) > 0){
        ct2.agg = aggregateRegionsByGroup(ct2, row_VARS)
    }else{
        rowData(ct2)$`__ALL__` = "all regions"
        ct2.agg = aggregateRegionsByGroup(ct2, "__ALL__")
        all_VARS = unique(c( "__ALL__", color_VAR, facet_VAR, extra_VARS))
        meta_VARS = setdiff(all_VARS, ct2@name_VAR)
    }
    agg_dt = getTidyProfile(ct2.agg, meta_VARS = meta_VARS)
    x_ = ct2@position_VAR
    x_ = ensym(x_)
    y_ = ct2@value_VAR
    y_ = ensym(y_)

    agg_dt$GROUP_ = apply(agg_dt[, c(all_VARS), with = FALSE], 1, function(x)paste(x, collapse = ","))
    group_ = "GROUP_"
    group_ = ensym(group_)
    if(!is.null(color_VAR))
        color_VAR = ensym(color_VAR)

    plot_subtitle = ifelse(is.null(group_VAR), "all regions", "")

    if(is.null(group_VAR)) group_VAR = "."
    if(is.null(facet_VAR)) facet_VAR = "."
    agg_dt = agg_dt[order(get(ct2@position_VAR))]
    if(moving_average_window > 1){
        agg_dt = seqsetvis::applyMovingAverage(agg_dt, n = moving_average_window, centered = TRUE, x_ = ct2@position_VAR, y_ = ct2@value_VAR, by_ = c(all_VARS, "GROUP_"))
    }
    if(n_splines > 1){
        agg_dt = seqsetvis::applySpline(agg_dt, n = n_splines, x_ = ct2@position_VAR, y_ = ct2@value_VAR, by_ = c(all_VARS, "GROUP_"))
    }

    if(return_data) return(agg_dt)

    ggplot(agg_dt,
           aes(x = !!x_, y = !!y_, color = !!color_VAR, group = !!group_)) +
        geom_path(linewidth = linewidth) +
        facet_grid(paste0(group_VAR, "~", facet_VAR), labeller = label_both) +
        labs(subtitle = plot_subtitle)
}

.plotSignalLinePlot_meta = function(ct2,
                                    group_VAR = ct2@region_VAR,
                                    color_VAR = NULL,
                                    facet_VAR = ct2@name_VAR,
                                    extra_VARS = character(),
                                    linewidth = 1.5,
                                    moving_average_window = 1,
                                    n_splines = 1,
                                    return_data = FALSE){
    args = get_args(to_ignore = NULL)
    do.call(.plotSignalLinePlot, args = args)
}


#' plotSignalLinePlot
#'
#' @param ct2 A ChIPtsne2 object
#' @param group_VAR A grouping variable used for vertical facets.
#' @param color_VAR A grouping variable used for coloring. Use scale_color_* from ggplot2 to control coloring.
#' @param facet_VAR A grouping variable used for horizontal facets.
#' @param extra_VARS Extra variables that should still be present in final ggplot data.
#' @param return_data If TRUE, return the data.table instead of creating a plot.
#' @param linewidth The pt width of lines plotted. Default is 1.5.
#' @param moving_average_window The size of windows for moving average applied profiles points. Default is 1 (disabled).
#' @param n_splines The number of splines to interpolate between profile points. Default is 1 (disabled).
#'
#' @return ggplot2 of averaged signal profiles, potentially faceted in interesting ways.
#' @export
#' @rdname plotSignalLinePlot
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF")
#' # smoothing profiles
#' plotSignalLinePlot(
#'   ct2,
#'   group_VAR = "peak_MCF10A_CTCF",
#'   n_splines = 5,
#'   moving_average_window = 5)
#' # map sample meta data to color
#' plotSignalLinePlot(
#'   ct2,
#'   group_VAR = "peak_MCF10A_CTCF",
#'   facet_VAR = NULL,
#'   color_VAR = "cell"
#' )
#' # map region meta data to facets
#' plotSignalLinePlot(
#'   ct2,
#'   color_VAR = "cell",
#'   facet_VAR = "peak_MCF10AT1_CTCF",
#'   group_VAR = "peak_MCF10A_CTCF"
#' )
#'
#' # more complicated examples using clustering
#' ct2 = groupRegionsBySignalCluster(ct2, group_VAR = "cluster")
#' ct2 = groupRegionsByOverlap(
#'   ct2,
#'   seqsetvis::CTCF_in_10a_narrowPeak_grs[1:2],
#'   group_VAR = "overlap"
#' )
#' plotSignalLinePlot(
#'   ct2,
#'   color_VAR = "cell",
#'   group_VAR = "cluster",
#'   facet_VAR = "peak_MCF10A_CTCF"
#' )
#'
#' plotSignalLinePlot(
#'   ct2,
#'   color_VAR = "peak_MCF10A_CTCF",
#'   group_VAR = "cluster",
#'   facet_VAR = "cell"
#' )
setGeneric("plotSignalLinePlot", function(
        ct2,
        group_VAR = NULL,
        color_VAR = NULL,
        facet_VAR = ct2@name_VAR,
        extra_VARS = character(),
        linewidth = 1.5,
        moving_average_window = 1,
        n_splines = 1,
        return_data = FALSE)
    standardGeneric("plotSignalLinePlot"),
    signature = "ct2")

#' @export
#' @rdname plotSignalLinePlot
setMethod("plotSignalLinePlot", c("ChIPtsne2"), .plotSignalLinePlot)

#' @export
#' @rdname plotSignalLinePlot
setMethod("plotSignalLinePlot", c("ChIPtsne2_no_rowRanges"), .plotSignalLinePlot_meta)


