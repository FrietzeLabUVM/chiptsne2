
#' Title
#'
#' @param ct2
#' @param group_VAR
#' @param new_meta_VAR
#'
#' @return
#' @export
#'
#' @examples
aggregateRegionsByGroup = function(ct2, group_VAR, new_meta_VAR = ifelse(length(group_VAR) == 1, group_VAR, "meta_id")){
    centroid = calculateGroupCentroid(ct2, group_VAR)

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
    colnames(df) = c(new_meta_VAR, ct2@position_VAR, ct2@value_VAR, ct2@name_VAR)
    df[[new_meta_VAR]] = factor(df[[new_meta_VAR]], levels = rownames(centroid))

    rd = unique(rowData(ct2)[, group_VAR, drop = FALSE])
    rownames(rd) = apply(rd, 1, paste, collapse = ",")
    rd[[new_meta_VAR]] = rownames(rd)

    ct2.meta = ChIPtsne2.from_tidy(
        df,
        query_gr = NULL,
        sample_metadata = colData(ct2),
        region_metadata = rd,
        position_VAR = ct2@position_VAR,
        name_VAR = ct2@name_VAR,
        value_VAR = ct2@value_VAR,
        region_VAR = new_meta_VAR)
    ct2.meta
}

#' Title
#'
#' @param ct2
#' @param group_VAR
#' @param new_meta_VAR
#'
#' @return
#' @export
#'
#' @examples
aggregateSamplesByGroup = function(ct2, group_VAR, new_meta_VAR = ifelse(length(group_VAR) == 1, group_VAR, "meta_name")){
    carried_VARS = unique(c(group_VAR, new_meta_VAR))
    if(group_VAR != new_meta_VAR){
        ct2@colData[[new_meta_VAR]] =  apply(colData(ct2)[, group_VAR, drop = FALSE], 1, paste, collapse = ",")
    }
    ct2.sp = split(ct2, new_meta_VAR)
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
    ct2.meta = setNameVariable(ct2.meta, new_meta_VAR)
    ct2.meta
}

#' Title
#'
#' @param ct2
#' @param group_VAR
#'
#' @return
#' @export
#'
#' @examples
aggregateByGroup = function(ct2, group_VAR){
    group_VAR.col = group_VAR[group_VAR %in% colnames(colData(ct2))]
    group_VAR.row = group_VAR[group_VAR %in% colnames(rowData(ct2))]
    if(length(group_VAR.col) > 0){
        ct2.meta = aggregateSamplesByGroup(ct2, group_VAR.col)
    }else{
        ct2.meta = ct2
    }
    if(length(group_VAR.row) > 0){
        ct2.meta = aggregateRegionsByGroup(ct2.meta, group_VAR.row)
    }

    ct2.meta
}

#' plotSignalLinePlot
#'
#' @param ct2 A ChIPtsne2 object
#' @param group_VAR A grouping variable used for vertical facets.
#' @param color_VAR A grouping variable used for coloring. Use scale_color_* from ggplot2 to control coloring.
#' @param facet_VAR A grouping variable used for horizontal facets.
#' @param extra_VARS Extra variables that should still be present in final ggplot data.
#' @param return_data If TRUE, return the data.table instead of creating a plot.
#'
#' @return ggplot2 of averaged signal profiles, potentially faceted in interesting ways.
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF")
#' plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF", moving_average_window = 5)
#' plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF", n_splines = 5)
#' plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF", moving_average_window = 5, n_splines = 5)
#' plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF", facet_VAR = NULL, color_VAR = "cell")
#' plotSignalLinePlot(ct2, color_VAR = "cell", facet_VAR = "peak_MCF10AT1_CTCF", group_VAR = "peak_MCF10A_CTCF")
#'
#' ct2 = groupRegionsBySignalCluster(ct2, group_VAR = "cluster")
#' ct2 = groupRegionsByOverlap(ct2, seqsetvis::CTCF_in_10a_narrowPeak_grs[1:2], group_VAR = "overlap")
#' plotSignalLinePlot(ct2, color_VAR = "cell", group_VAR = "cluster", facet_VAR = NULL)
#' plotSignalLinePlot(ct2, color_VAR = "cell", group_VAR = "cluster", facet_VAR = "peak_MCF10A_CTCF")
#'
#' plotSignalLinePlot(ct2, color_VAR = "peak_MCF10A_CTCF", group_VAR = "cluster", facet_VAR = "cell", extra_VARS = "mark")
.plotSignalLinePlot = function(ct2,
                               group_VAR = NULL,
                               color_VAR = NULL,
                               facet_VAR = ct2@name_VAR,
                               linewidth = 1.5,
                               extra_VARS = character(),
                               moving_average_window = 1,
                               n_splines = 1,
                               return_data = FALSE){
    all_VARS = unique(c(group_VAR, color_VAR, facet_VAR, extra_VARS))
    meta_VARS = setdiff(all_VARS, ct2@name_VAR)
    row_VARS = meta_VARS[meta_VARS %in% colnames(rowData(ct2))]
    if(length(row_VARS) > 0){
        ct2.agg = aggregateRegionsByGroup(ct2, row_VARS)
    }else{
        ct2.agg = ct2
    }
    agg_dt = getTidyProfile(ct2.agg, meta_VARS = meta_VARS)
    # prof_dt = getTidyProfile(ct2, meta_VARS = meta_VARS)
    # agg_dt = prof_dt[, .(VALUE_ = mean(get(ct2@value_VAR))), c(unique(c(all_VARS, ct2@position_VAR, ct2@name_VAR)))]
    # data.table::setnames(agg_dt, "VALUE_", ct2@value_VAR)
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
        agg_dt = seqsetvis:::applyMovingAverage(agg_dt, n = moving_average_window, centered = TRUE, x_ = ct2@position_VAR, y_ = ct2@value_VAR, by_ = c(all_VARS, "GROUP_"))
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
                                    linewidth = 1.5,
                                    extra_VARS = character(),
                                    moving_average_window = 1,
                                    n_splines = 1,
                                    return_data = FALSE){
    args = get_args(to_ignore = NULL)
    do.call(.plotSignalLinePlot, args = args)
}

#' @export
setGeneric("plotSignalLinePlot", function(
        ct2,
        group_VAR = NULL,
        color_VAR = NULL,
        facet_VAR = ct2@name_VAR,
        linewidth = 1.5,
        extra_VARS = character(),
        moving_average_window = 1,
        n_splines = 1,
        return_data = FALSE)
    standardGeneric("plotSignalLinePlot"),
    signature = "ct2")

#' @export
setMethod("plotSignalLinePlot", c("ChIPtsne2"), .plotSignalLinePlot)

#' @export
setMethod("plotSignalLinePlot", c("ChIPtsne2_no_rowRanges"), .plotSignalLinePlot_meta)


