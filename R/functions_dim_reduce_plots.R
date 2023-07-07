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
    if(!hasDimReduce(ct2)){
        stop("No dimensional reduction data present in this ChIPtsne2 object. Run dimReduceTSNE/PCA/UMAP first then try again.")
    }
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
setGeneric("plotDimReducePoints",
           function(ct2, color_VAR = NULL, point_size = NULL)
               standardGeneric("plotDimReducePoints"),
           signature = "ct2")

#' @export
setMethod("plotDimReducePoints", c("ChIPtsne2"), .plotDimReducePoints)


aggregate_signals = function(profile_dt,
                             agg_FUN = max,
                             y_ = "y",
                             yout_ = "y",
                             xmin = -Inf,
                             xmax = Inf,
                             by_ = c("sample")){
    agg_dt = profile_dt[x >= xmin & x <= xmax,
                        .(val_ = agg_FUN(get(y_))),
                        by = c("id", by_)]
    agg_dt[[yout_]] = agg_dt$val_
    agg_dt[, c(yout_, "id", by_), with = FALSE]
}

plot_binned_aggregates = function(agg_dt,
                                  xbins = 50,
                                  ybins = xbins,
                                  xrng = NULL,
                                  yrng = NULL,
                                  val = "y",
                                  bxval = "tx",
                                  byval = "ty",
                                  facet_ = "wide_var",
                                  extra_vars = character(),
                                  bin_met = mean,
                                  min_size = 1, return_data = FALSE){



    if(is.null(xrng)) xrng = range(agg_dt[[bxval]])
    if(is.null(yrng)) yrng = range(agg_dt[[byval]])
    agg_dt[bxval >= min(xrng) & bxval <= max(xrng) &
               byval >= min(yrng) & byval <= max(yrng)]
    agg_dt[, bx := bin_values(get(bxval), n_bins = xbins, xrng = xrng)]
    agg_dt[, by := bin_values(get(byval), n_bins = ybins, xrng = yrng)]

    bin_dt = agg_dt[, .(y = bin_met(get(val)), N = .N), c(unique(c(facet_, extra_vars, "bx", "by")))]
    bxvc = bin_values_centers(n_bins = xbins, xrng)
    w = diff(bxvc[1:2])
    byvc = bin_values_centers(n_bins = ybins, yrng)
    h = diff(byvc[1:2])
    bin_dt[, tx := bxvc[bx]]
    bin_dt[, ty := byvc[by]]
    if(return_data){
        return(bin_dt)
    }
    if(nrow(bin_dt[N >= min_size]) == 0){
        "All bins would have been removed based on min_size. Relaxing min_size to 0."
        min_size = 0
    }
    ggplot(bin_dt[N >= min_size], aes(x = tx, y = ty, fill = y)) +
        geom_tile(width = w, height = h) +
        facet_wrap(facet_) +
        scale_fill_viridis_c() +
        coord_cartesian(xlim = xrng, ylim = yrng)
}

bin_values = function(x, n_bins, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(rescale_capped(x, 0:1, xrng) * (n_bins-.00001))+1
}

rescale_capped = function(x, to = c(0,1), from = range(x, na.rm = TRUE, finite = TRUE)){
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}

bin_values_centers = function(n_bins, rng){
    if(length(rng) != 2)
        rng = range(rng)
    stopifnot(length(rng) == 2)
    xspc = diff(rng)/n_bins/2
    xs = seq(min(rng)+xspc, max(rng)-xspc, diff(rng)/(n_bins))
    xs
}

.plotDimReduceBins = function(ct2,
                              xmin = -Inf,
                              xmax = Inf,
                              agg_FUN = max,
                              xbins = 50,
                              ybins = xbins,
                              bg_color = "gray60",
                              min_size = 5,
                              extra_vars = character()){
    if(!hasDimReduce(ct2)){
        stop("No dimensional reduction data present in this ChIPtsne2 object. Run dimReduceTSNE/PCA/UMAP first then try again.")
    }
    prof_dt = getTidyProfile(ct2)
    tsne_dt = getRegionMetaData(ct2, c("tx", "ty"))
    #aggregate_signals is retrieves a single representative value for a genomic region per sample
    #by default this is the maximum anywhere in the region but this can be
    #overriden using xmin/xmax and agg_FUN
    agg_dt = aggregate_signals(
        prof_dt,
        y_ = ct2@value_VAR,
        yout_ = ct2@value_VAR,
        agg_FUN = agg_FUN,
        xmin = xmin,
        xmax = xmax,
        by_ = ct2@name_VAR)
    agg_dt = merge(agg_dt, tsne_dt, by = "id")
    # m_vars = setdiff(colnames(sts$signal_config$meta_data), colnames(agg_dt))
    # agg_dt = as.data.table(merge(sts$signal_config$meta_data[, m_vars], agg_dt, by.y = "wide_var", by.x = "name"))
    # facet_str = paste0(sts$signal_config$color_by, "~", sts$signal_config$run_by)
    facet_str = paste0(ct2@name_VAR, "~", ".")
    #TODO add extra vars
    #TODO add ssvQC win_size
    extra_vars =  unique(c(
        # sts$signal_config$run_by,
        # sts$signal_config$color_by,
        # extra_vars,
        ct2@name_VAR

    ))
    extra_vars = extra_vars[extra_vars %in% colnames(agg_dt)]
    plot_binned_aggregates(
        agg_dt = agg_dt,
        xbins = xbins,
        ybins = ybins,
        val = ct2@value_VAR,
        extra_vars = extra_vars,
        facet_ = ct2@name_VAR,
        min_size = min_size
    ) +
        facet_grid(facet_str) +
        theme(
            panel.background = element_rect(fill = bg_color),
            panel.grid = element_blank()
        )
}

