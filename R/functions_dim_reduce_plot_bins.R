make_facet_str = function(facet_rows = character(), facet_columns = character()){
    if(length(facet_rows) > 0){
        row_str = paste0("`", paste(facet_rows, collapse = "`+`"), "`")
    }else{
        row_str = "."
    }
    if(length(facet_columns) > 0){
        col_str = paste0("`", paste(facet_columns, collapse = "`+`"), "`")
    }else{
        col_str = "."
    }

    paste0(row_str, "~", col_str)
}

.plotDimReduceBins = function(ct2,
                              facet_rows = ct2@name_VAR,
                              facet_columns = NULL,
                              xmin = -Inf,
                              xmax = Inf,
                              agg_FUN = max,
                              x_bins = NULL,
                              y_bins = x_bins,
                              bg_color = "gray60",
                              min_size = 1,
                              bin_fill_limits = c(NA, NA),
                              has_symmetrical_limits = NULL,
                              bin_colors = NULL,
                              extra_VARS = character(),
                              return_data = FALSE){
    if(!hasDimReduce(ct2)){
        stop("No dimensional reduction data present in this ChIPtsne2 object. Run dimReduceTSNE/PCA/UMAP first then try again.")
    }
    if(is.null(x_bins)){
        x_bins = min(round(nrow(ct2)^.5), 10)
    }
    if(is.null(y_bins)){
        y_bins = x_bins
    }
    meta_VARS = c(facet_rows, facet_columns, "tx", "ty", extra_VARS)
    prof_dt = getTidyProfile(ct2, meta_VARS = meta_VARS)
    #aggregate_signals is retrieves a single representative value for a genomic region per sample
    #by default this is the maximum anywhere in the region but this can be
    #overridden using xmin/xmax and agg_FUN
    extra_VARS = union(c(facet_rows, facet_columns), extra_VARS)
    agg_dt = aggregate_signals(
        prof_dt,
        x_ = ct2@position_VAR,
        y_ = ct2@value_VAR,
        id_ = ct2@region_VAR,
        yout_ = ct2@value_VAR,
        agg_FUN = agg_FUN,
        xmin = xmin,
        xmax = xmax,
        by_ = meta_VARS
    )
    if(is.null(facet_rows)){
        facet_rows = "."
    }
    if(is.null(facet_columns)){
        facet_columns = "."
    }
    # facet_str = paste0(facet_rows, "~", facet_columns)
    facet_str = make_facet_str(facet_rows = facet_rows, facet_columns = facet_columns)
    bin_dt = bin_signals(
        agg_dt = agg_dt,
        x_bins = x_bins,
        y_bins = y_bins,
        val = ct2@value_VAR,
        extra_VARS = extra_VARS[-1],
        facet_ = extra_VARS[1],
        min_size = min_size
    )

    bin_colors = .prep_color_scale(bin_dt[[ct2@value_VAR]], color_scale = bin_colors)
    bin_fill_limits = .prep_symmetrical(bin_dt[[ct2@value_VAR]], has_symmetrical_limits = has_symmetrical_limits, scale_limits = bin_fill_limits)
    bin_dt[[ct2@value_VAR]] = .apply_limits(bin_dt[[ct2@value_VAR]], bin_fill_limits)

    if(return_data){
        return(bin_dt)
    }

    p_bin = plot_binned_aggregates(
        bin_dt = bin_dt,
        fill_VAR = ct2@value_VAR,
        facet_ = extra_VARS[1],
        min_size = min_size
    ) +
        facet_grid(facet_str) +
        theme(
            panel.background = ggplot2::element_rect(fill = bg_color),
            panel.grid = ggplot2::element_blank()
        ) +
        labs(caption = paste("Binned to", y_bins, "rows by", x_bins, "columns."))
    p_bin = .apply_scale(p_bin, bin_colors, bin_fill_limits)
    p_bin
}

generic_plotDimReduceBins = function(ct2,
                                     facet_rows = ct2@name_VAR,
                                     facet_columns = NULL,
                                     xmin = -Inf,
                                     xmax = Inf,
                                     agg_FUN = max,
                                     x_bins = 50,
                                     y_bins = x_bins,
                                     bg_color = "gray60",
                                     min_size = 1,
                                     bin_fill_limits = c(NA, NA),
                                     has_symmetrical_limits = NULL,
                                     bin_colors = NULL,
                                     extra_VARS = character(),
                                     return_data = FALSE){
    standardGeneric("plotDimReduceBins")
}


#' #' plotDimReduceBins
#'
#' @param ct2 `r doc_ct2()`
#' @param facet_rows The colData or rowData attribute for facetting row.
#' @param facet_columns The colData or rowData attribute for facetting column
#' @param xmin Minimum value of profile position allowed. Default is -Inf.
#' @param xmax Maximum value of profile position allowed. Default is Inf.
#' @param agg_FUN Function used to summarize each profile to a single value.
#'   Default is max.
#' @param x_bins Resolution in dim reduced x-axis. Defaults to change with sqrt
#'   of number of regions.
#' @param y_bins Resolution in dim reduced y-axis. Defaults to same as `x_bins`.
#' @param bg_color Background color for plot. Passed to fill in plot.background
#'   of ggplot2 theme.
#' @param min_size Bins with fewer items than this value will be omitted.
#'   Default is 1.
#' @param bin_fill_limits Fill color scale limits to apply. Default is to fit
#'   all data.
#' @param has_symmetrical_limits If TRUE, fill color scale limits will be
#'   symmetrical around 0. May override `bin_fill_limits`. Default is FALSE.
#' @param bin_colors Colors to interpolate for fill scale.
#' @param extra_VARS `r doc_extra_VARS()`
#' @param return_data `r doc_return_data()`
#'
#' @return A ggplot2 object summarizing profiles by tiles at a resolution
#'   determined by `x_bins` and `y_bins`.
#' @rdname plotDimReduceBins
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = dimReducePCA(ct2)
#' plotDimReduceBins(ct2)
#'
#' # built in support for facet_grid
#' plotDimReduceBins(ct2, facet_rows = "cell", facet_columns = "mark")
#' # alternatively you can use extra_VARS and control facetting yourself
#' plotDimReduceBins(ct2, extra_VARS = c("cell", "mark")) +
#'   facet_grid(cell~mark)
#'
#' # attributes from rowData can be used as well
#' plotDimReduceBins(
#'   ct2,
#'   facet_rows = "peak_MCF10AT1_CTCF",
#'   facet_columns = "mark"
#' )
#' plotDimReduceBins(ct2, extra_VARS = c("peak_MCF10AT1_CTCF", "mark")) +
#'   facet_grid(peak_MCF10AT1_CTCF~mark)
#'
#' plotDimReduceBins(ct2, xmin = -300, xmax = 0)
#' plotDimReduceBins(
#'   ct2,
#'   xmin =  0,
#'   xmax = 300,
#'   bin_fill_limits = c(0, 20),
#'   bin_colors = c("blue", "green")
#' )
setGeneric("plotDimReduceBins",
           generic_plotDimReduceBins,
           signature = "ct2")

#' @export
#' @rdname plotDimReduceBins
setMethod("plotDimReduceBins", c("ChIPtsne2_no_rowRanges"), .plotDimReduceBins)

aggregate_signals = function(profile_dt,
                             agg_FUN = max,
                             x_ = "x",
                             y_ = "y",
                             id_ = "id",
                             yout_ = "y",
                             xmin = -Inf,
                             xmax = Inf,
                             by_ = c("sample")){
    agg_dt = profile_dt[get(x_) >= xmin & get(x_) <= xmax,
                        list(val_ = agg_FUN(get(y_))),
                        by = c(id_, by_)]
    agg_dt[[yout_]] = agg_dt$val_
    agg_dt[, c(yout_, id_, by_), with = FALSE]
}

bin_signals = function(agg_dt,
                       x_bins = 50,
                       y_bins = x_bins,
                       xrng = NULL,
                       yrng = NULL,
                       val = "y",
                       bxval = "tx",
                       byval = "ty",
                       facet_ = "wide_var",
                       extra_VARS = character(),
                       bin_met = mean,
                       min_size = 1, return_data = FALSE){
    #visible binding NOTE
    bx = .N = `:=` = tx = ty = NULL
    if(is.null(xrng)) xrng = range(agg_dt[[bxval]])
    if(is.null(yrng)) yrng = range(agg_dt[[byval]])
    agg_dt[bxval >= min(xrng) & bxval <= max(xrng) &
               byval >= min(yrng) & byval <= max(yrng)]
    agg_dt[, bx := bin_values(get(bxval), n_bins = x_bins, xrng = xrng)]
    agg_dt[, by := bin_values(get(byval), n_bins = y_bins, xrng = yrng)]

    bin_dt = agg_dt[, list(y = bin_met(get(val)), N = .N),
                    c(unique(c(facet_, extra_VARS, "bx", "by")))]
    data.table::setnames(bin_dt, "y", val)
    bxvc = bin_values_centers(n_bins = x_bins, xrng)
    w = diff(bxvc[1:2])
    byvc = bin_values_centers(n_bins = y_bins, yrng)
    h = diff(byvc[1:2])
    bin_dt[, tx := bxvc[bx]]
    bin_dt[, ty := byvc[by]]

    bin_dt
}

plot_binned_aggregates = function(bin_dt,
                                  xrng = NULL,
                                  yrng = NULL,
                                  fill_VAR = "y",
                                  bxval = "tx",
                                  byval = "ty",
                                  facet_ = "wide_var",
                                  min_size = 1){
    #visible binding NOTE
    tx = ty = N = NULL
    if(is.null(xrng)) xrng = range(bin_dt[[bxval]])
    if(is.null(yrng)) yrng = range(bin_dt[[byval]])
    if(nrow(bin_dt[N >= min_size]) == 0){
        "All bins would have been removed based on min_size. Relaxing min_size to 0."
        min_size = 0
    }
    w = min(diff(unique(sort(bin_dt$tx))))
    h = min(diff(unique(sort(bin_dt$ty))))
    fill_ = fill_VAR
    fill_ = ensym(fill_)
    ggplot(bin_dt[N >= min_size], aes(x = tx, y = ty, fill = !!fill_)) +
        geom_tile(width = w, height = h) +
        facet_wrap(facet_) +
        coord_cartesian(xlim = xrng, ylim = yrng)
}

bin_values = function(x, n_bins, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(.rescale_capped(x, 0:1, xrng) * (n_bins-.00001))+1
}

bin_values_centers = function(n_bins, rng){
    if(length(rng) != 2)
        rng = range(rng)
    stopifnot(length(rng) == 2)
    xspc = diff(rng)/n_bins/2
    xs = seq(min(rng)+xspc, max(rng)-xspc, diff(rng)/(n_bins))
    xs
}
