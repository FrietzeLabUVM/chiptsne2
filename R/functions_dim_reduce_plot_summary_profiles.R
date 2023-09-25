#' plotDimReduceSummaryProfiles
#'
#' @param ct2
#' @param color_VAR
#' @param x_bins
#' @param y_bins
#' @param extra_vars
#' @param xrng
#' @param yrng
#' @param value_limits
#' @param ma_size
#' @param n_splines
#' @param p
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#' @param return_data
#'
#' @return
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = dimReducePCA(ct2)
#' plotDimReduceSummaryProfiles(ct2)
#' plotDimReduceSummaryProfiles(ct2, value_limits = c(0, 50))
#' plotDimReduceSummaryProfiles(ct2, value_limits = c(0, 200))
#' plotDimReduceSummaryProfiles(ct2, value_limits = c(0, 500))
.plotDimReduceSummaryProfiles = function(
        ct2,
        color_VAR = ct2@name_VAR,
        x_bins = min(nrow(ct2)^.5, 10),
        y_bins = x_bins,
        extra_vars = character(),
        xrng = NULL,
        yrng = NULL,
        value_limits = c(0, 1),
        ma_size = 2,
        n_splines = 10,
        p = NULL,
        N_floor = 0,
        N_ceiling = NULL,
        min_size = 0.3,
        return_data = FALSE){

    position_VAR = ct2@position_VAR
    value_VAR = ct2@value_VAR
    region_VAR = ct2@region_VAR

    profile_dt = getTidyProfile(ct2, meta_VARS = TRUE)
    cn = colnames(getSampleMetaData(ct2))
    profile_dt = profile_dt[order(x)]
    position_dt = data.table::as.data.table(getRegionMetaData(ct2))
    extra_vars = union(extra_vars, cn)

    stopifnot(length(value_limits) == 2)
    if(is.na(value_limits[1])){
        value_limits[1] = min(profile_dt[[ct2@value_VAR]])
    }
    if(is.na(value_limits[2])){
        value_limits[2] = max(profile_dt[[ct2@value_VAR]])
    }
    n_total = nrow(profile_dt)
    n_over = sum(profile_dt[[ct2@value_VAR]] > max(value_limits))
    n_under = sum(profile_dt[[ct2@value_VAR]] < min(value_limits))
    rng = range(profile_dt[[ct2@value_VAR]])
    top_wasted = max(max(value_limits) - max(rng), 0)
    bottom_wasted = min(min(value_limits) - min(rng), 0)
    limits_total = diff(range(value_limits))
    if((n_under + n_over) > .2*n_total){
        msg = paste0(round(100*(n_under + n_over) / n_total, 2), "% of signal values are outside of value_limits. Please verify.")
        warning(msg)
    }
    if((top_wasted + bottom_wasted)/limits_total > .2){
        msg = paste0(round(100*(top_wasted + bottom_wasted)/limits_total, 2), "% of value_limits spans no signal value. Please verify.")
        warning(msg)
    }

    plot_summary_profiles(
        profile_dt = profile_dt,
        position_dt = position_dt,
        position_VAR = position_VAR,
        value_VAR = value_VAR,
        region_VAR = region_VAR,
        color_VAR = color_VAR,
        x_bins = x_bins,
        y_bins = x_bins,
        extra_vars = extra_vars,
        xrng = xrng,
        yrng = yrng,
        value_limits = value_limits,
        ma_size = ma_size,
        n_splines = n_splines,
        p = p,
        N_floor = N_floor,
        N_ceiling = N_ceiling,
        min_size = min_size,
        return_data = return_data)
}

#' plotDimReduceSummaryProfiles
#'
#' @param ct2
#' @param color_VAR
#' @param x_bins
#' @param y_bins
#' @param extra_vars
#' @param xrng
#' @param yrng
#' @param value_limits
#' @param ma_size
#' @param n_splines
#' @param p
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#' @param return_data
#'
#' @return
#'
#' @examples
#' @export
setGeneric("plotDimReduceSummaryProfiles", function(
        ct2,
        color_VAR = ct2@name_VAR,
        x_bins = nrow(ct2)^.5,
        y_bins = x_bins,
        extra_vars = character(),
        xrng = NULL,
        yrng = NULL,
        value_limits = c(0, 1),
        ma_size = 2,
        n_splines = 10,
        p = NULL,
        N_floor = 0,
        N_ceiling = NULL,
        min_size = 0.3,
        return_data = FALSE)
        standardGeneric("plotDimReduceSummaryProfiles"),
        signature = "ct2")

#' @export
setMethod("plotDimReduceSummaryProfiles", c("ChIPtsne2"), .plotDimReduceSummaryProfiles)

plot_summary_profiles = function (profile_dt,
                                  position_dt,
                                  position_VAR,
                                  value_VAR,
                                  region_VAR,
                                  color_VAR,
                                  x_bins = 8,
                                  y_bins = x_bins,
                                  extra_vars = character(),
                                  xrng = NULL,
                                  yrng = NULL,
                                  value_limits = c(0, 1),
                                  ma_size = 2,
                                  n_splines = 10,
                                  p = NULL,
                                  N_floor = 0,
                                  N_ceiling = NULL,
                                  min_size = 0.3,
                                  return_data = FALSE) {



    if(is.null(xrng)){
        xrng = range(position_dt$tx)
    }
    if(is.null(yrng)){
        yrng = range(position_dt$ty)
    }
    summary_dt = prep_summary(profile_dt = profile_dt,
                              position_dt = position_dt,
                              x_bins = x_bins,
                              y_bins = y_bins,
                              xrng = xrng,
                              yrng = yrng,
                              ma_size = 2,
                              n_splines = 10,
                              facet_by = NULL,
                              position_VAR = position_VAR,
                              value_VAR = value_VAR,
                              region_VAR = region_VAR,
                              color_VAR = color_VAR,
                              extra_vars = extra_vars)
    plot_summary_glyph(
        summary_dt = summary_dt,
        p = p,
        xrng = xrng,
        yrng = yrng,
        x_bins = x_bins,
        y_bins = y_bins,
        value_limits = value_limits,
        N_floor = N_floor,
        N_ceiling = N_ceiling,
        min_size = min_size,
        return_data = return_data,
        position_VAR = position_VAR,
        value_VAR = value_VAR,
        color_VAR = color_VAR)
}

prep_summary = function (profile_dt,
                         position_dt,
                         x_bins,
                         y_bins = x_bins,
                         xrng = range(position_dt$tx),
                         yrng = range(position_dt$ty),
                         ma_size = 2,
                         n_splines = 10,
                         facet_by = NULL,
                         position_VAR = "x",
                         value_VAR = "y",
                         region_VAR = "id",
                         color_VAR = "name",
                         extra_vars = character()){
    position_dt = data.table::copy(position_dt[tx >= min(xrng) & tx <= max(xrng) &
                                                   ty >= min(yrng) & ty <= max(yrng)])
    position_dt = position_dt[get(region_VAR) %in% unique(profile_dt[[region_VAR]])]
    if (is.null(position_dt$bx))
        position_dt[, `:=`(bx, bin_values(tx, x_bins,
                                          xrng = xrng))]
    if (is.null(position_dt$by))
        position_dt[, `:=`(by, bin_values(ty, y_bins,
                                          xrng = yrng))]
    summary_dt = merge(profile_dt, position_dt[, c("bx", "by", region_VAR), with = FALSE],
                       allow.cartesian = TRUE,
                       by = intersect(colnames(profile_dt), c(region_VAR)))
    if (is.null(summary_dt[[color_VAR]]))
        summary_dt[[color_VAR]] = "signal"
    if (is.null(facet_by)) {
        summary_dt = summary_dt[, list(y_tmp_ = mean(get(value_VAR))), c(unique(c("bx", "by", position_VAR, color_VAR, extra_vars)))]
    }
    else {
        summary_dt = summary_dt[, list(y_tmp_ = mean(get(value_VAR))), c(unique(c("bx", "by", position_VAR, color_VAR, facet_by, extra_vars)))]
    }
    data.table::setnames(summary_dt, "y_tmp_", value_VAR)

    merge_var_names = c(unique(c("bx", "by", intersect(extra_vars, colnames(position_dt)))))
    N_dt = position_dt[, .(.N), by = merge_var_names]
    summary_dt = merge(summary_dt, N_dt, by = merge_var_names)
    summary_dt[, `:=`(plot_id, paste(bx, by, sep = "_"))]

    #
    summary_dt = summary_dt[order(get(position_VAR))][order(get(color_VAR))][order(plot_id)]

    if(ma_size > 1){
        summary_dt[, y_ := applyMovingAverage(x = get(value_VAR), ma_size), c("plot_id", color_VAR)]
        summary_dt[[value_VAR]] = NULL
        data.table::setnames(summary_dt, "y_", value_VAR)
    }
    if(n_splines > 1){
        suppressWarnings({
            summary_dt = seqsetvis::applySpline(
                summary_dt,
                n_splines,
                x_ = position_VAR,
                y_ = value_VAR,
                by_ = c("plot_id", color_VAR))
        })
    }

    summary_dt[]
}


#' plot_summary_glyph
#'
#' @param summary_dt results from prep_summary()
#' @param x_bins numeric.  number of grid points to use in x dimension.
#' @param y_bins numeric.  number of grid points to use in y dimension.
#' @param xrng view domain in x dimension.
#' @param yrng view domain in y dimension.
#' @param value_limits value_limitsits per glyph
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param N_floor The value of N to consider 0.  bins with N values <= N_floor
#'   will be ignored.
#' @param N_ceiling The value of N to consider 1.  bins with N values >=
#'   N_ceiling will have images drawn at full size.
#' @param min_size Numeric (0, 1]. The minimum size images to draw.  The default
#'   of .3 draws images for all bins with N values >= 30% of the way from
#'   N_floor to N_ceiling.
#' @param return_data if TRUE, data.table that would have been used to create
#'   ggplot is returned instead.
#'
#' @return a ggplot containing glyphs of local profile summaries arranged in
#'   t-sne space.
#' @importFrom GGally glyphs
#'
plot_summary_glyph = function (summary_dt,
                               x_bins,
                               y_bins = x_bins,
                               xrng = c(-0.5, 0.5),
                               yrng = c(-0.5, 0.5),
                               value_limits = NULL,
                               p = NULL,
                               N_floor = 0,
                               N_ceiling = NULL,
                               min_size = 0.3,
                               return_data = FALSE,
                               position_VAR = "x",
                               value_VAR = "y",
                               color_VAR = "name",
                               extra_vars = character())
{
    group_size = gx = gy = gid = NULL
    summary_dt = set_size(summary_dt, N_floor, N_ceiling, size.name = "group_size")
    if(nrow(summary_dt[group_size >= min_size]) == 0){
        "All summary groups would have been removed based on min_size. Relaxing min_size to 0."
        min_size = 0
    }
    summary_dt = summary_dt[group_size >= min_size]
    if (is.null(value_limits)) {
        value_limits = range(summary_dt[[value_VAR]])
    }
    value_limits = range(value_limits)
    data.table::set(summary_dt, i = which(summary_dt[[value_VAR]] < min(value_limits)), j = value_VAR, value = min(value_limits))
    data.table::set(summary_dt, i = which(summary_dt[[value_VAR]] > max(value_limits)), j = value_VAR, value = max(value_limits))
    data.table::set(summary_dt, j = position_VAR, value = summary_dt[[position_VAR]] * summary_dt$group_size)
    data.table::set(summary_dt, j = value_VAR, value = summary_dt[[value_VAR]] * summary_dt$group_size)
    xs = bin_values_centers(x_bins, rng = xrng)
    ys = bin_values_centers(y_bins, rng = yrng)
    summary_dt[, `:=`(tx, xs[bx])]
    summary_dt[, `:=`(ty, ys[by])]
    down_scale = max(summary_dt$group_size)

    #scale height of glyphs based on value_limits
    limits_size = diff(range(value_limits))
    values_size = diff(range(summary_dt$y))

    glyph_dt = data.table::as.data.table(
        GGally::glyphs(summary_dt,
                       x_major = "tx",
                       x_minor = position_VAR,
                       y_major = "ty",
                       y_minor = value_VAR,
                       width = diff(xrng)/x_bins * 0.95 * down_scale,
                       height = diff(yrng)/y_bins * 0.95 * down_scale*values_size/limits_size
        )
    )
    if (return_data) {
        return(glyph_dt)
    }
    if (is.null(p)) {
        p = ggplot()
    }
    glyph_groups = apply(glyph_dt[, c(unique(c("gid", color_VAR, extra_vars))), with = FALSE],
                         1,
                         function(x){
                             paste(x, collapse = " ")
                         }
    )
    data.table::set(glyph_dt, j = "glyph_group", value = glyph_groups)

    color_VAR = ensym(color_VAR)

    p = p + geom_path(data = glyph_dt,
                      aes(x = gx,
                          y = gy,
                          group = glyph_group,
                          color = !!color_VAR)) +
        labs(x = "tx", y = "ty") +
        coord_cartesian(xlim = xrng, ylim = yrng)
    p
}

#### helpers ####

bin_values = function(x, n_bins, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(rescale_capped(x, 0:1, xrng) * (n_bins-.00001))+1
}

bin_values_centers = function (n_bins, rng){
    if (length(rng) != 2)
        rng = range(rng)
    stopifnot(length(rng) == 2)
    xspc = diff(rng)/n_bins/2
    xs = seq(min(rng) + xspc, max(rng) - xspc, diff(rng)/(n_bins))
    xs
}

rescale_capped = function(x, to = c(0,1), from = range(x, na.rm = TRUE, finite = TRUE)){
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}


set_size = function (dt, N_floor, N_ceiling, size.name = "img_size"){
    tmp_var = NULL
    stopifnot("N" %in% colnames(dt))
    if (is.null(N_ceiling)) {
        N_ceiling = max(dt$N)
    }
    dt[, `:=`(tmp_var, N)]
    dt[tmp_var > N_ceiling, `:=`(tmp_var, N_ceiling)]
    dt[tmp_var < N_floor, `:=`(tmp_var, N_floor)]
    dt[, `:=`(tmp_var, tmp_var - N_floor)]
    dt[, `:=`(tmp_var, tmp_var/N_ceiling)]
    dt[[size.name]] = dt$tmp_var
    dt$tmp_var = NULL
    dt
}

applyMovingAverage = function(x, n = 1, centered = TRUE) {

    if (centered) {
        before <- floor((n - 1)/2)
        after <- ceiling((n - 1)/2)
    } else {
        before <- n - 1
        after <- 0
    }

    # Track the sum and count of number of non-NA items
    s <- rep(0, length(x))
    count <- rep(0, length(x))

    # Add the centered data
    new <- x
    # Add to count list wherever there isn't a
    count <- count + (!is.na(new))
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new

    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new <- c(rep(NA, i), x[seq_len(length(x) - i)])

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new <- c(x[(i + 1):length(x)], rep(NA, i))

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # return sum divided by count
    s/count
}
