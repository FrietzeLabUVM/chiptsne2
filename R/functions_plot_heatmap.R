
#' @importFrom cowplot get_plot_component plot_grid
#' @importFrom scales pretty_breaks
.plotSignalHeatmap = function(ct2,
                              group_VARS = NULL,
                              sort_VAR = NULL,
                              balance_VAR = NULL,
                              max_rows = 500,
                              sort_strategy =  c("hclust", "sort", "left", "right")[2],
                              heatmap_fill_limits = c(NA, NA),
                              has_symmetrical_limits = NULL,
                              heatmap_colors = NULL,
                              heatmap_format_FUN = NULL,
                              heatmap_theme = .heatmap_theme.no_y,
                              annotation_colors = NULL,
                              annotation_format_FUN = NULL,
                              annotation_theme = .annotation_theme,
                              annotation_text_size = 8,
                              name_FUN = .prep_names,
                              color_key_strategy = c("if_not_sorted", "except_sort_VAR", "all")[2],
                              n_legend_rows = 1,
                              relative_heatmap_width = .5,
                              relative_heatmap_height = .66,
                              return_data = FALSE
){
    if(!is.null(group_VARS) & is.null(sort_VAR)){
        sort_VAR = group_VARS[length(group_VARS)]
    }
    if(is.null(sort_VAR)){
        sort_VAR = FALSE
    }
    stopifnot(relative_heatmap_width > 0 & relative_heatmap_width < 1)
    stopifnot(relative_heatmap_height > 0 & relative_heatmap_height < 1)
    stopifnot(n_legend_rows >= 1)

    meta_dt = getRegionMetaData(ct2)
    req_vars = unique(c(group_VARS, sort_VAR, balance_VAR))
    req_vars = setdiff(req_vars, FALSE)
    if(!all(req_vars %in% colnames(meta_dt))){
        stop(paste(c("Missing variables from region metadata:",
                     setdiff(req_vars, colnames(meta_dt))), collapse = "\n"))
    }
    fake_VAR = "__FAKE_CLUSTER__"
    if(sort_VAR == FALSE){
        sort_VAR = fake_VAR
        meta_dt[[sort_VAR]] = 1
    }
    if(!is.null(balance_VAR)){
        group_dt = unique(meta_dt[, c(ct2@region_VAR, balance_VAR), with = FALSE])
        group_l = split(as.character(group_dt[[ct2@region_VAR]]), group_dt[[balance_VAR]])
        min_group = min(lengths(group_l))
        if((length(group_l)*min_group) > max_rows){
            min_group = floor(max_rows / length(group_l))
        }
        min_l = lapply(group_l, function(x){
            sample(x, min_group)
        })
        all_ids = unlist(min_l)
        names(all_ids) = NULL
    }else{
        all_ids = unique(meta_dt[[ct2@region_VAR]])
    }
    if(!is.infinite(max_rows)){
        if(max_rows < length(all_ids)){
            all_ids = sample(all_ids, max_rows)
        }
    }
    meta_dt = dplyr::filter(meta_dt, get(ct2@region_VAR) %in% all_ids)
    if(is.factor(meta_dt[[ct2@region_VAR]])){
        meta_dt[[ct2@region_VAR]] = droplevels(meta_dt[[ct2@region_VAR]])
    }
    if(!is.factor(meta_dt[[sort_VAR]])){
        meta_dt[[sort_VAR]] = factor(meta_dt[[sort_VAR]])
    }
    #### Fetch profile ####
    tidy_vars = unique(c(group_VARS, sort_VAR, balance_VAR))
    tidy_vars = setdiff(tidy_vars, '__FAKE_CLUSTER__')

    prof_dt = getTidyProfile(ct2[meta_dt[[ct2@region_VAR]],], tidy_vars)

    heatmap_colors = .prep_color_scale(values = prof_dt[[ct2@value_VAR]], has_symmetrical_limits = has_symmetrical_limits, color_scale = heatmap_colors)
    heatmap_fill_limits = .prep_symmetrical(values = prof_dt[[ct2@value_VAR]], has_symmetrical_limits, heatmap_fill_limits)

    if(sort_VAR == fake_VAR){
        prof_dt[[sort_VAR]] = 1
    }
    if(is.numeric(prof_dt[[sort_VAR]])){
        id_lev = prof_dt[order(prof_dt[[sort_VAR]]),][[ct2@region_VAR]] %>%
            as.character %>%
            unique
        clust_dt = data.table::copy(prof_dt)
        clust_dt[[ct2@region_VAR]] = factor(clust_dt[[ct2@region_VAR]], levels = id_lev)
    }else{
        clust_dt = seqsetvis::within_clust_sort(
            prof_dt,
            row_ = ct2@region_VAR,
            fill_ = ct2@value_VAR,
            column_ = ct2@position_VAR,
            facet_ = ct2@name_VAR,
            cluster_ = sort_VAR, dcast_fill = 0,
            within_order_strategy = sort_strategy)
    }


    x_ = ct2@position_VAR
    x_ = ensym(x_)
    y_ = ct2@region_VAR
    y_ = ensym(y_)
    fill_ = ct2@value_VAR
    fill_ = ensym(fill_)

    assign_dt = clust_dt %>% dplyr::select(dplyr::all_of(c(ct2@region_VAR, sort_VAR))) %>%
        unique

    clust_dt[[ct2@region_VAR]] = factor(clust_dt[[ct2@region_VAR]], levels = rev(levels(clust_dt[[ct2@region_VAR]])))
    clust_dt[[ct2@name_VAR]] = name_FUN(clust_dt[[ct2@name_VAR]])
    # apply limits
    if(!is.na(heatmap_fill_limits[1])){
        #TODO
        clust_dt[[ct2@value_VAR]] = ifelse(
            clust_dt[[ct2@value_VAR]] < heatmap_fill_limits[1],
            heatmap_fill_limits[1],
            clust_dt[[ct2@value_VAR]]
        )
    }
    if(!is.na(heatmap_fill_limits[2])){
        clust_dt[[ct2@value_VAR]] = ifelse(
            clust_dt[[ct2@value_VAR]] > heatmap_fill_limits[2],
            heatmap_fill_limits[2],
            clust_dt[[ct2@value_VAR]]
        )
    }
    if(return_data) return(clust_dt)
    p_heat = ggplot(clust_dt, aes(x = !!x_, y = !!y_, fill = !!fill_)) +
        geom_raster() +
        scale_x_continuous(expand = c(0,0), breaks = scales::pretty_breaks(n = 3)) +
        heatmap_theme +
        facet_grid(paste0(".~", ct2@name_VAR)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
    p_heat = .apply_scale(p_heat, heatmap_colors, heatmap_fill_limits)
    if(!is.null(heatmap_format_FUN)){
        p_heat = heatmap_format_FUN(p_heat)
    }
    # cowplot::get_legend() now returning warning
    p_heat.leg = cowplot::get_plot_component(p_heat, "guide-box", return_all = TRUE)[[1]]
    p_heat = p_heat + guides(fill = "none")

    #### annotation ####
    anno_VARS = group_VARS
    anno_VARS = anno_VARS[!anno_VARS %in% fake_VAR]
    anno_df = getRegionMetaData(ct2, anno_VARS)[, c(ct2@region_VAR, anno_VARS), drop = FALSE]
    # anno_df = anno_df[anno_df[[ct2@region_VAR]] %in% assign_dt[[ct2@region_VAR]], , drop = FALSE]
    anno_df = dplyr::filter(anno_df, get(ct2@region_VAR) %in% assign_dt[[ct2@region_VAR]])
    anno_df[[ct2@region_VAR]] = factor(anno_df[[ct2@region_VAR]],
                                       levels = levels(assign_dt[[ct2@region_VAR]]))
    anno_df = anno_df[order(anno_df[[ct2@region_VAR]]),, drop = FALSE]

    if(is.null(annotation_format_FUN)){
        annotation_format_FUN = function(p)p
    }
    if(!is.list(annotation_format_FUN)){
        annotation_format_FUN = lapply(seq_along(anno_VARS), function(i){
            annotation_format_FUN
        })
    }
    stopifnot(length(annotation_format_FUN) == length(anno_VARS))
    if(!is.list(annotation_colors)){
        annotation_colors = lapply(seq_along(anno_VARS), function(i){
            annotation_colors
        })
    }
    stopifnot(length(annotation_colors) == length(anno_VARS))

    anno_plots = list()
    legend_plots = list()
    for(i in seq_along(anno_VARS)){
        var = anno_VARS[i]
        anno_rle = rle(as.character(anno_df[[var]]))
        is_sorted = length(anno_rle$values) == length(unique(anno_rle$values))
        if(color_key_strategy  == COLOR_KEY_STRAT$if_not_sorted){
            # if(assume_sorted_is_clustered){
            is_clustered = is_sorted
        }else if(color_key_strategy == COLOR_KEY_STRAT$except_sort_VAR){
            #only the final instance of duplicate var should be considered clustered
            is_clustered = (var == sort_VAR) & (i == max(which(var == anno_VARS)))
        }else if(color_key_strategy == COLOR_KEY_STRAT$all){
            is_clustered = FALSE
        }else{
            stop("Unrecognized color_key_strategy: ", color_key_strategy,
                 "\nAllowed values are:\n\"",
                 paste(as.character(COLOR_KEY_STRAT), collapse = "\"\n\""), "\"")
        }
        if(is_clustered){
            p_anno = add_cluster_annotation(
                anno_df,
                row_ = ct2@region_VAR,
                cluster_ = var,
                rect_colors = annotation_colors[[i]],
                annotation_theme = annotation_theme,
                text_size = annotation_text_size,
            )
            p_anno = annotation_format_FUN[[i]](p_anno)
            if(is.numeric(anno_df[[var]])){
                p_leg = add_cluster_annotation.legend(
                    anno_df,
                    row_ = ct2@region_VAR,
                    cluster_ = var,
                    rect_colors = annotation_colors[[i]],
                    plot_format_FUN = annotation_format_FUN[[i]],
                    annotation_theme = annotation_theme
                )
                legend_plots[[length(legend_plots) + 1]] = p_leg
            }
        }else{
            p_anno = add_group_annotation(
                anno_df,
                row_ = ct2@region_VAR,
                cluster_ = var,
                rect_colors = annotation_colors[[i]],
                annotation_theme = annotation_theme
            )
            p_anno = annotation_format_FUN[[i]](p_anno)
            p_leg = add_group_annotation.legend(
                anno_df,
                row_ = ct2@region_VAR,
                cluster_ = var,
                rect_colors = annotation_colors[[i]],
                plot_format_FUN = annotation_format_FUN[[i]],
                annotation_theme = annotation_theme
            )
            legend_plots[[length(legend_plots) + 1]] = p_leg
        }
        anno_plots[[length(anno_plots) + 1]] = p_anno
    }
    legend_plots = c(legend_plots, list(p_heat.leg))
    rel_widths = c(1 - relative_heatmap_width, relative_heatmap_width)
    rel_widths = c(rep(rel_widths[1] / length(anno_plots), length(anno_plots)), rel_widths[2])
    row1 = seqsetvis::assemble_heatmap_cluster_bars(
        c(anno_plots, list(p_heat)),
        rel_widths = rel_widths
    )

    leg_grps = ceiling((seq_along(legend_plots) / length(legend_plots)) * n_legend_rows)
    legend_plots.sp = split(legend_plots, leg_grps)
    leg_rows = lapply(legend_plots.sp, function(x){
        leg_rel_widths = get_rel_widths(x, sync_width = TRUE)
        row = cowplot::plot_grid(
            plotlist = x, nrow = 1, rel_widths = leg_rel_widths
        )
        row
    })
    row2 = cowplot::plot_grid(plotlist = leg_rows, ncol = 1)
    rel_heights = c(relative_heatmap_height, 1 - relative_heatmap_height)
    cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = rel_heights)
}


#' plotSignalHeatmap
#'
#' @param ct2 A ChIPtsne2 object
#' @param group_VARS Optional region metadata variables
#' @param sort_VAR Optional region metadata variable to use for sorting data.
#'   Default to last entry in group_VARS.
#' @param sort_strategy Strategy to use for sorting within groups. Valid choices
#'   are: 1) "sort", which sorts decreasing top to bottom, 2) "hclust" which
#'   uses hierarchical clustering, 3) "left" which puts most left tiled profiles
#'   at top, and 4) "right" which puts most right tilted profiles at top.
#' @param heatmap_fill_limits Passed to limits of scale_fill_gradientn. Default
#'   of c(NA, NA) uses natural range of the data.
#' @param heatmap_colors Either a scale_fill_* or vector of R colors passed to
#'   scale_fill_gradientn.
#' @param heatmap_format_FUN A function that applies any additional formatting
#'   to the heatmap ggplot.
#' @param annotation_colors Colors to use for annotation, should be a character
#'   vector or list of character vectors mirroring group_VARS for fine control.
#' @param annotation_format_FUN A function that applies any additional
#'   formatting to the annotation ggplots. May be a single function to apply to
#'   every annotation plot or list mirroring group_VARS for finer control.
#' @param annotation_text_size Font size for text appearing in clustered group
#'   annotations.
#' @param name_FUN A function to apply to facet labels (heatmap and annotation
#'   columns). By default, underscores are replaced with newlines.
#' @param color_key_strategy Strategy to use for selection of annotation method.
#'   The 2 annotation methods are arbitrary group annotation (has legend and no
#'   labels) and clustered group annotation (no legend, labels on boxes).  All
#'   members of a clustered group must be contiguous. Valid choices are: 1)
#'   "if_not_sorted" applies clustered method whenever groups are sorted 2)
#'   "except_sort_VAR" only uses clustered method for final instance of sort_VAR
#'   3) "all" treats all as arbitrary groupings.
#' @param n_legend_rows How many rows to split the region legends into. Default
#'   is 1.
#' @param relative_heatmap_width Fraction of final plot width dedicated to the
#'   heatmap.
#' @param relative_heatmap_height Fraction of final plot height dedicated to the
#'   heatmap.
#' @param return_data If TRUE, return the data.table instead of creating a plot.
#' @param balance_VAR When set to a categorical region metadata variable,
#'   heatmap grouping will not evenly select among all regions but instead try
#'   to pick comparable numbers of regions between groups defined by
#'   `balance_VAR`.
#' @param max_rows The maximum number of rows in the heatmap. Heatmaps should not have rows allocated less than 1 pixel. Default is 500.
#' @param has_symmetrical_limits If TRUE color scale limits will extend to equal magnitude in positive and negative direction. Default is TRUE when negative values are present and FALSE otherwise.
#' @param heatmap_theme Theme applied to heatmap. Set to NULL if you wish to disable themeing. This is useful if you set the theme prior to plotting.
#' @param annotation_theme Theme applied to annotation plots. Set to NULL if you wish to disable themeing. This is useful if you set the theme prior to plotting.
#'
#' @return A grob of ggplots assembled using cowplot::plot_grid
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = groupRegionsBySignalCluster(ct2, group_VAR = "cluster")
#' ct2 = groupRegionsByOverlap(ct2, seqsetvis::CTCF_in_10a_narrowPeak_grs[1:2], group_VAR = "overlap")
#'
#' meta_df = getRegionMetaData(ct2)
#' meta_df = meta_df %>% dplyr::mutate(overlap_num = as.numeric(overlap))
#' meta_df$chr_num = as.numeric(GenomicRanges::seqnames(rowRanges(ct2)))
#' ct2 = setRegionMetaData(ct2, meta_df)
#'
#' plotSignalHeatmap(ct2)
#'
#' plotSignalHeatmap(ct2, group_VARS = c("cluster", "overlap", "chr_num"))
#' plotSignalHeatmap(ct2, group_VARS = c("cluster", "chr_num", "overlap"), sort_VAR = "chr_num")
#' plotSignalHeatmap(ct2, group_VARS = c("overlap", "cluster"))
#' plotSignalHeatmap(
#'     ct2,
#'     group_VARS = c(
#'         "overlap",
#'         "peak_MCF10A_CTCF",
#'         "peak_MCF10AT1_CTCF",
#'         "peak_MCF10CA1_CTCF",
#'         "cluster",
#'         "overlap"
#'     ),
#'     color_key_strategy = "if_not_sorted"
#' )
#' ct2_diff = subsetCol(ct2, cell == "MCF10AT1") - subsetCol(ct2, cell == "MCF10A")
#' plotSignalHeatmap(ct2_diff, group_VARS = c("overlap", "cluster"))
#'
#' plotSignalHeatmap(
#'     ct2,
#'     group_VARS = c(
#'         "cluster",
#'         "peak_MCF10A_CTCF",
#'         "peak_MCF10AT1_CTCF",
#'         "peak_MCF10CA1_CTCF",
#'         "overlap",
#'         "cluster",
#'         "cluster"
#'     ),
#'     color_key_strategy = "all"
#' )
#'
#' plotSignalHeatmap(
#'     ct2,
#'     group_VARS = c(
#'         "cluster",
#'         "peak_MCF10A_CTCF",
#'         "peak_MCF10AT1_CTCF",
#'         "peak_MCF10CA1_CTCF",
#'         "overlap",
#'         "cluster",
#'         "cluster"
#'     ), n_legend_rows = 2, relative_heatmap_height = .5,
#'     color_key_strategy = "except_sort_VAR"
#' )
#'
#' plotSignalHeatmap(ct2, group_VARS = c("overlap",
#'                                        "peak_MCF10A_CTCF",
#'                                        "peak_MCF10AT1_CTCF",
#'                                        "peak_MCF10CA1_CTCF",
#'                                        "cluster",
#'                                        "overlap"), sort_VAR = FALSE, n_legend_rows = 2)
#'
#' plotSignalHeatmap(ct2, group_VARS = c("cluster", "overlap"))
#' plotSignalHeatmap(ct2, group_VARS = c("overlap", "cluster"))
#' plotSignalHeatmap(
#'     ct2,
#'     group_VARS = c(
#'         "overlap",
#'         "peak_MCF10A_CTCF",
#'         "peak_MCF10AT1_CTCF",
#'         "peak_MCF10CA1_CTCF",
#'         "cluster",
#'         "overlap"
#'     ),
#'     color_key_strategy = "if_not_sorted"
#' )
#'
#' plotSignalHeatmap(
#'     ct2,
#'     group_VARS = c(
#'         "cluster",
#'         "peak_MCF10A_CTCF",
#'         "peak_MCF10AT1_CTCF",
#'         "peak_MCF10CA1_CTCF",
#'         "overlap",
#'         "cluster",
#'         "cluster"
#'     ),
#'     color_key_strategy = "all"
#' )
#'
#' plotSignalHeatmap(
#'     ct2,
#'     group_VARS = c(
#'         "cluster",
#'         "peak_MCF10A_CTCF",
#'         "peak_MCF10AT1_CTCF",
#'         "peak_MCF10CA1_CTCF",
#'         "overlap",
#'         "cluster",
#'         "cluster"
#'     ),
#'     color_key_strategy = "except_sort_VAR"
#' )
#'
#' plotSignalHeatmap(ct2, group_VARS = c("overlap",
#'                                        "peak_MCF10A_CTCF",
#'                                        "peak_MCF10AT1_CTCF",
#'                                        "peak_MCF10CA1_CTCF",
#'                                        "cluster",
#'                                        "overlap"), sort_VAR = FALSE, n_legend_rows = 2)
setGeneric("plotSignalHeatmap", function(
        ct2,
        group_VARS = NULL,
        sort_VAR = NULL,
        balance_VAR = NULL,
        max_rows = 500,
        sort_strategy =  c("hclust", "sort", "left", "right")[2],
        heatmap_fill_limits = c(NA, NA),
        has_symmetrical_limits = NULL,
        heatmap_colors = NULL,
        heatmap_format_FUN = NULL,
        heatmap_theme = .heatmap_theme.no_y,
        annotation_colors = NULL,
        annotation_format_FUN = NULL,
        annotation_theme = .annotation_theme,
        annotation_text_size = 8,
        name_FUN = .prep_names,
        color_key_strategy = c("if_not_sorted", "except_sort_VAR", "all")[2],
        n_legend_rows = 1,
        relative_heatmap_width = .5,
        relative_heatmap_height = .66,
        return_data = FALSE)
        standardGeneric("plotSignalHeatmap"),
        signature = "ct2")

#' @export
setMethod("plotSignalHeatmap", c("ChIPtsne2_no_rowRanges"), .plotSignalHeatmap)
