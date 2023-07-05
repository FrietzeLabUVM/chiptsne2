.heatmap_theme = theme(panel.background = element_blank(), panel.grid = element_blank())
.heatmap_theme.no_y = .heatmap_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

.annotation_theme = theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank())

.prep_names = function(name){
    name = gsub("_", "\n", name)
    name
}

.prep_ids = function(ids, row_, cluster_){
    if(is.data.frame(ids)){
        assign_dt = unique(ids[, c(row_, cluster_), drop = FALSE])
        assign_dt = assign_dt[order(assign_dt[[row_]]),]
        ids = assign_dt[[cluster_]]

    }
    if(is.character(ids)) ids = factor(ids, levels = unique(ids))
    if(is.numeric(ids)) ids = factor(ids)
    if(is.logical(ids)) ids = factor(ids, levels = c("TRUE", "FALSE"))

    stopifnot(is.factor(ids))

    ids
}


add_group_annotation = function(group_ids,
                                xleft = 0,
                                xright = 1,
                                rect_colors = NULL,
                                text_colors = NULL,
                                label_angle = 0,
                                row_ = "id",
                                cluster_ = "group_id",
                                annotation_theme = .annotation_theme,
                                name_FUN = .prep_names,
                                show_legend = FALSE){
    if(is.null(rect_colors)){
        rect_colors = RColorBrewer::brewer.pal(8, "Dark2")

    }
    if(is.null(text_colors)){
        text_colors = rep("black", length(rect_colors))
    }
    group_ids = .prep_ids(group_ids, row_, cluster_)
    n_grps = length(unique(group_ids))
    if(n_grps > length(rect_colors)){
        stop("Not enough rect_colors provided for number of unique group_ids. Defautl supports up to 8 groups.")
    }
    if(is.null(names(rect_colors))){
        rect_colors = rect_colors[seq(n_grps)]
        names(rect_colors) = unique(group_ids)
    }
    if(is.null(names(text_colors))){
        text_colors = text_colors[seq(n_grps)]
        names(text_colors) = unique(group_ids)
    }
    stopifnot(all(group_ids %in% names(rect_colors)))
    stopifnot(all(group_ids %in% names(text_colors)))

    group_ids = rev(group_ids)

    ends = c(which(!(group_ids[-1] == group_ids[-length(group_ids)])), length(group_ids))
    # ends = cumsum(rev(table(group_ids)))
    starts = c(1, ends[-length(ends)] + 1)
    grps = group_ids[starts]
    starts = starts-.5
    ends = ends+.5

    df_rects = data.frame(xmin = xleft,
                          xmax = xright,
                          ymin = starts,
                          ymax = ends,
                          grp = grps)
    if(setequal(names(rect_colors), df_rects$grp)){
        df_rects$fill = rect_colors[df_rects$grp]
    }else{
        stop("invalid rect_colors, check names match group ids.")
    }
    if(setequal(names(text_colors), df_rects$grp)){
        df_rects$color = text_colors[df_rects$grp]
    }else{
        stop("invalid text_colors, check names match group ids")
    }

    df_rects[[cluster_]] = name_FUN(cluster_)
    p = ggplot(df_rects) +
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(group_ids))+.5, expand = FALSE) +
        facet_grid(paste0(".~", cluster_)) +
        labs(fill = cluster_)

    p = p + geom_rect(data = df_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = grp)) +
        scale_fill_manual(values = rect_colors)
    if(!show_legend){
        p = p + guides(fill = "none")
    }
    p + annotation_theme
}
add_group_annotation.legend = function(group_ids,
                                       rect_colors = NULL,
                                       row_ = "id",
                                       cluster_ = "group_id",
                                       plot_format_FUN = NULL,
                                       annotation_theme = .annotation_theme){
    if(is.null(rect_colors)){
        rect_colors = RColorBrewer::brewer.pal(8, "Dark2")
    }
    xleft = 0
    xright = 1
    text_colors = rep("black", length(rect_colors))
    label_angle = 0
    p = add_group_annotation(
        group_ids = group_ids,
        xleft = xleft,
        xright = xright,
        rect_colors = rect_colors,
        text_colors = text_colors,
        label_angle = label_angle,
        row_ = row_,
        cluster_ = cluster_,
        show_legend = TRUE,
        annotation_theme = annotation_theme
    )
    if(!is.null(plot_format_FUN)){
        p = plot_format_FUN(p)
    }
    cowplot::get_legend(p)
}

add_cluster_annotation = function(cluster_ids,
                                  xleft = 0, xright = 1,
                                  rect_colors = NULL,
                                  text_colors = NULL,
                                  text_size = 8,
                                  show_labels = TRUE,
                                  label_angle = 0,
                                  row_ = "id",
                                  cluster_ = "cluster_id",
                                  annotation_theme = .annotation_theme,
                                  name_FUN = .prep_names,
                                  show_legend = FALSE){
    if(is.null(rect_colors)){
        rect_colors = c("black", "gray")
    }
    if(is.null(text_colors)){
        text_colors = rev(rect_colors)
    }
    cluster_ids = .prep_ids(cluster_ids, row_, cluster_)

    ends = cumsum(rev(table(cluster_ids)))
    starts = c(1, ends[-length(ends)] + 1)
    starts = starts - .5
    names(starts) = names(ends)
    ends = ends + .5

    df_rects = data.frame(xmin = xleft,
                          xmax = xright,
                          ymin = starts,
                          ymax = ends,
                          grp = cluster_)
    if(setequal(names(rect_colors), rownames(df_rects))){
    }else{
        rect_colors = rect_colors[seq_len(nrow(df_rects))%%length(rect_colors)+1]
        names(rect_colors) = rownames(df_rects)
    }
    if(setequal(names(text_colors), rownames(df_rects))){
    }else{
        text_colors = text_colors[seq_len(nrow(df_rects))%%length(text_colors)+1]
        names(text_colors) = rownames(df_rects)
    }

    df_rects = df_rects[rev(seq_len(nrow(df_rects))),]
    cluster_labels = levels(cluster_ids)
    df_rects[[cluster_]] = name_FUN(cluster_)
    df_rects[["grp"]] = rownames(df_rects)

    p = ggplot(df_rects) +
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(cluster_ids))+.5, expand = FALSE) +
        facet_grid(.~grp)
    p = p + geom_rect(aes(
        xmin = xmin,
        xmax = xmax,
        ymin= ymin,
        ymax = ymax,
        fill = grp
    )) +
        scale_fill_manual(values = rect_colors) +
        facet_grid(paste0(".~", cluster_)) +
        labs(fill = cluster_)
    if(show_labels){
        for(i in seq_len(nrow(df_rects))){
            p = p + annotate("text",
                             x = mean(c(df_rects$xmin[i], df_rects$xmax[i])),
                             y = mean(c(df_rects$ymin[i], df_rects$ymax[i])),
                             label = cluster_labels[i],
                             color = text_colors[rownames(df_rects)[i]],
                             angle = label_angle,
                             size = text_size/ggplot2::.pt)
        }
    }
    if(!show_legend){
        p = p + guides(fill = "none")
    }
    p + annotation_theme
}

add_cluster_annotation.legend = function(cluster_ids,
                                         rect_colors = NULL,
                                         row_ = "id",
                                         cluster_ = "cluster_id",
                                         annotation_theme = .annotation_theme){
    if(is.null(rect_colors)){
        rect_colors = c("black", "gray")
    }
    xleft = 0
    xright = 1
    text_colors = rep("black", length(rect_colors))
    label_angle = 0
    p = add_cluster_annotation(
        cluster_ids = cluster_ids,
        xleft = xleft,
        xright = xright,
        rect_colors = rect_colors,
        text_colors = text_colors,
        label_angle = label_angle,
        row_ = row_,
        cluster_ = cluster_,
        show_legend = TRUE,
        annotation_theme = annotation_theme
    )
    cowplot::get_legend(p)
}

get_rel_widths = function(my_plots, sync_width = FALSE){
    stopifnot(is(my_plots, "list"))
    is_ok = vapply(my_plots, function(x){
        "ggplot" %in% class(x) | "grob" %in% class(x)
    }, FUN.VALUE = TRUE)
    stopifnot(all(is_ok))
    my_grobs = lapply(my_plots, function(x){
        if(grid::is.grob(x)){
            x
        }else{
            ggplotGrob(x)
        }
    })

    if(sync_width){
        val = "widths"
    }else{
        val = "heights"
    }
    dim_sizes = lapply(my_grobs, function(gt){
        gt[[val]]
    })
    rel_widths = sapply(dim_sizes, function(x){
        sum(sapply(x, function(xx){
            as.numeric(xx)
        }))
    })
    rel_widths
}

COLOR_KEY_STRAT = list(
    if_not_sorted = "if_not_sorted",
    except_sort_VAR = "except_sort_VAR",
    all = "all"
)

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
#' @param heatmap_colors Either a scale_fill_* or vector of R colors passed to
#'   scale_fill_gradientn.
#' @param heatmap_format_FUN A function that applies any additional formatting
#'   to the heatmap ggplot.
#' @param annotation_colors Colors to use for annotation, should be a character
#'   vector or list of character vectors mirroring group_VARS for fine control.
#' @param annotation_format_FUN A function that applies any additional
#'   formatting to the annotation ggplots. May be a single function to apply to
#'   every annotation plot or list mirroring group_VARS for finer control.
#' @param annotation_text_size Font size for text appearing in clustered group annotations.
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
#' @param relative_heatmap_width Fraction of final plot width dedicated to the heatmap.
#' @param relative_heatmap_height Fraction of final plot height dedicated to the heatmap.
#' @param return_data If TRUE, return the data.table instead of creating a plot.
#'
#' @return A grob of ggplots assembled using cowplot::plot_grid
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = groupRegionsBySignalCluster(ct2, group_VAR = "cluster")
#' ct2 = groupRegionsByOverlap(ct2, seqsetvis::CTCF_in_10a_narrowPeak_grs[1:2], group_VAR = "overlap")
#' meta_df = getRegionMetaData(ct2)
#' meta_df = meta_df %>% dplyr::mutate(overlap_num = as.numeric(overlap))
#' ct2 = setRegionMetaData(ct2, meta_df)
#' plotSignalHeatmap(ct2, group_VARS = c("cluster", "overlap", "overlap_num"))
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
.plotSignalHeatmap = function(ct2,
                              group_VARS = NULL,
                              sort_VAR = NULL,
                              sort_strategy =  c("hclust", "sort", "left", "right")[2],
                              heatmap_colors = scale_fill_viridis_c(option = "magma"),
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
                              return_data = FALSE){
    if(!is.null(group_VARS) & is.null(sort_VAR)){
        sort_VAR = group_VARS[length(group_VARS)]
    }
    stopifnot(relative_heatmap_width > 0 & relative_heatmap_width < 1)
    stopifnot(relative_heatmap_height > 0 & relative_heatmap_height < 1)
    stopifnot(n_legend_rows >= 1)

    prof_dt = .getTidyProfile(ct2, group_VARS)

    fake_VAR = "__FAKE_CLUSTER__"
    if(sort_VAR == FALSE){
        sort_VAR = fake_VAR
        prof_dt[[sort_VAR]] = 1
    }
    clust_dt = seqsetvis::within_clust_sort(
        prof_dt,
        row_ = ct2@region_VAR,
        fill_ = ct2@value_VAR,
        column_ = ct2@position_VAR,
        facet_ = ct2@name_VAR,
        cluster_ = sort_VAR, dcast_fill = 0,
        within_order_strategy = sort_strategy)

    x_ = ct2@position_VAR
    x_ = ensym(x_)
    y_ = ct2@region_VAR
    y_ = ensym(y_)
    fill_ = ct2@value_VAR
    fill_ = ensym(fill_)

    assign_dt = clust_dt %>% dplyr::select(all_of(c(ct2@region_VAR, sort_VAR))) %>%
        unique
    assign_dt[order(assign_dt$peak_regions)]

    clust_dt[[ct2@region_VAR]] = factor(clust_dt[[ct2@region_VAR]], levels = rev(levels(clust_dt[[ct2@region_VAR]])))
    clust_dt[[ct2@name_VAR]] = name_FUN(clust_dt[[ct2@name_VAR]])
    if(return_data) return(clust_dt)
    p_heat = ggplot(clust_dt, aes(x = !!x_, y = !!y_, fill = !!fill_)) +
        geom_raster() +
        scale_x_discrete(expand = c(0,0)) +
        heatmap_theme +
        facet_grid(paste0(".~", ct2@name_VAR))
    if(is(heatmap_colors, "Scale")){
        p_heat = p_heat + heatmap_colors
    }else{
        p_heat = p_heat + scale_fill_gradientn(colours = heatmap_colors)
    }
    if(!is.null(heatmap_format_FUN)){
        p_heat = heatmap_format_FUN(p_heat)
    }
    p_heat.leg = cowplot::get_legend(p_heat)
    p_heat = p_heat + guides(fill = "none")

    #### annotation ####
    anno_VARS = group_VARS
    anno_VARS = anno_VARS[!anno_VARS %in% fake_VAR]
    anno_df = getRegionMetaData(ct2, anno_VARS)[, c(ct2@region_VAR, anno_VARS)]
    anno_df[[ct2@region_VAR]] = factor(anno_df[[ct2@region_VAR]],
                                       levels = levels(assign_dt[[ct2@region_VAR]]))
    anno_df = anno_df[order(anno_df[[ct2@region_VAR]]),]

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
        }else{
            p_anno = add_group_annotation(
                anno_df,
                row_ = ct2@region_VAR,
                cluster_ = var,
                rect_colors = annotation_colors[[i]]
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



#' @export
setGeneric("plotSignalHeatmap", function(
        ct2,
        group_VARS = NULL,
        sort_VAR = NULL,
        sort_strategy =  c("hclust", "sort", "left", "right")[2],
        heatmap_colors = scale_fill_viridis_c(option = "magma"),
        heatmap_format_FUN = NULL,
        heatmap_theme = .heatmap_theme.no_y,
        annotation_colors = seqsetvis::safeBrew(n = 8),
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
setMethod("plotSignalHeatmap", c("ChIPtsne2"), .plotSignalHeatmap)

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
                               extra_VARS = character(),
                               return_data = FALSE){
    all_VARS = unique(c(group_VAR, color_VAR, facet_VAR, extra_VARS))
    meta_VARS = setdiff(all_VARS, ct2@name_VAR)
    prof_dt = getTidyProfile(ct2, meta_VARS = meta_VARS)
    agg_dt = prof_dt[, .(VALUE_ = mean(get(ct2@value_VAR))), c(unique(c(all_VARS, ct2@position_VAR, ct2@name_VAR)))]
    data.table::setnames(agg_dt, "VALUE_", ct2@value_VAR)
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
    if(return_data) return(agg_dt)

    ggplot(agg_dt,
           aes(x = !!x_, y = !!y_, color = !!color_VAR, group = !!group_)) +
        geom_path() +
        facet_grid(paste0(group_VAR, "~", facet_VAR), labeller = label_both) +
        labs(subtitle = plot_subtitle)
}

#' @export
setGeneric("plotSignalLinePlot", function(
        ct2,
        group_VAR = NULL,
        color_VAR = NULL,
        facet_VAR = ct2@name_VAR,
        extra_VARS = character(),
        return_data = FALSE)
    standardGeneric("plotSignalLinePlot"),
    signature = "ct2")

#' @export
setMethod("plotSignalLinePlot", c("ChIPtsne2"), .plotSignalLinePlot)


