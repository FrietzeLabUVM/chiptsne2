
library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()
ct2@position_VAR = "bp"
ct2@value_VAR = "reads"
ct2@region_VAR = "peak_regions"
group_VAR = "cluster"
n_to_plot = 500
.heatmap_theme = theme(panel.background = element_blank(), panel.grid = element_blank())
.heatmap_theme.no_y = .heatmap_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

.annotation_theme = theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank())

plot_theme = .heatmap_theme.no_y
# undebug(groupRegionsBySignalCluster, "ChIPtsne2")
ct2 = groupRegionsBySignalCluster(ct2, group_VAR = "cluster")
ct2 = groupRegionsByOverlap(ct2, seqsetvis::CTCF_in_10a_narrowPeak_grs[1:2], group_VAR = "overlap")
colData(ct2)
rowRanges(ct2)
prof_dt = getTidyProfile.with_meta(ct2, region_meta_VARS = group_VAR)
prof_dt = seqsetvis::within_clust_sort(prof_dt, cluster_ = group_VAR, fill_ = ct2@value_VAR, column_ = ct2@position_VAR, row_ = ct2@region_VAR)

x_ = ct2@position_VAR
x_ = ensym(x_)
y_ = ct2@region_VAR
y_ = ensym(y_)
fill_ = ct2@value_VAR
fill_ = ensym(fill_)

plot_dt = data.table::copy(prof_dt)
plot_dt[[ct2@region_VAR]] = factor(plot_dt[[ct2@region_VAR]], levels = rev(levels(prof_dt[[ct2@region_VAR]])))
p_heat = ggplot(plot_dt, aes(x = !!x_, y = !!y_, fill = !!fill_)) +
    geom_raster() +
    scale_x_discrete(expand = c(0,0)) +
    plot_theme +
    facet_grid(paste0(".~", ct2@name_VAR))
p_heat

x_per_strata = .1
anno_VARS = c("cluster", "overlap")
anno_df = getRegionMetaData(ct2)[, c(ct2@region_VAR, anno_VARS)]
anno_df[[ct2@region_VAR]] = factor(anno_df[[ct2@region_VAR]], levels = levels(prof_dt[[ct2@region_VAR]]))
anno_df = anno_df[order(anno_df[[ct2@region_VAR]]),]
group_ids = anno_df$overlap
cluster_ = "overlap"
row_ = "peak_regions"

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
                                rect_colors = RColorBrewer::brewer.pal(8, "Dark2"),
                                text_colors = rep("black", length(rect_colors)),
                                label_angle = 0,
                                row_ = "id",
                                cluster_ = "group_id",
                                annotation_theme = .annotation_theme,
                                name_FUN = .prep_names,
                                show_legend = FALSE){
    group_ids = .prep_ids(group_ids, row_, cluster_)
    # if(is.data.frame(group_ids)){
    #     assign_dt = unique(group_ids[, c(row_, cluster_), drop = FALSE])[order(assign_dt[[row_]]),]
    #     group_ids = assign_dt[[cluster_]]
    #
    # }
    # if(is.character(group_ids)) group_ids = factor(group_ids, levels = unique(group_ids))
    # if(is.numeric(group_ids)) group_ids = factor(group_ids)
    # if(is.logical(group_ids)) group_ids = factor(group_ids, levels = c("TRUE", "FALSE"))
    #
    # stopifnot(is.factor(group_ids))
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
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(group_ids)), expand = FALSE) +
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
                                       rect_colors = RColorBrewer::brewer.pal(8, "Dark2"),
                                       row_ = "id",
                                       cluster_ = "group_id",
                                       plot_format_FUN = NULL){
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
        show_legend = TRUE
    )
    if(!is.null(plot_format_FUN)){
        p = plot_format_FUN(p)
    }
    cowplot::get_legend(p)
}

add_group_annotation(group_ids, cluster_ = "overlap")
add_group_annotation(group_ids, annotation_theme = NULL)

p_overlap = add_group_annotation(anno_df, cluster_ = "overlap", row_ = ct2@region_VAR)
p_overlap.leg = add_group_annotation.legend(anno_df, cluster_ = "overlap", row_ = ct2@region_VAR)
p_overlap + cowplot::draw_grob(p_overlap.leg, x = 0, y = 50)

p_cluster_bar2 = seqsetvis::add_cluster_annotation(anno_df$cluster)
p_cluster_bar2 = seqsetvis::add_cluster_annotation(anno_df$overlap)

assign_dt = anno_df
rect_colors = seqsetvis::safeBrew(anno_df$cluster)
text_colors = "black"
show_labels = TRUE
label_angle = 0

add_cluster_annotation = function(cluster_ids,
                                  xleft = 0, xright = 1,
                                  rect_colors = c("black", "gray"),
                                  text_colors = rev(rect_colors),
                                  show_labels = TRUE,
                                  label_angle = 0,
                                  row_ = "id",
                                  cluster_ = "cluster_id",
                                  annotation_theme = .annotation_theme,
                                  name_FUN = .prep_names,
                                  show_legend = FALSE){
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
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(cluster_ids)), expand = FALSE) +
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
                             color = text_colors[rev(rownames(df_rects))[i]],
                             angle = label_angle)
        }
    }
    if(!show_legend){
        p = p + guides(fill = "none")
    }
    p + annotation_theme
}

add_cluster_annotation.legend = function(cluster_ids,
                                         rect_colors = c("black", "gray"),
                                         row_ = "id",
                                         cluster_ = "cluster_id"){
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
        show_legend = TRUE
    )
    cowplot::get_legend(p)
}

p_cluster_bar = add_cluster_annotation(
    cluster_ids = anno_df[[group_VAR]],
    # p = p_cluster_bar,
    xleft = 0,
    xright = 1,
    rect_colors = rect_colors,
    text_colors = "black",
    show_labels = show_labels,
    label_angle = label_angle)
p_cluster_bar

p_cluster_bar.leg = add_cluster_annotation.legend(
    cluster_ids = anno_df[[group_VAR]],
    rect_colors = rect_colors)

p_cluster_bar2 = add_cluster_annotation(
    cluster_ids = anno_df[[group_VAR]],
)
p_cluster_bar2

p_cluster_bar2.leg = add_cluster_annotation.legend(
    anno_df[[group_VAR]],
    row_ = ct2@region_VAR)

p_cluster_bar2 + cowplot::draw_grob(p_cluster_bar2.leg, y = 50)



overlap_cols = c("gray", "MCF10A_CTCF" = "red", "MCF10AT1_CTCF" = "blue", "MCF10A_CTCF & MCF10AT1_CTCF" = "purple")
p_overlap = add_group_annotation(
    anno_df,
    cluster_ = "overlap",
    row_ = ct2@region_VAR,
    rect_colors = overlap_cols)
p_overlap.leg = add_group_annotation.legend(
    anno_df,
    cluster_ = "overlap",
    row_ = ct2@region_VAR,
    rect_colors = overlap_cols)


p_heat.leg = cowplot::get_legend(p_heat)
p_heat = p_heat + guides(fill = "none")

row1 = seqsetvis::assemble_heatmap_cluster_bars(list(p_overlap, p_cluster_bar, p_cluster_bar2, p_heat), rel_widths = c(1, 1, 1, 5))

plist = list(p_overlap.leg, p_cluster_bar.leg, p_cluster_bar2.leg, p_heat.leg)
my_plots = plist
sync_width = TRUE
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
leg_rel_widths = get_rel_widths(plist, sync_width = TRUE)

row2 = cowplot::plot_grid(
    p_overlap.leg, p_cluster_bar.leg, p_cluster_bar2.leg, p_heat.leg, nrow = 1, rel_widths = leg_rel_widths
)
cowplot::plot_grid(row1, row2, ncol = 1, rel_heights = c(2, 1))

.getTidyProfile = chiptsne2:::.getTidyProfile


COLOR_KEY_STRAT = list(
    if_not_sorted = "if_not_sorted",
    except_sort_VAR = "except_sort_VAR",
    all = "all"
)

.plotSignalHeatmap = function(ct2,
                              group_VARS = NULL,
                              sort_VAR = NULL,
                              sort_strategy =  c("hclust", "sort", "left", "right")[2],
                              heatmap_colors = scale_fill_viridis_c(option = "magma"),
                              heatmap_format_FUN = NULL,
                              annotation_colors = seqsetvis::safeBrew(n = 8),
                              annotation_format_FUN = NULL,
                              name_FUN = .prep_names,
                              color_key_strategy = c("if_not_sorted", "except_sort_VAR", "all")[2],
                              # assume_sorted_is_clustered = TRUE,
                              n_legend_rows = 1){
    if(!is.null(group_VARS) & is.null(sort_VAR)){
        sort_VAR = group_VARS[length(group_VARS)]
    }

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
    p_heat = ggplot(clust_dt, aes(x = !!x_, y = !!y_, fill = !!fill_)) +
        geom_raster() +
        scale_x_discrete(expand = c(0,0)) +
        plot_theme +
        facet_grid(paste0(".~", ct2@name_VAR))
    if(is(heatmap_colors, "Scale")){
        p_heat = p_heat + heatmap_colors
    }else{
        p_heat = p_heat + scale_color_gradientn(colours = heatmap_colors)
    }
    if(!is.null(heatmap_format_FUN)){
        p_heat = heatmap_format_FUN(p_heat)
    }
    p_heat.leg = cowplot::get_legend(p_heat)
    p_heat = p_heat + guides(fill = "none")

    #### annotation ####
    # browser()
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
                 "\nAllowed values are:\n\"", paste(as.character(COLOR_KEY_STRAT), collapse = "\"\n\""), "\"")
        }
        if(is_clustered){
            p_anno = add_cluster_annotation(anno_df, row_ = ct2@region_VAR, cluster_ = var, rect_colors = annotation_colors[[i]])
            p_anno = annotation_format_FUN[[i]](p_anno)
        }else{
            p_anno = add_group_annotation(anno_df, row_ = ct2@region_VAR, cluster_ = var, rect_colors = annotation_colors[[i]])
            p_anno = annotation_format_FUN[[i]](p_anno)
            p_leg = add_group_annotation.legend(anno_df, row_ = ct2@region_VAR, cluster_ = var, rect_colors = annotation_colors[[i]], plot_format_FUN = annotation_format_FUN[[i]])
            legend_plots[[length(legend_plots) + 1]] = p_leg
        }
        anno_plots[[length(anno_plots) + 1]] = p_anno
    }
    legend_plots = c(legend_plots, list(p_heat.leg))
    row1 = seqsetvis::assemble_heatmap_cluster_bars(
        c(anno_plots, list(p_heat)),
        rel_widths = c(rep(1, length(anno_plots)), 5))

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
    cowplot::plot_grid(row1, row2, ncol = 1)
}


group_VARS = "cluster"
sort_VAR = NULL
sort_strategy =  c("hclust", "sort", "left", "right")[2]

head(getRegionMetaData(ct2))
# debug(add_group_annotation)


.plotSignalHeatmap(ct2, group_VARS = c("cluster", "overlap"))
.plotSignalHeatmap(ct2, group_VARS = c("overlap", "cluster"))
.plotSignalHeatmap(
    ct2,
    group_VARS = c(
        "overlap",
        "peak_MCF10A_CTCF",
        "peak_MCF10AT1_CTCF",
        "peak_MCF10CA1_CTCF",
        "cluster",
        "overlap"
    ),
    color_key_strategy = "if_not_sorted"
)

.plotSignalHeatmap(
    ct2,
    group_VARS = c(
        "cluster",
        "peak_MCF10A_CTCF",
        "peak_MCF10AT1_CTCF",
        "peak_MCF10CA1_CTCF",
        "overlap",
        "cluster",
        "cluster"
    ),
    color_key_strategy = "all"
)

.plotSignalHeatmap(
    ct2,
    group_VARS = c(
        "cluster",
        "peak_MCF10A_CTCF",
        "peak_MCF10AT1_CTCF",
        "peak_MCF10CA1_CTCF",
        "overlap",
        "cluster",
        "cluster"
    ),
    color_key_strategy = "except_sort_VAR"
)

.plotSignalHeatmap(ct2, group_VARS = c("overlap",
                                       "peak_MCF10A_CTCF",
                                       "peak_MCF10AT1_CTCF",
                                       "peak_MCF10CA1_CTCF",
                                       "cluster",
                                       "overlap"), sort_VAR = FALSE, n_legend_rows = 2)

p_heat

.plotSignalLinePlot = function(ct2, group_VAR = NULL){

}
