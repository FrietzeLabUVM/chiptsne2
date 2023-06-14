library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()
ct2@position_VAR = "bp"
ct2@value_VAR = "reads"
ct2@region_VAR = "peak_regions"
group_VAR = "cluster"
n_to_plot = 500
heatmap_theme = theme(panel.background = element_blank(), panel.grid = element_blank())
hetamap_theme.no_y = heatmap_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
plot_theme = heatmap_theme
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

add_group_annotation = function(group_ids,
                                xleft = 0,
                                xright = 1,
                                rect_colors = RColorBrewer::brewer.pal(8, "Dark2"),
                                text_colors = rep("black", length(rect_colors)),
                                label_angle = 0,
                                row_ = "id",
                                cluster_ = "group_id",
                                show_legend = FALSE){
    if(data.table::is.data.table(group_ids)){
        assign_dt = unique(group_ids[, c(row_, cluster_), with = FALSE])[order(get(row_))]
        group_ids = assign_dt[[cluster_]]

    }else if(is.data.frame(group_ids)){
        assign_dt = group_ids %>% dplyr::select(dplyr::all_of(c(row_, cluster_))) %>% dplyr::arrange(get(row_))
        group_ids = assign_dt[[cluster_]]
    }
    if(is.character(group_ids)) group_ids = factor(group_ids, levels = unique(group_ids))
    if(is.numeric(group_ids)) group_ids = factor(group_ids)

    stopifnot(is.factor(group_ids))
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

    p = ggplot(df_rects) +
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(group_ids)), expand = FALSE) +
        facet_grid(.~grp)

    p = p + geom_rect(data = df_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = grp)) +
        scale_fill_manual(values = rect_colors)
    if(!show_legend){
        p = p + guides(fill = "none")
    }
    p
}
add_group_annotation.legend = function(group_ids,
                                       rect_colors = RColorBrewer::brewer.pal(8, "Dark2"),
                                       row_ = "id",
                                       cluster_ = "group_id"){
    xleft = 0
    xright = 1
    text_colors = rep("black", length(rect_colors))
    label_angle = 0
    p = add_group_annotation(
        group_ids,
        xleft,
        xright,
        rect_colors,
        text_colors,
        label_angle,
        row_,
        cluster_,
        show_legend = TRUE
    )
    cowplot::get_legend(p)
}

add_group_annotation(group_ids)
debug(add_group_annotation)
p = add_group_annotation(anno_df, cluster_ = "overlap", row_ = ct2@region_VAR)
leg = add_group_annotation.legend(anno_df, cluster_ = "overlap", row_ = ct2@region_VAR)
p + cowplot::draw_grob(leg, x = 0, y = 50)

p_cluster_bar2 = seqsetvis::add_cluster_annotation(anno_df$cluster)
p_cluster_bar2 = seqsetvis::add_cluster_annotation(anno_df$overlap)

assign_dt = anno_df
rect_colors = seqsetvis::safeBrew(anno_df$cluster)
text_colors = "black"
show_labels = TRUE
label_angle = 0

colData(ct2)



add_cluster_annotation = function(cluster_ids,
                                  xleft = 0, xright = 1,
                                  rect_colors = c("black", "gray"),
                                  text_colors = rev(rect_colors),
                                  show_labels = TRUE,
                                  label_angle = 0,
                                  row_ = "id",
                                  cluster_ = "cluster_id"){
    if(data.table::is.data.table(cluster_ids)){
        assign_dt = unique(cluster_ids[, c(row_, cluster_), with = FALSE])[order(get(row_))]
        cluster_ids = assign_dt[[cluster_]]

    }
    if(is.character(cluster_ids)) cluster_ids = factor(cluster_ids, levels = unique(cluster_ids))
    if(is.numeric(cluster_ids)) cluster_ids = factor(cluster_ids)
    stopifnot(is.factor(cluster_ids))



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
        df_rects$fill = rect_colors[rownames(df_rects)]
    }else{
        df_rects$fill = rect_colors[seq_len(nrow(df_rects))%%length(rect_colors)+1]
    }
    if(setequal(names(text_colors), rownames(df_rects))){
        df_rects$color = text_colors[rownames(df_rects)]
    }else{
        df_rects$color = text_colors[seq_len(nrow(df_rects))%%length(text_colors)+1]
    }


    df_rects = df_rects[rev(seq_len(nrow(df_rects))),]
    cluster_labels = levels(cluster_ids)

    p = ggplot(df_rects) +
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(cluster_ids)), expand = FALSE) +
        facet_grid(.~facet)
    p = p + geom_rect(aes(
        xmin = xmin,
        xmax = xmax,
        ymin= ymin,
        ymax = ymax,
        fill = fill
    )) +
        scale_fill_identity() +
        facet_grid(.~grp)
    if(show_labels){
        for(i in seq_len(nrow(df_rects))){

            p = p + annotate("text",
                             x = mean(c(df_rects$xmin[i], df_rects$xmax[i])),
                             y = mean(c(df_rects$ymin[i], df_rects$ymax[i])),
                             label = cluster_labels[i],
                             color = df_rects$color[i],
                             angle = label_angle)
        }
    }
    p
}

p_cluster_bar = add_cluster_annotation(
    cluster_ids = anno_df[[group_VAR]],
    # p = p_cluster_bar,
    xleft = 0,
    xright = 1,
    rect_colors = rect_colors,
    text_colors = text_colors,
    show_labels = show_labels,
    label_angle = label_angle)

overlap_cols = c("gray", "MCF10A_CTCF" = "red", "MCF10AT1_CTCF" = "blue", "MCF10A_CTCF & MCF10AT1_CTCF" = "purple")
p_overlap = add_group_annotation(anno_df, cluster_ = "overlap", row_ = ct2@region_VAR, rect_colors = overlap_cols)
p_overlap
p_heat

seqsetvis::assemble_heatmap_cluster_bars(list(p_overlap + labs(subtitle = "overlap"), p_cluster_bar, p_heat), rel_widths = c(1, 1, 7))
seqsetvis::ssvSignalHeatmap.ClusterBars()

group_dt = prof_dt[order(get(ct2@region_VAR)), c(ct2@region_VAR, group_VAR), with = FALSE] %>% unique
group_dt = group_dt[, .N, c(group_VAR)]
group_dt[, ]
group_tab =table(group_dt[[group_VAR]])

.plotSignalHeatmap = function(ct2, group_VARS = NULL, sort_VAR = NULL){
    if(!is.null(group_VARS) & is.null(sort_VAR)){
        sort_VAR = group_VARS[1]
    }
    prof_dt = .getTidyProfile(ct2, group_VAR)
    if(is.null(group))

}

.plotSignalLinePlot = function(ct2, group_VAR = NULL){

}
