COLOR_KEY_STRAT = list(
    if_not_sorted = "if_not_sorted",
    except_sort_VAR = "except_sort_VAR",
    all = "all"
)

.heatmap_theme = theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.direction = "horizontal"
)
.heatmap_theme.no_y = .heatmap_theme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

.annotation_theme = theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
)

.prep_names = function(name){
    if(is.factor(name)){
        levels(name) = gsub("_", "\n", levels(name))
    }else{
        name = gsub("_", "\n", name)
    }
    name
}

.prep_ids = function(ids, row_, cluster_, force_factor = TRUE){
    if(is.data.frame(ids)){
        assign_dt = unique(ids[, c(row_, cluster_), drop = FALSE])
        assign_dt = assign_dt[order(assign_dt[[row_]]),]
        ids = assign_dt[[cluster_]]

    }
    if(force_factor){
        if(!is.factor(ids)){
            ids = factor(ids, levels = unique(ids))
        }
        stopifnot(is.factor(ids))
    }
    ids
}


add_group_annotation = function(anno_ids,
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
    if(is.data.frame(anno_ids)){#yes this could be replaced with S4 method dispatch
        if(is.numeric(anno_ids[[cluster_]])){
            argg <- as.list(environment())
            return(do.call(add_cluster_annotation.numeric, argg))
        }
    }
    if(is.null(rect_colors)){
        rect_colors = RColorBrewer::brewer.pal(8, "Dark2")

    }
    if(is.null(text_colors)){
        text_colors = rep("black", length(rect_colors))
    }
    anno_ids = .prep_ids(anno_ids, row_, cluster_)
    n_grps = length(unique(anno_ids))
    if(n_grps > length(rect_colors)){
        stop("Not enough rect_colors provided for number of unique anno_ids. Default supports up to 8 groups.")
    }
    if(is.null(names(rect_colors))){
        rect_colors = rect_colors[seq(n_grps)]
        names(rect_colors) = unique(anno_ids)
    }
    if(is.null(names(text_colors))){
        text_colors = text_colors[seq(n_grps)]
        names(text_colors) = unique(anno_ids)
    }
    stopifnot(all(anno_ids %in% names(rect_colors)))
    stopifnot(all(anno_ids %in% names(text_colors)))

    anno_ids = rev(anno_ids)

    ends = c(which(!(anno_ids[-1] == anno_ids[-length(anno_ids)])), length(anno_ids))
    # ends = cumsum(rev(table(anno_ids)))
    starts = c(1, ends[-length(ends)] + 1)
    grps = anno_ids[starts]
    starts = starts-.5
    ends = ends+.5

    df_rects = data.frame(xmin = xleft,
                          xmax = xright,
                          ymin = starts,
                          ymax = ends,
                          grp = grps)
    if(all(unique(df_rects$grp) %in% names(rect_colors))){
        df_rects$fill = rect_colors[df_rects$grp]
    }else{
        stop("invalid rect_colors, check names match group ids.")
    }
    if(all(unique(df_rects$grp) %in% names(text_colors))){
        df_rects$color = text_colors[df_rects$grp]
    }else{
        stop("invalid text_colors, check names match group ids")
    }

    df_rects[[cluster_]] = name_FUN(cluster_)
    p = ggplot(df_rects) +
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(anno_ids))+.5, expand = FALSE) +
        facet_grid(paste0(".~", cluster_)) +
        labs(fill = cluster_)

    p = p + geom_rect(data = df_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = grp)) +
        scale_fill_manual(values = rect_colors)
    if(!show_legend){
        p = p + guides(fill = "none")
    }
    p + annotation_theme
}
add_group_annotation.legend = function(anno_ids,
                                       rect_colors = NULL,
                                       row_ = "id",
                                       cluster_ = "group_id",
                                       plot_format_FUN = NULL,
                                       annotation_theme = .annotation_theme){
    # if(is.null(rect_colors)){
    #     rect_colors = RColorBrewer::brewer.pal(8, "Dark2")
    # }
    xleft = 0
    xright = 1
    text_colors = rep("black", length(rect_colors))
    label_angle = 0
    p = add_group_annotation(
        anno_ids = anno_ids,
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

add_cluster_annotation.numeric = function(anno_ids,
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
        rect_colors = c("#132B43",
                        "#56B1F7")
    }
    anno_ids = .prep_ids(anno_ids, row_, cluster_, force_factor = FALSE)

    anno_rle = rle(anno_ids)
    ends = cumsum(rev(anno_rle$lengths))
    starts = c(1, ends[-length(ends)] + 1)
    starts = starts - .5
    names(starts) = names(ends)
    ends = ends + .5

    df_rects = data.frame(xmin = xleft,
                          xmax = xright,
                          ymin = starts,
                          ymax = ends,
                          grp = cluster_)
    df_rects = df_rects[rev(seq_len(nrow(df_rects))),]
    df_rects[[cluster_]] = name_FUN(cluster_)
    df_rects[["grp"]] = anno_rle$values
    p = ggplot(df_rects) +
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(anno_ids))+.5, expand = FALSE) +
        facet_grid(.~grp)
    p = p + geom_rect(aes(
        xmin = xmin,
        xmax = xmax,
        ymin= ymin,
        ymax = ymax,
        fill = grp
    )) +
        scale_fill_gradientn(colours = rect_colors) +
        facet_grid(paste0(".~", cluster_)) +
        labs(fill = cluster_)
    if(!show_legend){
        p = p + guides(fill = "none")
    }
    p + annotation_theme

}

add_cluster_annotation = function(anno_ids,
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
    if(is.data.frame(anno_ids)){#yes this could be replaced with S4 method dispatch
        if(is.numeric(anno_ids[[cluster_]])){
            argg <- as.list(environment())
            return(do.call(add_cluster_annotation.numeric, argg))
        }
    }
    if(is.null(rect_colors)){
        rect_colors = c("black", "gray")
    }
    if(is.null(text_colors)){
        text_colors = rev(rect_colors)
    }
    anno_ids = .prep_ids(anno_ids, row_, cluster_)

    ends = cumsum(rev(table(anno_ids)))
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
    cluster_labels = levels(anno_ids)
    df_rects[[cluster_]] = name_FUN(cluster_)
    df_rects[["grp"]] = rownames(df_rects)

    p = ggplot(df_rects) +
        coord_cartesian(xlim = c(xleft, xright), ylim = c(0, length(anno_ids))+.5, expand = FALSE) +
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

add_cluster_annotation.legend = function(anno_ids,
                                         rect_colors = NULL,
                                         row_ = "id",
                                         cluster_ = "cluster_id",
                                         plot_format_FUN = NULL,
                                         annotation_theme = .annotation_theme){
    # if(is.null(rect_colors)){
    #     rect_colors = c("black", "gray")
    # }
    xleft = 0
    xright = 1
    text_colors = rep("black", length(rect_colors))
    label_angle = 0
    p = add_cluster_annotation(
        anno_ids = anno_ids,
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

