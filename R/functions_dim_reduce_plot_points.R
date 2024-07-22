
.add_labels = function(p, xy_df, label_VAR, label_FUN, label_size, map_label_colors){
    #visible binding for global variable
    group = value = group_value = tx = ty = NULL
    if(label_VAR == TMP_group_VAR){
        # label_VAR = c("group_value")
        xy_df = dplyr::mutate(xy_df, group_value = paste(group, value))
        label_VAR = TMP_value_VAR
        label_ = label_VAR
        label_ = ensym(label_)
        lab_df = xy_df %>% dplyr::group_by(group_value) %>% dplyr::summarise(tx = mean(tx), ty = mean(ty), value = unique(value), group = unique(group))
    }else{
        label_ = label_VAR
        label_ = ensym(label_)
        lab_df = xy_df %>% dplyr::group_by(!!label_) %>% dplyr::summarise(tx = mean(tx), ty = mean(ty))
    }
    if(map_label_colors){
        p = p + label_FUN(data = lab_df,
                          mapping = aes(label = !!label_, color = !!label_),
                          show.legend = FALSE,
                          size = label_size / ggplot2::.pt)
    }else{
        p = p + label_FUN(data = lab_df,
                          mapping = aes(label = !!label_),
                          show.legend = FALSE,
                          size = label_size / ggplot2::.pt)
    }

    p
}

#' .enforce_extra_VARS
#'
#' @param ct2 chiptsne2 object
#' @param df current data.frame to add extra_VARS to
#' @param extra_VARS extra variables that must be in metadata or in df already
#' @param expected_missing it's ok if these are in extra_VARS and not present in
#'   df or meta data
#'
#' @return df with extra_VARS added from meta data.
#'
.enforce_extra_VARS = function(ct2, df, extra_VARS, expected_missing = NULL){
    extra_VARS_missed = setdiff(extra_VARS, colnames(df))
    all_meta_cn = c(
        colnames(getSampleMetaData(ct2)),
        colnames(getRegionMetaData(ct2))
    )
    .validate_allowed_input(
        input = extra_VARS_missed,
        allowed = union(all_meta_cn, expected_missing),
        msg_prefix = "Some VAR are missing from metadata:")

    if(length(extra_VARS_missed) > 0){
        if(ct2@region_VAR %in% colnames(df) & ct2@name_VAR %in% colnames(df)){
            full_region_cn = colnames(getRegionMetaData(ct2))
            extra_region_cn = intersect(extra_VARS, full_region_cn)
            if(!is.null(extra_region_cn)){
                region_df = getRegionMetaData(ct2, select_VARS = extra_region_cn)
                df = merge(df, region_df, by = ct2@region_VAR)
            }
            full_sample_cn = colnames(getSampleMetaData(ct2))
            extra_sample_cn = intersect(extra_VARS, full_sample_cn)
            if(!is.null(extra_sample_cn)){
                sample_df = getSampleMetaData(ct2, select_VARS = extra_sample_cn)
                df = merge(df, sample_df, by = ct2@name_VAR)
            }
        }else if(ct2@region_VAR %in% colnames(df)){
            full_region_cn = colnames(getRegionMetaData(ct2))
            extra_region_cn = intersect(extra_VARS, full_region_cn)
            if(!is.null(extra_region_cn)){
                region_df = getRegionMetaData(ct2, select_VARS = extra_region_cn)
                df = merge(df, region_df, by = ct2@region_VAR)
            }
        }else if(ct2@name_VAR %in% colnames(df)){
            full_sample_cn = colnames(getSampleMetaData(ct2))
            extra_sample_cn = intersect(extra_VARS, full_sample_cn)
            if(!is.null(extra_sample_cn)){
                sample_df = getSampleMetaData(ct2, select_VARS = extra_sample_cn)
                df = merge(df, sample_df, by = ct2@name_VAR)
            }
        }else{
            stop("confusing")
        }
    }
    df
}
.background_FUN = function(p, xy_df, point_size, background_annotation_color){
    if(!is.null(background_annotation_color)){
        bg_df = unique(xy_df[, c("tx", "ty")])
        p = p + annotate("point", x = bg_df$tx, y = bg_df$ty, color = background_annotation_color, size = .7*point_size)
    }
    p
}

TMP_group_VAR = "TMP___group"
TMP_value_VAR = "TMP___value"

.plotDimReducePoints = function(ct2,
                                color_VAR = NULL,
                                label_VAR = NULL,
                                label_FUN = geom_label,
                                label_size = 10,
                                point_size = NULL,
                                point_color_limits = c(NA, NA),
                                has_symmetrical_limits = NULL,
                                point_colors = NULL,
                                extra_VARS = NULL,
                                background_annotation_color = NULL,
                                underlayer_FUN = function(p, ...)p,
                                return_data = FALSE){
    #visible binding NOTE
    tx = ty = value = NULL
    if(!hasDimReduce(ct2)){
        stop("No dimensional reduction data present in this ChIPtsne2 object. Run dimReduceTSNE/PCA/UMAP first then try again.")
    }
    xy_df = getRegionMetaData(ct2, select_VARS = c("tx", "ty"))
    if(is.null(point_size)){
        nr = nrow(xy_df)
        point_size = 1/nr*100
        if(point_size < .05) point_size = .05
        if(point_size > 1) point_size = 1
    }
    if(is.null(color_VAR)){
        color_VAR = colnames(ct2)
    }
    extra_VARS = union(extra_VARS, color_VAR)
    if(!is.null(label_VAR)){
        if(label_VAR == TRUE){
            label_VAR = TMP_value_VAR
        }
        if(label_VAR == FALSE){
            label_VAR = NULL
        }
        extra_VARS = union(extra_VARS, label_VAR)
    }
    map_label_colors = FALSE
    if(!is.null(label_VAR)){
        if(label_VAR %in% color_VAR | label_VAR == TMP_value_VAR){
            map_label_colors = TRUE
        }
    }

    if(all(is.na(color_VAR))){
        # no color
        xy_df = .enforce_extra_VARS(ct2, xy_df, extra_VARS, expected_missing = NA)
        if(return_data){
            return(xy_df)
        }
        p = ggplot(xy_df, aes(x = tx, y = ty))
        p = underlayer_FUN(p, xy_df, point_size, background_annotation_color)
        p = .background_FUN(p, xy_df, point_size, background_annotation_color)
        p = p +
            geom_point(size = point_size)

    }else if(all(color_VAR %in% colnames(getRegionMetaData(ct2)))){
        # color with region variable
        xy_df = getRegionMetaData(ct2) %>%
            dplyr::select(dplyr::all_of(c("tx", "ty", ct2@region_VAR, color_VAR)))
        #color_VAR must (error) all be same class and should (warning) have items in common to make sense
        xy_df[, color_VAR, drop = FALSE]
        if(length(color_VAR) > 1){
            color_classes = sapply(color_VAR, function(cv){
                class(xy_df[[cv]])
            })
            color_is_num = sapply(color_VAR, function(cv){
                is.numeric(xy_df[[cv]])
            })
            if(length(unique(color_classes)) != 1){
                msg = .message_list(split(names(color_classes), color_classes))
                stop("Classes of all color_VAR items must match.\n", msg)
            }
            if(!all(color_is_num)){
                color_values = lapply(color_VAR, function(cv){
                    as.character(unique(xy_df[[cv]]))
                })
                for(i in seq_len(length(color_values) - 1)){
                    for(j in seq(i + 1, length(color_values))){
                        in_common = intersect(color_values[[i]], color_values[[j]])
                        if(length(in_common) == 0){
                            warning("There are entries in color_VAR with no items in common. The same scale may be inappropriate.:\n",
                                    names(color_values)[i], " and ", names(color_values)[j])
                        }
                    }
                }
            }
        }
        xy_df = tidyr::pivot_longer(xy_df, setdiff(colnames(xy_df), c(ct2@region_VAR, "tx", "ty")), names_to = TMP_group_VAR, values_to = TMP_value_VAR)
        xy_df = .enforce_extra_VARS(ct2, xy_df, extra_VARS)
        if(return_data){
            return(xy_df)
        }

        point_colors = .prep_color_scale(xy_df[[TMP_value_VAR]], color_scale = point_colors)
        p = ggplot(xy_df, aes(x = tx, y = ty))
        p = .apply_scale(p, point_colors, point_color_limits, fill = FALSE)
        p = underlayer_FUN(p, xy_df, point_size, background_annotation_color)
        p = .background_FUN(p, xy_df, point_size, background_annotation_color)
        p = p +
            geom_point(aes(color = !!ensym(TMP_value_VAR)), size = point_size) +
            labs(color = NULL) +
            facet_wrap(paste0("~", TMP_group_VAR))
    }else if(all(color_VAR %in% colnames(ct2))){
        # color by max signal
        signal_df = SummarizedExperiment::assay(ct2, "max") %>%
            as.data.frame
        signal_df = signal_df[, color_VAR, drop = FALSE]
        signal_df[[ct2@region_VAR]] = rownames(signal_df)
        xy_df = merge(xy_df, signal_df, by = ct2@region_VAR)
        xy_df = tidyr::pivot_longer(xy_df, setdiff(colnames(xy_df), c(ct2@region_VAR, "tx", "ty")), names_to = ct2@name_VAR, values_to = "max")
        xy_df[[ct2@name_VAR]] = factor(xy_df[[ct2@name_VAR]], levels = colnames(ct2))
        point_colors = .prep_color_scale(xy_df$max, has_symmetrical_limits, point_colors)
        point_color_limits = .prep_symmetrical(xy_df$max, has_symmetrical_limits, point_color_limits)
        xy_df$max = .apply_limits(xy_df$max, point_color_limits)
        xy_df = .enforce_extra_VARS(ct2, xy_df, extra_VARS, expected_missing = colnames(ct2))
        if(return_data){
            return(xy_df)
        }
        p = ggplot(xy_df, aes(x = tx, y = ty))
        p = underlayer_FUN(p, xy_df, point_size, background_annotation_color)
        p = .background_FUN(p, xy_df, point_size, background_annotation_color)
        p = p +
            geom_point(aes(color = max), size = point_size) +
            facet_wrap(paste0("~", ct2@name_VAR)) +
            labs(color = paste("max", ct2@value_VAR, "\nper", ct2@region_VAR))
        p = .apply_scale(p, point_colors, point_color_limits, fill = FALSE)
    }else{
        stop("color_VAR: \"", color_VAR, "\" was not recognized. Check vs colnames of ct2 object or colnames of rowData(ct2).")
    }
    if(!is.null(label_VAR)){
        p = .add_labels(
            p = p,
            xy_df = xy_df,
            label_VAR = label_VAR,
            label_FUN = label_FUN,
            label_size = label_size,
            map_label_colors
        )
    }
    p
}

generic_plotDimReducePoints = function(ct2,
                                       color_VAR = NULL,
                                       label_VAR = NULL,
                                       label_FUN = geom_label,
                                       label_size = 10,
                                       point_size = NULL,
                                       point_color_limits = c(NA, NA),
                                       has_symmetrical_limits = NULL,
                                       point_colors = NULL,
                                       extra_VARS = NULL,
                                       background_annotation_color = NULL,
                                       underlayer_FUN = function(p, ...)p,
                                       return_data = FALSE){
    standardGeneric("plotDimReducePoints")
}


#' plotDimReducePoints
#'
#' @param ct2 valid ChIPtsne2 after dimReduce has been run
#' @param color_VAR Control color assignment in plot. Can match entries in
#'   *either* sample metadata (colnames) or region metadata (rowRanges). Default
#'   of NULL will plot max signal for all sample profiles. NA will perform no
#'   color mapping.
#' @param label_VAR Categorical variable to label the mean position of. Use with `label_FUN` to control geom used.
#' @param label_FUN Function to add labels to plot. Only used when `label_VAR` is specified. Should be equivalent to geom_label: geom_text, ggrepel::geom_text_repel, or ggrepel::geom_label_repel. Must accept parameters, data, mapping, and show.legend.
#' @param label_size Font size of label.
#' @param point_size Size of points in plot.
#' @param point_color_limits color scale limits for continuous color_VAR.
#' @param has_symmetrical_limits If TRUE color scale limits will extend to equal
#'   magnitude in positive and negative direction. Default is TRUE when negative
#'   values are present and FALSE otherwise.
#' @param point_colors Either a vectors of colors to pass to
#'   [ggplot2::scale_colour_gradientn] for continuous data, or a named vector of
#'   colors for categorical data.
#' @param extra_VARS `r doc_extra_VARS()`
#' @param background_annotation_color Color to use for points not in facet. Default of NULL will not draw any background points.
#' @param underlayer_FUN Function to add to layer below background annotation.
#' @param return_data `r doc_return_data()`
#'
#' @return ggplot
#' @export
#' @rdname plotDimReducePoints
#'
#' @examples
#' library(ggplot2)
#' ct2 = exampleChIPtsne2.with_meta() %>%
#'    dimReduceUMAP() %>%
#'    groupRegionsByDimReduceCluster(group_VAR = "umap_cluster", nearest_neighbors = 20) %>%
#'    groupRegionsBySignalCluster(group_VAR = "signal_cluster")
#'
#' # default is max signal value per sample
#' plotDimReducePoints(ct2)
#' plotDimReducePoints(ct2, point_colors = c("gray", "red"))
#'
#' #NA disable color
#' plotDimReducePoints(ct2, color_VAR = NA)
#'
#' # color scale is different when negative values are present
#' ct2_diff = subsetSamples(ct2, cell == "MCF10A") -
#'   subsetSamples(ct2, cell == "MCF10AT1")
#' plotDimReducePoints(ct2_diff)
#'
#' # categorical variables
#' plotDimReducePoints(ct2, "umap_cluster")
#' plotDimReducePoints(ct2, "umap_cluster", label_VAR = "umap_cluster")
#' #a named vector of colors for point_colors
#' plotDimReducePoints(ct2, "umap_cluster",
#'   point_colors = seqsetvis::safeBrew(as.character(1:4), "paired"))
#' plotDimReducePoints(ct2, c("umap_cluster", "signal_cluster"))
#' plotDimReducePoints(ct2, c("MCF10A_CTCF", "MCF10AT1_CTCF"))
#'
#' # layer plot elements beneath the final plot with a function like this:
#' base_plot = function(p, xy_df, point_size, background_annotation_color){
#'   p +
#'     annotate("rect",
#'              xmin = 0, xmax = .3,
#'              ymin = -.05, ymax = .13,
#'              fill = "lightblue", color = "red")
#' }
#'
#' plotDimReducePoints(ct2, extra_VARS = "peak_MCF10CA1_CTCF",
#'   background_annotation_color = "gray50",
#'   underlayer_FUN = base_plot) +
#'     facet_grid(peak_MCF10CA1_CTCF~sample) +
#'     labs(title = "signal per sample facetted by peak_MCF10CA1_CTCF")
setGeneric("plotDimReducePoints",
           generic_plotDimReducePoints,
           signature = "ct2")

#' @export
#' @rdname plotDimReducePoints
setMethod("plotDimReducePoints", c("ChIPtsne2_no_rowRanges"), .plotDimReducePoints)
