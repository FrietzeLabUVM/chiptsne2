
check_FetchConfig = function(object){
    errors <- character()
    #is_null can't be invalid
    if(object@is_null){
        return(TRUE)
    }
    #read_mode sensitive checks.
    if(grepl("bam", object@read_mode)){#bam checks

    }else{#bigwig checks
        if(object@cluster_value == "RPM"){
            errors = c(errors, "cluster_value cannot be RPM for bigwig files.")
        }
        if(object@sort_value == "RPM"){
            errors = c(errors, "sort_value cannot be RPM for bigwig files.")
        }
    }
    #check attributes are present
    if(is.null(object$meta_data[["name"]])){
        msg = "'name' attribute must be present in meta_data."
        errors = c(errors, msg)
    }
    if(is.null(object$meta_data[["name_split"]])){
        msg = "'name_split' attribute must be present in meta_data."
        errors = c(errors, msg)
    }
    #check attributes are valid
    if(any(duplicated(object$meta_data[["name"]]))){
        msg = paste0("'name' values must be unique. Following have been duplicated: ", paste(unique(object$meta_data$name[duplicated(object$meta_data$name)]), collapse = ", "))
        errors = c(errors, msg)
    }
    if(any(duplicated(object$meta_data[["name_split"]]))){
        msg = paste0("'name_split' values must be unique. Following have been duplicated: ", paste(unique(object$meta_data$name_split[duplicated(object$meta_data$name_split)]), collapse = ", "))
        errors = c(errors, msg)
    }
    if(!is.factor(object$meta_data[["name"]])){
        msg = "'name' attribute must be a factor."
        errors = c(errors, msg)
    }
    if(!is.factor(object$meta_data[["name_split"]])){
        msg = "'name_split' attribute must be a factor."
        errors = c(errors, msg)
    }
    #check to_run and reference
    if(!object@center_signal_at_max %in% c(FALSE, TRUE)){
        msg = paste0("center_signal_at_max must be TRUE or FALSE, was: ", object@center_signal_at_max, ".")
        errors = c(errors, msg)
    }
    if(!object@flip_signal_mode %in% flip_signal_modes){
        msg = paste0("flip_signal_mode must be one of: ", paste(flip_signal_modes, collapse = ", "), ". But it was : ", object@flip_signal_mode, ".")
        errors = c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
}

#' FetchConfig
#'
#' @slot view_size Consistent size to use when viewing assessment regions. Uses
#'   CT_VIEW_SIZE option or 3kb as default.
#' @slot read_mode Read mode of signal data, one of bam_SE, bam_PE, or bigwig.
#'   Use CT_READ_MODES$.
#' @slot fetch_options Named list of additional arguments to pass to signal
#'   fetch function.
#' @slot cluster_value Value in CT_SIGNAL_VALUES$ to use for clustering. RPM
#'   values are not valid if read_mode is "bigwig".
#' @slot linearQuantile_cutoff Quantile to use for linearQuantile normalization
#'   procedure. Values above this cutoff are treated as outliers.
#' @slot sort_value Value in CT_SIGNAL_VALUES$ to use for sorting. RPM values
#'   are not valid if read_mode is "bigwig".
#' @slot sort_method One of two available method to use when sorting within
#'   clusters. If "hclust", hierarchical clustering is applied. If the default
#'   of "sort", regions are sorting by decreasing signal.
#' @slot plot_value Value in CT_SIGNAL_VALUES$ to represent in plots. RPM
#'   values are not valid if read_mode is "bigwig".
#' @slot heatmap_limit_values Color scale limits for heatmaps. Default is 0 to
#'   10.
#' @slot lineplot_free_limits If TRUE (default), lineplot facets per cluster
#'   will have free axis. If FALSE a consistnet y-axis is used for all clusters.
#' @slot center_signal_at_max If TRUE, signal is centered at local maxima prior
#'   to any clustering. The default is FALSE. See details for explanation or
#'   interaction with assessment features.
#' @slot flip_signal_mode  Value is CT_FLIP_SIGNAL_MODES$.  If not "none"
#'   (Default) signal profiles are flipped so that highest signal is on one side
#'   or the other.  See details for explanation or interaction with assessment
#'   features.
#' @slot n_clusters The number of k-means clusters for the heatmap.
#'
#' @rdname FetchConfig
#' @export
#'
setClass("FetchConfig", contains = "Config",
         representation = list(
             view_size = "numeric",
             win_size = "numeric",
             read_mode = "character",
             fetch_options = "list",
             cluster_value = "character",
             linearQuantile_cutoff = "numeric",
             sort_value = "character",
             sort_method = "character",
             plot_value = "character",
             heatmap_limit_values = "ANY",
             lineplot_free_limits = "logical",
             center_signal_at_max = "logical",
             flip_signal_mode = "character",
             n_clusters = "numeric"
         ), validity = check_FetchConfig)

setMethod("initialize","FetchConfig", function(.Object,...){
    .Object <- callNextMethod()
    .Object
})

#' @export
setMethod("plot", "FetchConfig", definition = function(x).plot_Config(x))

setMethod("names", "FetchConfig",
          function(x)
          {
              c(
                  "view_size",
                  "window_size",
                  "read_mode",
                  "fetch_options",
                  "cluster_value",
                  "linearQuantile_cutoff",
                  "sort_value",
                  "sort_method",
                  "plot_value",
                  "heatmap_limit_values",
                  "lineplot_free_limits",
                  "meta_data",
                  "color_by",
                  "color_mapping",
                  "center_signal_at_max",
                  "flip_signal_mode",
                  "n_clusters"
              )

          })


setMethod("$", "FetchConfig",
          function(x, name)
          {
              switch (name,
                      view_size = x@view_size,
                      window_size = x@win_size,
                      read_mode = x@read_mode,
                      fetch_options = x@fetch_options,
                      cluster_value = x@cluster_value,
                      linearQuantile_cutoff = x@linearQuantile_cutoff,
                      sort_value = x@sort_value,
                      sort_method = x@sort_method,
                      plot_value = x@plot_value,
                      heatmap_limit_values = x@heatmap_limit_values,
                      lineplot_free_limits = x@lineplot_free_limits,
                      meta_data = x@meta_data,
                      color_by = x@color_by,
                      color_mapping = as.list(x@color_mapping),
                      center_signal_at_max = x@center_signal_at_max,
                      flip_signal_mode = x@flip_signal_mode,
                      n_clusters = x@n_clusters
              )
          })

setReplaceMethod("$", "FetchConfig",
                 function(x, name, value)
                 {
                     warn_msg = "This assignment is not supported.  No effect."
                     switch (name,
                             view_size = {
                                 x@view_size = value
                             },
                             window_size = {
                                 x@win_size = value
                             },
                             read_mode = {
                                 stopifnot(value %in% c("bam_SE", "bam_PE", "bigwig", "null"))
                                 x@read_mode = value
                             },
                             fetch_options = {
                                 x@fetch_options = value
                             },
                             cluster_value = {
                                 if(!value %in% signal_vars){
                                     stop("cluster_value must be one of: ", paste(signal_vars, collapse = ", "))
                                 }
                                 x@cluster_value = value
                             },
                             linearQuantile_cutoff = {
                                 stopifnot(value > 0 & value <= 1)
                                 x@linearQuantile_cutoff = value
                             },
                             sort_value = {
                                 if(!value %in% signal_vars){
                                     stop("sort_value must be one of: ", paste(signal_vars, collapse = ", "))
                                 }
                                 x@sort_value = value
                             },
                             sort_method = {
                                 if(!value %in% c("hclust", "sort")){
                                     stop("sort_method must be one of: ", paste(c("hclust", "sort"), collapse = ", "))
                                 }
                                 x@sort_method = value
                             },
                             plot_value = {
                                 if(!value %in% signal_vars){
                                     stop("plot_value must be one of: ", paste(signal_vars, collapse = ", "))
                                 }
                                 x@plot_value = value
                             },
                             heatmap_limit_values = {
                                 if(is.null(value)){

                                 }else if(is.numeric(value)){
                                     if(length(value) != 2){
                                         stop("2 values required for limits, may use NA for one to use min/max.")
                                     }
                                 }else if(is.function(value)){

                                 }else{
                                     stop("heatmap_limit_values must be defined in a way compatible with ggplot2::scale_continuous limits.  NULL, numeric, or function allowed.")
                                 }
                                 x@heatmap_limit_values = value
                             },
                             lineplot_free_limits = {
                                 if(!value %in% c(FALSE, TRUE)){
                                     stop("lineplot_free_limits must be one of FALSE or TRUE")
                                 }
                                 x@lineplot_free_limits = value
                             },
                             meta_data = {
                                 value = .enforce_file_var(value)
                                 value = .enforce_name_var(value)
                                 value = .enforce_found_order(value)
                                 x@meta_data = value
                             },
                             color_by = {
                                 x@color_by = value
                                 message("Applying option CT_COLORS for updated color_mapping.")
                                 x$color_mapping = getOption("CT_COLORS")
                             },
                             color_mapping = {
                                 value = unlist(value)
                                 col_lev = unique(x@meta_data[[x@color_by]])
                                 if(is.null(names(value))){
                                     if(length(value) >= length(col_lev)){
                                         value = value[seq_along(col_lev)]
                                     }else{
                                         stop("Insufficient colors supplied. Mapping requires at least ", length(col_lev), ".")
                                     }
                                     names(value) = col_lev
                                 }else{
                                     if(!all(col_lev %in% names(value))){
                                         stop(paste(collapse = "\n",
                                                    c("Missing name values from color mapping. Required:",
                                                      setdiff(col_lev, names(value)))))
                                     }
                                 }
                                 x@color_mapping = value
                             },
                             center_signal_at_max = {
                                 x@center_signal_at_max = value
                             },
                             flip_signal_mode = {
                                 x@flip_signal_mode = value
                             },
                             n_clusters = {
                                 x@n_clusters = value
                             },
                             warning(warn_msg)
                     )
                     validObject(x)
                     x
                 })

#' FetchConfig
#'
#' @param config_df A data.frame containing configuration information for signal
#'   (bam or bigwig) files. Should contain a "file" attribute and entries for
#'   and color_by.
#' @param color_by Name of the attribute specify how signal data should be
#'   colored in relevant plots
#' @param color_mapping Name character vector where names are values of color_by
#'   and values are valid R colors (color names or hex values).
#' @param read_mode Read mode of signal data, one of bam_SE, bam_PE, or bigwig.
#'   Use CT_READ_MODES$.
#' @param view_size Consistent size to use when viewing assessment regions. Uses
#'   CT_VIEW_SIZE option or 3kb as default.
#' @param window_size The window size used when fetching signal. Lower values
#'   increase resolution but also RAM usage. Default is 200 bp.
#' @param fetch_options Named list of additional arguments to pass to signal
#'   fetch function.
#' @param cluster_value Value in CT_SIGNAL_VALUES$ to use for clustering. RPM
#'   values are not valid if read_mode is "bigwig".
#' @param linearQuantile_cutoff Quantile to use for linearQuantile normalization
#'   procedure. Values above this cutoff are treated as outliers.
#' @param sort_value Value in CT_SIGNAL_VALUES$ to use for sorting. RPM values
#'   are not valid if read_mode is "bigwig".
#' @param sort_method One of two available method to use when sorting within
#'   clusters. If "hclust", hierarchical clustering is applied. If the default
#'   of "sort", regions are sorting by decreasing signal.
#' @param plot_value Value in CT_SIGNAL_VALUES$ to represent in plots. RPM
#'   values are not valid if read_mode is "bigwig".
#' @param heatmap_limit_values Color scale limits for heatmaps. Default is 0 to
#'   10.
#' @param lineplot_free_limits If TRUE (default), lineplot facets per cluster
#'   will have free axis. If FALSE a consistnet y-axis is used for all clusters.
#' @param center_signal_at_max If TRUE, signal is centered at local maxima prior
#'   to any clustering. The default is FALSE. See details for explanation or
#'   interaction with assessment features.
#' @param flip_signal_mode  Value is CT_FLIP_SIGNAL_MODES$.  If not "none"
#'   (Default) signal profiles are flipped so that highest signal is on one side
#'   or the other.  See details for explanation or interaction with assessment
#'   features.
#' @param is_null If TRUE, this FetchConfig is considered null/empty.
#' @param n_clusters The number of k-means clusters for the heatmap.
#'
#'   Since center_signal_at_max and flip_signal_mode have the potential to
#'   modify the assessment regions based on signal run groups, the modified
#'   assessment feature set is store in signal_data for each run group. They can
#'   be accessed like this: sqc$signal_data$FEATURE_NAME$SIGNAL_NAME$query_gr
#'
#' @return A FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config_df = .parse_config_body(bam_config_file)
#' sig_conf = FetchConfig(bam_config_df)
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config_df = .parse_config_body(bigwig_config_file)
#' sig_conf.bw = FetchConfig(bigwig_config_df)
FetchConfig = function(config_df,
                          color_by = "file",
                          color_mapping = NULL,
                          read_mode = NULL,
                          view_size = getOption("CT_VIEW_SIZE", 3e3),
                          window_size = 200,
                          fetch_options = list(),
                          cluster_value = NULL,
                          linearQuantile_cutoff = .98,
                          sort_value = NULL,
                          sort_method = c("hclust", "sort")[2],
                          plot_value = NULL,
                          heatmap_limit_values = c(0, 10),
                          lineplot_free_limits = TRUE,
                          center_signal_at_max = FALSE,
                          flip_signal_mode = flip_signal_modes$none,
                          is_null = FALSE,
                          n_clusters = 6){
    config_df = .enforce_file_var(config_df)
    config_df = .enforce_name_var(config_df)
    if(!color_by %in% colnames(config_df)){
        stop("color_by ", color_by, " was not in column names.")
    }
    if(linearQuantile_cutoff <= 0 | linearQuantile_cutoff > 1){
        stop("linearQuantile_cutoff must be between 0 and 1. Was ", linearQuantile_cutoff)
    }
    if(!sort_method %in% c("hclust", "sort")){
        stop("sort_method of ", sort_method, " was not one of : ", paste(c("hclust", "sort"), collapse = ", "))
    }

    if(!is.null(color_mapping)){
        if(!is.null(names(color_mapping))){
            color_names = names(color_mapping)
        }else if(is.factor(config_df[[color_by]])){
            color_names = levels(config_df[[color_by]])
        }else{
            color_names = unique(config_df[[color_by]])
        }
        stopifnot(length(color_names) == length(color_mapping))
        names(color_mapping) = color_names

    }else{
        if(is.factor(config_df[[color_by]])){
            color_names = levels(config_df[[color_by]])
        }else{
            color_names = unique(config_df[[color_by]])
        }
        color_mapping = seqsetvis::safeBrew(length(color_names))
        names(color_mapping) = color_names
    }

    stopifnot(config_df[[color_by]] %in% names(color_mapping))

    #Guess read mode
    if(is.null(read_mode)){
        read_mode = guess_read_mode(config_df$file[1])
    }
    #Value defaults
    if(is.null(cluster_value)){
        cluster_value = get_default_signal_var(read_mode)
    }
    if(is.null(sort_value)){
        sort_value = get_default_signal_var(read_mode)
    }
    if(is.null(plot_value)){
        plot_value = get_default_signal_var(read_mode)
    }
    #Value checks
    if(!cluster_value %in% signal_vars){
        stop("cluster_value of ", cluster_value, " was not one of : ", paste(signal_vars, collapse = ", "))
    }
    if(!sort_value %in% signal_vars){
        stop("sort_value of ", sort_value, " was not one of : ", paste(signal_vars, collapse = ", "))
    }
    if(!plot_value %in% signal_vars){
        stop("plot_value of ", plot_value, " was not one of : ", paste(signal_vars, collapse = ", "))
    }

    stopifnot(read_mode %in% sqc_read_modes)

    new("FetchConfig",
        meta_data =  config_df,
        color_by = color_by,
        color_mapping = color_mapping,
        read_mode = read_mode,
        view_size = view_size,
        win_size = window_size,
        fetch_options = fetch_options,
        cluster_value = cluster_value,
        linearQuantile_cutoff = linearQuantile_cutoff,
        sort_value = sort_value,
        sort_method = sort_method,
        plot_value = plot_value,
        heatmap_limit_values = heatmap_limit_values,
        lineplot_free_limits = lineplot_free_limits,
        center_signal_at_max = center_signal_at_max,
        flip_signal_mode = flip_signal_mode,
        is_null = is_null,
        n_clusters = n_clusters)
}

#' FetchConfig null placeholder
#'
#' @return A null/empty FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' FetchConfig.null()
FetchConfig.null = function(){
    qc = suppressWarnings({FetchConfig(data.frame(file = "null", stringsAsFactors = FALSE), is_null = TRUE)})
    qc
}

#' @param signal_config_file Configuration file for signal data.
#'
#' @return A FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' FetchConfig.parse(bam_config_file)
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' FetchConfig.parse(bigwig_config_file)
FetchConfig.parse = function(signal_config_file){
    signal_config_dt = .parse_config_body(signal_config_file)
    valid_signal_var = c("main_dir",
                         "view_size",
                         "window_size",
                         "read_mode",
                         "color_by",
                         "color_mapping",
                         "fetch_options",
                         "is_null",
                         "cluster_value",
                         "linearQuantile_cutoff",
                         "sort_value",
                         "plot_value",
                         "sort_method",
                         "heatmap_limit_values",
                         "lineplot_free_limits",
                         "center_signal_at_max",
                         "flip_signal_mode",
                         "n_clusters"

    )
    cfg_vals = .parse_config_header(signal_config_file, valid_signal_var)

    if(!is.null(cfg_vals[["main_dir"]])){
        choose_file_path = function(main_dir, files){
            abs_path = files
            rel_path = file.path(main_dir, files)
            ifelse(file.exists(rel_path), rel_path, abs_path)
        }
        signal_config_dt$file = choose_file_path(cfg_vals[["main_dir"]], signal_config_dt$file)
        cfg_vals[["main_dir"]] = NULL
    }
    if(!all(file.exists(signal_config_dt$file))){
        stop(paste(c("Files specified in config do not exist:",
                     signal_config_dt$file[!file.exists(signal_config_dt$file)]), collapse = "\n  "))
    }

    tfun = function(config_dt,
                    main_dir = NULL,
                    read_mode = NULL,
                    view_size = getOption("CT_VIEW_SIZE", 3e3),
                    window_size = getOption("CT_WINDOW_SIZE", 200),
                    color_by = NULL, color_mapping = NULL,
                    cluster_value = NULL,
                    linearQuantile_cutoff = .98,
                    sort_value = NULL,
                    sort_method = c("hclust", "sort")[2],
                    plot_value = NULL,
                    fetch_options = list(),
                    heatmap_limit_values = c(0, 10),
                    lineplot_free_limits = TRUE,
                    center_signal_at_max = FALSE,
                    flip_signal_mode = flip_signal_modes$none,
                    is_null = FALSE,
                    n_clusters = 6){
        FetchConfig(config_df = config_dt,
                       color_by = color_by,
                       color_mapping = color_mapping,
                       read_mode = read_mode,
                       view_size = view_size,
                       window_size = window_size,
                       fetch_options = fetch_options,
                       cluster_value = cluster_value,
                       linearQuantile_cutoff = linearQuantile_cutoff,
                       sort_value = sort_value,
                       sort_method = sort_method,
                       plot_value = plot_value,
                       heatmap_limit_values = heatmap_limit_values,
                       lineplot_free_limits = lineplot_free_limits,
                       center_signal_at_max = center_signal_at_max,
                       flip_signal_mode = flip_signal_mode,
                       is_null = is_null,
                       n_clusters = n_clusters)
    }
    do.call(tfun, c(list(config_dt = signal_config_dt), cfg_vals))
}

#
#' FetchConfig.files
#'
#' @param file_paths character paths to files
#' @param groups numeric vector of group assignments. 1 is first item in group_names, 2 is second, etc. Default is seq_along(file_path)
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#' @param view_size view size in bp to apply. Defaults to 3000.
#'
#' @return a FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_files = dir(system.file(package = "ssvQC", "extdata"), pattern = "CTCF.+bam$", full.names = TRUE)
#' object = FetchConfig.files(bam_files)
#' plot(object)
#'
#' object2 = FetchConfig.files(bam_files,
#'   sample_names = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1a_CTCF"),
#'   group_names = c("10A", "AT1", "CA1"),
#'   group_colors = c("firebrick", "slategray2", "forestgreen")
#' )
#' plot(object2)
FetchConfig.files = function(file_paths,
                                file_paths.input = character(),
                                sample_names = NULL,
                                sample_names.split = NULL,
                                group_names = NULL,
                                group_name.input = "input",
                                group_colors = NULL,
                                view_size = getOption("CT_VIEW_SIZE", 3e3),
                                window_size = getOption("CT_WINDOW_SIZE", 200),
                                read_mode = NULL,
                                cluster_value = NULL,
                                linearQuantile_cutoff = .98,
                                sort_value = NULL,
                                plot_value = NULL,
                                sort_method = c("hclust", "sort")[2],
                                center_signal_at_max = FALSE,
                                flip_signal_mode = flip_signal_modes$none,
                                n_clusters = 6
){
    if(is.null(group_names)){
        group_names = paste(seq_along(file_paths), basename(file_paths))
    }
    if(is.null(group_colors)){
        group_colors = get_group_colors(group_names)
    }
    if(is.null(names(group_colors))){
        names(group_colors) = group_names
    }

    config_df = data.frame(file = c(as.character(file_paths), file_paths.input),
                           group = c(group_names, rep(group_name.input, length(file_paths.input))),
                           All = "All",
                           stringsAsFactors = FALSE)

    if(is.null(sample_names)){
        config_df$name = basename(config_df$file)
    }else{
        config_df$name = sample_names
    }
    if(is.null(sample_names.split)){
        config_df$name_split = gsub("[_\\. ]", "\n", config_df$name)
    }else{
        config_df$name_split = sample_names.split
    }

    FetchConfig(config_df,
                   color_by = "group",
                   color_mapping = group_colors,
                   view_size = view_size,
                   window_size = window_size,
                   cluster_value = cluster_value,
                   linearQuantile_cutoff = linearQuantile_cutoff,
                   sort_value = sort_value,
                   sort_method = sort_method,
                   plot_value = plot_value,
                   center_signal_at_max = center_signal_at_max,
                   flip_signal_mode = flip_signal_mode,
                   n_clusters = n_clusters)
}

get_fetch_fun = function(read_mode){
    stopifnot(read_mode %in% c("bam_SE", "bam_PE", "bigwig"))
    switch(read_mode,
           bam_SE = {
               seqsetvis::ssvFetchBam
           },
           bam_PE = {
               seqsetvis::ssvFetchBamPE
           },
           bigwig = {
               seqsetvis::ssvFetchBigwig
           })
}

#' @param qc_signal A FetchConfig object
#' @param query_gr A GRanges to fetch data for
#'
#' @return A list of 2 items prof_dt and query_gr.  prof_dt is a tidy data.table
#'   of signal profiles.  query_gr is a GRanges that may have been modified from
#'   input query_gr if signal profiles are flipped or centered according to
#'   center_signal_at_max or flip_signal_mode in the signal config.
#' @rdname FetchConfig
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' qc_signal = FetchConfig.parse(bam_config_file)
#'
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' fetch_signal_at_features(qc_signal, query_gr)
fetch_signal_at_features = function(qc_signal, query_gr, bfc = new_cache()){
    extra_args = qc_signal@fetch_options
    ### JRB commenting out for now. user provided fragLens should be used.
    # if(!is.null(qc_signal@meta_data$fragLens)){ # fragLens is in meta data
    #   if(!is.null(extra_args$fragLens)){ # fragLens is in extra_args
    #     if(extra_args$fragLens != "auto"){
    #       warning("Overwriting configured fragLens with detected fragLens for fetch call.")
    #     }
    #   }
    #   extra_args$fragLens = qc_signal@meta_data$fragLens
    # }
    if(is.null(extra_args$fragLens)){ # no fragLens in extra_args
        if(!is.null(qc_signal@meta_data$fragLens)){ # fragLens is in meta_data
            extra_args$fragLens = qc_signal@meta_data$fragLens # apply fragLens to extra_args from meta_data
        }
    }
    if("win_size" %in% names(extra_args)){
        warning("win_size found in configured fetch_options. ignored. Please use FetchConfig$window_size.")
        extra_args[["win_size"]] = NULL
    }
    call_args = c(list(
        file_paths = qc_signal@meta_data,
        win_size = qc_signal@win_size,
        qgr = query_gr,
        return_data.table = TRUE,
        names_variable = "name"),
        extra_args)
    fetch_FUN = get_fetch_fun(qc_signal@read_mode)
    prof_dt = bfcif(bfc, digest(list(fetch_FUN, call_args)), function(){
        do.call(fetch_FUN, call_args)
    })
    #### apply center_signal_at_max ####
    if(qc_signal@center_signal_at_max == TRUE){
        query_gr.center = centerGRangesAtMax(prof_dt = prof_dt, qgr = query_gr, width = width(query_gr))
        message("Centering signal profiles...")
        call_args.center = call_args
        call_args.center$qgr = query_gr.center
        prof_dt = bfcif(bfc, digest(list(fetch_FUN, call_args.center)), function(){
            do.call(fetch_FUN, call_args.center)
        })
        query_gr = query_gr.center
    }
    #### apply flip_signal_mode ####
    if(qc_signal@flip_signal_mode != flip_signal_modes$none){
        #need to bring query_gr with flip info out
        balance_dt = prof_dt[, list(right_sum = sum(y[x > 0]),
                                    left_sum = sum(y[x < 0])),
                             by = list(name, id)]
        if(qc_signal@flip_signal_mode == flip_signal_modes$high_on_right){
            balance_dt = balance_dt[, list(needs_flip = left_sum > right_sum,
                                           name,
                                           id)]
        }else{
            balance_dt = balance_dt[, list(needs_flip = left_sum < right_sum,
                                           name,
                                           id)]
        }

        most_flipped = balance_dt[,
                                  list(fraction_flipped = sum(needs_flip) / .N),
                                  by = list(id)]
        most_flipped[, flip_strand := fraction_flipped > .5]
        GenomicRanges::strand(query_gr) = "+"
        GenomicRanges::strand(query_gr)[most_flipped$flip_strand] = "-"
        prof_dt = merge(prof_dt, balance_dt, by = c("id", "name"))
        remove(balance_dt)
        prof_dt[needs_flip == TRUE, x := -x]
        prof_dt$needs_flip = NULL
    }

    list(prof_dt = prof_dt, query_gr = query_gr)
}

setMethod("split", signature = c("FetchConfig", "factor"), definition = function(x, f){
    split(x, f, FALSE)
})
setMethod("split", signature = c("FetchConfig", "factor", "logical"), definition = function(x, f, drop = FALSE){
    config_df = x@meta_data
    meta_split = split(config_df, f)
    refs = config_df[config_df[[x@run_by]] %in% x@to_run_reference,]

    meta_split = meta_split[x@to_run]
    meta_split = lapply(meta_split, function(sel_meta_df)rbind(sel_meta_df, refs))

    lapply(meta_split, function(sel_meta_df){
        new("FetchConfig",
            meta_data =  sel_meta_df,
            color_by = x@color_by,
            color_mapping = x@color_mapping,
            read_mode = x@read_mode,
            view_size = x@view_size,
            win_size = x@win_size,
            cluster_value = x@cluster_value,
            linearQuantile_cutoff = x@linearQuantile_cutoff,
            sort_value = x@sort_value,
            sort_method = x@sort_method,
            plot_value = x@plot_value,
            fetch_options = x@fetch_options,
            heatmap_limit_values = x@heatmap_limit_values,
            lineplot_free_limits = x@lineplot_free_limits,
            center_signal_at_max = x@center_signal_at_max,
            flip_signal_mode = x@flip_signal_mode,
            is_null = FALSE,
            n_clusters = x@n_clusters
        )
    })
})
setMethod("split", signature = c("FetchConfig", "character"), definition = function(x, f){
    split(x, f, FALSE)
})
setMethod("split", signature = c("FetchConfig", "character", "logical"), definition = function(x, f, drop = FALSE){
    split(x, factor(f, levels = unique(f)), FALSE)
})

#' FetchConfig.save_config
#'
#' @return Invisibly returns path to saved config file.
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = FetchConfig.parse(bam_config_file)
#' #FetchConfig.save_config(bam_config, "bam_config.csv")
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config = FetchConfig.parse(bigwig_config_file)
#' #FetchConfig.save_config(bigwig_config, "bigwig_config.csv")
FetchConfig.save_config = function(object, file){
    slots_to_save = c(
        "view_size",
        "win_size",
        "read_mode",
        "color_by",
        "is_null",
        "cluster_value",
        "linearQuantile_cutoff",
        "sort_value",
        "sort_method",
        "plot_value",
        "heatmap_limit_values",
        "lineplot_free_limits",
        "center_signal_at_max",
        "flip_signal_mode",
        "n_clusters"

    )
    # key value pair slots are saved/loaded differently
    kvp_slots = c("color_mapping", "fetch_options")
    # FetchConfig.parse(file)
    .save_config(object, file, slots_to_save, kvp_slots)
}

