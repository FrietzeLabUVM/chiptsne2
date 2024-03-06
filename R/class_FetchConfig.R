
#### Validity ####

.check_FetchConfig = function(object){
    errors <- character()
    #is_null can't be invalid
    if(object@is_null){
        return(TRUE)
    }
    name_VAR = object@name_VAR
    #read_mode sensitive checks.
    if(grepl("bam", object@read_mode)){#bam checks
    }else{#bigwig checks
    }
    #check attributes are present
    if(is.null(object$meta_data[[name_VAR]])){
        msg = paste("'name_VAR':", name_VAR, " attribute must be present in meta_data.")
        errors = c(errors, msg)
    }
    #check attributes are valid
    if(any(duplicated(object$meta_data[[name_VAR]]))){
        msg = paste0("'name_VAR': ", name_VAR, " values must be unique. Following have been duplicated: ", paste(unique(object$meta_data$name[duplicated(object$meta_data$name)]), collapse = ", "))
        errors = c(errors, msg)
    }
    if(!is.factor(object$meta_data[[name_VAR]])){
        msg = paste("'name_VAR':", name_VAR, " attribute must be a factor.")
        errors = c(errors, msg)
    }

    if(!all(file.exists(object$meta_data$file))){
        stop(paste(c("Files specified in config do not exist:",
                     object$meta_data$file[!file.exists(object$meta_data$file)]), collapse = "\n  "))
    }

    #check to_run and reference
    if (length(errors) == 0) TRUE else errors
}

#' @importFrom S4Vectors setValidity2
#' @importFrom BiocGenerics NCOL NROW
S4Vectors::setValidity2("FetchConfig", .check_FetchConfig)

setMethod("initialize","FetchConfig", function(.Object,...){
    .Object <- callNextMethod()
    .Object
})

setMethod("names", "FetchConfig",
          function(x)
          {
              c(
                  "view_size",
                  "window_size",
                  "read_mode",
                  "fetch_options",
                  "meta_data"
              )

          })


setMethod("$", "FetchConfig",
          function(x, name)
          {
              switch (name,
                      view_size = x@view_size,
                      window_size = x@window_size,
                      read_mode = x@read_mode,
                      fetch_options = x@fetch_options,
                      meta_data = x@meta_data
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
                                 x@window_size = value
                             },
                             read_mode = {
                                 stopifnot(value %in% c("bam_SE", "bam_PE", "bigwig", "null"))
                                 x@read_mode = value
                             },
                             fetch_options = {
                                 x@fetch_options = value
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
#' @param read_mode Read mode of signal data, one of bam_SE, bam_PE, or bigwig.
#'   Use CT_READ_MODES$.
#' @param view_size Consistent size to use when viewing assessment regions. Uses
#'   CT_VIEW_SIZE option or 3kb as default.
#' @param window_size The window size used when fetching signal. Lower values
#'   increase resolution but also RAM usage. Default is 200 bp.
#' @param fetch_options Named list of additional arguments to pass to signal
#'   fetch function.
#' @param is_null If TRUE, this FetchConfig is considered null/empty.
#'
#' @return A FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config_df = chiptsne2:::.parse_config_body(bam_config_file)
#' sig_conf = FetchConfig(bam_config_df)
#'
#' bigwig_config_file = system.file(
#'   package = "ssvQC",
#'   "extdata/ssvQC_bigwig_config.csv"
#' )
#' bigwig_config_df = chiptsne2:::.parse_config_body(bigwig_config_file)
#' sig_conf.bw = FetchConfig(bigwig_config_df)
FetchConfig = function(config_df,
                       read_mode = NULL,
                       view_size = getOption("CT_VIEW_SIZE", 3e3),
                       window_size = 200,
                       fetch_options = list(),
                       is_null = FALSE,
                       name_VAR = "name"){
    config_df = .enforce_file_var(config_df)
    config_df = .enforce_name_var(config_df, name_VAR = name_VAR)

    #Guess read mode
    if(is.null(read_mode)){
        read_mode = guess_read_mode(config_df$file[1])
    }

    stopifnot(read_mode %in% sqc_read_modes)

    .FetchConfig(
        meta_data =  config_df,
        read_mode = read_mode,
        view_size = view_size,
        window_size = window_size,
        fetch_options = fetch_options,
        name_VAR = name_VAR,
        is_null = is_null)
}

#' FetchConfig null placeholder
#'
#' @return A null/empty FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' FetchConfig.null()
FetchConfig.null = function(){
    qc = suppressWarnings({FetchConfig(data.frame(file = "null", name = "null", name_split = "null", stringsAsFactors = FALSE), is_null = TRUE)})
    qc
}

#' isFetchConfigNull
#'
#' @param cfg A FetchConfig object
#'
#' @return TRUE if object is null placeholder
#' @export
#' @rdname FetchConfig
#'
#' @examples
#' cfg.null = FetchConfig.null()
#' isFetchConfigNull(cfg.null)
isFetchConfigNull = function(fetch_config){
    fetch_config@is_null
}

#' @param signal_config_file Configuration file for signal data.
#'
#' @return A FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_config_file = system.file(package = "chiptsne2", "extdata/bam_config.csv", mustWork = TRUE)
#' FetchConfig.load_config(bam_config_file)
#'
#' bigwig_config_file = system.file(
#'   package = "chiptsne2",
#'   "extdata/bigwig_config.csv", mustWork = TRUE
#' )
#' FetchConfig.load_config(bigwig_config_file)
FetchConfig.load_config = function(signal_config_file, name_VAR = NULL){
    cfg_vals = .parse_config_header(signal_config_file)
    if(!is.null(cfg_vals$name_VAR)){
        name_VAR = cfg_vals$name_VAR
    }else{
        name_VAR = "name"
    }
    signal_config_dt = .parse_config_body(signal_config_file, name_VAR = name_VAR)
    if(any(c("main_dir", "data_dir", "file_prefix") %in% names(cfg_vals))){
        #ADD PREFIX TO FILE AND REMOVE VAR
        path_VAR = intersect(c("main_dir", "data_dir", "file_prefix"), names(cfg_vals))
        if(length(path_VAR) > 1){
            stop("only one of following allowed: ", paste(path_VAR, collapse = ", "))
        }
        path_val = cfg_vals[[path_VAR]]
        if(path_val == "$SSV_DATA"){ #special value indicating included package data
            path_val = system.file("extdata", package = "seqsetvis", mustWork = TRUE)
        }
        if(path_val == "$PACKAGE_DATA"){ #special value indicating included package data
            path_val = system.file("extdata", package = "chiptsne2", mustWork = TRUE)
        }
        signal_config_dt$file = file.path(path_val, signal_config_dt$file)
        cfg_vals[[path_VAR]] = NULL
    }

    tfun = function(config_dt,
                    read_mode = NULL,
                    view_size = getOption("CT_VIEW_SIZE", 3e3),
                    window_size = getOption("CT_WINDOW_SIZE", 200),
                    fetch_options = list(),
                    is_null = FALSE,
                    name_VAR = "name"){
        FetchConfig(config_df = config_dt,
                    read_mode = read_mode,
                    view_size = view_size,
                    window_size = window_size,
                    fetch_options = fetch_options,
                    is_null = is_null,
                    name_VAR = name_VAR
        )
    }
    do.call(tfun, c(list(config_dt = signal_config_dt), cfg_vals))
}

#
#' FetchConfig.from_files
#'
#' @param file_paths character paths to files
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#' @param view_size view size in bp to apply. Defaults to 3000.
#'
#' @return a FetchConfig object
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_files = dir(
#'   system.file(
#'     package = "ssvQC",
#'     "extdata"),
#'   pattern = "CTCF.+bam$",
#'   full.names = TRUE
#' )
#' object = FetchConfig.from_files(bam_files)
#'
#' object2 = FetchConfig.from_files(bam_files,
#'   group_names = c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1a_CTCF")
#' )
FetchConfig.from_files = function(file_paths,
                                  group_names = NULL,
                                  name_VAR = "name",
                                  view_size = getOption("CT_VIEW_SIZE", 3e3),
                                  window_size = getOption("CT_WINDOW_SIZE", 200),
                                  read_mode = NULL,
                                  fetch_options = list()
){
    if(is.null(group_names)){
        if(is.null(names(file_paths))){
            group_names = basename(file_paths)
        }else{
            group_names = names(file_paths)
        }
    }
    config_df = data.frame(file = as.character(file_paths),
                           group = group_names,
                           stringsAsFactors = FALSE)
    colnames(config_df)[2] = name_VAR

    FetchConfig(config_df,
                view_size = view_size,
                window_size = window_size,
                read_mode =  read_mode,
                name_VAR = name_VAR,
                fetch_options = fetch_options
    )
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

#' @param fetch_config A FetchConfig object
#' @param query_gr A GRanges to fetch data for
#'
#' @return A list of 2 items prof_dt and query_gr.  prof_dt is a tidy data.table
#'   of signal profiles.  query_gr is a GRanges that may have been modified from
#'   input query_gr if signal profiles are flipped or centered according to
#'   center_signal_at_max or flip_signal_mode in the signal config.
#' @rdname FetchConfig
#' @export
#' @examples
#' bam_config_file = system.file(package = "chiptsne2", "extdata/bam_config.csv")
#' fetch_config = FetchConfig.load_config(bam_config_file)
#'
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' chiptsne2:::fetch_signal_at_features(fetch_config, query_gr)
fetch_signal_at_features = function(fetch_config, query_gr, bfc = NULL){
    extra_args = fetch_config@fetch_options
    ### JRB commenting out for now. user provided fragLens should be used.
    # if(!is.null(fetch_config@meta_data$fragLens)){ # fragLens is in meta data
    #   if(!is.null(extra_args$fragLens)){ # fragLens is in extra_args
    #     if(extra_args$fragLens != "auto"){
    #       warning("Overwriting configured fragLens with detected fragLens for fetch call.")
    #     }
    #   }
    #   extra_args$fragLens = fetch_config@meta_data$fragLens
    # }
    if(is.null(extra_args$fragLens)){ # no fragLens in extra_args
        if(!is.null(fetch_config@meta_data$fragLens)){ # fragLens is in meta_data
            extra_args$fragLens = fetch_config@meta_data$fragLens # apply fragLens to extra_args from meta_data
        }
    }
    if("window_size" %in% names(extra_args)){
        warning("window_size found in configured fetch_options. ignored. Please use FetchConfig$window_size.")
        extra_args[["window_size"]] = NULL
    }
    .n_region_splits = 50
    if("n_region_splits" %in% names(extra_args)){
        .n_region_splits = extra_args[["n_region_splits"]]
        extra_args[["n_region_splits"]] = NULL
    }

    call_args = c(list(
        file_paths = fetch_config@meta_data,
        win_size = fetch_config@window_size,
        qgr = query_gr,
        return_data.table = TRUE,
        n_region_splits = .n_region_splits,
        names_variable = fetch_config@name_VAR),
        extra_args)
    fetch_FUN = get_fetch_fun(fetch_config@read_mode)
    prof_dt = bfcif(
        FUN = function(){
            do.call(fetch_FUN, call_args)
        },
        rname = digest::digest(list(fetch_FUN, call_args)),
        bfc = bfc)
    list(prof_dt = prof_dt, query_gr = query_gr)
}

#' FetchConfig.save_config
#'
#' @return Invisibly returns path to saved config file.
#' @export
#' @rdname FetchConfig
#' @examples
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = FetchConfig.load_config(bam_config_file)
#' #FetchConfig.save_config(bam_config, "bam_config.csv")
#'
#' bigwig_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bigwig_config.csv")
#' bigwig_config = FetchConfig.load_config(bigwig_config_file)
#' #FetchConfig.save_config(bigwig_config, "bigwig_config.csv")
FetchConfig.save_config = function(object, file){
    slots_to_save = c(
        "view_size",
        "window_size",
        "read_mode",
        "is_null"
    )
    # key value pair slots are saved/loaded differently
    kvp_slots = c("fetch_options")
    .save_config(object, file, slots_to_save, kvp_slots)
}

#### show ####
#### show config ####

#' .show_Config
#'
#' used by show
#'
#' @param qc A FetchConfig object
#'
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' qc = ConfigFeatures.parse(feature_config_file)
#' qc
.show_Config = function(qc){
    if(qc@is_null){
        msg = "This FetchConfig is a NULL placeholder."
    }else{
        msg = paste(sep = "\n",
                    paste("Configuration for", nrow(qc@meta_data), "items.")
        )
    }
    message(msg)
}

#' @param FetchConfig
#'
#' @export
#' @rdname FetchConfig
setMethod("show", "FetchConfig", definition = function(object).show_Config(object))
