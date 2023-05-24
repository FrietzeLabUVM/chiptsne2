
#' Config
#'
#' @slot file_paths character.
#' @slot groups numeric.
#' @slot group_names character.
#' @slot group_colors character.
#' @rdname Config
#' @export
setClass("Config",
         representation = list(
             meta_data = "data.frame",
             color_by = "character",
             color_mapping = "character",
             is_null = "logical"
         ))

setMethod("initialize","Config", function(.Object,...){
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})

#' Config
#'
#' @param file_paths character paths to files
#' @param groups numeric vector of group assignments. 1 is first item in group_names, 2 is second, etc. Default is seq_along(file_path)
#' @param group_names vector of group names to assign from according to groups
#' @param group_colors vector of colors to use per group
#'
#' @return A Config object
#' @export
#' @importFrom methods new
#' @rdname Config
#' @examples
#' Config.files(c("A", "B"))
Config.files = function(file_paths,
                          groups = NULL,
                          group_names = NULL,
                          group_colors = NULL){
    if(is.null(groups)){
        groups = seq_along(file_paths)
    }
    if(is.null(group_names)){
        group_names = LETTERS[seq_along(unique(groups))]
    }
    if(is.null(group_colors)){
        group_colors = seqsetvis::safeBrew(length(group_names))
    }
    if(is.null(names(group_colors))){
        names(group_colors) = group_names
    }
    config_df = data.frame(file = as.character(file_paths),  color = group_names[groups], group = "All", stringsAsFactors = FALSE)
    new("Config",
        meta_data =  config_df,
        color_by = "color",
        color_mapping = group_colors
    )
}

#### color_mapping ####

#' color_mapping
#'
#' @param object a Config object
#'
#' @return a named vector of colors
#' @export
#' @rdname Config
setGeneric("color_mapping", function(object){standardGeneric("color_mapping")})

#' color_mapping for Config
#'
#' @export
#' @rdname Config
#' @aliases color_mapping,Config-method
#' @examples
#' color_mapping(Config(c("A", "B")))
setMethod("color_mapping", c("Config"), function(object){
    cols = object@group_colors
    names(cols) = as.character(object@group_names)
    cols
})

#### scale_color_config ####

#' scale_color_config
#'
#' @param object a Config object
#'
#' @return a ggplot2 scale_fill_manual
#' @export
#' @rdname scale_color_config-methods
setGeneric("scale_color_config", function(object){standardGeneric("scale_color_config")})

#' scale_color_config for Config
#'
#' @export
#' @rdname Config
#' @aliases scale_color_config,Config-method
#' @examples
#' my_df = data.frame(group = c("A", "B"))
#' my_df$x = 1:2
#' my_df$y = 1
#' library(ggplot2)
#' ggplot(my_df, aes(x = x, y = y, color = group)) +
#'   geom_point(size = 20) +
#'   expand_limits(x = c(0, 3)) +
#'   scale_color_config(Config(my_df$group))
setMethod("scale_color_config", c("Config"), function(object){
    cols = color_mapping(object)
    ggplot2::scale_color_manual(values = cols)
})

#### scale_fill_config ####

#' scale_fill_config
#'
#' @param object a Config object
#'
#' @return a ggplot2 scale_fill_manual
#' @export
#' @rdname Config
setGeneric("scale_fill_config", function(object){standardGeneric("scale_fill_config")})

#' scale_fill_config for Config
#'
#' @export
#' @aliases scale_fill_config,Config-method
#' @rdname Config
#' @examples
#' my_df = data.frame(group = c("A", "B"))
#' my_df$x = 1:2
#' my_df$y = 1
#' library(ggplot2)
#' ggplot(my_df, aes(x = x, y = y, fill = group)) +
#'   geom_point(size = 20, pch = 22) +
#'   expand_limits(x = c(0, 3)) +
#'   scale_fill_config(Config(my_df$group))
setMethod("scale_fill_config", c("Config"), function(object){
    cols = object@group_colors
    names(cols) = as.character(object@group_names)
    ggplot2::scale_fill_manual(values = cols)
})

#### show config ####

#' .show_Config
#'
#' used by show
#'
#' @param qc A Config object
#'
#' @examples
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' qc = ConfigFeatures.parse(feature_config_file)
#' qc
.show_Config = function(qc){
    if(qc@is_null){
        msg = "This Config is a NULL placeholder."
    }else{
        msg = paste(sep = "\n",
                    paste("Configuration for", nrow(qc@meta_data), "items."),
                    paste0("Use plot() to view color mapping for '", qc@color_by , "'.")
        )
    }
    message(msg)
}

.plot_Config = function(qc){
    if(qc@is_null){
        msg = "This Config is a NULL placeholder."
        ggplot() + labs(title = "This Config is a NULL placeholder.")
    }else{
        plot_dt = as.data.table(qc@meta_data)
        browser()
        ggplot(plot_dt, aes_string(x = qc@run_by, y = "y", fill = qc@color_by, label = "name_split")) +
            geom_label() +
            scale_fill_manual(values = qc@color_mapping) +
            theme(panel.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            labs(y = "") +
            guides(
                fill = guide_legend(
                    override.aes = aes(label = "")
                ))
    }

}

#' @export
setMethod("plot", "Config", definition = function(x).plot_Config(x))

#' @param Config
#'
#' @export
#' @rdname Config
setMethod("show", "Config", definition = function(object).show_Config(object))

#### save config ####

#' Config.save_config
#'
#' @param object Config object to save to text file.
#'
#' @return invisibly returns path to saved config file.
#' @export
#' @rdname Config
#' @examples
#' cfg_file = system.file("extdata/ssvQC_peak_config.csv", package = "ssvQC")
#' qc_features = ConfigFeatures.parse(cfg_file)
#' f = tempfile()
#' Config.save_config(qc_features, f)
Config.save_config = function(object, file){
    slots_to_save = c(
        "color_by",
        "is_null"
    )
    kvp_slots = c("color_mapping")
    # ConfigSignal.parse(file)
    .save_config(object, file, slots_to_save, kvp_slots)
}
