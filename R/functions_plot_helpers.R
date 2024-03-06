.prep_color_scale = function(values, has_symmetrical_limits = NULL, color_scale = NULL){
    if(is.numeric(values)){
        if(is.null(color_scale)){
            apply_symm = any(values < 0)
            if(!is.null(has_symmetrical_limits)){
                apply_symm = has_symmetrical_limits
            }
            if(apply_symm){
                #for negative values, default to different colors and symmetrical limits
                color_scale = c("orange", "gray80", "purple")
            }else{
                # this does not show up well on light color background
                # color_scale = c("#000004FF", "#51127CFF", "#B63679FF", "#FB8861FF", "#FCFDBFFF")
                color_scale = c("gray80", "yellow", "orange", "red")
            }
        }
    }else{
        if(!all(unique(values) %in% names(color_scale))){
            auto_point_colors = seqsetvis::safeBrew(setdiff(values, names(color_scale)))
            color_scale = c(auto_point_colors, color_scale)
        }
    }
    color_scale
}

.prep_symmetrical = function(values, has_symmetrical_limits, scale_limits){
    if(any(values < 0)){
        #for negative values, default to different colors and symmetrical limits
        if(is.null(has_symmetrical_limits)){
            has_symmetrical_limits = TRUE
        }
    }else{
        if(is.null(has_symmetrical_limits)){
            has_symmetrical_limits = FALSE
        }
    }
    if(has_symmetrical_limits){
        if(is.na(scale_limits[1])){
            scale_limits[1] = min(values)
        }
        if(is.na(scale_limits[2])){
            scale_limits[2] = max(values)
        }
        lim_max = max(abs(scale_limits))
        scale_limits = c(-lim_max, lim_max)
    }
    scale_limits
}

.apply_limits = function(values, scale_limits){
    values[values > max(scale_limits)] = max(scale_limits)
    values[values < min(scale_limits)] = min(scale_limits)
    values
}

#' @import ggplot2
.apply_scale = function(p, color_scale, scale_limits, fill = TRUE){
    if(is(color_scale, "Scale")){
        p = p + color_scale
    }else{
        if(fill){
            p = p + ggplot2::scale_fill_gradientn(colours = color_scale, limits = scale_limits)
        }else{
            p = p + ggplot2::scale_color_gradientn(colours = color_scale, limits = scale_limits)
        }

    }
    p
}
