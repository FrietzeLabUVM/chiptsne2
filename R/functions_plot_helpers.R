.prep_color_scale = function(values, color_scale){
    if(is.null(color_scale)){
        if(any(values < 0)){
            #for negative values, default to different colors and symmetrical limits
            color_scale = c("orange", "white", "purple")
        }else{
            color_scale = c("#000004FF", "#51127CFF", "#B63679FF", "#FB8861FF", "#FCFDBFFF")
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

.apply_scale = function(p, color_scale, scale_limits, fill = TRUE){
    if(is(color_scale, "Scale")){
        p = p + color_scale
    }else{
        if(fill){
            p = p + scale_fill_gradientn(colours = color_scale, limits = scale_limits)
        }else{
            p = p + scale_color_gradientn(colours = color_scale, limits = scale_limits)
        }

    }
    p
}
