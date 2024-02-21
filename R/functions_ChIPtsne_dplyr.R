#longterm should use:
#https://dplyr.tidyverse.org/reference/dplyr_extending.html

#' @export
subsetRegions = function(ct2, expr){
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, eval(ps))
}
#' @export
subsetSamples = function(ct2, expr){
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, TRUE, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, TRUE, eval(ps))
}

#' @export
subsetRow = function(ct2, expr){
    .Deprecated("subsetRegions")
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, eval(ps))
}
#' @export
subsetCol = function(ct2, expr){
    .Deprecated("subsetSamples")
    #https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
    # ssubset = deparse(substitute(expr))
    # subset(ct2, TRUE, eval(parse(text = ssubset)))
    ps <- substitute(expr)
    subset(ct2, TRUE, eval(ps))
}

#' @export
filter.ChIPtsne2 = function(.data, ...){
    subset(.data, ...)
}

