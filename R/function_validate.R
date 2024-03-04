.validate_allowed_input = function(input, allowed, msg_prefix = "Invalid values detected:"){
    if(!all(input %in% allowed)){
        missed = setdiff(input, allowed)
        missed = paste0('"', missed, '"')
        stop(paste(c(msg_prefix, missed), collapse = "\n"))
    }
}

.validate_names_match = function(args, dim_FUN, str){
    ref = args[[1]]
    for(test in args[-1]){
        is_match = dim_FUN(ref) == dim_FUN(test)
        if(!all(is_match)){
            a = dim_FUN(ref)[!is_match]
            b = dim_FUN(test)[!is_match]
            stop(paste(c(paste0(str, " names must be identical for all ChIPtsne2_no_rowRanges objects. Example mismatches: "), head(paste(a, b, sep = " != "))), collapse = "\n"))
        }
    }
}

.validate_names_unique = function(args, dim_FUN, str){
    cns = unname(unlist(lapply(args, dim_FUN)))
    cn_dupes = duplicated(cns)
    if(any(cn_dupes)){
        stop(paste0("Duplicated ", str, " names are not allowed when combining ChIPtsne2_no_rowRanges objects. You may need to use setNameVariable to differentiate names between ChIPtsne2_no_rowRanges objects. Duplicated examples:\n"),
             paste(head(unique(cns[cn_dupes])), collapse = "\n"))
    }
}
