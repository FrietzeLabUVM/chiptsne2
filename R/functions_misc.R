
#' guess_read_mode
#'
#' @param signal_file A bam of bigwig file path.
#'
#' @return A character describing format of signal_file.
#' @export
#'
#' @examples
#' bam_file = dir(system.file("extdata", package = "seqsetvis"), pattern = "bam$", full.names = TRUE)[1]
#' bw_file = dir(system.file("extdata", package = "seqsetvis"), pattern = "bw$", full.names = TRUE)[1]
#' guess_read_mode(bam_file)
#' guess_read_mode(bw_file)
guess_read_mode = function(signal_file){
    if(signal_file[1] == "null"){
        return("null")
    }
    if(grepl(".bam$", signal_file[1])){
        mode = "bam_SE"
    }else{
        mode = "bigwig"
    }
    message("read_mode has been guessed as ", mode)
    if(mode == "bam_SE"){
        message("Currently chiptsne2 cannot guess whether a bam file is SE or PE.  Please manually specify bam_PE if appropriate.")
    }
    mode
}

#' sampleCap
#'
#' A safer version of sample with additional code so that if a sample size of n is greater than length of x, n is reduced and a reordered x is returned.
#'
#' @inheritParams base::sample
#'
#' @return A vector of length size or original length of x if size is too larger, with elements drawn randomly from x.
#' @export
#'
#' @examples
#' x = LETTERS
#' #this would cause an error is normal sample
#' sampleCap(x, 50)
#' #this is equivalent to sample
#' sampleCap(x, 5)
sampleCap = function(x, size = 500){
    size = min(size, length(unique(x)))
    out = sample(unique(x), size)
    # if(is.factor(out)) out = as.character(out) #not sure why this would be necessary
    out
}

#' get_mapped_reads
#'
#' @param bam_file A bam file.  Matching .bai file must exist.
#'
#' @return The number of mapped reads in bam file.
#' @export
#'
#' @examples
#' data_dir = system.file("extdata", package = "seqsetvis", mustWork = TRUE)
#' bam_file = dir(data_dir, pattern = "bam$", full.names = TRUE)[1]
#'
#' get_mapped_reads(bam_file)
get_mapped_reads = function(bam_file){
    stopifnot(file.exists(bam_file))
    stopifnot(file.exists(paste0(bam_file, ".bai")))
    stats = Rsamtools::idxstatsBam(bam_file)
    sum(stats[,3])
}

.enforce_file_var = function(my_df){
    if(!"file" %in% colnames(my_df)){
        file_var = colnames(my_df)[1]
        message("Guessing file paths are in first column, ", file_var)
        colnames(my_df)[colnames(my_df) == file_var] = "file"
    }
    my_df
}

.enforce_name_var = function(my_df, name_VAR = "name"){
    if(!name_VAR %in% colnames(my_df)){
        my_df[[name_VAR]] = basename(my_df$file)
        message("Using \"file\" for sample names. Include ", name_VAR, " in config to specify more meaningful names.")
        # stop("Config must specify \"", name_VAR, "\" for use in chiptsne2")
    }
    if(any(duplicated(my_df[[name_VAR]]))){
        stop("Duplicate entries in 'name', all values must be unique.")
    }
    if(!is.factor(my_df[[name_VAR]]))
        my_df[[name_VAR]] = factor(my_df[[name_VAR]], levels = my_df[[name_VAR]])
    my_df
}

.enforce_found_order = function(my_df){
    for(var in colnames(my_df)){
        if(is.character(my_df[[var]])){
            if(var != "file")
                my_df[[var]] = factor(my_df[[var]], levels = unique(my_df[[var]]))
        }
    }
    my_df
}

#' @importFrom utils read.table
.parse_config_body = function(f, name_VAR = "name"){
    config_dt = utils::read.table(f, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    config_dt = .enforce_file_var(config_dt)
    config_dt = .enforce_name_var(config_dt, name_VAR = name_VAR)
    config_dt = .enforce_found_order(config_dt)
    #move file to first column
    config_dt = config_dt[, colnames(config_dt)[order(colnames(config_dt) != "file")]]
    config_dt
}


#' .parse_fetch_options
#'
#' @param fop character of options
#'
#' @return named list of options
#'
#' @examples
#' fop = "win_size:10,win_method:\"summary\",summary_FUN:mean"
#' fop = strsplit(fop, ",")[[1]]
#' chiptsne2:::.parse_fetch_options(fop)
.parse_fetch_options = function(fop){
    if(is.null(fop)) return(list())
    if(any(is.na(fop))) return(list())
    # fop = strsplit(fop, ",")[[1]]
    fop = strsplit(fop, ":")
    fop_names = sapply(fop, function(x)x[1])
    fop_opts = sapply(fop, function(x)x[2])
    if(any(is.na(fop_names))){
        stop("problem parsing fetch_options, verify fetch_options=key1:val1,key2:val2,... syntax")
    }
    if(any(is.na(fop_opts))){
        stop("problem parsing fetch_options, verify fetch_options=key1:val1,key2:val2,... syntax")
    }
    if("summary_FUN" %in% fop_names){
        message("summary_FUN is not yet fully supported.  Definitions of summary_FUN should work but its value will not be saved via FetchConfig.save_config.")
    }
    names(fop_opts) = fop_names
    lapply(fop_opts, function(x){
        #check if number
        if(suppressWarnings({!is.na(as.numeric(x))})){
            x = as.numeric(x)
        }else if(grepl('"', x)){
            x = gsub('"', "", x)
        }else{
            x = get(x)
        }
        x
    })
}

#' .parse_config_header
#'
#' @param f a config file
#' @param valid_config_VARS supported variables to attempt to extract
#'
#' @return A named list containing configuration options mapped to values.
#' @importFrom dplyr filter
#' @examples
#'
#' bam_config_file = system.file(package = "chiptsne2", "extdata/bam_config.csv", mustWork = TRUE)
#' chiptsne2:::.parse_config_header(bam_config_file)
.parse_config_header = function(f,
                                valid_config_VARS = c(
                                    "main_dir",
                                    "data_dir",
                                    "file_prefix",
                                    "view_size",
                                    "window_size",
                                    "read_mode",
                                    "fetch_options",
                                    "is_null",
                                    "name_VAR"

                                ), allowed_deprecated_VARS = c(
                                    "color_by",
                                    "color_mapping",
                                    "run_by",
                                    "to_run",
                                    "to_run_reference"
                                )
){
    cfg_txt = read.table(f, sep = "\n", header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "") %>%
        dplyr::filter(grepl("^#CFG", V1))
    if(nrow(cfg_txt) == 0){
        return(list())
    }
    cfg_txt = paste(sub("#CFG ?", "", cfg_txt$V1), collapse = " ")
    sp = strsplit(cfg_txt, " +")[[1]]
    sp = strsplit(sp, "=")
    sp
    cfg_names = sapply(sp, function(x){
        x[1]
    })
    cfg_vals = sapply(sp, function(x){
        x[2]
    })
    names(cfg_vals) = cfg_names
    cfg_vals = strsplit(cfg_vals, ",")

    bad_var = setdiff(cfg_names, valid_config_VARS)
    bad_var = setdiff(bad_var, allowed_deprecated_VARS)
    if(length(bad_var) > 0){
        stop("Unrecogized variables in config file: ", paste(bad_var, collapse = ", "))
    }

    ### special cases
    #numeric: consensus_n, consensus_fraction, n_peaks, view_size
    for(var in c("consensus_n", "consensus_fraction", "n_peaks", "view_size", "overlap_extension", "heatmap_limit_values", "linearQuantile_cutoff", "window_size")){
        if(!is.null(cfg_vals[[var]])){
            cfg_vals[[var]] = as.numeric(cfg_vals[[var]])
            if(!is.numeric(cfg_vals[[var]])){
                stop("The variable, '", var, "' is only supported as a numeric currently.")
            }
        }
    }
    #logical: is_null
    for(var in c("is_null", "lineplot_free_limits", "balance_groups")){
        if(!is.null(cfg_vals[[var]])){
            cfg_vals[[var]] = as.logical(cfg_vals[[var]])
        }
    }
    #read mode synonym
    if(!is.null(cfg_vals[["read_mode"]])){
        read_mode = cfg_vals[["read_mode"]]
        if(read_mode == "SE") read_mode = "bam_SE"
        if(read_mode == "PE") read_mode = "bam_PE"
        if(!read_mode %in% c("bam_SE", "bam_PE", "bigwig")){
            stop("read_mode '", read_mode, "' is not recognized as a valid choice from 'bam_SE', 'bam_PE', or 'bigwig'")
        }
        cfg_vals[["read_mode"]] = read_mode
    }
    #fetch_options
    if(!is.null(cfg_vals[["fetch_options"]])){
        cfg_vals[["fetch_options"]] = .parse_fetch_options(cfg_vals[["fetch_options"]])
    }
    cfg_vals = cfg_vals[setdiff(names(cfg_vals), allowed_deprecated_VARS)]
    k_char = sapply(cfg_vals, is.character)
    #strip quotes
    cfg_vals[k_char] = gsub('"', "", cfg_vals[k_char])
    cfg_vals
}

.test_suff = function(files, suff){
    sapply(files, function(f){
        any(sapply(suff, function(s){
            grepl(paste0(s, "$"), f)
        }))
    })
}


is_signal_file = function(files, suff = getOption("SQC_SIGNAL_FILE_SUFF", c("bam", "bigwig", "bw", "bigWig", "BigWig"))){
    .test_suff(files, suff)
}

#' internal function used by FetchConfig.save_config FetchConfigSignal.save_config and FetchConfigFeatures.save_config
#' @importFrom data.table fwrite
.save_config = function(object, file, slots_to_save, kvp_slots, toss_names = "summary_FUN"){
    hdr1 = sapply(slots_to_save, function(x){
        val = slot(object, x)
        ifelse(length(val) > 0,
               paste0("#CFG ", x, "=", paste(val, collapse = ",")),
               character())
    })
    hdr1 = hdr1[!is.na(hdr1)]

    hdr2 = sapply(kvp_slots, function(x){
        val = slot(object, x)
        if(any(toss_names %in% names(val))){
            warning("Could not save all FetchConfig slots: ", paste(intersect(toss_names, names(val)), collapse = ", "))
            val = val[!names(val) %in% toss_names]
        }
        #characters must be protected with quotes
        val = lapply(val, function(x){
            if(is.character(x)){
                x = paste0('"', x, '"')
            }else{
                x
            }
        })

        val = paste(names(val), val, sep = ":", collapse = ",")
        ifelse(length(val) > 0,
               paste0("#CFG ", x, "=", val),
               character())
    })
    hdr2 = hdr2[!is.na(hdr2)]

    hdr3 = paste(colnames(object@meta_data), collapse = ",")


    hdr = c(hdr1, hdr2, hdr3)
    names(hdr) = NULL
    writeLines(hdr, file)
    data.table::fwrite(object@meta_data, file, append = TRUE)
    invisible(file)
}

get_group_colors = function(group_names){
    cols = getOption("SQC_COLORS", seqsetvis::safeBrew(group_names, "Dark2"))
    cols
}

#' ssv_mclapply
#'
#' @return result of either pblapply or pbmclapply
#'
#' @importFrom pbapply pblapply
#' @importFrom pbmcapply pbmclapply
ssv_mclapply = function(X, FUN, mc.cores = getOption("mc.cores", 1), ...){
    if(.Platform$OS.type == "windows" || mc.cores == 1) {
        pbapply::pblapply(X = X, FUN = FUN, ...)

    } else {
        pbmcapply::pbmclapply(X = X, FUN = FUN, mc.cores = mc.cores, ...)
    }
}


#' get_args
#'
#' returns parameters of calling function as a named list.
#'
#' @param env Environment to collect variable from. By default this is the
#'   calling environment, i.e. the function containing a call to get_args. This
#'   may also be a function that returns an environment.
#' @param to_ignore character names of variables to ignore.
#' @param ... Additional variables that should be considered as part of
#'   environment.
#'
#' @return parameters of calling function as a named list.
#' @export
#'
#' @examples
#' #The most common usage is to simply collect all local variables in a function
#' test_fun = function(x = 1, y = 2){
#'   get_args()
#' }
#' test_fun()
#'
#' #Specified variables may be ignored
#' test_fun2 = function(x = 1, y = 2){
#'   get_args(to_ignore = "x")
#' }
#' test_fun2()
#'
#' #Additional variables can also be added from higher environments
#' global_z = 3
#' test_fun3 = function(x = 1, y = 2){
#'   get_args(env = parent.frame, to_ignore = character(), z = global_z)
#' }
#' test_fun3()
get_args = function(env = parent.frame(), to_ignore = "ct2", ...){
    if(is.function(env)) env = env()
    stopifnot(is.environment(env))
    args = c(as.list(env), list(...))
    args = args[!names(args) %in% to_ignore]
    args[order(names(args))]
}

.add_region_metadata = function(query_gr, region_metadata, region_VAR, overwrite = FALSE){
    if(is.null(region_metadata[[region_VAR]])){
        region_metadata[[region_VAR]] = rownames(region_metadata)
    }
    stopifnot(setequal(
        names(query_gr),
        region_metadata[[region_VAR]]
    ))
    region_metadata[[region_VAR]] = factor(region_metadata[[region_VAR]], levels = names(query_gr))
    region_metadata = region_metadata %>% dplyr::arrange(get(region_VAR))
    conflicting_cn = intersect(colnames(region_metadata), colnames(mcols(query_gr)))
    old_mcols = GenomicRanges::mcols(query_gr)
    if(length(conflicting_cn) > 0){
     if(overwrite){
         old_mcols = old_mcols[, setdiff(colnames(old_mcols), conflicting_cn)]
     }else{
         stop(paste(c("Conflicting colnames in region_metadata already present in query_gr:", conflicting_cn), collapse = "\n"))
     }
    }
    new_mcols = cbind(
        old_mcols,
        as.data.frame(region_metadata %>% dplyr::select(!dplyr::all_of(c(region_VAR))))
    )
    GenomicRanges::mcols(query_gr) = NULL
    GenomicRanges::mcols(query_gr) = new_mcols
    query_gr
}
