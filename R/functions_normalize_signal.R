.form_merge_df = function(ct2, optional_data, optional_data_name, data_VAR, data_VAR_name){
    if(is.null(optional_data)){
        if(is.null(colData(ct2)[[data_VAR]])){
            stop(optional_data_name, " was not supplied and ", data_VAR_name, ": ", data_VAR, " is not in colData")
        }
        optional_data = getSampleMetaData(ct2)[, c(ct2@name_VAR, data_VAR)]
    }
    #handle numeric input and reform to data.frame
    if(is.numeric(optional_data)){
        if(is.null(names(optional_data))){
            stop("When ", optional_data_name, " is supplied, names must be set.")
        }else if(!setequal(rownames(colData(ct2)), names(optional_data))){
            stop("names of ", optional_data_name, " is not setequal to name_VAR: ", ct2@name_VAR, " in colData.")
        }
        optional_data = data.frame(names(optional_data), optional_data)
        colnames(optional_data) = c(ct2@name_VAR, data_VAR)
        optional_data
    }
    if(is.null(optional_data[[data_VAR]])){
        stop(optional_data_name, " does not contain ", data_VAR_name, ": ", data_VAR, "")
    }
    optional_data
}

#### RPM ####
.normalizeSignalRPM = function(ct2, mapped_reads_data = NULL, mapped_reads_VAR = "mapped_reads"){
    args = get_args()
    mapped_reads_data = .form_merge_df(ct2, mapped_reads_data, "mapped_reads_data", mapped_reads_VAR, "mapped_reads_VAR")
    #merge mapped_reads_data to tidy profile
    new_prof_dt = getTidyProfile(ct2)
    new_prof_dt = merge(new_prof_dt, mapped_reads_data, by = ct2@name_VAR)
    new_prof_dt = new_prof_dt %>%
        dplyr::mutate(NEW_SIGNAL_ = get(ct2@value_VAR) / get(mapped_reads_VAR) * 1e6)
    new_prof_dt[[ct2@value_VAR]] = NULL
    new_prof_dt[[mapped_reads_VAR]] = NULL
    data.table::setnames(new_prof_dt, "NEW_SIGNAL_", ct2@value_VAR)

    history_item = list(normalizeSignalRPM = list(FUN = .normalizeSignalRPM, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        prof_dt = new_prof_dt,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("normalizeSignalRPM",
           function(ct2, mapped_reads_data = NULL, mapped_reads_VAR = "mapped_reads")
               standardGeneric("normalizeSignalRPM"),
           signature = "ct2")

#' @export
setMethod("normalizeSignalRPM", c("ChIPtsne2"), .normalizeSignalRPM)

#### cap normalizeSignalCapValue #####
.normalizeSignalCapValue = function(ct2, signal_cap_data = NULL, signal_cap_VAR = "cap_value", norm_to_1 = TRUE, trim_values_to_cap = TRUE, minimum_ceiling = NULL){
    args = get_args()
    if(is.null(signal_cap_data)){
        signal_cap_data = .form_merge_df(ct2, signal_cap_data, "signal_cap_data", signal_cap_VAR, "signal_cap_VAR")
    }else{
        if(is.numeric(signal_cap_data)){
            if(is.null(names(signal_cap_data))){
                stop("vector signal_cap_data must be named.")
            }
            tmp = data.frame(names(signal_cap_data), signal_cap_data)
            colnames(tmp) = c(ct2@name_VAR, signal_cap_VAR)
            signal_cap_data = tmp
        }
        if(!signal_cap_VAR %in% colnames(signal_cap_data)){
            stop("Supplied signal_cap_data does not contain ", signal_cap_VAR)
        }
    }
    if(!is.null(minimum_ceiling)){
        k = signal_cap_data[[signal_cap_VAR]] < minimum_ceiling
        signal_cap_data[[signal_cap_VAR]][k] = minimum_ceiling
    }

    new_prof_dt = getTidyProfile(ct2)

    new_prof_dt = seqsetvis::append_ynorm(
        new_prof_dt,
        cap_dt = signal_cap_data,
        cap_value_ = signal_cap_VAR,
        norm_value_ = "NEW_SIGNAL_",
        by1 = ct2@region_VAR,
        by2 = ct2@name_VAR,
        do_not_cap = !trim_values_to_cap,
        do_not_scaleTo1 = !norm_to_1,
    )

    new_prof_dt[[ct2@value_VAR]] = NULL
    new_prof_dt[[signal_cap_VAR]] = NULL
    data.table::setnames(new_prof_dt, "NEW_SIGNAL_", ct2@value_VAR)

    history_item = list(normalizeSignalCapValue = list(FUN = .normalizeSignalCapValue, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        prof_dt = new_prof_dt,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("normalizeSignalCapValue",
           function(ct2, signal_cap_data = NULL, signal_cap_VAR = "cap_value", norm_to_1 = TRUE, trim_values_to_cap = TRUE, minimum_ceiling = NULL)
               standardGeneric("normalizeSignalCapValue"),
           signature = "ct2")

#' @export
setMethod("normalizeSignalCapValue", c("ChIPtsne2"), .normalizeSignalCapValue)

#### calculateSignalCapValue ####

.calculateSignalCapValue = function(ct2, signal_cap_VAR = "cap_value", cap_quantile = .95){
    args = get_args()
    prof_dt = getTidyProfile(ct2)
    cap_dt = seqsetvis::calc_norm_factors(prof_dt,
                                 value_ = ct2@value_VAR,
                                 cap_value_ = signal_cap_VAR,
                                 by1 = ct2@region_VAR,
                                 by2 = ct2@name_VAR,
                                 aggFUN2 = function(x)quantile(x, cap_quantile))
    new_meta_dt = merge(getSampleMetaData(ct2), cap_dt, by = ct2@name_VAR)
    new_meta_dt = new_meta_dt[order(new_meta_dt[[ct2@name_VAR]]), ]

    history_item = list(calculateSignalCapValue = list(FUN = .calculateSignalCapValue, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        sample_metadata = new_meta_dt,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("calculateSignalCapValue",
           function(ct2, signal_cap_VAR = "cap_value", cap_quantile = .95)
               standardGeneric("calculateSignalCapValue"),
           signature = "ct2")

#' @export
setMethod("calculateSignalCapValue", c("ChIPtsne2"), .calculateSignalCapValue)
