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
    prof_dt = getTidyProfile(ct2)
    prof_dt = merge(prof_dt, mapped_reads_data, by = ct2@name_VAR)
    prof_dt = prof_dt %>%
        dplyr::mutate(NEW_SIGNAL_ = get(ct2@value_VAR) / get(mapped_reads_VAR) * 1e6)
    prof_dt[[ct2@value_VAR]] = NULL
    prof_dt[[mapped_reads_VAR]] = NULL
    data.table::setnames(prof_dt, "NEW_SIGNAL_", ct2@value_VAR)

    history_item = list(normalizeSignalRPM = list(FUN = .normalizeSignalRPM, ARG = args))
    ChIPtsne2.from_tidy(prof_dt,
                        rowRanges(ct2),
                        sample_metadata = colData(ct2),
                        name_VAR = ct2@name_VAR,
                        position_VAR = ct2@position_VAR,
                        value_VAR = ct2@value_VAR,
                        region_VAR = ct2@region_VAR,
                        obj_history = c(ChIPtsne2.history(ct2), history_item),
                        init = FALSE
    )
}

#' @export
setGeneric("normalizeSignalRPM",
           function(ct2, mapped_reads_data = NULL, mapped_reads_VAR = "mapped_reads")
               standardGeneric("normalizeSignalRPM"),
           signature = "ct2")

#' @export
setMethod("normalizeSignalRPM", c("ChIPtsne2"), .normalizeSignalRPM)

#### cap value #####
.normalizeSignalCapValue = function(ct2, signal_cap_data = NULL, signal_cap_VAR = "cap_value", norm_to_1 = TRUE, trim_values_to_cap = TRUE){
    args = get_args()
    signal_cap_data = .form_merge_df(ct2, signal_cap_data, "signal_cap_data", signal_cap_VAR, "signal_cap_VAR")

    prof_dt = getTidyProfile(ct2)

    prof_dt = seqsetvis::append_ynorm(
        prof_dt,
        cap_dt = signal_cap_data,
        cap_value_ = signal_cap_VAR,
        norm_value_ = "NEW_SIGNAL_",
        by1 = ct2@region_VAR,
        by2 = ct2@name_VAR,
        do_not_cap = !norm_to_1,
        do_not_scaleTo1 = !trim_values_to_cap
    )

    prof_dt[[ct2@value_VAR]] = NULL
    prof_dt[[signal_cap_VAR]] = NULL
    data.table::setnames(prof_dt, "NEW_SIGNAL_", ct2@value_VAR)

    history_item = list(normalizeSignalCapValue = list(FUN = .normalizeSignalCapValue, ARG = args))
    ChIPtsne2.from_tidy(prof_dt,
                        rowRanges(ct2),
                        sample_metadata = colData(ct2),
                        name_VAR = ct2@name_VAR,
                        position_VAR = ct2@position_VAR,
                        value_VAR = ct2@value_VAR,
                        region_VAR = ct2@region_VAR,
                        obj_history = c(ChIPtsne2.history(ct2), history_item),
                        init = FALSE
    )
}

#' @export
setGeneric("normalizeSignalCapValue",
           function(ct2, signal_cap_data = NULL, signal_cap_VAR = "cap_value", norm_to_1 = TRUE, trim_values_to_cap = TRUE)
               standardGeneric("normalizeSignalCapValue"),
           signature = "ct2")

#' @export
setMethod("normalizeSignalCapValue", c("ChIPtsne2"), .normalizeSignalCapValue)
