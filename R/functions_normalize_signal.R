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
    message("normalizeSignalRPM ...")
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
        new_prof_dt = new_prof_dt,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' normalizeSignalRPM
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param mapped_reads_data Optionally, a data.frame containing `ct2`'s name
#'   variable and specified `mapped_reads_VAR`. Otherwise colData/sample
#'   metadata will be used.
#' @param mapped_reads_VAR Variable name in colData/sample metadata or
#'   `mapped_reads_data` if provided.
#'
#' @return `r doc_ct2_nrr()` where value_VAR has been normalized according to
#'   `mapped_reads_VAR`
#' @export
#' @rdname ct2-normrpm
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' # we're going to spoof in mapped reads here, but if you have bam files you can use [get_mapped_reads].
#' colData(ct2)$mapped_reads = c(3e6, 2e6, 1e6)
#' plotSignalLinePlot(ct2)
#' ct2.rpm = normalizeSignalRPM(ct2)
#' # RPM normalization corrects for depth, reversing magnitude trend in this exampl data.
#' plotSignalLinePlot(ct2.rpm)
setGeneric("normalizeSignalRPM",
           function(ct2, mapped_reads_data = NULL, mapped_reads_VAR = "mapped_reads")
               standardGeneric("normalizeSignalRPM"),
           signature = "ct2")

#' @export
#' @rdname ct2-normrpm
setMethod("normalizeSignalRPM", c("ChIPtsne2_no_rowRanges"), .normalizeSignalRPM)

#### cap normalizeSignalCapValue #####
.normalizeSignalCapValue = function(ct2, signal_cap_data = NULL, signal_cap_VAR = "cap_value", norm_to_1 = TRUE, trim_values_to_cap = TRUE, cap_floor = NULL){
    message("normalizeSignalCapValue ...")
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
    if(!is.null(cap_floor)){
        k = signal_cap_data[[signal_cap_VAR]] < cap_floor
        signal_cap_data[[signal_cap_VAR]][k] = cap_floor
    }

    if(!signal_cap_VAR %in% colnames(signal_cap_data)){
        stop("signal_cap_VAR: \"", signal_cap_VAR, "\" expected in colnames of signal_cap_data but was not found!")
    }
    if(!ct2@name_VAR %in% colnames(signal_cap_data)){
        stop("ChIPtsne2 name variable: \"", ct2@name_VAR, "\" expected in colnames of signal_cap_data but was not found!")
    }
    if(!all(colnames(ct2) %in% signal_cap_data[[ct2@name_VAR]])){
        err_extra = paste(setdiff(colnames(ct2), signal_cap_data[[ct2@name_VAR]]), collapse = "\n")
        stop(paste0("Some or all values of ChIPtsne2 name variable are missing from provided signal_cap_data.\n", err_extra))

    }

    new_prof_dt = getTidyProfile(ct2)

    new_prof_dt = seqsetvis::append_ynorm(
        new_prof_dt,
        cap_dt = signal_cap_data,
        value_ = ct2@value_VAR,
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
        new_prof_dt = new_prof_dt,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' normalizeSignalCapValue
#'
#' Applies a cap value for each sample. Meant to be run after
#' [calculateSignalCapValue] or after creating a data.frame to supply to
#' `signal_cap_data`.
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param signal_cap_data An optional data.frame containing name_VAR in `ct2`
#'   and `signal_cap_VAR`. Useful if you have calculated cap values outside of
#'   ChIPtsne2 functions or are reapplying previously calculated values. Does
#'   not require [calculateSignalCapValue] to be run.
#' @param signal_cap_VAR The attribute name in colData/sample metadata where cap
#'   data results have been added. Should be added with
#'   [calculateSignalCapValue]
#' @param norm_to_1 If TRUE, all values will be divided by cap_VAR values.
#'   Default is TRUE.
#' @param trim_values_to_cap If TRUE, all values will be capped at cap value
#'   prior to dividing. Default is TRUE.
#' @param cap_floor Optional floor for cap values. Useful when input or control
#'   samples are present that should not have true signal. Their calculated cap
#'   values would be so low they'd only magnify background noise. Default of
#'   NULL does nothing.
#' @rdname ct2-normcap
#'
#' @return `r doc_ct2_nrr()` where value_VAR has been normalized according to
#'   `signal_cap_VAR`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = calculateSignalCapValue(ct2, signal_cap_VAR = "higher_cap_value", cap_quantile = .98)
#' ct2.norm1 = normalizeSignalCapValue(
#'   ct2,
#'   signal_cap_VAR = "higher_cap_value"
#' )
#'
#' # we can extract sample meta data and modify cap value
#' # then manually supply using the signal_cap_data parameter
#' cap_df = getSampleMetaData(ct2)
#' # by lowering the cap value we'll increase the final normalized values across the board.
#' cap_df$higher_cap_value = cap_df$higher_cap_value / 2
#' ct2.norm2 = normalizeSignalCapValue(
#'   ct2,
#'   signal_cap_VAR = "higher_cap_value",
#'   signal_cap_data = cap_df
#' )
#'
#' # see the impact of each normalization
#' plotSignalLinePlot(ct2)
#' plotSignalLinePlot(ct2.norm1)
#' plotSignalLinePlot(ct2.norm2)
setGeneric("normalizeSignalCapValue",
           function(ct2, signal_cap_data = NULL, signal_cap_VAR = "cap_value", norm_to_1 = TRUE, trim_values_to_cap = TRUE, cap_floor = NULL)
               standardGeneric("normalizeSignalCapValue"),
           signature = "ct2")

#' @export
#' @rdname ct2-normcap
setMethod("normalizeSignalCapValue", c("ChIPtsne2_no_rowRanges"), .normalizeSignalCapValue)

#### calculateSignalCapValue ####

#' @importFrom stats quantile
.calculateSignalCapValue = function(ct2, signal_cap_VAR = "cap_value", cap_quantile = .95){
    message("calculateSignalCapValue ...")
    args = get_args()
    prof_dt = getTidyProfile(ct2)
    cap_dt = seqsetvis::calc_norm_factors(prof_dt,
                                 value_ = ct2@value_VAR,
                                 cap_value_ = signal_cap_VAR,
                                 by1 = ct2@region_VAR,
                                 by2 = ct2@name_VAR,
                                 aggFUN2 = function(x)stats::quantile(x, cap_quantile))

    # remove signal_cap_VAR if present to overwrite
    old_meta_dt = getSampleMetaData(ct2)
    old_meta_dt = old_meta_dt[, setdiff(colnames(old_meta_dt), signal_cap_VAR), drop = FALSE]

    new_meta_dt = merge(old_meta_dt, cap_dt, by = ct2@name_VAR)
    new_meta_dt = new_meta_dt[order(new_meta_dt[[ct2@name_VAR]]), ]

    history_item = list(calculateSignalCapValue = list(FUN = .calculateSignalCapValue, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        new_sample_metadata = new_meta_dt,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' calculateSignalCapValue
#'
#' Calculates but does not apply a cap value for each sample. Meant to be run before [normalizeSignalCapValue].
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param signal_cap_VAR The attribute name to store results in colData/sample metadata. Will be overwritten if present.
#' @param cap_quantile The quantile at which to set cap value. Default is .95.
#'
#' @export
#' @return `r doc_ct2_nrr()` where colData/sample metadata has had `signal_cap_VAR` added.
#' @rdname ct2-calccap
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = calculateSignalCapValue(ct2)
#' ct2 = calculateSignalCapValue(ct2, signal_cap_VAR = "higher_cap_value", cap_quantile = .98)
#' colData(ct2)
setGeneric("calculateSignalCapValue",
           function(ct2, signal_cap_VAR = "cap_value", cap_quantile = .95)
               standardGeneric("calculateSignalCapValue"),
           signature = "ct2")

#' @export
#' @rdname ct2-calccap
setMethod("calculateSignalCapValue", c("ChIPtsne2_no_rowRanges"), .calculateSignalCapValue)
