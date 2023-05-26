.normalizeSignalRPM = function(ct2, mapped_reads_data = NULL, mapped_reads_VAR = "mapped_reads"){
    args = get_args()
    if(is.null(mapped_reads_data)){
        if(is.null(colData(ct2)[[mapped_reads_VAR]])){
            stop("mapped_reads_data was not supplied and mapped_reads_VAR: ", mapped_reads_VAR, " is not in colData")
        }
        mapped_reads_data = getSampleMetaData(ct2)[, c(ct2@name_VAR, mapped_reads_VAR)]
    }
    #handle numeric input and reform to data.frame
    if(is.numeric(mapped_reads_data)){
        if(is.null(names(mapped_reads_data))){
            stop("When mapped_reads_data is supplied, names must be set.")
        }else if(!setequal(colData(ct2)[[ct2@name_VAR]], names(mapped_reads_data))){
            stop("names of mapped_reads_data is not setequal to name_VAR: ", ct2@name_VAR, " in colData.")
        }
        mapped_reads_data = data.frame(names(mapped_reads_data), mapped_reads_data)
        colnames(mapped_reads_data) = c(ct2@name_VAR, mapped_reads_VAR)
    }
    #merge mapped_reads_data to tidy profile
    browser()
    prof_dt = getTidyProfile(ct2)
    merge(prof_dt, mapped_reads_data, by = ct2@name_VAR)

    history_item = list(normalizeSignalRPM = list(FUN = .normalizeSignalRPM, ARG = args))
    ChIPtsne2.from_tidy(new_prof_dt,
                        new_query_gr,
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

