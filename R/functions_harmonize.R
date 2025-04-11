#
# ct2 = exampleChIPtsne2.with_meta()
# rowData(ct2)
# ct2.by_peak = split(ct2, c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
# rbind(ct2.by_peak$`FALSE FALSE`, ct2.by_peak$`FALSE TRUE`)
#
#
# ct2_ref = ct2.by_peak$`FALSE FALSE`
# ct2_harm = ct2.by_peak$`FALSE TRUE`
# #adding an inconsistent colData column normally causes an rbind error
# colData(ct2_harm)$col_test = "test"
# colData(ct2_ref)$col_test = "test2"
# rowData(ct2_harm)$row_test = "test"
# ct2_harm[seq(nrow(ct2_harm), 1),]
#
# ct2_new = rbind(ct2_ref, ct2_harm[seq(nrow(ct2_harm), 1),])
# ct2_new2 = rbind(ct2_harm[seq(nrow(ct2_harm), 1),], ct2_ref)
#
# colData(ct2_new)
#
#
# # ct2_harm2 = .harmonize_rbind(ct2_harm, ct2_ref)
#
# ct2_new = cbind(ct2_ref, ct2_harm2[seq(nrow(ct2_harm), 1),])
# ct2_new2 = cbind(ct2_harm2[seq(nrow(ct2_harm), 1),], ct2_ref)
#
# rowData(ct2_new)
# rowData(ct2_new2)
# .harmonize_rbind = function(to_harmonize, reference){
#     stop("NYI")
# }

.harmonize_cbind_list = function(ct2.l, allow_subset_match = FALSE, allow_missing_colData = TRUE){
    #simply return input if single item
    if(length(ct2.l) == 1) return(ct2.l)
    common_rn = unique(sapply(ct2.l, rownames))
    common_rn.start = common_rn
    for(i in seq(2, length(ct2.l))){
        common_rn = intersect(common_rn, rownames(ct2.l[[i]]))
    }
    if(length(common_rn) == 0){
        stop("No common rownames cound between ChIPtsne objects.")
    }
    if(!setequal(common_rn, common_rn.start)){
        message("Some subsetting of rows has occured to resolve mismatches.")
        ct2.l = lapply(ct2.l, function(x){
            x[common_rn,]
        })
    }
    ct2_ref = ct2.l[[1]]
    ct2.l = lapply(ct2.l, function(to_harmonize){
        .harmonize_cbind(to_harmonize, reference = ct2_ref, allow_subset_match = allow_subset_match, allow_missing_colData = allow_missing_colData)
    })

}

#' .harmonize_cbind
#'
#' @param to_harmonize ChIPtsne2 object to harmonize in preparation to cbind.
#' @param reference ChIPtsne2 or rowData to use as a reference when harmonizing.
#' @param allow_subset_match If TRUE, rows of to_harmonize may be a subset/superset of rows of reference. Default is FALSE.
#' @param allow_missing_colData If TRUE, colData may have additional attributes
#'
#' @return ChIPtsne2 object derived from to_harmonize that will be compatible to cbind with reference. If reference must be a modified further, a warning will be generated.
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2.by_cell = split(ct2, "cell")
#'
#' ct2_ref = ct2.by_cell$MCF10A
#' ct2_harm = ct2.by_cell$MCF10AT1
#' colData(ct2_harm)$col_test = "test"
#' rowData(ct2_harm)$row_test = "test"
#' rowData(ct2_harm)$peak_MCF10CA1_CTCF = TRUE
#' ct2_harm[seq(nrow(ct2_harm), 1),]
#'
#' ct2_harm2 = .harmonize_cbind(ct2_harm, ct2_ref)
#'
#' ct2_new = cbind(ct2_ref, ct2_harm2[seq(nrow(ct2_harm), 1),])
#' ct2_new2 = cbind(ct2_harm2[seq(nrow(ct2_harm), 1),], ct2_ref)
#'
#' rowData(ct2_new)
#' rowData(ct2_new2)
.harmonize_cbind = function(to_harmonize, reference, allow_subset_match = FALSE, allow_missing_colData = TRUE){
    args = get_args(to_ignore = "to_harmonize")
    to_harmonize = to_harmonize[rownames(reference),]
    if(is(reference, "ChIPtsne2_no_rowRanges")){
        df_ref = rowData(reference)
        args$reference = df_ref
    }else if(is(reference, "ChIPtsne2_no_rowRanges")){
        df_ref = reference
    }

    #verify rownames
    if(!setequal(rownames(to_harmonize), rownames(reference))){
        shared_rn = intersect(rownames(to_harmonize), rownames(reference))
        if(length(shared_rn) == 0){
            stop("There are no shared rownames with reference. Names have either been modified or objects are not compatible due to representing different region sets.")
        }else{
            if(!allow_subset_match){
                frac_harm = length(shared_rn) / nrow(to_harmonize)
                frac_ref = length(shared_rn) / nrow(reference)
                stop(
                    "There is an incomplete match with rownames of reference.\n",
                    paste(round(100*frac_harm, digits = 2), "% of object to harmonize matches.\n"),
                    paste(round(100*frac_ref, digits = 2), "% of reference matches.\n"),
                    "To allow this mismatch, enable with allow_subset_match = TRUE."
                )
            }else{
                frac_harm = length(shared_rn) / nrow(to_harmonize)
                frac_ref = length(shared_rn) / nrow(reference)
                message(
                    "With allow_subset_match = TRUE some rows have been dropped to match reference.\n",
                    paste(round(100*(1-frac_harm), digits = 2), "% of rows dropped..\n")
                )
                if(frac_ref < 1){
                    warning("You MUST call this function with input reference as to_harmonize swapped and allow_subset_match = TRUE.")
                }
            }
        }
    }

    common_rn = rownames(reference)[rownames(reference) %in% rownames(to_harmonize)]

    to_harmonize = to_harmonize[common_rn,]
    df_ref = df_ref[common_rn,]

    df_harm = rowData(to_harmonize)
    shared_cn = intersect(colnames(df_harm), colnames(df_ref))

    #cbind requires that values be identical for any shared colnames, this defaults to reference to resolve conflict
    for(cn in shared_cn){
        if(any(df_harm[[cn]] != df_ref[[cn]])){
            message("rowData attribute `", cn, "` forced equal to reference values.")
            df_harm[[cn]] = df_ref[[cn]]
        }
    }
    shared_colData = intersect(colnames(colData(to_harmonize)), colnames(colData(reference)))
    dropped_colData = setdiff(colnames(colData(to_harmonize)), shared_colData)
    if(length(dropped_colData) > 0){
        message("Dropping incompatible colData attributes:\n",
                paste(paste0("  `", dropped_colData, "`"), collapse = "\n"))
    }
    new_colData = colData(to_harmonize)[, shared_colData]
    cloneChIPtsne2(ct2 = to_harmonize[common_rn,], new_rowData = df_harm, new_colData = new_colData)
}

harmonizeStrand = function(to_harmonize, reference, track_history = TRUE){
    args = get_args(to_ignore = "to_harmonize")
    to_harmonize = to_harmonize[rownames(reference),]
    if(is(reference, "ChIPtsne2")){
        gr_ref = rowRanges(reference)
        args$reference = gr_ref
    }else if(is(reference, "GRanges")){
        gr_ref = reference
    }

    gr_harm = rowRanges(to_harmonize)
    stopifnot(all(start(gr_ref) == start(gr_harm)))
    stopifnot(all(end(gr_ref) == end(gr_harm)))

    ref_is_flipped = as.character(strand(gr_ref)) == "-"
    harm_is_flipped = as.character(strand(gr_harm)) == "-"

    needs_flipped = ref_is_flipped != harm_is_flipped
    ct2_out = .flip_regions(to_harmonize, needs_flipped)
    if(track_history){
        history_item = list(harmonizeStrand = list(FUN = harmonizeStrand, ARG = args))
        ct2_out@metadata = c(ct2_out@metadata, history_item)
    }
    ct2_out
}


