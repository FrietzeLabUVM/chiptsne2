#
# harmonizeColData = function(to_harmonize, reference = ct2_ref){
#
# }
#
# harmonizeRowData = function(to_harmonize, reference = ct2_ref){
#
# }
#
# harmonizeStrand = function(to_harmonize, reference = ct2_ref){
#     args = get_args(to_ignore = "to_harmonize")
#     to_harmonize = to_harmonize[rownames(reference),]
#     if(is(reference, "ChIPtsne2")){
#         gr_ref = rowRanges(reference)
#         args$reference = gr_ref
#     }else if(is(reference, "GRanges")){
#         gr_ref = reference
#     }
#
#     gr_harm = rowRanges(to_harmonize)
#     stopifnot(all(start(gr_ref) == start(gr_harm)))
#     stopifnot(all(end(gr_ref) == end(gr_harm)))
#
#     ref_is_flipped = as.character(strand(gr_ref)) == "-"
#     harm_is_flipped = as.character(strand(gr_harm)) == "-"
#
#     needs_flipped = ref_is_flipped != harm_is_flipped
#     ct2_out = .flip_regions(to_harmonize, needs_flipped)
#     history_item = list(harmonizeStrand = list(FUN = harmonizeStrand, ARG = args))
#     ct2_out@metadata = c(ct2_out@metadata, history_item)
#     ct2_out
# }
#
#
