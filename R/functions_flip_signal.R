#### flipProfilesToMatch ####

.flipProfilesToMatch = function(ct2, highest_on_right = FALSE){
    #visible binding NOTE
    left_sum = right_sum = needs_flip = `:=` = .N = fraction_flipped = TMP_POS_I__ = NEW_POS_X__ = NULL
    args = get_args()
    message("flipProfilesToMatch ...")
    prof_dt = getTidyProfile(ct2)
    new_rowRanges = rowRanges(ct2)

    balance_dt = prof_dt[, list(right_sum = sum(get(ct2@value_VAR)[get(ct2@position_VAR) > 0]),
                                left_sum = sum(get(ct2@value_VAR)[get(ct2@position_VAR) < 0])),
                         by = c(ct2@region_VAR, ct2@name_VAR)]
    balance_dt = balance_dt[, list(needs_flip = left_sum > right_sum),
                            c(ct2@region_VAR, ct2@name_VAR)]
    most_flipped = balance_dt[,
                              list(fraction_flipped = sum(needs_flip) / .N),
                              by = c(ct2@region_VAR)]
    most_flipped[, needs_flip := fraction_flipped > .5]
    if(!highest_on_right){
        most_flipped$needs_flip = !most_flipped$needs_flip
    }
    most_flipped$fraction_flipped = NULL
    GenomicRanges::strand(new_rowRanges) = "+"
    GenomicRanges::strand(new_rowRanges)[most_flipped$needs_flip] = "-"
    prof_dt = merge(prof_dt, most_flipped, by = c(ct2@region_VAR))
    prof_dt = prof_dt[order(get(ct2@position_VAR))]
    x_vals = unique(prof_dt[[ct2@position_VAR]])
    remove(balance_dt)
    flip_i = which(prof_dt$needs_flip)
    #it is not enough to simply negate position values when win_size is odd
    # data.table::set(prof_dt, i = flip_i, j = ct2@position_VAR, value = -prof_dt[[ct2@position_VAR]][flip_i])
    prof_dt[, TMP_POS_I__ := seq(.N), c(ct2@region_VAR, ct2@name_VAR)]
    prof_dt[needs_flip == TRUE, TMP_POS_I__ := seq(.N, 1), c(ct2@region_VAR, ct2@name_VAR)]
    prof_dt[, NEW_POS_X__ := x_vals[TMP_POS_I__]]
    prof_dt[[ct2@position_VAR]] = NULL
    prof_dt$TMP_POS_I__ = NULL
    prof_dt$needs_flip = NULL
    data.table::setnames(prof_dt, "NEW_POS_X__", ct2@position_VAR)
    prof_dt = prof_dt[order(get(ct2@position_VAR))]

    history_item = list(.flipProfilesToMatch = list(FUN = .flipProfilesToMatch, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2,
        new_prof_dt = prof_dt,
        new_rowRanges = new_rowRanges,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item))
}


#' flipProfilesToMatch
#'
#' For unstranded data, a peak with higher signal on one side is not
#' meaningfully different from its mirror image. This procedure identifies
#' "tilted" peaks and flips those with more signal right-of-center. By flipping,
#' rowRanges becomes stranded, flipped regions get assigned (-) strand with
#' uaffected regions getting (+) strand
#'
#' @param ct2 A ChIPtsne2 object
#' @param highest_on_right
#'
#' @return A chiptsne2 object updated such that signal "tilts" in the same way.
#'   When signal is flipped, strand of rowRanges is set to negative.
#'
#' @importFrom seqsetvis centerGRangesAtMax
#' @export
#' @rdname ct2-flip
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' #add a metadata column
#' colData(ct2)$flip = "none"
#' ct2_left = flipProfilesToMatch(ct2)
#' colData(ct2_left)$flip = "left"
#' #colnames will need to be different
#' colnames(ct2_left) = paste0(colnames(ct2_left), "_left")
#' ct2_right = flipProfilesToMatch(ct2, highest_on_right = TRUE)
#' colData(ct2_right)$flip = "right"
#' colnames(ct2_right) = paste0(colnames(ct2_right), "_right")
#'
#' # flipping imposes strandedness on regions so rowRanges needs to be removed
#' rowRanges(ct2) = NULL
#' rowRanges(ct2_left) = NULL
#' rowRanges(ct2_right) = NULL
#' ct2.combined = cbind(ct2, ct2_left, ct2_right)
#'
#' debug(plotSignalLinePlot)
#' plotSignalLinePlot(ct2.combined, facet_VAR = "cell", color_VAR = "flip")
#'
#' prof_original = getTidyProfile(ct2)
#' prof_original$group = "original"
#'
#' prof_left = getTidyProfile(ct2_left)
#' prof_left$group = "flip_left"
#'
#' prof_right = getTidyProfile(ct2_right)
#' prof_right$group = "flip_right"
#'
#' prof_dt = rbind(prof_left, prof_right, prof_original)
#' prof_dt = prof_dt[, list(y = mean(y)), .(x, group)]
#'
#' library(ggplot2)
#' ggplot(prof_dt, aes(x = x, y = y, color = group)) + geom_path()
setGeneric("flipProfilesToMatch", function(ct2, highest_on_right = TRUE) standardGeneric("flipProfilesToMatch"))

#' @export
#' @rdname ct2-flip
setMethod("flipProfilesToMatch", c("ChIPtsne2_no_rowRanges"), .flipProfilesToMatch)

