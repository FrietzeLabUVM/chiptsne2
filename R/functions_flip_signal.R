#### flipProfilesToMatch ####

.flipProfilesToMatch = function(ct2, highest_on_right = FALSE){
    #visible binding NOTE
    left_sum = right_sum = needs_flip = `:=` = .N = fraction_flipped = TMP_POS_I__ = NEW_POS_X__ = NULL
    args = get_args()
    message("flipProfilesToMatch ...")

    # to_flip =
    c2rmc = colToRowMatCols(ct2)
    r2rm = rowToRowMat(ct2)
    left_cn = character()
    right_cn = character()
    for(cn in c2rmc){
        left_cn = c(left_cn, cn[seq(1, length(cn)/2)])
        right_cn = c(right_cn, setdiff(cn, left_cn))
    }
    left_sums = rowSums(r2rm[, left_cn, drop = FALSE])
    right_sums = rowSums(r2rm[, right_cn, drop = FALSE])

    to_flip = right_sums > left_sums
    if(highest_on_right) to_flip = !to_flip
    ct2.flipped = .flip_regions(ct2, to_flip)


    # prof_dt = getTidyProfile(ct2)
    # new_rowRanges = rowRanges(ct2)
    #
    # balance_dt = prof_dt[, list(right_sum = sum(get(ct2@value_VAR)[get(ct2@position_VAR) > 0]),
    #                             left_sum = sum(get(ct2@value_VAR)[get(ct2@position_VAR) < 0])),
    #                      by = c(ct2@region_VAR, ct2@name_VAR)]
    # balance_dt = balance_dt[, list(needs_flip = sum(left_sum) > sum(right_sum)),
    #                         c(ct2@region_VAR)]
    # # most_flipped = balance_dt[,
    # #                           list(fraction_flipped = sum(needs_flip) / .N),
    # #                           by = c(ct2@region_VAR)]
    # # most_flipped[, needs_flip := fraction_flipped > .5]
    # if(!highest_on_right){
    #     balance_dt$needs_flip = !balance_dt$needs_flip
    # }
    # # most_flipped$fraction_flipped = NULL
    # GenomicRanges::strand(new_rowRanges) = "+"
    # GenomicRanges::strand(new_rowRanges)[balance_dt$needs_flip] = "-"
    # prof_dt = merge(prof_dt, balance_dt, by = c(ct2@region_VAR))
    # prof_dt = prof_dt[order(get(ct2@position_VAR))]
    # x_vals = unique(prof_dt[[ct2@position_VAR]])
    # remove(balance_dt)
    # flip_i = which(prof_dt$needs_flip)
    # #it is not enough to simply negate position values when win_size is odd
    # # data.table::set(prof_dt, i = flip_i, j = ct2@position_VAR, value = -prof_dt[[ct2@position_VAR]][flip_i])
    # prof_dt[, TMP_POS_I__ := seq(.N), c(ct2@region_VAR, ct2@name_VAR)]
    # prof_dt[needs_flip == TRUE, TMP_POS_I__ := seq(.N, 1), c(ct2@region_VAR, ct2@name_VAR)]
    # prof_dt[, NEW_POS_X__ := x_vals[TMP_POS_I__]]
    # prof_dt[[ct2@position_VAR]] = NULL
    # prof_dt$TMP_POS_I__ = NULL
    # prof_dt$needs_flip = NULL
    # data.table::setnames(prof_dt, "NEW_POS_X__", ct2@position_VAR)
    # prof_dt = prof_dt[order(get(ct2@position_VAR))]

    # cloneChIPtsne2_fromTidy(
    #     ct2,
    #     new_prof_dt = prof_dt,
    #     new_rowRanges = new_rowRanges,
    #     new_obj_history = c(ChIPtsne2.history(ct2), history_item))

    history_item = list(flipProfilesToMatch = list(FUN = .flipProfilesToMatch, ARG = args))

    ct2.flipped@metadata = c(ct2.flipped@metadata, history_item)
    ct2.flipped

}


#' flipProfilesToMatch
#'
#' For unstranded data, a peak with higher signal on one side is not
#' meaningfully different from its mirror image. This procedure identifies
#' "tilted" peaks and flips those with more signal right-of-center. By flipping,
#' rowRanges becomes stranded, flipped regions get assigned (-) strand with
#' unaffected regions getting (+) strand
#'
#' @param ct2 A ChIPtsne2 object
#' @param highest_on_right If TRUE, majority of signal will be on the right half rather than left. Default is FALSE.
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
#' prof_dt = prof_dt[, list(value = mean(value)), .(position, group)]
#'
#' library(ggplot2)
#' ggplot(prof_dt, aes(x = position, y = value, color = group)) + geom_path()
setGeneric("flipProfilesToMatch", function(ct2, highest_on_right = FALSE) standardGeneric("flipProfilesToMatch"))

#' @export
#' @rdname ct2-flip
setMethod("flipProfilesToMatch", c("ChIPtsne2_no_rowRanges"), .flipProfilesToMatch)


#' .flip_regions
#'
#' @param ct2 valid chiptsne2 object
#' @param to_flip selector for rownames in ct2. accepts logical or character vector.
#'
#' @return ChIPtsne2 object with regions flipped according to to_flip.
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = groupRegionsBySignalCluster(ct2)
#' ct2 = sortRegions(ct2, sort_strategy = "left", group_VAR = "cluster_id")
#' p1 = plotSignalHeatmap(ct2, sort_strategy = "none", group_VARS = "cluster_id")
#'
#' to_flip = rowData(ct2)$cluster_id == 1
#' ct2_flipped = .flip_regions(ct2, to_flip)
#' p2 = plotSignalHeatmap(ct2_flipped, sort_strategy = "none", group_VARS = "cluster_id")
#' cowplot::plot_grid(p1, p2)
#'
#' ct2_flipped2 = .flip_regions(ct2, which(to_flip))
#' p3 = plotSignalHeatmap(ct2_flipped2, sort_strategy = "none", group_VARS = "cluster_id")
#' cowplot::plot_grid(p1, p3)
#'
#' ct2_flipped3 = .flip_regions(ct2, rownames(ct2)[to_flip])
#' p4 = plotSignalHeatmap(ct2_flipped3, sort_strategy = "none", group_VARS = "cluster_id")
#' cowplot::plot_grid(p1, p4)
.flip_regions = function(ct2, to_flip){

    reg_var = getRegionVariable(ct2)
    pos_var = getPositionVariable(ct2)
    nam_var = getNameVariable(ct2)

    #convert to logical for easier downstream handling
    if(is.numeric(to_flip)){
        to_flip = seq_len(nrow(ct2)) %in% to_flip
    }
    if(is.character(to_flip)){
        to_flip = rownames(ct2) %in% to_flip
    }
    if(is.factor(to_flip)){
        to_flip = rownames(ct2) %in% as.character(to_flip)
    }
    stopifnot(is.logical(to_flip))

    ct2.to_flip = ct2[to_flip,]
    new_rowRanges = rowRanges(ct2.to_flip)
    GenomicRanges::strand(new_rowRanges)[GenomicRanges::strand(new_rowRanges) == "*"] = "+"
    new_rowRanges = GenomicRanges::invertStrand(new_rowRanges)

    # iterate through each sample and reverse columns
    c2rmc = colToRowMatCols(ct2.to_flip)
    new_r2rm = rowToRowMat(ct2.to_flip)
    for(cn in c2rmc){
        tmp = new_r2rm[, rev(cn), drop = FALSE]
        colnames(tmp) = cn
        new_r2rm[, cn] = tmp
    }
    # rowToRowMat(ct2.to_flip) = new_r2rm





    # prof_dt = getTidyProfile(ct2.to_flip)
    # prof_dt = prof_dt[order(get(pos_var))]
    # x_vals = unique(prof_dt[[pos_var]])
    # #it is not enough to simply negate position values when win_size is odd
    # # data.table::set(prof_dt, i = flip_i, j = pos_var, value = -prof_dt[[pos_var]][flip_i])
    # prof_dt[, TMP_POS_I__ := seq(.N), c(reg_var, nam_var)]
    # prof_dt[, TMP_POS_I__ := seq(.N, 1), c(reg_var, nam_var)]
    # prof_dt[, NEW_POS_X__ := x_vals[TMP_POS_I__]]
    # prof_dt[[pos_var]] = NULL
    # prof_dt$TMP_POS_I__ = NULL
    # data.table::setnames(prof_dt, "NEW_POS_X__", pos_var)
    # prof_dt = prof_dt[order(get(pos_var))]

    ct2.flipped = cloneChIPtsne2(
        ct2.to_flip,
        new_rowToRowMat = new_r2rm,
        new_rowRanges = new_rowRanges
    )
    ct2_out = rbind(
        ct2[!to_flip,],
        ct2.flipped
    )
    ct2_out = ct2_out[rownames(ct2),]
    #if input data was unstranded, default to + strand.
    GenomicRanges::strand(rowRanges(ct2_out))[GenomicRanges::strand(rowRanges(ct2_out)) == "*"] = "+"
    #I'm not sure why this doesn't drop unused * level on assignment
    #GenomicRanges::strand(rowRanges(ct2_out)) = GenomicRanges::droplevels(GenomicRanges::strand(rowRanges(ct2_out)))
    ct2_out
}


