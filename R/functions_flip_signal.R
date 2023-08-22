#### flipProfilesToMatch ####

#' .centerSignal
#'
#' @param ct2 A ChIPtsne2 object
#'
#' @return A chiptsne2 object updated to reflect centering procedure. Width will be the same as original but this requires a second fetch.
#'
#' @importFrom seqsetvis centerGRangesAtMax
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2_left = flipProfilesToMatch(ct2)
#' ct2_right = flipProfilesToMatch(ct2, highest_on_right = TRUE)
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
#' prof_dt = prof_dt[, .(y = mean(y)), .(x, group)]
#' ggplot(prof_dt, aes(x = x, y = y, color = group)) + geom_path()
.flipProfilesToMatch = function(ct2, highest_on_right = FALSE){
    args = get_args()
    prof_dt = getTidyProfile(ct2)
    new_query_gr = rowRanges(ct2)

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
        GenomicRanges::strand(new_query_gr) = "+"
        GenomicRanges::strand(new_query_gr)[most_flipped$needs_flip] = "-"
        prof_dt = merge(prof_dt, most_flipped, by = c(ct2@region_VAR))
        remove(balance_dt)
        # prof_dt[needs_flip == TRUE, x := -x]
        flip_i = which(prof_dt$needs_flip)
        data.table::set(prof_dt, i = flip_i, j = ct2@position_VAR, value = -prof_dt[[ct2@position_VAR]][flip_i])
        prof_dt$needs_flip = NULL

    history_item = list(.flipProfilesToMatch = list(FUN = .flipProfilesToMatch, ARG = args))

    ChIPtsne2.from_tidy(prof_dt = prof_dt,
                        new_query_gr,
                        obj_history = c(ChIPtsne2.history(ct2), history_item),
                        init = FALSE)
}

#' @export
setGeneric("flipProfilesToMatch", function(ct2, highest_on_right = TRUE) standardGeneric("flipProfilesToMatch"))

#' @export
setMethod("flipProfilesToMatch", c("ChIPtsne2"), .flipProfilesToMatch)

