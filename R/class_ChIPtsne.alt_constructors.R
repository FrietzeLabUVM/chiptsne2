#' ChIPtsne2.from_tidy
#'
#' @param prof_dt Profile data.table, as returned from seqsetvis::ssvFetch* functions
#' @param query_gr The query GRanges object used to fetch prof_dt.
#' @param meta_dt Additional information for each name_VAR entry in prof_dt.
#' @param name_VAR Variable name that contains sample ids/names. Links prof_dt to meta_dt. Default is "sample".
#' @param position_VAR Variable name that contains positional information in prof_dt. Default is "x".
#' @param value_VAR Variable name that contains signal value information in prof_dt. Default is "y".
#' @param region_VAR Variable name that contains region ID information in prof_dt. Default is "id".
#'
#' @return ChIPtsne2 object with tsne results
#' @export
#'
#' @examples
#' query_gr = CTCF_in_10a_overlaps_gr
#' prof_dt = CTCF_in_10a_profiles_dt
#' ChIPtsne2.from_tidy(prof_dt, query_gr)
ChIPtsne2.from_tidy = function(prof_dt,
                               query_gr,
                               meta_dt = NULL,
                               name_VAR = "sample",
                               position_VAR = "x",
                               value_VAR = "y",
                               region_VAR = "id"){

    prof_dt[, setdiff(colnames(prof_dt), c("seqnames", "start", "end", "width", "strand")), with = FALSE]

    #create wide profile matrix
    tmp_wide = tidyr::pivot_wider(prof_dt, names_from = all_of(c(name_VAR, position_VAR)), values_from = value_VAR, id_cols = region_VAR)
    prof_mat = as.matrix(tmp_wide[, -1])
    rownames(prof_mat) = tmp_wide$id

    #create max assay
    prof_max = prof_dt[, .(y = max(y)), .(id, sample)] %>%
        tidyr::pivot_wider(names_from = name_VAR, id_cols = region_VAR, values_from = value_VAR)
    prof_max_mat = as.matrix(prof_max[, -1])
    rownames(prof_max_mat) = prof_max$id

    xy_dt = tsne_from_profile_mat(prof_mat)

    if(is.null(meta_dt)){
        meta_dt = unique(prof_dt[, c(name_VAR), with = FALSE])
    }


    map_dt = unique(prof_dt[, c(position_VAR, name_VAR), with = FALSE])
    map_dt[, cn := paste(name_VAR, position_VAR, sep = "_")]
    map_dt = dplyr::mutate(map_dt, cn = paste(get(name_VAR), get(position_VAR), sep = "_"))
    cn = map_dt$cn
    stopifnot(cn == colnames(prof_mat))
    map_dt = dplyr::mutate(map_dt, nr = seq(nrow(map_dt)))
    map_list = split(map_dt$n, map_dt[[name_VAR]])

    ChIPtsne2(assay = list(max = prof_max_mat),
              rowRanges = query_gr,
              rowToRowMat = prof_mat,
              colToRowMatCols = map_list,
              colData = meta_dt,
              metadata = list(time = date()))
}
