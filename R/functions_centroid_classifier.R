
calculate_centroid_per_group = function(ct2, group_VARS){
    if(length(group_VARS) > 1){
        tmp_df = data.frame(TMP_GROUP__ = apply(GenomicRanges::mcols(rowRanges(ct2))[,group_VARS], 1, paste, collapse = ","))
        ct2 = chiptsne2::setRegionMetaData(ct2, tmp_df)
        rowRanges(ct2)
        group_VARS = "TMP_GROUP__"
    }
    ct2.r_sp = split(ct2, group_VARS)
    centroids = t(sapply(ct2.r_sp, function(x){
        colMeans(rowToRowMat(x))
    }))
    centroids
}

calculate_distance_to_centroids = function(ct2, centroids, name_FUN = function(x){rownames(x)}){
    r2rm = rowToRowMat(ct2)
    euclidean_distance <- function(p,q){
        sqrt(sum((p - q)^2))
    }
    dist_mat = as.matrix(pdist::pdist(r2rm, centroids))
    rownames(dist_mat) = rownames(r2rm)
    colnames(dist_mat) = name_FUN(centroids)
    dist_mat
}


classify_by_centroid_distances = function(distances, centroids,
                                          ambiguous_value = "ambiguous_match", ambiguous_distance = 0,
                                          no_match_value = "no_match", match_distance = Inf){
    apply(distances,1 , FUN = function(x){
        rnk = order(x)

        if(diff(x[rnk][1:2]) <= ambiguous_distance){
            ambiguous_value
        }else if(x[rnk][1] < match_distance){
            no_match_value
        }else{
            rownames(centroids)[rnk][1]
        }
    })
}

#' calculateGroupCentroid
#'
#' @param ct2 `R doc_ct2()`
#' @param group_VARS A least 1 categorical attribute in rowData of `ct2`. When
#'   multiple attributes are specified, centroids will be calculated for all
#'   combinations of groups.
#'
#' @return A matrix of centroid profiles per group. Input to
#'   [groupRegionsByCentroidDistance].
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' cents = calculateGroupCentroid(ct2, c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
#' cents
#' ct2 = groupRegionsByCentroidDistance(ct2, cents, group_VAR = "peak_overlap_profile")
#' rowData(ct2)
calculateGroupCentroid = function(ct2, group_VARS){
    cent = calculate_centroid_per_group(ct2, group_VARS)
    cent
}

#' groupRegionsByCentroidDistance
#'
#' @param ct2 `r doc_ct2()`
#' @param centroid A matrix of centroid profiles. Use [calculateGroupCentroid]
#'   to generate.
#' @param group_VAR `r doc_group_VAR()`
#' @param ambiguous_value Value used to indicate regions/rows that don't
#'   unambiguously match centroid profiles.
#' @param ambiguous_distance Increase to make the definition of ambiguous
#'   mathces more lenient. When a row/region matches multiple centroids with
#'   distance values within `ambiguous_distance`, it is considered ambiguous and
#'   assigned `ambiguous_value`.
#' @param no_match_value  Value used to indicate regions/rows that don't closely
#'   match any centroid profiles.
#' @param match_distance Maximum allowed distance for a row/region to be
#'   considered matching a centroid profile. If distances to all centroid
#'   profiles is over `match_distance`, it is assigned the value of
#'   `no_match_value`.
#'
#'   With default values of ambiguous_distance = 0 and match_distance = Inf, all
#'   rows/regions will be assigned to a single centroid.
#'
#' @return `r doc_return_group()`
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' cents = calculateGroupCentroid(ct2, c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
#' ct2 = groupRegionsByCentroidDistance(ct2, cents, group_VAR = "peak_overlap_profile")
#' rowData(ct2)
groupRegionsByCentroidDistance = function(ct2,
                                          centroid,
                                          group_VAR,
                                          ambiguous_value = "ambiguous_match",
                                          ambiguous_distance = 0,
                                          no_match_value = "no_match",
                                          match_distance = Inf){
    message("groupRegionsByCentroidDistance ...")
    cent_dist = calculate_distance_to_centroids(ct2, centroid)
    new_grps = classify_by_centroid_distances(
        cent_dist,
        centroid,
        ambiguous_value = ambiguous_value,
        ambiguous_distance = ambiguous_distance,
        no_match_value = no_match_value,
        match_distance = match_distance)
    df_grps = data.frame(names(new_grps), new_grps)
    colnames(df_grps) = c(ct2@region_VAR, group_VAR)
    setRegionMetaData.no_history(ct2, df_grps)
}
