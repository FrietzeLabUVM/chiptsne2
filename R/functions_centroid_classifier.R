
calculate_centroid_per_group = function(ct2, group_VAR){
    ct2.r_sp = split(ct2, group_VAR)
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

    # install.packages("pdist")
    # library(pdist)
    dist_mat = as.matrix(pdist::pdist(r2rm, centroids))
    # colnames(dist_mat) = paste0(group_VAR, ":", rownames(centroids))
    rownames(dist_mat) = rownames(r2rm)
    colnames(dist_mat) = name_FUN(centroids)
    dist_mat
}

#' classify_by_centroid
#'
#' @param ct2
#' @param centroids
#'
#' @return
#' @import pdist
#'
#' @examples
classify_by_centroid_distances = function(distances, centroids, ambiguous_value = "ambiguous", tolerance = 0){
    apply(distances,1 , FUN = function(x){
        rnk = order(x)

        if(diff(x[rnk][1:2]) <= tolerance){
            ambiguous_value
        }else{
            rownames(centroids)[rnk][1]
        }
    })
}

#' calculateGroupCentroid
#'
#' @param ct2
#' @param group_VAR
#'
#' @return
#' @export
#'
#' @examples
calculateGroupCentroid = function(ct2, group_VAR){
    cent = calculate_centroid_per_group(ct2, group_VAR)
    cent
}

#' groupRegionsByCentroidDistance
#'
#' @param ct2
#' @param centroid
#' @param group_VAR
#' @param ambiguous_value
#' @param tolerance
#'
#' @return
#' @export
#'
#' @examples
groupRegionsByCentroidDistance = function(ct2, centroid, group_VAR, ambiguous_value = "ambiguous", tolerance = 0){
    message("groupRegionsByCentroidDistance ...")
    cent_dist = calculate_distance_to_centroids(ct2, centroid)
    new_grps = classify_by_centroid_distances(cent_dist, centroid, ambiguous_value = ambiguous_value, tolerance = tolerance)
    df_grps = data.frame(names(new_grps), new_grps)
    colnames(df_grps) = c(ct2@region_VAR, group_VAR)
    chiptsne2:::isolated_setRegionMetaData(ct2, df_grps)
}
