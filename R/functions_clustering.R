
.groupRegionByMembershipTable = function(ct)
.groupRegionManually
.groupRegionByDimReduceCluster

.groupRegionBySignalCluster = function(ct2, group_VAR = "cluster_id", n_clusters = 6){
    args = get_args()
    prof_dt = getTidyProfile(ct2)
    clust_dt = seqsetvis::ssvSignalClustering(prof_dt,
                                            nclust = n_clusters,
                                            fill_ =ct2@value_VAR,
                                            facet_ = ct2@name_VAR,
                                            column_ = ct2@position_VAR,
                                            row_ = ct2@region_VAR)
    assign_dt = clust_dt %>%
        dplyr::select(all_of(c(group_VAR, ct2@region_VAR))) %>%
        unique
    new_meta_dt = merge(getSampleMetaData(ct2), assign_dt, by = ct2@name_VAR)

    history_item = list(groupRegionBySignalCluster = list(FUN = .groupRegionBySignalCluster, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = new_meta_dt,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("groupRegionBySignalCluster",
           function(ct2, signal_cap_VAR = "cap_value", cap_quantile = .95)
               standardGeneric("groupRegionBySignalCluster"),
           signature = "ct2")

#' @export
setMethod("groupRegionBySignalCluster", c("ChIPtsne2"), .groupRegionBySignalCluster)
