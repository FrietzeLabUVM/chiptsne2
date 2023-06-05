#### signal clustering ####
.groupRegionBySignalCluster = function(ct2, group_VAR = "cluster_id", n_clusters = 6){
    args = get_args()
    prof_dt = getTidyProfile(ct2)
    clust_dt = seqsetvis::ssvSignalClustering(prof_dt,
                                              cluster_ = group_VAR,
                                              nclust = n_clusters,
                                              fill_ =ct2@value_VAR,
                                              facet_ = ct2@name_VAR,
                                              column_ = ct2@position_VAR,
                                              row_ = ct2@region_VAR)
    assign_dt = clust_dt %>%
        dplyr::select(all_of(c(group_VAR, ct2@region_VAR))) %>%
        unique
    # new_query_gr = .add_region_metadata(rowRanges(ct2), region_metadata, region_VAR)


    history_item = list(groupRegionBySignalCluster = list(FUN = .groupRegionBySignalCluster, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = assign_dt,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("groupRegionBySignalCluster",
           function(ct2, group_VAR = "cluster_id", n_clusters = 6)
               standardGeneric("groupRegionBySignalCluster"),
           signature = "ct2")

#' @export
setMethod("groupRegionBySignalCluster", c("ChIPtsne2"), .groupRegionBySignalCluster)



#### dim clustering ####
.knn_clustering = function(xy_df,
                           nn = 100,
                           id_var = "id"){
    valid_vars = c(id_var, "tx", "ty")
    stopifnot(valid_vars %in% colnames(xy_df))

    if(nn > .2 * nrow(xy_df)){
        nn = floor(.2 * nrow(xy_df))
        message("Decreasing nearest-neighbors to ", nn, ".  Original value was too high for dataset.")
    }
    mat = t(as.matrix(xy_df[, c("tx", "ty")]))
    colnames(mat) = xy_df[[id_var]]
    knn.info <- RANN::nn2(t(mat), k = nn)
    knn <- knn.info$nn.idx
    colnames(knn) = c("tid", paste0("V", seq(nn - 1)))
    knn = data.table::as.data.table(knn)
    mknn = reshape2::melt(knn, id.vars = "tid")
    ADJ = Matrix::Matrix(0, ncol(mat), ncol(mat))
    ADJ[cbind(mknn$tid, mknn$value)] = 1
    rownames(ADJ) = colnames(mat)
    colnames(ADJ) = colnames(mat)
    g <- igraph::graph.adjacency(ADJ, mode = "undirected")
    g <- igraph::simplify(g)
    km <- igraph::cluster_walktrap(g)
    com <- km$membership
    names(com) <- km$names
    com_dt = data.table::data.table(tid = names(com), cluster_id = com)
    data.table::setnames(com_dt, "tid", id_var)
    p_dt = merge(xy_df, com_dt, by = id_var)
    p_dt$cluster_id = factor(p_dt$cluster_id)

    p_dt
}

.groupRegionByDimReduceCluster = function(ct2, group_VAR = "knn_id", nearest_neighbors = 100){
    args = get_args()
    xy_df = GenomicRanges::mcols(rowRanges(ct2)) %>%
        as.data.frame() %>%
        dplyr::select(tx, ty)
    xy_df[[ct2@region_VAR]] = rownames(xy_df)

    knn_res = .knn_clustering(xy_df,
                    nn = nearest_neighbors,
                    id_var = ct2@region_VAR)
    colnames(knn_res)[4] = group_VAR
    knn_res = knn_res[, c(ct2@region_VAR, group_VAR)]

    history_item = list(groupRegionByDimReduceCluster = list(FUN = .groupRegionByDimReduceCluster, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = knn_res,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}


#' @export
setGeneric("groupRegionByDimReduceCluster",
           function(ct2, group_VAR = "knn_id", nearest_neighbors = 100)
               standardGeneric("groupRegionByDimReduceCluster"),
           signature = "ct2")

#' @export
setMethod("groupRegionByDimReduceCluster", c("ChIPtsne2"), .groupRegionByDimReduceCluster)

#### memb table grouping ####

.groupRegionByMembershipTable = function(ct2, membership, group_VAR = "membership_id"){
    args = get_args()
    memb_df = seqsetvis::ssvMakeMembTable(membership)
    group_df = seqsetvis::ssvFactorizeMembTable(memb_df)
    colnames(group_df) = c(ct2@region_VAR, group_VAR)
    history_item = list(groupRegionByMembershipTable = list(FUN = .groupRegionByMembershipTable, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = group_df,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}


#' @export
setGeneric("groupRegionByMembershipTable",
           function(ct2, membership, group_VAR = "membership_id")
               standardGeneric("groupRegionByMembershipTable"),
           signature = "ct2")

#' @export
setMethod("groupRegionByMembershipTable", c("ChIPtsne2"), .groupRegionByMembershipTable)

#### manual grouping ####
.groupRegionManually = function(ct2, assignment, group_VAR = "group_id"){
    args = get_args()
    if(is.data.frame(assignment)){
        if(!ct2@region_VAR %in% colnames(assignment)){
            if(is.null(rownames(assignment))){
                stop("assignment must have rownames or region_VAR:", ct2@region_VAR, " in colnames")
            }else{
                assignment[[ct2@region_VAR]] = rownames(assignment)
            }
        }
        if(!group_VAR %in% colnames(assignment)){
            stop("group_VAR:", group_VAR, " must be in colnames of assignment.")
        }
        assignment = assignment[, c(ct2@region_VAR, group_VAR)]
    }else{
        if(is.null(names(assignment))){
            stop("assignment must be a data.frame or named vector.")
        }else{
            df = data.frame(names(assignment), assignment)
            colnames(df) = c(ct2@region_VAR, group_VAR)
            assignment = df
        }
    }
    history_item = list(groupRegionManually = list(FUN = .groupRegionManually, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = assignment,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("groupRegionManually",
           function(ct2, assignment, group_VAR = "group_id")
               standardGeneric("groupRegionManually"),
           signature = "ct2")

#' @export
setMethod("groupRegionManually", c("ChIPtsne2"), .groupRegionManually)
