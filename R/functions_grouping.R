#### signal clustering ####
#' groupRegionsBySignalCluster
#'
#' @param ct2
#' @param group_VAR
#' @param n_clusters
#'
#' @return
#' @rdname groupRegionsBySignalCluster
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = groupRegionsBySignalCluster(ct2)
#' ct2 = groupRegionsBySignalCluster(ct2)
.groupRegionsBySignalCluster = function(ct2, group_VAR = "cluster_id", n_clusters = 6){
    message("groupRegionsBySignalCluster ...")
    args = get_args()
    prof_dt = getTidyProfile(ct2)
    clust_dt = suppressMessages({
        seqsetvis::ssvSignalClustering(
            prof_dt,
            cluster_ = group_VAR,
            nclust = n_clusters,
            fill_ =ct2@value_VAR,
            facet_ = ct2@name_VAR,
            column_ = ct2@position_VAR,
            row_ = ct2@region_VAR,
            max_rows = Inf,
            max_cols = Inf)
    })
    assign_dt = clust_dt %>%
        dplyr::select(all_of(c(group_VAR, ct2@region_VAR))) %>%
        unique

    history_item = list(groupRegionsBySignalCluster = list(FUN = .groupRegionsBySignalCluster, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = assign_dt,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("groupRegionsBySignalCluster",
           function(ct2, group_VAR = "cluster_id", n_clusters = 6)
               standardGeneric("groupRegionsBySignalCluster"),
           signature = "ct2")

#' @export
setMethod("groupRegionsBySignalCluster", c("ChIPtsne2"), .groupRegionsBySignalCluster)



#### dim clustering ####
#' .knn_clustering
#'
#' @param xy_df data.frame contain "tx", "ty", and id_var
#' @param nn number of nearest neighbors, passed as k to RANN::nn2
#' @param id_var must be in xy_df, default is "id"
#'
#' @return cluster assignment table for xy_df based on tx and ty coordinates
#' @importFrom Matrix Matrix
#' @importFrom igraph graph.adjacency simplify cluster_walktrap
#' @importFrom RANN nn2
#' @importFrom reshape2 melt
#'
#' @examples
#' xy_df = data.frame(tx = runif(100), ty = runif(100))
#' xy_df$id = paste0("id_", seq_len(nrow(xy_df)))
#' .knn_clustering(xy_df, 20)
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

.groupRegionsByDimReduceCluster = function(ct2, group_VAR = "knn_id", nearest_neighbors = 100){
    message("groupRegionsByDimReduceCluster ...")
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

    history_item = list(groupRegionsByDimReduceCluster = list(FUN = .groupRegionsByDimReduceCluster, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = knn_res,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}


#' @export
setGeneric("groupRegionsByDimReduceCluster",
           function(ct2, group_VAR = "knn_id", nearest_neighbors = 100)
               standardGeneric("groupRegionsByDimReduceCluster"),
           signature = "ct2")

#' @export
setMethod("groupRegionsByDimReduceCluster", c("ChIPtsne2"), .groupRegionsByDimReduceCluster)

#### region overlap ####

#' Title
#'
#' @param ct2
#' @param gr_list
#' @param group_VAR
#'
#' @return
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' peak_grs = seqsetvis::CTCF_in_10a_narrowPeak_grs
#' ct2.olap = groupRegionsByOverlap(ct2, peak_grs[1:2], group_VAR = "10A_AT1_overlap")
.groupRegionsByOverlap = function(ct2, gr_list, group_VAR = "overlap_id"){
    message("groupRegionsByOverlap ...")
    args = get_args()
    gr = rowRanges(ct2)
    GenomicRanges::mcols(gr) = NULL
    memb_gr = seqsetvis::ssvOverlapIntervalSets(c(list(TMP__ = gr), as.list(gr_list)))
    memb_gr$TMP__ = NULL
    memb_df = as.data.frame(GenomicRanges::mcols(memb_gr))
    group_df = seqsetvis::ssvFactorizeMembTable(memb_df)
    colnames(group_df) = c(ct2@region_VAR, group_VAR)
    history_item = list(groupRegionsByMembershipTable = list(FUN = .groupRegionsByMembershipTable, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = group_df,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}


#' @export
setGeneric("groupRegionsByOverlap",
           function(ct2, gr_list, group_VAR = "overlap_id")
               standardGeneric("groupRegionsByOverlap"),
           signature = "ct2")

#' @export
setMethod("groupRegionsByOverlap", c("ChIPtsne2"), .groupRegionsByOverlap)

#### memb table grouping ####

.groupRegionsByMembershipTable = function(ct2, membership, group_VAR = "membership_id"){
    message("groupRegionsByMembershipTable ...")
    args = get_args()
    memb_df = seqsetvis::ssvMakeMembTable(membership)
    group_df = seqsetvis::ssvFactorizeMembTable(memb_df)
    colnames(group_df) = c(ct2@region_VAR, group_VAR)
    history_item = list(groupRegionsByMembershipTable = list(FUN = .groupRegionsByMembershipTable, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = group_df,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}


#' @export
setGeneric("groupRegionsByMembershipTable",
           function(ct2, membership, group_VAR = "membership_id")
               standardGeneric("groupRegionsByMembershipTable"),
           signature = "ct2")

#' @export
setMethod("groupRegionsByMembershipTable", c("ChIPtsne2"), .groupRegionsByMembershipTable)

#### manual grouping ####
.groupRegionsManually = function(ct2, assignment, group_VAR = "group_id"){
    message("groupRegionsManually ...")
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
    history_item = list(groupRegionsManually = list(FUN = .groupRegionsManually, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        region_metadata = assignment,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("groupRegionsManually",
           function(ct2, assignment, group_VAR = "group_id")
               standardGeneric("groupRegionsManually"),
           signature = "ct2")

#' @export
setMethod("groupRegionsManually", c("ChIPtsne2"), .groupRegionsManually)

#### sort regions ####
.sortRegions = function(ct2, group_VAR = NULL){
    message("sortRegions ...")
    args = get_args()
    prof_dt = getTidyProfile(ct2, meta_VARS = group_VAR)
    clust_dt = seqsetvis::within_clust_sort(
        prof_dt,
        row_ = ct2@region_VAR,
        column_ = ct2@position_VAR,
        fill_ = ct2@value_VAR,
        facet_ = ct2@name_VAR,
        cluster_ = group_VAR,
        dcast_fill = 0)
    region_lev = clust_dt[[ct2@region_VAR]]

    new_query_gr = rowRanges(ct2)[region_lev]

    history_item = list(sortRegions  = list(FUN = .sortRegions , ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        query_gr = new_query_gr,
        obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' @export
setGeneric("sortRegions",
           function(ct2, assignment, group_VAR = NULL)
               standardGeneric("sortRegions"),
           signature = "ct2")

#' @export
setMethod("sortRegions", c("ChIPtsne2"), .sortRegions)
