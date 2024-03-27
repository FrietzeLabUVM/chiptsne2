#### signal clustering ####

.groupRegionsBySignalCluster = function(ct2, group_VAR = "cluster_id", n_clusters = 6, iter.max = 30){
    message("groupRegionsBySignalCluster ...")
    args = get_args()
    # prof_dt = getTidyProfile(ct2)
    # clust_dt = suppressMessages({
    #     seqsetvis::ssvSignalClustering(
    #         prof_dt,
    #         cluster_ = group_VAR,
    #         nclust = n_clusters,
    #         fill_ =ct2@value_VAR,
    #         facet_ = ct2@name_VAR,
    #         column_ = ct2@position_VAR,
    #         row_ = ct2@region_VAR,
    #         max_rows = Inf,
    #         max_cols = Inf)
    # })
    # assign_dt = clust_dt %>%
    #     dplyr::select(dplyr::all_of(c(group_VAR, ct2@region_VAR))) %>%
    #     unique
    assign_dt = seqsetvis::clusteringKmeansNestedHclust(ct2@rowToRowMat,
                                             nclust = n_clusters,
                                             iter.max = iter.max)
    data.table::setnames(assign_dt, c(ct2@region_VAR, group_VAR))


    history_item = list(groupRegionsBySignalCluster = list(FUN = .groupRegionsBySignalCluster, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        new_region_metadata = assign_dt,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}


#' groupRegionsBySignalCluster
#'
#' @param ct2 `r doc_ct2()`
#' @param group_VAR `r doc_group_VAR()`
#' @param n_clusters Number of clusters specified for k-means.
#'
#' @return `r doc_return_group()`
#' @export
#' @rdname ct2-group-regions-signal
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = groupRegionsBySignalCluster(ct2)
#' ct2 = groupRegionsBySignalCluster(ct2)
setGeneric("groupRegionsBySignalCluster",
           function(ct2, group_VAR = "cluster_id", n_clusters = 6, iter.max = 30)
               standardGeneric("groupRegionsBySignalCluster"),
           signature = "ct2")

#' @export
#' @rdname ct2-group-regions-signal
setMethod("groupRegionsBySignalCluster", c("ChIPtsne2_no_rowRanges"), .groupRegionsBySignalCluster)



#### dim clustering ####
#' .knn_clustering
#'
#' @param xy_df data.frame contain "tx", "ty", and id_var
#' @param nn number of nearest neighbors, passed as k to RANN::nn2
#' @param id_var must be in xy_df, default is "id"
#' @param ... passed to RANN::nn2
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
#' chiptsne2:::.knn_clustering(xy_df, 20)
.knn_clustering = function(xy_df,
                           nn = 100,
                           id_var = "id",
                           ...){
    valid_vars = c(id_var, "tx", "ty")
    stopifnot(valid_vars %in% colnames(xy_df))

    if(nn > .2 * nrow(xy_df)){
        nn = floor(.2 * nrow(xy_df))
        message("Decreasing nearest-neighbors to ", nn, ".  Original value was too high for dataset.")
    }
    mat = t(as.matrix(xy_df[, c("tx", "ty")]))
    colnames(mat) = xy_df[[id_var]]
    knn.info <- RANN::nn2(t(mat), k = nn, ...)
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

.groupRegionsByDimReduceCluster = function(ct2, group_VAR = "knn_id", nearest_neighbors = 100, ...){
    #visible binding NOTE
    tx = ty = NULL
    message("groupRegionsByDimReduceCluster ...")
    if(!hasDimReduce(ct2)){
        stop("No dimensional reduction data present in this ChIPtsne2 object. Run dimReduceTSNE/PCA/UMAP first then try again.")
    }
    args = get_args()
    xy_df = rowData(ct2) %>%
        as.data.frame() %>%
        dplyr::select(tx, ty)
    xy_df[[ct2@region_VAR]] = rownames(xy_df)

    knn_res = .knn_clustering(xy_df,
                              nn = nearest_neighbors,
                              id_var = ct2@region_VAR,
                              ...)
    colnames(knn_res)[4] = group_VAR
    knn_res = knn_res[, c(ct2@region_VAR, group_VAR)]

    history_item = list(groupRegionsByDimReduceCluster = list(FUN = .groupRegionsByDimReduceCluster, ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2,
        new_region_metadata = knn_res,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}



#' groupRegionsByDimReduceCluster
#'
#' @param ct2 `r doc_ct2()`
#' @param group_VAR `r doc_group_VAR()`
#' @param nearest_neighbors The number of nearest neighbors to use when clustering. Higher numbers result in fewer clusters. Will be automatically reduced if set too high. See documentation of k in [RANN::nn2()].
#' @param ... passed to RANN::nn2
#'
#' @return `r doc_return_group()`
#' @export
#' @rdname ct2-group-regions-dr
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' ct2 = dimReduceTSNE(ct2, perplexity = 20)
#' ct2 = groupRegionsByDimReduceCluster(ct2, group_VAR = "knn_id", nearest_neighbors = 20)
#' rowData(ct2)
#' plotDimReducePoints(ct2, color_VAR = "knn_id")
setGeneric("groupRegionsByDimReduceCluster",
           function(ct2, group_VAR = "knn_id", nearest_neighbors = 100, ...)
               standardGeneric("groupRegionsByDimReduceCluster"),
           signature = "ct2")

#' @export
#' @rdname ct2-group-regions-dr
setMethod("groupRegionsByDimReduceCluster", c("ChIPtsne2_no_rowRanges"), .groupRegionsByDimReduceCluster)

#### region overlap ####


.groupRegionsByOverlap = function(ct2, query, group_VAR = "overlap_id", use_priority = FALSE, ...){
    message("groupRegionsByOverlap ...")
    args = get_args()
    gr = rowRanges(ct2)
    GenomicRanges::mcols(gr) = NULL
    if(is(query, "GRangesList")){
        query = as.list(query)
    }
    if(!is.list(query)){#GRanges input
        if(!is(query, "GRanges")){
            stop("query must be a named list or a GRanges object.")
        }
        if(use_priority){
            stop("use_priority = TRUE may not be used (does not make sense) with non-list input.")
        }
        to_overlap = list(query)
        names(to_overlap) = group_VAR
        memb_gr = seqsetvis::ssvOverlapIntervalSets(c(list(TMP__ = gr), to_overlap), use_first = TRUE, ...)
        names(memb_gr) = names(gr)
        memb_df = data.frame(GenomicRanges::mcols(memb_gr), check.names = FALSE)
        group_df = memb_df[, group_VAR, drop = FALSE]
    }else{# list input
        if(is.null(names(query))){
            stop("query must have names.")
        }
        if(any(duplicated(names(query)))){
            stop("All names of query must be unique.")
        }
        if(use_priority){
            group_df = data.frame(id = rownames(ct2), group = "no_hit")
            for(i in rev(seq_along(query))){
                qgr = query[[i]]
                olaps = GenomicRanges::findOverlaps(gr, qgr, ...)
                group_df[S4Vectors::queryHits(olaps),]$group = names(query)[i]
            }
        }else{
            memb_gr = seqsetvis::ssvOverlapIntervalSets(c(list(TMP__ = gr), as.list(query)), use_first = TRUE, ...)
            names(memb_gr) = names(gr)
            memb_df = data.frame(GenomicRanges::mcols(memb_gr), check.names = FALSE)
            group_df = seqsetvis::ssvFactorizeMembTable(memb_df[, names(query), drop = FALSE])

        }
        rownames(group_df) = group_df$id
        group_df$id = NULL
        colnames(group_df) = group_VAR
    }
    ct2 = setRegionMetaData(ct2, group_df, silent = TRUE)
    ct2 = .add_history_entry(ct2, "groupRegionsByMembershipTable", FUN = .groupRegionsByOverlap, ARG = args)
    ct2
}

#' groupRegionsByOverlap
#'
#' @param ct2 `r doc_ct2()`
#' @param query Either a single GRanges or named list of GRanges.
#' @param group_VAR `r doc_group_VAR()`
#' @param ... arguments passed to [IRanges::findOverlaps()], i.e. maxgap, minoverlap, type, select, invert.
#'
#' @return `r doc_return_group()`
#' @export
#' @rdname ct2-group-regions-olap
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' peak_grs = seqsetvis::CTCF_in_10a_narrowPeak_grs
#' ct2.olap = ct2
#' ct2.olap = groupRegionsByOverlap(
#'   ct2.olap,
#'   peak_grs[1:2],
#'   group_VAR = "10A_AT1_combos"
#' )
#' ct2.olap = groupRegionsByOverlap(
#'   ct2.olap,
#'   peak_grs$MCF10A_CTCF,
#'   group_VAR = "10A_peaks"
#' )
#' ct2.olap = groupRegionsByOverlap(
#'   ct2.olap,
#'   peak_grs$MCF10AT1_CTCF,
#'   group_VAR = "AT1_peaks"
#' )
#' ct2.olap = groupRegionsByOverlap(
#'   ct2.olap,
#'   peak_grs[1:2],
#'   group_VAR = "10A_AT1_priority",
#'   use_priority = TRUE
#' )
#' rowData(ct2.olap)
#'
#' # a GRangesList is allowed too
#' gr_ls = GenomicRanges::GRangesList(peak_grs)
#' ct2.olap = groupRegionsByOverlap(ct2.olap, gr_ls, group_VAR = "gr_ls_combos")
setGeneric("groupRegionsByOverlap",
           function(ct2, query, group_VAR = "overlap_id", use_priority = FALSE, ...)
               standardGeneric("groupRegionsByOverlap"),
           signature = "ct2")

#' @export
#' @rdname ct2-group-regions-olap
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
        new_region_metadata = group_df,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}


#' groupRegionsByMembershipTable
#'
#' Adds a single grouping variable to rowData/region metadata based on a region overlap membership table.
#'
#' This membership table can be one of 2 different main formats.
#'
#' 1) a data.frame as returned from [seqsetvis::ssvMakeMembTable], where rownames are the same as rownames in `ct2` and columns are all logical.
#'
#' 2) a GenomicRanges object as returned from [seqsetvis::ssvOverlapIntervalSets] or [seqsetvis::ssvConsensusIntervalSets], where names are the same as rownames in `ct2` and mcols are all logical.
#'
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param membership Either a data.frame or GenomicRanges. See details.
#' @param group_VAR `r doc_group_VAR()`
#'
#' @return `r doc_return_group()`
#' @rdname ct2-group-regions-memb
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' rowData(ct2)
#' np_grs = seqsetvis::CTCF_in_10a_narrowPeak_grs
#'
#' # membership data.frame
#' memb_df = seqsetvis::ssvMakeMembTable(np_grs)
#' ct2.with_memb1 = chiptsne2:::.groupRegionsByMembershipTable(ct2, memb_df)
#' rowData(ct2.with_memb1)
#'
#' # membership GRanges
#' memb_grs = seqsetvis::ssvOverlapIntervalSets(np_grs)
#' ct2.with_memb2 = chiptsne2:::.groupRegionsByMembershipTable(ct2, memb_grs)
#' rowData(ct2.with_memb2)
#'
#' # rowData of example ct2 is actually a membership table
#' ct2.with_memb3 = chiptsne2:::.groupRegionsByMembershipTable(ct2, rowData(ct2))
#' rowData(ct2.with_memb3)
setGeneric("groupRegionsByMembershipTable",
           function(ct2, membership, group_VAR = "membership_id")
               standardGeneric("groupRegionsByMembershipTable"),
           signature = "ct2")

#' @export
#' @rdname ct2-group-regions-memb
setMethod("groupRegionsByMembershipTable", c("ChIPtsne2_no_rowRanges"), .groupRegionsByMembershipTable)

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
        if(!all(group_VAR %in% colnames(assignment))){
            stop("group_VAR:", group_VAR, " must be in colnames of assignment.")
        }
        assignment = assignment[, c(ct2@region_VAR, group_VAR)]
    }else{
        if(length(group_VAR) > 1){
            stop("When assignment is not a data.frame only 1 group_VAR is allowed.")
        }
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
        new_region_metadata = assignment,
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' groupRegionsManually
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param assignment A data.frame or named character vector. See details.
#' @param group_VAR `r doc_group_VAR()`
#'
#' Assignment may be either a data.frame or named character vector.
#'
#' If data.frame, must have rownames setequal to rownames of `ct2` or contain the same region variable as `ct2`. Multiple `group_VAR` may be specified and all must be present in data.frame.
#'
#' If a named character vector, names must be setequal to rownames of `ct2`. Only a signle `group_VAR` is allowed.
#'
#' @return `r doc_return_group()`
#' @rdname ct2-group-regions-manual
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#'
#' # data.frame assignment
#' np_grs = seqsetvis::CTCF_in_10a_narrowPeak_grs
#' memb_grs = seqsetvis::ssvOverlapIntervalSets(np_grs)
#' assign_df = seqsetvis::ssvFactorizeMembTable(memb_grs)
#' assign_df$not_used = "A"
#' assign_df$group2 = "B"
#' ct2.grouped = groupRegionsManually(ct2, assign_df, group_VAR = c("group", "group2"))
#' rowData(ct2.grouped)
#'
#' # character vector assignment
#' manual_groups = assign_df$group
#' names(manual_groups) = assign_df$id
#' ct2.grouped2 = groupRegionsManually(ct2, manual_groups, group_VAR = c("group_by_vector"))
#' rowData(ct2.grouped2)
setGeneric("groupRegionsManually",
           function(ct2, assignment, group_VAR = "group_id")
               standardGeneric("groupRegionsManually"),
           signature = "ct2")

#' @export
#' @rdname ct2-group-regions-manual
setMethod("groupRegionsManually", c("ChIPtsne2_no_rowRanges"), .groupRegionsManually)

#### sort regions ####

.sortRegions = function(ct2, group_VAR = NULL){
    message("sortRegions ...")
    args = get_args()
    prof_dt = getTidyProfile(ct2, meta_VARS = group_VAR)
    if(length(group_VAR) > 1){
        prof_dt$TMP_GROUP__ = apply(prof_dt[, group_VAR, with = FALSE], 1, paste, collapse = "_")
        group_VAR = "TMP_GROUP__"
    }
    clust_dt = seqsetvis::within_clust_sort(
        prof_dt,
        row_ = ct2@region_VAR,
        column_ = ct2@position_VAR,
        fill_ = ct2@value_VAR,
        facet_ = ct2@name_VAR,
        cluster_ = group_VAR,
        dcast_fill = 0)
    region_lev = levels(clust_dt[[ct2@region_VAR]])
    history_item = list(sortRegions  = list(FUN = .sortRegions , ARG = args))
    cloneChIPtsne2_fromTidy(
        ct2 = ct2[region_lev,],
        new_obj_history = c(ChIPtsne2.history(ct2), history_item)
    )
}

#' sortRegions
#'
#' @param ct2 `r doc_ct2_nrr()`
#' @param group_VAR `r doc_group_VAR()`
#'
#' @return `r doc_ct2_nrr()` with rows sorted by `group_VAR` and within groups by signal decreasing.
#' @rdname ct2-sort
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' rowData(ct2)
#' ct2.sorted = sortRegions(ct2, group_VAR = "peak_MCF10CA1_CTCF")
#' rowData(ct2.sorted)
#'
#' ct2.sorted2 = sortRegions(ct2, group_VAR = c("peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))
#' rowData(ct2.sorted2)
setGeneric("sortRegions",
           function(ct2, assignment, group_VAR = NULL)
               standardGeneric("sortRegions"),
           signature = "ct2")

#' @export
#' @rdname ct2-sort
setMethod("sortRegions", c("ChIPtsne2_no_rowRanges"), .sortRegions)
