#### tSNE ####

#' @export
setGeneric("dimReduceTSNE", function(ct2, perplexity = 50, ...) standardGeneric("dimReduceTSNE"), signature = "ct2")

#' @export
setMethod("dimReduceTSNE", c("ChIPtsne2"), .tsne_from_ct2)

#' .tsne_from_profile_mat
#'
#' @param prof_mat matrix of profile data
#' @param perplexity perplexity, will be reduced to 1/4 of nrow of prof_mat if too large.
#' @param ... passed to Rtsne::Rtsne
#'
#' @return data.frame of reduced dimensions by tSNE
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' prof_mat = rowToRowMat(ct2)
#' ct2.tsne = dimReduceTSNE(ct2, perplexity = 10)
#' rowRanges(ct2.tsne)
.tsne_from_profile_mat = function(prof_mat, perplexity = 50, ...){
    if(perplexity > nrow(prof_mat)/4){
        perplexity = nrow(prof_mat)/4
        warning("autoreduced perplexity to ", perplexity, " due to small number of regions.")
    }
    tsne_res = Rtsne::Rtsne(prof_mat, check_duplicates = FALSE, perplexity = perplexity, ...)
    xy_df = as.data.frame(tsne_res$Y)
    colnames(xy_df) =  c("tx", "ty")
    rownames(xy_df) = rownames(prof_mat)
    .post_dim_reduce(xy_df[])
}

.tsne_from_ct2 = function(ct2, perplexity = 50, ...){
    args = get_args()
    xy_df = .tsne_from_profile_mat(rowToRowMat(ct2), perplexity = perplexity, ...)
    xy_df[[ct2@region_VAR]] = rownames(xy_df)

    history_item = list(centerSignalProfile = list(FUN = .tsne_from_ct2, ARG = args))
    ChIPtsne2.from_tidy(getTidyProfile(ct2),
                        rowRanges(ct2),
                        region_metadata = xy_df,
                        sample_metadata = colData(ct2),
                        name_VAR = ct2@name_VAR,
                        position_VAR = ct2@position_VAR,
                        value_VAR = ct2@value_VAR,
                        region_VAR = ct2@region_VAR,
                        obj_history = c(ChIPtsne2.history(ct2), history_item),
                        init = FALSE
    )
}

#### UMAP ####
#' @export
setGeneric("dimReduceUMAP", function(ct2, config = umap::umap.defaults, ...) standardGeneric("dimReduceUMAP"), signature = "ct2")

#' @export
setMethod("dimReduceUMAP", c("ChIPtsne2"), .umap_from_ct2)

#' .umap_from_profile_mat
#'
#' @param prof_mat matrix of profile data
#' @param config umap.config object, uses umap::umap.defaults by default
#' @param ... passed to umap::umap
#'
#' @return data.frame of reduced dimensions by UMAP
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' prof_mat = rowToRowMat(ct2)
#' ct2.umap = dimReduceUMAP(ct2)
.umap_from_profile_mat = function(prof_mat, config = umap::umap.defaults, ...){
    tsne_res = umap::umap(prof_mat, config = config, ...)
    xy_df = as.data.frame(tsne_res$layout)
    colnames(xy_df) = c("tx", "ty")
    rownames(xy_df) = rownames(prof_mat)
    .post_dim_reduce(xy_df[])
}

.umap_from_ct2 = function(ct2, config = umap::umap.defaults, ...){
    args = get_args()
    xy_df = .umap_from_profile_mat(rowToRowMat(ct2), config = config, ...)
    xy_df[[ct2@region_VAR]] = rownames(xy_df)

    history_item = list(centerSignalProfile = list(FUN = .umap_from_ct2, ARG = args))
    ChIPtsne2.from_tidy(getTidyProfile(ct2),
                        rowRanges(ct2),
                        region_metadata = xy_df,
                        sample_metadata = colData(ct2),
                        name_VAR = ct2@name_VAR,
                        position_VAR = ct2@position_VAR,
                        value_VAR = ct2@value_VAR,
                        region_VAR = ct2@region_VAR,
                        obj_history = c(ChIPtsne2.history(ct2), history_item),
                        init = FALSE
    )
}

#### PCA ####
#' @export
setGeneric("dimReducePCA", function(ct2) standardGeneric("dimReducePCA"), signature = "ct2")

#' @export
setMethod("dimReducePCA", c("ChIPtsne2"), .pca_from_ct2)

#' .pca_from_profile_mat
#'
#' @param prof_mat matrix of profile data
#'
#' @return data.frame of reduced dimensions by PCA
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' dimReducePCA(ct2)
.pca_from_profile_mat = function(prof_mat){
    pca_res = prcomp(prof_mat)
    xy_df = as.data.frame(pca_res$x[, c(1, 2)])
    colnames(xy_df) = c("tx", "ty")
    rownames(xy_df) = rownames(prof_mat)
    .post_dim_reduce(xy_df[])
}

.pca_from_ct2 = function(ct2){
    args = get_args()
    xy_df = .pca_from_profile_mat(rowToRowMat(ct2))
    xy_df[[ct2@region_VAR]] = rownames(xy_df)

    history_item = list(centerSignalProfile = list(FUN = .pca_from_ct2, ARG = args))
    ChIPtsne2.from_tidy(getTidyProfile(ct2),
                        rowRanges(ct2),
                        region_metadata = xy_df,
                        sample_metadata = colData(ct2),
                        name_VAR = ct2@name_VAR,
                        position_VAR = ct2@position_VAR,
                        value_VAR = ct2@value_VAR,
                        region_VAR = ct2@region_VAR,
                        obj_history = c(ChIPtsne2.history(ct2), history_item),
                        init = FALSE
    )
}

#### helpers ####

.rescale_capped = function(x, to = c(0,1), from = range(x, na.rm = TRUE, finite = TRUE)){
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}

.post_dim_reduce = function(xy_df, norm1 = TRUE, high_topright = TRUE){
    if (norm1) {
        xy_df$tx = .rescale_capped(xy_df$tx) - 0.5
        xy_df$ty = .rescale_capped(xy_df$ty) - 0.5
    }
    if (high_topright) {#flip tx/ty if needed so that
        rs = rowSums(prof_mat)
        xy_df$rs = rs[rownames(xy_df)]
        x_cutoff = mean(range(xy_df$tx))
        x_flip = sum(xy_df[xy_df$tx > x_cutoff,]$rs) < sum(xy_df[xy_df$tx < x_cutoff,]$rs)
        if (x_flip) {
            xy_df = xy_df %>% dplyr::mutate(tx = max(tx) - tx + min(tx))
        }
        y_cutoff = mean(range(xy_df$ty))
        right_sum = sum((xy_df %>% dplyr::filter(ty > y_cutoff))$rs)
        left_sum = sum((xy_df %>% dplyr::filter(ty <= y_cutoff))$rs)
        y_flip = right_sum < left_sum
        if (y_flip) {
            xy_df = xy_df %>% dplyr::mutate(ty = max(ty) - ty + min(ty))
        }
        xy_df$rs = NULL
    }
    xy_df[]
}
