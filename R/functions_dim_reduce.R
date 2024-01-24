#### tSNE ####


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
    tsne_res = Rtsne::Rtsne(prof_mat,
                            check_duplicates = FALSE,
                            perplexity = perplexity,
                            num_threads = getOption("mc.cores", 1),
                            ...)
    xy_df = as.data.frame(tsne_res$Y)
    colnames(xy_df) =  c("tx", "ty")
    rownames(xy_df) = rownames(prof_mat)
    ggplot(xy_df, aes(x = tx, y = ty)) + geom_point()
    .post_dim_reduce(xy_df, prof_mat)
}

.dimReduceTSNE = function(ct2, perplexity = 50, ...){
    message("dimReduceTSNE ...")
    args = get_args()
    xy_df = .tsne_from_profile_mat(rowToRowMat(ct2), perplexity = perplexity, ...)
    xy_df[[ct2@region_VAR]] = rownames(xy_df)
    ct2 = setRegionMetaData(ct2, xy_df, silent = TRUE)
    ct2 = .add_history_entry(ct2 = ct2, entry_name = "dimReduceTSNE", FUN = .dimReduceTSNE, ARG = args)
    ct2
}

#' @export
setGeneric("dimReduceTSNE", function(ct2, perplexity = 50, ...) standardGeneric("dimReduceTSNE"), signature = "ct2")

#' @export
setMethod("dimReduceTSNE", c("ChIPtsne2_no_rowRanges"), .dimReduceTSNE)

#### UMAP ####

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
#' rowRanges(ct2.umap)
.umap_from_profile_mat = function(prof_mat, config = umap::umap.defaults, ...){
    tsne_res = umap::umap(prof_mat, config = config, ...)
    xy_df = as.data.frame(tsne_res$layout)
    colnames(xy_df) = c("tx", "ty")
    rownames(xy_df) = rownames(prof_mat)
    .post_dim_reduce(xy_df, prof_mat)
}

.dimReduceUMAP = function(ct2, config = umap::umap.defaults, ...){
    message("dimReduceUMAP ...")
    args = get_args()
    xy_df = .umap_from_profile_mat(rowToRowMat(ct2), config = config, ...)
    xy_df[[ct2@region_VAR]] = rownames(xy_df)
    ct2 = setRegionMetaData(ct2, xy_df, silent = TRUE)
    ct2 = .add_history_entry(ct2 = ct2, entry_name = "dimReduceUMAP", FUN = .dimReduceUMAP, ARG = args)
    ct2
}

#' @export
setGeneric("dimReduceUMAP", function(ct2, config = umap::umap.defaults, ...) standardGeneric("dimReduceUMAP"), signature = "ct2")

#' @export
setMethod("dimReduceUMAP", c("ChIPtsne2_no_rowRanges"), .dimReduceUMAP)

#### PCA ####

#' .pca_from_profile_mat
#'
#' @param prof_mat matrix of profile data
#'
#' @return data.frame of reduced dimensions by PCA
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' ct2.pca = dimReducePCA(ct2)
#' rowRanges(ct2.pca)
.pca_from_profile_mat = function(prof_mat){
    pca_res = prcomp(prof_mat)
    xy_df = as.data.frame(pca_res$x[, c(1, 2)])
    colnames(xy_df) = c("tx", "ty")
    rownames(xy_df) = rownames(prof_mat)
    .post_dim_reduce(xy_df, prof_mat)
}

.dimReducePCA = function(ct2){
    message("dimReducePCA ...")
    args = get_args()
    xy_df = .pca_from_profile_mat(rowToRowMat(ct2))
    xy_df[[ct2@region_VAR]] = rownames(xy_df)
    ct2 = setRegionMetaData(ct2, xy_df, silent = TRUE)
    ct2 = .add_history_entry(ct2 = ct2, entry_name = "dimReducePCA", FUN = .dimReducePCA, ARG = args)
    ct2
}

#' @export
setGeneric("dimReducePCA", function(ct2) standardGeneric("dimReducePCA"), signature = "ct2")

#' @export
setMethod("dimReducePCA", c("ChIPtsne2_no_rowRanges"), .dimReducePCA)
#### helpers ####

.rescale_capped = function(x, to = c(0,1), from = range(x, na.rm = TRUE, finite = TRUE)){
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}

.post_dim_reduce = function(xy_df, prof_mat, norm1 = TRUE, high_topright = TRUE){
    if (norm1) {
        xy_df$tx = .rescale_capped(xy_df$tx) - 0.5
        xy_df$ty = .rescale_capped(xy_df$ty) - 0.5
    }
    if (high_topright) {#flip tx/ty if needed so that strongest signal is at topright of plot
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

#' hasDimReduce
#'
#' @param ct2
#'
#' @return TRUE if input ct2 has coordinates from dimensional reduction data
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' hasDimReduce(ct2)
#' ct2 = dimReduceTSNE(ct2, perplexity = 25)
#' hasDimReduce(ct2)
hasDimReduce = function(ct2){
    meta_dt = getRegionMetaData(ct2)
    all(c("tx", "ty") %in% colnames(meta_dt))
}

#' reportDimReduceMethod
#'
#' @param ct2
#'
#' @return Character describing active dimReduction method.
#' @export
#'
#' @examples
#' ct2 = exampleChIPtsne2()
#' reportDimReduceMethod(ct2)
#' ct2 = dimReduceTSNE(ct2, perplexity = 25)
#' reportDimReduceMethod(ct2)
reportDimReduceMethod = function(ct2){
    obj_history = ChIPtsne2.history(ct2)
    is_dim_red = which(grepl("dimReduce", names(obj_history)))
    if(length(is_dim_red) < 1){
        out = "none"
    }else{
        most_recent = names(obj_history[max(is_dim_red)])
        out = sub("dimReduce", "", most_recent)
        out = sub("TSNE", "t-SNE", out)
    }
    out
}
