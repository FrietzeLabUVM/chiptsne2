# debug(plotSignalHeatmap)
clust_cols = seqsetvis::safeBrew(getRegionMetaData(ct2)$cluster)
bin_cols = c("FALSE" = "white", "TRUE" = "lightblue")
olap_cols = seqsetvis::safeBrew(meta_df$overlap, pal = "set1")

plotSignalHeatmap(
    ct2,
    group_VARS = c(
        # "cluster",
        "peak_MCF10A_CTCF",
        "peak_MCF10AT1_CTCF",
        "peak_MCF10CA1_CTCF",
        "overlap",
        "cluster"
    ),
    # balance_VAR = "cluster",
    annotation_colors = list(
     # cluster = clust_cols,
     peak_MCF10A_CTCF = bin_cols,
     peak_MCF10AT1_CTCF = bin_cols,
     peak_MCF10CA1_CTCF = bin_cols,
     overlap = olap_cols
    ),
    annotation_legends_to_show = c(1, 4),
    # n_legend_rows = 2,
    relative_heatmap_height = .7,
    relative_heatmap_width = .7
)

####
ct2 = exampleChIPtsne2.with_meta()
ct2 = groupRegionsBySignalCluster(ct2, group_VAR = "cluster")
ct2 = groupRegionsByOverlap(ct2, seqsetvis::CTCF_in_10a_narrowPeak_grs[1:2], group_VAR = "overlap")

meta_df = getRegionMetaData(ct2)
meta_df = meta_df %>% dplyr::mutate(overlap_num = as.numeric(overlap))
meta_df$chr_num = as.numeric(GenomicRanges::seqnames(rowRanges(ct2)))
ct2 = setRegionMetaData(ct2, meta_df)

plotSignalHeatmap(ct2)

plotSignalHeatmap(ct2, group_VARS = c("cluster", "overlap", "chr_num"), sort_VAR = "cluster", annotation_colors = list("chr_num" = c("blue", "red")))
