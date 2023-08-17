ct2 = exampleChIPtsne2.with_meta()
ct2 = groupRegionsBySignalCluster(ct2)
ct2 = dimReduceTSNE(ct2)
plotDimReduceBins(ct2)
plotDimReduceBins(ct2, facet_columns = "sample", facet_rows = NULL)
colData(ct2)
plotDimReduceBins(ct2, facet_columns = "cell", facet_rows = "mark", xbins = 8, ybins = 8, min_size = 1) +
    coord_fixed()

plotDimReduceBins(ct2, facet_columns = "cell", facet_rows = "mark", xbins = 10, ybins = 10, min_size = 1) +
    coord_fixed() +
    labs(caption = "10x10 bins")

plotDimReducePoints(ct2)

# debug(chiptsne2:::add_cluster_annotation)
# debug(chiptsne2:::.prep_ids)
ct2@rowRanges$peak_MCF10AT1_CTCF = factor(ct2@rowRanges$peak_MCF10AT1_CTCF, levels = c("FALSE", "TRUE"))
ct2@rowRanges$peak_MCF10AT1_CTCF = factor(ct2@rowRanges$peak_MCF10AT1_CTCF, levels = c("TRUE", "FALSE"))

chiptsne2::plotSignalHeatmap(ct2, group_VARS = "peak_MCF10AT1_CTCF", sort_strategy = "right")


chiptsne2::plotSignalHeatmap(ct2, group_VARS = "cluster_id")
chiptsne2::plotSignalLinePlot(ct2, group_VAR = "peak_MCF10AT1_CTCF")
chiptsne2::plotSignalLinePlot(ct2, group_VAR = "cluster_id")

round(nrow(ct2)^.5)

chiptsne2::setRegionMetaData
