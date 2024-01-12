library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()
ct2 = groupRegionsBySignalCluster(ct2)
ct2.sp = chiptsne2::split(ct2, ct2$cell)
ct2.sp$MCF10A - ct2.sp$MCF10AT1


meta_df = getRegionMetaData(ct2)
head(meta_df)
meta_df$peak_MCF10A_CTCF
rowToRowMat(ct2)

# cent = calculate_centroid_per_group(ct2, "cluster_id")
# cent_dist = calculate_distance_to_centroids(ct2, cent)
# classify_by_centroid_distances(cent_dist, cent)

rowRanges(ct2)
calculateGroupCentroid(ct2, c("peak_MCF10CA1_CTCF", "cluster_id"))
calculateGroupCentroid(ct2, c("peak_MCF10CA1_CTCF", "cluster_id"))

group_VAR = c("peak_MCF10CA1_CTCF", "cluster_id")
#
# aggregateRegionsByGroup = function(ct2, group_VAR, new_meta_VAR = ifelse(length(group_VAR) == 1, group_VAR, "meta_id")){
#     centroid = calculateGroupCentroid(ct2, group_VAR)
#
#     df = do.call(rbind,
#                  lapply(names(ct2@colToRowMatCols), function(nam){
#                      i = ct2@colToRowMatCols[[nam]]
#                      df = reshape2::melt(centroid[, i])
#                      df$Var2 = factor(df$Var2)
#                      levels(df$Var2) = as.numeric(sub(paste0(nam, "_"), "", levels(df$Var2), fixed = TRUE))
#                      df[[ct2@name_VAR]] = nam
#                      df
#                  })
#     )
#     colnames(df) = c(new_meta_VAR, ct2@position_VAR, ct2@value_VAR, ct2@name_VAR)
#     df[[new_meta_VAR]] = factor(df[[new_meta_VAR]], levels = rownames(centroid))
#
#     # rd_sp = split(rowData(ct2), rowData(ct2)[["TMP_GROUP__"]])
#     rd = unique(rowData(ct2)[, group_VAR])
#     rownames(rd) = apply(rd, 1, paste, collapse = ",")
#     rd[[new_meta_VAR]] = rownames(rd)
#
#     ct2.meta = ChIPtsne2.from_tidy(
#         df,
#         query_gr = NULL,
#         sample_metadata = colData(ct2),
#         region_metadata = rd,
#         position_VAR = ct2@position_VAR,
#         name_VAR = ct2@name_VAR,
#         value_VAR = ct2@value_VAR,
#         region_VAR = new_meta_VAR)
#     ct2.meta
# }
# debug(aggregateRegionsByGroup)

function(ct2, group_VAR, new_meta_VAR = ifelse(length(group_VAR) == 1, group_VAR, "meta_id")){
    centroid = calculateGroupCentroid(ct2, group_VAR)

    df = do.call(rbind,
                 lapply(names(ct2@colToRowMatCols), function(nam){
                     i = ct2@colToRowMatCols[[nam]]
                     df = reshape2::melt(centroid[, i])
                     df$Var2 = factor(df$Var2)
                     levels(df$Var2) = as.numeric(sub(paste0(nam, "_"), "", levels(df$Var2), fixed = TRUE))
                     df[[ct2@name_VAR]] = nam
                     df
                 })
    )
    colnames(df) = c(new_meta_VAR, ct2@position_VAR, ct2@value_VAR, ct2@name_VAR)
    df[[new_meta_VAR]] = factor(df[[new_meta_VAR]], levels = rownames(centroid))

    # rd_sp = split(rowData(ct2), rowData(ct2)[["TMP_GROUP__"]])
    rd = unique(rowData(ct2)[, group_VAR, drop = FALSE])
    rownames(rd) = apply(rd, 1, paste, collapse = ",")
    rd[[new_meta_VAR]] = rownames(rd)

    ct2.meta = ChIPtsne2.from_tidy(
        df,
        query_gr = NULL,
        sample_metadata = colData(ct2),
        region_metadata = rd,
        position_VAR = ct2@position_VAR,
        name_VAR = ct2@name_VAR,
        value_VAR = ct2@value_VAR,
        region_VAR = new_meta_VAR)
    ct2.meta
}

debug(aggregateRegionsByGroup)
ct2.meta = aggregateRegionsByGroup(ct2, "cluster_id")

plotSignalLinePlot(ct2)
cowplot::plot_grid(


    plotSignalLinePlot(ct2, group_VAR = "cluster_id"),
    plotSignalLinePlot(ct2.meta)
)


getSampleMetaData(ct2)
getRegionMetaData(ct2)

getSampleMetaData(ct2.meta)
getRegionMetaData(ct2.meta)

getTidyProfile(ct2.meta)
getTidyProfile(ct2)


rowData(ct2.meta)
plotSignalLinePlot(ct2)
plotSignalLinePlot(ct2.meta)
rownames(ct2.meta)
df
group_VAR = "cent_class"

centroid = calculateGroupCentroid(ct2, "cluster_id")


ct2 = dimReduceTSNE(ct2)

ct2 = groupRegionsByCentroidDistance(ct2, centroid, "new_group")
plotDimReducePoints(ct2, color_VAR = "new_group", point_colors = c("ambiguous" = "gray"))

ct2 = groupRegionsByCentroidDistance(ct2, centroid, "new_group", tolerance = 10)
plotDimReducePoints(ct2, color_VAR = "new_group", point_colors = c("ambiguous" = "gray"))

ct2.by_cell = split(ct2, "cell")
cbind(ct2.by_cell$MCF10A, ct2.by_cell$MCF10AT1)
