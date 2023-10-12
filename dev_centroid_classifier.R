library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()
ct2 = groupRegionsBySignalCluster(ct2)
ct2.sp = chiptsne2::split(ct2, ct2$cell)
ct2.sp$MCF10A - ct2.sp$MCF10AT1


meta_df = getRegionMetaData(ct2)
head(meta_df)
meta_df$peak_MCF10A_CTCF
rowToRowMat(ct2)

cent = calculate_centroid_per_group(ct2, "cluster_id")
cent_dist = calculate_distance_to_centroids(ct2, cent)
classify_by_centroid_distances(cent_dist, cent)



centroid = calculateGroupCentroid(ct2, "cluster_id")

group_VAR = "cent_class"


ct2 = groupRegionsByCentroidDistance(ct2, centroid, "new_group")
ct2 = dimReduceTSNE(ct2)
ct2 = groupRegionsByCentroidDistance(ct2, centroid, "new_group", tolerance = 10)
plotDimReducePoints(ct2, color_VAR = "new_group", point_colors = c("ambiguous" = "gray"))

ct2.by_cell = split(ct2, "cell")
cbind(ct2.by_cell$MCF10A, ct2.by_cell$MCF10AT1)
