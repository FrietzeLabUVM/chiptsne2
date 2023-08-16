ct2 = exampleChIPtsne2.with_meta()
ct2 = dimReduceTSNE(ct2)
hasDimReduce = chiptsne2:::hasDimReduce
.plotDimReduceBins(ct2)
# meta_dt = getSampleMetaData(ct2)
# meta_dt = meta_dt %>% separate(sample, c("cell", "mark"), "_", remove = FALSE)
# ct2 = setSampleMetaData(ct2, meta_dt)
.plotDimReduceBins(ct2, facet_columns = "sample", facet_rows = NULL)
colData(ct2)
debug(.plotDimReduceBins)
.plotDimReduceBins(ct2, facet_columns = "cell", facet_rows = "mark", xbins = 8, ybins = 8, min_size = 1) +
    coord_fixed()
