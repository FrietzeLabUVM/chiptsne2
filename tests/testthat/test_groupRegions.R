testthat::context("ChIPtsne_class")
# flipping viewGranges
library(chiptsne2)
library(testthat)

ct2 = exampleChIPtsne2.with_meta()
peak_grs = seqsetvis::CTCF_in_10a_narrowPeak_grs
ct2.olap = groupRegionsByOverlap(ct2, peak_grs[1:2], group_VAR = "10A_AT1_overlap")

test_that("groupRegionsByOverlap list", {
    expect_equal(colnames(rowData(ct2.olap)), c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF", "10A_AT1_overlap"))
    expect_equal(levels(rowData(ct2.olap)[["10A_AT1_overlap"]]), c("MCF10A_CTCF & MCF10AT1_CTCF", "MCF10A_CTCF", "MCF10AT1_CTCF",""))
})

ct2.olap = groupRegionsByOverlap(ct2.olap, peak_grs$MCF10A_CTCF, group_VAR = "10A_peaks")

test_that("groupRegionsByOverlap GRanges", {
    expect_equal(colnames(rowData(ct2.olap)), c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF", "10A_AT1_overlap", "10A_peaks"))
    expect_setequal(rowData(ct2.olap)[["10A_peaks"]], c(FALSE, TRUE))
})

test_that("groupRegionsByOverlap errors", {
    expect_error(groupRegionsByOverlap(ct2.olap, list(peak_grs$MCF10A_CTCF), group_VAR = "10A_AT1_overlap"), regexp = "query must have names")
    expect_error(groupRegionsByOverlap(ct2.olap, list(a = peak_grs$MCF10A_CTCF, a = peak_grs$MCF10A_CTCF), group_VAR = "10A_AT1_overlap"), regexp = "All names of query must be unique.")
})


ct2 = exampleChIPtsne2.with_meta()
#note that you can groupRegions using either a list of GRanges or a single GRanges
ct2 = ct2 %>%
    groupRegionsByOverlap(peak_grs$MCF10A_CTCF, "10A_peaks") %>%
    groupRegionsByOverlap(peak_grs$MCF10AT1_CTCF, "AT1_peaks") %>%
    groupRegionsByOverlap(peak_grs$MCF10CA1_CTCF, "ca1a_peaks") %>%
    groupRegionsByOverlap(peak_grs, "peak_overlap")

ct2 = dimReduceTSNE(ct2, perplexity = 25)


theme_update(panel.background = element_rect(fill = "gray80"), panel.grid = element_blank())

plotDimReducePoints(ct2, "peak_overlap")
plotDimReducePoints(ct2, extra_VARS = "peak_overlap") +
    facet_grid(sample~`peak_overlap`)


plotDimReducePoints(ct2, "AT1_peaks")
plotDimReducePoints(ct2, extra_VARS = "AT1_peaks") +
    facet_grid(sample~`AT1_peaks`)


plotDimReducePoints(ct2, c("10A_peaks", "AT1_peaks", "ca1a_peaks"))


plotDimReducePoints(ct2, "10A_peaks")
plotDimReducePoints(ct2, extra_VARS = "10A_peaks") +
    facet_grid(sample~`10A_peaks`)

plotDimReducePoints(subsetSamples(ct2, cell == "MCF10A"), extra_VARS = c("AT1_peaks", "ca1a_peaks")) +
    facet_grid(AT1_peaks~ca1a_peaks)

plotSignalHeatmap(ct2, c("AT1_peaks", "ca1a_peaks", "peak_overlap"))

plotSignalHeatmap(ct2, c("10A_peaks", "AT1_peaks", "ca1a_peaks", "peak_overlap"))

plotDimReduceBins(ct2)

plotDimReduceBins(ct2, facet_columns = c("peak_overlap"))
plotDimReduceBins(ct2, facet_columns = c("10A_peaks"))

plotDimReduceBins(ct2, facet_columns = c("10A_peaks", "AT1_peaks", "ca1a_peaks", "peak_overlap"))
