testthat::context("plotDimReducePoints details")
library(testthat)
library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta() %>%
   dimReduceUMAP() %>%
   groupRegionsByDimReduceCluster(group_VAR = "umap_cluster") %>%
   groupRegionsBySignalCluster(group_VAR = "signal_cluster")

expect_s3_class(class = "ggplot",
                plotDimReducePoints(ct2, color_VAR = "MCF10A_CTCF", label_size = 14))

expect_s3_class(class = "ggplot",
                plotDimReducePoints(ct2, color_VAR = c("umap_cluster"),
                                    label_VAR = "umap_cluster", label_size = 14))

expect_s3_class(class = "ggplot",
                plotDimReducePoints(ct2, color_VAR = c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"),
                                    label_VAR = "umap_cluster", label_size = 14))



expect_error(
    plotDimReducePoints(ct2, color_VAR = c("umap_cluster", "peak_MCF10A_CTCF"),
                        label_VAR = "umap_cluster"),
    "Classes of all color_VAR items must match.")




expect_error(
    plotDimReducePoints(ct2, color_VAR = c("umap_cluster"),
                        label_VAR = "groupasdf"),
    "Some VAR are missing from metadata")

expect_s3_class(class = "ggplot",
                plotDimReducePoints(ct2, color_VAR = c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"), label_VAR = "group", label_size = 14)
)
# plotDimReducePoints(ct2, color_VAR = c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"), label_VAR = "value", label_size = 14)

expect_s3_class(class = "ggplot",
                plotDimReducePoints(ct2, label_VAR = "umap_cluster", label_size = 14)
)


expect_s3_class(class = "ggplot",
                plotDimReducePoints(ct2,
                                    label_VAR = "umap_cluster",
                                    label_size = 14,
                                    extra_VARS = "umap_cluster",
                                    background_annotation_color = "gray") +
                    ggplot2::facet_wrap(~umap_cluster) +
                    ggplot2::theme(
                        panel.background = ggplot2::element_rect(fill = "gray50"),
                        panel.grid = ggplot2::element_blank())
)
