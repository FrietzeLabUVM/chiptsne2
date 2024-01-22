testthat::context("all_functions")
# flipping viewGranges
library(chiptsne2)
library(testthat)

ct2_objects = list(
    ct2_noMeta = exampleChIPtsne2(),
    ct2_wMeta = exampleChIPtsne2.with_meta()
)

ct2_objects.no_rowRanges = lapply(ct2_objects, function(x){
    rowRanges(x) = NULL
    x
})
names(ct2_objects.no_rowRanges) = paste0(names(ct2_objects), ".no_rowRanges")

ct2_objects.all = c(ct2_objects, ct2_objects.no_rowRanges)
obs_classes = sapply(ct2_objects.all, class)
names(obs_classes) = NULL

test_that("classes", {
    expect_equal(obs_classes, c("ChIPtsne2", "ChIPtsne2", "ChIPtsne2_no_rowRanges", "ChIPtsne2_no_rowRanges"))
})


ct2_objects.all.wPCA = lapply(ct2_objects.all, dimReducePCA)
ct2_objects.all.wTSNE = lapply(ct2_objects.all, dimReduceTSNE, perplexity = 25)
ct2_objects.all.wUMAP = lapply(ct2_objects.all, dimReduceUMAP)

sapply(ct2_objects.all, hasDimReduce)
sapply(ct2_objects.all.wPCA, hasDimReduce)
sapply(ct2_objects.all.wTSNE, hasDimReduce)
sapply(ct2_objects.all.wUMAP, hasDimReduce)

names(ct2_objects.all.wTSNE$ct2_noMeta@metadata)
ct2_objects.all.wPCA$ct2_noMeta@metadata$dimReducePCA$ARG
ct2_objects.all.wTSNE$ct2_noMeta@metadata$dimReduceTSNE$ARG
ct2_objects.all.wUMAP$ct2_noMeta@metadata$dimReduceUMAP$ARG

plot
