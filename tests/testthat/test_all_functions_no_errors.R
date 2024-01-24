testthat::context("all_functions")
# flipping viewGranges
library(chiptsne2)
library(testthat)

ct2_objects = list(
    ct2_noMeta = exampleChIPtsne2(),
    ct2_wMeta = exampleChIPtsne2.with_meta()
)

#all functions should adjust to renamed variables
ct2_objects = lapply(ct2_objects, function(x){
    x = setPositionVariable(x, "new_pos")
    x = setNameVariable(x, "new_name")
    x = setRegionVariable(x, "new_reg")
    x = setValueVariable(x, "new_value")
    x
})


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

test_that("reportDimReduceMethod", {
    expect_true(all(!sapply(ct2_objects.all, hasDimReduce)))
    expect_true(all(sapply(ct2_objects.all.wPCA, hasDimReduce)))
    expect_true(all(sapply(ct2_objects.all.wTSNE, hasDimReduce)))
    expect_true(all(sapply(ct2_objects.all.wUMAP, hasDimReduce)))

    expect_true(all(sapply(ct2_objects.all, reportDimReduceMethod) == "none"))
    expect_true(all(sapply(ct2_objects.all.wPCA, reportDimReduceMethod) == "PCA"))
    expect_true(all(sapply(ct2_objects.all.wTSNE, reportDimReduceMethod) == "t-SNE"))
    expect_true(all(sapply(ct2_objects.all.wUMAP, reportDimReduceMethod) == "UMAP"))
})

FUN_to_test = list(
    groupRegionsByDimReduceCluster = groupRegionsByDimReduceCluster,
    plotDimReduceBins = plotDimReduceBins,
    plotDimReducePoints = plotDimReducePoints,
    plotDimReduceSummaryProfiles = plotDimReduceSummaryProfiles,
    plotSignalHeatmap = plotSignalHeatmap,
    plotSignalLinePlot = plotSignalLinePlot
)

ct2_objects.combined = list(
    none = ct2_objects.all,
    PCA = ct2_objects.all.wPCA,
    TSNE = ct2_objects.all.wTSNE,
    UMAP = ct2_objects.all.wUMAP
)

all_res = list()
for(dim_method in names(ct2_objects.combined)){
    ct2_objects.sel = ct2_objects.combined[[dim_method]]
    message(dim_method)
    for(ct2_obj in names(ct2_objects.sel)){
        ct2.sel = ct2_objects.sel[[ct2_obj]]
        message("..", ct2_obj)
        for(test_fun in names(FUN_to_test)){
            message("....", test_fun)
            fun.sel = FUN_to_test[[test_fun]]
            res = tryCatch(
                expr = {
                    suppressWarnings({
                        fun.sel(ct2.sel)
                    })
                },
                error = function(e){
                    paste0("error: ", as.character(e))
                })
            if(!is.character(res)){
                res = class(res)[1]
            }
            new_df = data.frame(method = dim_method, obj_type = ct2_obj, FUN_name = test_fun, result = res)
            if(ncol(new_df) == 3) browser()
            all_res = c(all_res, list(new_df))

        }
    }
}
all_df = do.call(rbind, all_res)
dim(all_df)
all_df$simple_result = sub(":.+", "", all_df$result)
table(all_df$simple_result)


caught_error = grepl("No dimensional reduction data present in this ChIPtsne2 object.", all_df$result)
sum(caught_error)
all_df$simple_result[caught_error] = "caught"

fail_df = subset(all_df, simple_result == "error")
nrow(fail_df)

if(nrow(fail_df) > 0){
    message(paste(fail_df$method[1], fail_df$obj_type[1], fail_df$FUN_name[1]))
    ct2.sel = ct2_objects.combined[[fail_df$method[1]]][[fail_df$obj_type[1]]]
    fun.sel = FUN_to_test[[fail_df$FUN_name[1]]]
    fun.sel(ct2.sel)
}

test_that("all test FUN", {
    expect_equal(nrow(fail_df), 0)
})
