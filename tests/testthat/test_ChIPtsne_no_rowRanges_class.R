testthat::context("ChIPtsne2_no_rowRanges_class")
# flipping viewGranges
library(chiptsne2)
library(testthat)

ct2 = exampleChIPtsne2.with_meta()

test_that("valid ChIPtsne2_no_rowRanges from constructor", {
    ct2.nrr = ChIPtsne2_no_rowRanges(assays = assays(ct2),
                                     rowData = rowData(ct2),
                                     colData = colData(ct2),
                                     rowToRowMat = ct2@rowToRowMat,
                                     colToRowMatCols = ct2@colToRowMatCols,
                                     name_VAR = ct2@name_VAR,
                                     position_VAR = ct2@position_VAR, value_VAR = ct2@value_VAR, region_VAR = ct2@region_VAR)
    expect_equal(as.character(class(ct2.nrr)), "ChIPtsne2_no_rowRanges")
    expect_equal(colnames(rowData(ct2.nrr)), c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))
    expect_equal(rownames(rowData(ct2.nrr))[1:5], c("1", "2", "3", "4", "5"))
    reg_meta = getRegionMetaData(ct2.nrr)

    expect_equal(colnames(reg_meta), c("id", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))
    expect_equal(rownames(reg_meta)[1:5], c("1", "2", "3", "4", "5"))
})

test_that("valid ChIPtsne2_no_rowRanges from nullify rowRanges", {
    # nullify rowRanges
    ct2.nrr = ct2
    rowRanges(ct2.nrr) = NULL

    expect_equal(as.character(class(ct2.nrr)), "ChIPtsne2_no_rowRanges")
    expect_equal(colnames(rowData(ct2.nrr)), c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))
    expect_equal(rownames(rowData(ct2.nrr))[1:5], c("1", "2", "3", "4", "5"))
    reg_meta = getRegionMetaData(ct2.nrr)

    expect_equal(colnames(reg_meta), c("id", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))
    expect_equal(rownames(reg_meta)[1:5], c("1", "2", "3", "4", "5"))
})

test_that("valid ChIPtsne2_no_rowRanges from tidy", {
    # nullify rowRanges
    prof_dt = getTidyProfile(ct2)
    ct2.nrr = ChIPtsne2.from_tidy(prof_dt, query_gr = NULL, region_metadata = getRegionMetaData(ct2))

    expect_equal(as.character(class(ct2.nrr)), "ChIPtsne2_no_rowRanges")
    expect_equal(colnames(rowData(ct2.nrr)), c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))
    expect_equal(rownames(rowData(ct2.nrr))[1:5], c("1", "2", "3", "4", "5"))
    reg_meta = getRegionMetaData(ct2.nrr)

    expect_equal(colnames(reg_meta), c("id", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))
    expect_equal(rownames(reg_meta)[1:5], c("1", "2", "3", "4", "5"))
})

test_that("valid ChIPtsne2_no_rowRanges by aggregate", {
    ct2.nrr = aggregateRegionsByGroup(ct2, c("peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
    rowData(ct2.nrr)
    colData(ct2.nrr)
    getRegionMetaData(ct2.nrr)

    expect_equal(as.character(class(ct2.nrr)), "ChIPtsne2_no_rowRanges")
    expect_equal(rownames(rowData(ct2.nrr)), c("FALSE,FALSE", "FALSE,TRUE", "TRUE,FALSE", "TRUE,TRUE"))
    reg_meta = getRegionMetaData(ct2.nrr)

    expect_equal(colnames(reg_meta), c("meta_id", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF"))
    expect_equal(rownames(reg_meta), c("FALSE,FALSE", "FALSE,TRUE", "TRUE,FALSE", "TRUE,TRUE"))
})

test_that("aggregateSamplesByGroup mark", {
    ct2.agg = aggregateSamplesByGroup(ct2, "mark")
    expect_equal(
        rowData(ct2.agg),
        rowData(ct2)
    )

    colData(ct2)

    getSampleMetaData(ct2.agg)
    colData(ct2.agg)

})

test_that("aggregateByGroup variable names", {
    rowData(ct2)$all = "all"

    ct2.1col = aggregateSamplesByGroup(ct2, "mark")
    expect_equal(getRegionVariable(ct2.1col), "id")
    expect_equal(getNameVariable(ct2.1col), "mark")


    ct2.1row = aggregateRegionsByGroup(ct2, "all")
    expect_equal(getRegionVariable(ct2.1row), "all")
    expect_equal(getNameVariable(ct2.1row), "sample")

    ct2.1x1 = aggregateByGroup(ct2, c("all", "mark"))
    expect_equal(getRegionVariable(ct2.1x1), "all")
    expect_equal(getNameVariable(ct2.1x1), "mark")

    p.1row = plotSignalLinePlot(ct2.1row, group_VAR = "all")
    p.1col = plotSignalLinePlot(ct2.1col, group_VAR = "all")
    p.1x1 = plotSignalLinePlot(ct2.1x1, group_VAR = "all")

    # plotSignalHeatmap(ct2.1row)
    # plotSignalHeatmap(ct2.1col)
    # plotSignalHeatmap(ct2.1x1)

    expect_equal(colnames(getRegionMetaData(ct2.1row)), c("all"))
    expect_equal(colnames(getRegionMetaData(ct2.1col)), c("id", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF", "all"))
    expect_equal(colnames(getRegionMetaData(ct2.1x1)), c("all"))

    expect_equal(nrow(getRegionMetaData(ct2.1row)), 1)
    expect_equal(nrow(getRegionMetaData(ct2.1col)), 100)
    expect_equal(nrow(getRegionMetaData(ct2.1x1)), 1)
})



