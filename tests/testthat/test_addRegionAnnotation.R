testthat::context("center signal")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(ggplot2)

ct2 = exampleChIPtsne2.with_meta()
test_gr = head(rowRanges(ct2), n = 20)
memb_df = seqsetvis::ssvFactorizeMembTable(test_gr)
test_gr$group = memb_df[names(test_gr),]$group

ct2 = addRegionAnnotation(ct2, test_gr, anno_VAR = "group")
# debug(addRegionAnnotation)
ct2 = addRegionAnnotation(ct2, test_gr, anno_VAR = "group", anno_VAR_renames = "group2", no_overlap_value = factor("no_hit"))
rowData(ct2)

table(rowData(ct2)$group)
table(rowData(ct2)$group2)

test_that("addRegionAnnotation preserves factors", {
    expect_is(rowData(ct2)$group, "factor")
    expect_is(rowData(ct2)$group2, "factor")
    expect_equal(rowData(ct2)$group, rowData(ct2)$group2)
    expect_equal(table(rowData(ct2)$group)["no_hit"], c("no_hit" = 80))
    expect_equal(table(rowData(ct2)$group)["peak_MCF10A_CTCF"], c("peak_MCF10A_CTCF" = 5))

})

test_that("addRegionAnnotation preserves factors", {
    expect_is(rowData(ct2)$group, "factor")
})
