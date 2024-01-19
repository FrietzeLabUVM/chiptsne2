testthat::context("colnames rownames dimnames")
library(chiptsne2)
library(testthat)

test_that("update_rownames", {
    update_ct2_rownames = chiptsne2:::.update_ct2_rownames
    ct2 = exampleChIPtsne2.with_meta()
    new_rn = paste0("row_", rownames(ct2))
    ct2.new_rn = update_ct2_rownames(ct2, new_names = new_rn)
    expect_equal(ct2.new_rn@region_VAR, "id")
    expect_equal(rownames(ct2.new_rn), new_rn)
    expect_setequal(getTidyProfile(ct2.new_rn)$id, new_rn)
    expect_equal(rownames(ct2.new_rn@assays@data$max), new_rn)
    expect_equal(rownames(ct2.new_rn@rowToRowMat), new_rn)
    expect_equal(names(rowRanges(ct2.new_rn)), new_rn)
})

test_that("update_rownames old_name_VAR and new_name_VAR", {
    update_ct2_rownames = chiptsne2:::.update_ct2_rownames
    ct2 = exampleChIPtsne2.with_meta()
    rowData(ct2)$new_id = paste0("new_", rownames(ct2))
    ct2.new_rn = update_ct2_rownames(ct2, old_name_VAR = "id", new_name_VAR = "new_id")


    new_rn = paste0("new_", rownames(ct2))

    expect_equal(ct2.new_rn@region_VAR, "new_id")
    expect_equal(rownames(ct2.new_rn), new_rn)
    expect_setequal(getTidyProfile(ct2.new_rn)$new_id, new_rn)
    expect_equal(rownames(ct2.new_rn@assays@data$max), new_rn)
    expect_equal(rownames(ct2.new_rn@rowToRowMat), new_rn)
    expect_equal(names(rowRanges(ct2.new_rn)), new_rn)
})

test_that("update_colnames", {
    update_ct2_colnames = chiptsne2:::.update_ct2_colnames
    ct2 = exampleChIPtsne2.with_meta()
    getTidyProfile(ct2)

    new_cn = sub("_CTCF", "", colnames(ct2))

    ct2.new_cn = update_ct2_colnames(ct2, new_cn)
    expect_equal(ct2.new_cn@name_VAR, "sample")
    expect_setequal(getTidyProfile(ct2.new_cn)$sample, c("MCF10A", "MCF10AT1", "MCF10CA1"))
    expect_setequal(colnames(getTidyProfile(ct2.new_cn)), c("id", "x", "y", "sample"))

    expect_equal(rownames(ct2.new_cn@colData), c("MCF10A", "MCF10AT1", "MCF10CA1"))
    expect_equal(rownames(colData(ct2.new_cn)), c("MCF10A", "MCF10AT1", "MCF10CA1"))
    expect_equal(colnames(ct2.new_cn@assays@data$max), c("MCF10A", "MCF10AT1", "MCF10CA1"))

    expect_equal(colnames(ct2.new_cn@rowToRowMat)[1:3], c("MCF10A_-325", "MCF10A_-275", "MCF10A_-225"))
    expect_equal(colnames(ct2.new_cn@rowToRowMat)[40:42], c("MCF10CA1_225", "MCF10CA1_275", "MCF10CA1_325"))
    expect_equal(names(ct2.new_cn@colToRowMatCols), c("MCF10A", "MCF10AT1", "MCF10CA1"))
    cns = unlist(ct2.new_cn@colToRowMatCols)
    names(cns) = NULL
    expect_equal(cns[1:3], c("MCF10A_-325", "MCF10A_-275", "MCF10A_-225"))
    expect_equal(cns[40:42], c("MCF10CA1_225", "MCF10CA1_275", "MCF10CA1_325"))

})

test_that("update_colnames old_name_VAR and new_name_VAR", {
    update_ct2_colnames = chiptsne2:::.update_ct2_colnames
    ct2 = exampleChIPtsne2.with_meta()

    ct2.new_cn = update_ct2_colnames(ct2, old_name_VAR = "sample", new_name_VAR = "cell")
    expect_equal(ct2.new_cn@name_VAR, "cell")
    expect_setequal(getTidyProfile(ct2.new_cn)$cell, c("MCF10A", "MCF10AT1", "MCF10CA1"))
    expect_setequal(colnames(getTidyProfile(ct2.new_cn)), c("id", "x", "y", "cell"))

    expect_equal(rownames(ct2.new_cn@colData), c("MCF10A", "MCF10AT1", "MCF10CA1"))
    expect_equal(rownames(colData(ct2.new_cn)), c("MCF10A", "MCF10AT1", "MCF10CA1"))
    expect_equal(colnames(ct2.new_cn@assays@data$max), c("MCF10A", "MCF10AT1", "MCF10CA1"))

    expect_equal(colnames(ct2.new_cn@rowToRowMat)[1:3], c("MCF10A_-325", "MCF10A_-275", "MCF10A_-225"))
    expect_equal(colnames(ct2.new_cn@rowToRowMat)[40:42], c("MCF10CA1_225", "MCF10CA1_275", "MCF10CA1_325"))
    expect_equal(names(ct2.new_cn@colToRowMatCols), c("MCF10A", "MCF10AT1", "MCF10CA1"))
    cns = unlist(ct2.new_cn@colToRowMatCols)
    names(cns) = NULL
    expect_equal(cns[1:3], c("MCF10A_-325", "MCF10A_-275", "MCF10A_-225"))
    expect_equal(cns[40:42], c("MCF10CA1_225", "MCF10CA1_275", "MCF10CA1_325"))

})
