testthat::context("ChIPtsne_class")
library(chiptsne2)
library(testthat)

ct2 = exampleChIPtsne2.with_meta()

test_that("NameVariable", {
    expect_equal(getNameVariable(ct2), "sample")
    ct2.new = setNameVariable(ct2, "Name")
    expect_equal(getNameVariable(ct2.new), "Name")
})

test_that("RegionVariable", {
    expect_equal(getRegionVariable(ct2), "id")
    ct2.new = setRegionVariable(ct2, "Region")
    expect_equal(getRegionVariable(ct2.new), "Region")
})

test_that("PositionVariable", {
    expect_equal(getPositionVariable(ct2), "x")
    ct2.new = setPositionVariable(ct2, "Position")
    expect_equal(getPositionVariable(ct2.new), "Position")
})

test_that("ValueVariable", {
    expect_equal(getValueVariable(ct2), "y")
    ct2.new = setValueVariable(ct2, "Value")
    expect_equal(getValueVariable(ct2.new), "Value")
})

test_that("Set Secondary Impacts", {
    ct2.new = ct2
    ct2.new = setNameVariable(ct2.new, "Name")
    ct2.new = setRegionVariable(ct2.new, "Region")
    ct2.new = setPositionVariable(ct2.new, "Position")
    ct2.new = setValueVariable(ct2.new, "Value")

    expect_equal(getNameVariable(ct2.new), "Name")
    expect_equal(getRegionVariable(ct2.new), "Region")
    expect_equal(getPositionVariable(ct2.new), "Position")
    expect_equal(getValueVariable(ct2.new), "Value")

    hist1 = ChIPtsne2.history(ct2.new)
    expect_setequal(names(hist1),
        c("birthday", "session_info", "chiptsne2_version", "setNameVariable", "setRegionVariable", "setPositionVariable", "setValueVariable"))

    ct2.rerun = rerun_history(ct2, history_source = ct2.new)
    expect_equal(getNameVariable(ct2.rerun), "Name")
    expect_equal(getRegionVariable(ct2.rerun), "Region")
    expect_equal(getPositionVariable(ct2.rerun), "Position")
    expect_equal(getValueVariable(ct2.rerun), "Value")

    hist2 = ChIPtsne2.history(ct2.rerun)
    expect_setequal(names(hist2),
                    c("birthday", "session_info", "chiptsne2_version", "setNameVariable", "setRegionVariable", "setPositionVariable", "setValueVariable"))

    prof_dt = getTidyProfile(ct2.new)
    expect_setequal(colnames(prof_dt), c("Region", "Position", "Value", "Name"))

    row_meta_df = getRegionMetaData(ct2.new)
    expect_setequal(colnames(row_meta_df), c("Region", "peak_MCF10A_CTCF", "peak_MCF10AT1_CTCF", "peak_MCF10CA1_CTCF"))

    col_meta_df = getSampleMetaData(ct2.new)
    expect_setequal(colnames(col_meta_df), c("cell", "mark", "Name"))

    expect_setequal(colnames(ct2.new), colnames(ct2))

    expect_setequal(rownames(ct2.new), rownames(ct2))
})

#setNameVariable has some special checks and situations

test_that("Error when existing not unique", {
    expect_error(setNameVariable(ct2, "mark"), regexp = "New name variable is already present in colData\\(ct2\\). This is only allowed when all values are unique. Offending values: CTCF",)
})

test_that("Switch name variable to smaller set", {
    ct2.by_cell = split(ct2, "cell")
    ct2.10a = setNameVariable(ct2.by_cell$MCF10A, new_name_VAR = "mark")
    expect_equal(getNameVariable(ct2.10a), "mark")
    expect_equal(getSampleMetaData(ct2.10a)$mark, factor("CTCF"))
    expect_equal(getSampleMetaData(ct2.10a)$sample, "MCF10A_CTCF")
    expect_equal(colnames(rowToRowMat(ct2.10a))[1], "CTCF_-325")
})

test_that("When names match for operator", {
    ct2.by_cell = split(ct2, "cell")
    ct2.by_cell = lapply(ct2.by_cell, setNameVariable, new_name_VAR = "mark")
    ct2.diff = ct2.by_cell$MCF10A - ct2.by_cell$MCF10AT1
    expect_equal(colnames(rowToRowMat(ct2.diff))[1], "CTCF_-325")
    ct2.diff1 = setNameVariable(ct2.diff, new_name_VAR = "cell")
    ct2.diff2 = setNameVariable(ct2.diff, new_name_VAR = "sample")
    expect_equal(colnames(ct2.diff1), "MCF10A - MCF10AT1")
    expect_equal(colnames(rowToRowMat(ct2.diff1))[1], "MCF10A - MCF10AT1_-325")
    expect_equal(colnames(ct2.diff2), "MCF10A_CTCF - MCF10AT1_CTCF")
    expect_equal(colnames(rowToRowMat(ct2.diff2))[1], "MCF10A_CTCF - MCF10AT1_CTCF_-325")
})

test_that("Switch name variable same set", {
    ct2.name_cell = setNameVariable(ct2, "cell")
    expect_equal(getNameVariable(ct2.name_cell), "cell")
    expect_equal(getSampleMetaData(ct2.name_cell)$mark, rep("CTCF", 3))
    expect_equal(colnames(rowToRowMat(ct2.name_cell))[1], "MCF10A_-325")
    expect_equal(getSampleMetaData(ct2.name_cell)$sample, c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF"))
})


