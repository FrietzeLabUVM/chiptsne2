testthat::context("expression")
# flipping viewGranges
library(chiptsne2)
library(testthat)

#several chiptsne2 functions allow the use of expression as arguments
#I figure these will break in similar ways so test them together here

ct2 = exampleChIPtsne2.with_meta()
rowData(ct2)

test_that("subset - works", {
    #SummarizedExperiment provides subset
    ct2_1 = subset(ct2, peak_MCF10AT1_CTCF == TRUE, cell == "MCF10A")
    expect_equal(dim(ct2_1), c(74, 1))
    sel_cell = "MCF10A"
    ct2_2 = subset(ct2, peak_MCF10AT1_CTCF == TRUE, cell == sel_cell)
    expect_equal(dim(ct2_2), c(74, 1))

})

test_that("subsetRegions - works", {
    ct2_1 = subsetRegions(ct2, peak_MCF10AT1_CTCF == TRUE)
    expect_equal(dim(ct2_1), c(74, 3))

    ct2_2 = subsetRegions(ct2, id == "1")
    expect_equal(dim(ct2_2), c(1, 3))

    ct2_history = ChIPtsne2.history(ct2_1)
    expect_equal(names(ct2_history)[[length(ct2_history)]], "subsetRegions")
})


test_that("subsetSamples - works", {
    ct2_1 = subsetSamples(ct2, cell == "MCF10A")
    ct2_2 = subsetSamples(ct2, sample == "MCF10A_CTCF")
    expect_equal(dim(ct2_1), c(100, 1))
    expect_equal(dim(ct2_2), c(100, 1))

    ct2_history = ChIPtsne2.history(ct2_1)
    expect_equal(names(ct2_history)[[length(ct2_history)]], "subsetSamples")
})

test_that("subsetValues - works", {
    ct2_1 = subsetValues(ct2, MCF10A_CTCF > 30)
    expect_equal(dim(ct2_1), c(68, 3))

    ct2_history = ChIPtsne2.history(ct2_1)
    expect_equal(names(ct2_history)[[length(ct2_history)]], "subsetValues")
})

test_that("mutateRegions - works", {
    ct2_1 = mutateRegions(ct2, "either_10a_or_at1", peak_MCF10A_CTCF | peak_MCF10AT1_CTCF)
    expect_equal(dim(ct2_1), c(100, 3))
    expect_equal(sum(rowData(ct2_1)$either_10a_or_at1), 99)

    ct2_2 = mutateRegions(ct2, "silly", paste(id, peak_MCF10A_CTCF))
    expect_equal(dim(ct2_1), c(100, 3))
    expect_equal(sum(rowData(ct2_1)$either_10a_or_at1), 99)


    ct2_history = ChIPtsne2.history(ct2_1)
    expect_equal(names(ct2_history)[[length(ct2_history)]], "mutateRegions")
})

test_that("mutateSamples - works", {
    colData(ct2)$cell
    ct2_1 = mutateSamples(ct2, "cell_short", sub("MCF", "", cell))
    ct2_1 = mutateSamples(ct2_1, "cell_lower", paste0("mcf", "", cell_short))
    expect_equal(dim(ct2_1), c(100, 3))
    expect_equal(colData(ct2_1)$cell_short, c("10A", "10AT1", "10CA1"))
    expect_equal(colData(ct2_1)$cell_lower, c("mcf10A", "mcf10AT1", "mcf10CA1"))


    ct2_history = ChIPtsne2.history(ct2_1)
    expect_equal(names(ct2_history)[[length(ct2_history)]], "mutateSamples")
})

test_that("separateRegions - works", {
    head(getRegionMetaData(ct2))
    ct2_1 = mutateRegions(ct2, "silly", paste(id, as.character(peak_MCF10A_CTCF)))
    rowData(ct2_1)
    ct2_1 = separateRegions(ct2_1, "silly", sep = " ", into = c("id2", "peak2"))
    rowData(ct2_1)
    expect_equal(dim(ct2_1), c(100, 3))

    ct2_history = ChIPtsne2.history(ct2_1)
    expect_equal(names(ct2_history)[[length(ct2_history)]], "separateRegions")
})

test_that("separateSamples - works", {
    ct2_1 = separateSamples(ct2, col = "sample", into = c("cell2", "mark2"))
    expect_equal(dim(ct2_1), c(100, 3))
    colData(ct2_1)
    expect_equal(colData(ct2_1)$cell, c("MCF10A", "MCF10AT1", "MCF10CA1"))


    ct2_history = ChIPtsne2.history(ct2_1)
    expect_equal(names(ct2_history)[[length(ct2_history)]], "separateSamples")
})
