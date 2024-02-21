testthat::context("split cbind rbind List")
library(chiptsne2)
library(testthat)
ct2 = exampleChIPtsne2.with_meta()
suppressWarnings({
    ct2 = dimReduceTSNE(ct2)
})


ct2.by_cell = split(ct2, "cell")
ct2.by_peak = split(ct2, "peak_MCF10CA1_CTCF")

test_that("split returns ChIPtsne2List", {
    expect_is(ct2.by_cell, "ChIPtsne2List")
    expect_is(ct2.by_peak, "ChIPtsne2List")
})

test_that("ChIPtsne2List has names", {
    expect_equal(names(ct2.by_cell), c("MCF10A", "MCF10AT1", "MCF10CA1"))
    expect_equal(names(as.list(ct2.by_cell)), c("MCF10A", "MCF10AT1", "MCF10CA1"))

    expect_equal(names(ct2.by_peak), c("FALSE", "TRUE"))
    expect_equal(names(as.list(ct2.by_peak)), c("FALSE", "TRUE"))
})

test_that("cbind ChIPtsne2", {
    ct2.cbind = cbind(ct2.by_cell$MCF10A, ct2.by_cell$MCF10AT1, ct2.by_cell$MCF10CA1)
    expect_equal(ncol(ct2.cbind), 3)
    expect_equal(colnames(ct2.cbind), c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF"))
    ct2.cbind2 = cbind(ct2.by_cell$MCF10AT1, ct2.by_cell$MCF10CA1, ct2.by_cell$MCF10A)
    expect_equal(colnames(ct2.cbind2), c("MCF10AT1_CTCF", "MCF10CA1_CTCF", "MCF10A_CTCF"))

    expect_error(cbind(ct2.by_peak$`FALSE`, ct2.by_peak$`TRUE`), regexp = "Duplicated Column names are not allowed")
})

test_that("rbind ChIPtsne2", {
    ct2.rbind = rbind(ct2.by_peak$`FALSE`, ct2.by_peak$`TRUE`)
    ct2.rbind2 = rbind(ct2.by_peak$`TRUE`, ct2.by_peak$`FALSE`)
    expect_setequal(rownames(ct2.rbind), rownames(ct2.rbind2))
    expect_setequal(rownames(ct2.rbind), c(rownames(ct2.by_peak$`FALSE`), rownames(ct2.by_peak$`TRUE`)))
    expect_setequal(rownames(ct2.rbind2), c(rownames(ct2.by_peak$`TRUE`), rownames(ct2.by_peak$`FALSE`)))
    expect_error(rbind(ct2.by_cell$MCF10A, ct2.by_cell$MCF10AT1, ct2.by_cell$MCF10CA1), regexp = "Duplicated Row names are not allowed ")
})

test_that("cbind ChIPtsne2List", {
    ct2.cbind = cbind(ct2.by_cell)
    expect_equal(ncol(ct2.cbind), 3)
    expect_equal(colnames(ct2.cbind), c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF"))
    ct2.cbind2 = cbind(ct2.by_cell[c(2, 3, 1)])
    expect_equal(colnames(ct2.cbind2), c("MCF10AT1_CTCF", "MCF10CA1_CTCF", "MCF10A_CTCF"))

    expect_error(cbind(ct2.by_peak), regexp = "Duplicated Column names are not allowed")
})

test_that("rbind ChIPtsne2List", {
    ct2.rbind = rbind(ct2.by_peak)
    ct2.rbind2 = rbind(ct2.by_peak[c(2, 1)])
    expect_setequal(rownames(ct2.rbind), rownames(ct2.rbind2))
    expect_setequal(rownames(ct2.rbind), c(rownames(ct2.by_peak$`FALSE`), rownames(ct2.by_peak$`TRUE`)))
    expect_setequal(rownames(ct2.rbind2), c(rownames(ct2.by_peak$`TRUE`), rownames(ct2.by_peak$`FALSE`)))
    expect_error(rbind(ct2.by_cell$MCF10A, ct2.by_cell$MCF10AT1, ct2.by_cell$MCF10CA1), regexp = "Duplicated Row names are not allowed ")
})


