testthat::context("transform")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(ggplot2)

ct2 = exampleChIPtsne2.with_meta()
# plotSignalLinePlot(ct2) + expand_limits(y = 30)
# a transformation function takes a matrix and returns the same size matrix

test_that("transform - works", {
    trans_fun_1 = function(mat, ...){
        mat/2
    }

    ct2.t1 = transformSignal(ct2, trans_fun_1)
    expect_equal(colMeans(rowToRowMat(ct2)), 2*colMeans(rowToRowMat(ct2.t1)))
    # plotSignalLinePlot(ct2.t1) + expand_limits(y = 30)

    # as a second argument, sample metadata is passed in per sample
    trans_fun_2 = function(mat, sample_meta){
        if(sample_meta$cell == "MCF10CA1"){
            mat / 4
        }else{
            mat
        }
    }
    ct2.t2 = transformSignal(ct2, trans_fun_2)
    expect_equal(
        colMeans(rowToRowMat(subsetSamples(ct2, cell == "MCF10CA1"))),
        4*colMeans(rowToRowMat(subsetSamples(ct2.t2, cell == "MCF10CA1")))
    )
    expect_equal(
        colMeans(rowToRowMat(subsetSamples(ct2, cell != "MCF10CA1"))),
        colMeans(rowToRowMat(subsetSamples(ct2.t2, cell != "MCF10CA1")))
    )
    # plotSignalLinePlot(ct2.t2) + expand_limits(y = 30)
})

test_that("transform - errors", {
    trans_fun_bad_class = function(mat, ...){
        as.data.frame(mat)
    }
    expect_error(transformSignal(ct2, trans_fun_bad_class), "Output of transformation must be a matrix. Was: data.frame")

    trans_fun_bad_dim = function(mat, ...){
        mat[-1,]
    }
    expect_error(transformSignal(ct2, trans_fun_bad_dim), "Dimensions of transformation result \\(99x14\\) not equal to input matrix \\(100x14\\).")
})
