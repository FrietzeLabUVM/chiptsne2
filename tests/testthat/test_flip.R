testthat::context("flip")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(ggplot2)

ct2 = exampleChIPtsne2.with_meta()
ct2 = ct2[, "MCF10A_CTCF"]
#add a metadata column
colData(ct2)$flip = "none"
ct2_left = flipProfilesToMatch(ct2)
colData(ct2_left)$flip = "left"
#colnames will need to be different
colnames(ct2_left) = paste0(colnames(ct2_left), "_left")
ct2_right = flipProfilesToMatch(ct2, highest_on_right = TRUE)
colData(ct2_right)$flip = "right"
colnames(ct2_right) = paste0(colnames(ct2_right), "_right")

#compare profile matrixes
mat = rowToRowMat(ct2)
nc = ncol(mat)
left_i = seq(1, nc/2)
right_i = seq(nc/2+1, nc)
mat_left = rowToRowMat(ct2_left)
mat_right = rowToRowMat(ct2_right)

m_raw_left = mean(mat[, left_i])
m_raw_right = mean(mat[, right_i])

m_left_left = mean(mat_left[, left_i])
m_left_right = mean(mat_left[, right_i])

m_right_left = mean(mat_right[, left_i])
m_right_right = mean(mat_right[, right_i])

test_that("flipProfilesToMatch changes profile positions", {
    #flip for left means higher on left
    expect_gt(m_left_left, m_raw_left)
    #flip for left means lower on right
    expect_lt(m_left_right, m_raw_right)

    #flip for right means low on left
    expect_lt(m_right_left, m_raw_left)
    #flip for right means high on right
    expect_gt(m_right_right, m_raw_right)

    #flip for right right equals flip for left left
    expect_equal(m_right_right, m_left_left)
    #flip for right left equals flip for left right
    expect_equal(m_left_right, m_right_left)
})

test_that("flipProfilesToMatch adds strand orientation", {
    not_flipped_no_strand = all(GenomicRanges::strand(rowRanges(ct2)) == "*")
    expect_true(not_flipped_no_strand)

    flipped_has_strand = all(as.character(GenomicRanges::strand(rowRanges(ct2_left))) %in% c("-", "+"))
    expect_true(flipped_has_strand)

    left_and_right_inverted = all(rowRanges(ct2_left) == GenomicRanges::invertStrand(rowRanges(ct2_right)))
    expect_true(left_and_right_inverted)
})
