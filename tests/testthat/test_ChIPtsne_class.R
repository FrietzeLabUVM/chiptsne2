testthat::context("ChIPtsne_class")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(SummarizedExperiment)

query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
prof_dt = seqsetvis::CTCF_in_10a_profiles_dt

meta_dt = unique(prof_dt[, .(sample)])
meta_dt[, c("cell", "mark") := data.table::tstrsplit(sample, "_")]

map_dt = unique(prof_dt[, .(x, sample)])
map_dt[, cn := paste(sample, x, sep = "_")]
cn = map_dt$cn
setequal(cn, colnames(prof_mat))
cn == colnames(prof_mat)
map_dt$cn
map_dt[, nr := seq(.N)]
map_list = split(map_dt$n, map_dt$sample)

tmp_wide = tidyr::pivot_wider(prof_dt, names_from = c("sample", "x"), values_from = "y", id_cols = "id")
prof_mat = as.matrix(tmp_wide[, -1])
rownames(prof_mat) = tmp_wide$id


prof_max = prof_dt[, .(y = max(y)), .(id, sample)] %>%
    tidyr::pivot_wider(names_from = "sample", id_cols = "id", values_from = "y")
prof_max_mat = as.matrix(prof_max[, -1])
rownames(prof_max_mat) = prof_max$id

ct = ChIPtsne2(assay = list(max = prof_max_mat),
               rowRanges = query_gr,
               rowToRowMat = prof_mat,
               colToRowMatCols = map_list,
               colData = meta_dt,
               metadata = list(time = date()))

# debug(ChIPtsne2.from_tidy)
debug(ChIPtsne2.from_tidy)
ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr)

test_that("Constructors - valid", {
    expect_true(validObject(ct2))
    expect_true(validObject(.ChIPtsne2())) # internal
    expect_true(validObject(ChIPtsne2())) # exported

    se = as(ct2, "SummarizedExperiment")
    expect_true(validObject(se))
    # conv <- as(se, "ExampleClass")
    # expect_true(validObject(conv))
})

test_that("Constructors - invalid", {
    expect_error(ChIPtsne2(rowToRowMat=rbind(1)), "nrow\\(rowToRowMat\\)")
    expect_error(ChIPtsne2(colToRowMatCols=list(1)), "length\\(colToRowMatCols\\)")
})

test_that("Gettters", {

    expect_identical(rowToRowMat(ct2, withDimnames=FALSE), prof_mat)
    expect_identical(rownames(rowToRowMat(ct2)), rownames(ct2))

    expect_identical(names(colToRowMatCols(ct2)), colnames(ct2))


})
