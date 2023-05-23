testthat::context("ChIPtsne_class")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(SummarizedExperiment)
library(data.table)

query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
prof_dt = seqsetvis::CTCF_in_10a_profiles_dt

meta_dt = prof_dt %>% dplyr::select(sample) %>% unique
meta_dt = meta_dt %>% tidyr::separate(sample, c("cell", "mark"), sep = "_", remove = FALSE)

map_dt = prof_dt %>% dplyr::select(sample, x) %>% unique %>%
    dplyr::mutate(cn = paste(sample, x, sep = "_")) %>%
    dplyr::mutate(nr = seq_along(x))
map_list = split(map_dt$n, map_dt$sample)

tmp_wide = tidyr::pivot_wider(prof_dt, names_from = c("sample", "x"), values_from = "y", id_cols = "id")
prof_mat = as.matrix(tmp_wide[, -1])
rownames(prof_mat) = tmp_wide$id


prof_max = prof_dt %>%
    dplyr::group_by(id, sample) %>%
    dplyr::summarise(y = max(y)) %>%
    tidyr::pivot_wider(names_from = "sample", id_cols = "id", values_from = "y")
prof_max_mat = as.matrix(prof_max[, -1])
rownames(prof_max_mat) = prof_max$id

ct = ChIPtsne2(assay = list(max = prof_max_mat[names(query_gr),]),
               rowRanges = query_gr,
               rowToRowMat = prof_mat,
               colToRowMatCols = map_list,
               colData = meta_dt,
               metadata = list(time = date()))


prof_dt = seqsetvis::ssvSignalClustering(prof_dt, nclust = 4)
region_meta_dt = prof_dt %>% dplyr::select(id, cluster_id) %>% unique

# debug(ChIPtsne2.from_tidy)
ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, region_meta_dt = region_meta_dt)

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
    expect_identical(rowToRowMat(ct2, withDimnames=FALSE), prof_mat[levels(prof_dt$id),])
    expect_identical(rownames(rowToRowMat(ct2)), rownames(ct2))
    expect_identical(names(colToRowMatCols(ct2)), colnames(ct2))
})
