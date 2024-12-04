testthat::context("ChIPtsne_class")
# flipping viewGranges
library(chiptsne2)
library(testthat)

query_gr = exampleQueryGR()
prof_dt = exampleProfDT()

metadata = prof_dt %>% dplyr::select(name) %>% unique
metadata = metadata %>% tidyr::separate(name, c("cell", "mark"), sep = "_", remove = FALSE)

map_dt = prof_dt %>% dplyr::select(name, position) %>% unique %>%
    dplyr::mutate(cn = paste(name, position, sep = "_")) %>%
    dplyr::mutate(nr = seq_along(position))

map_list = split(map_dt$nr, map_dt$name)

tmp_wide = tidyr::pivot_wider(prof_dt, names_from = c("name", "position"), values_from = "value", id_cols = "region")
prof_mat = as.matrix(tmp_wide[, -1])
rownames(prof_mat) = tmp_wide$region


prof_max = prof_dt %>%
    dplyr::group_by(region, name) %>%
    dplyr::summarise(value = max(value)) %>%
    tidyr::pivot_wider(names_from = "name", id_cols = "region", values_from = "value")
prof_max_mat = as.matrix(prof_max[, -1])
rownames(prof_max_mat) = prof_max$region

ct = ChIPtsne2(assay = list(max = prof_max_mat[names(query_gr),]),
               rowRanges = query_gr,
               rowToRowMat = prof_mat,
               colToRowMatCols = map_list,
               colData = metadata,
               metadata = list(time = date()))


clust_dt = seqsetvis::ssvSignalClustering(prof_dt, nclust = 4, facet_ = "name", row_ = "region", column_ = "position", fill_ = "value")
# prof_dt = translateSSVtoCT2(prof_dt)
# clust_dt = translateSSVtoCT2(clust_dt)

region_metadata = clust_dt %>% dplyr::select(region, cluster_id) %>% unique


ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr, region_metadata = region_metadata)

test_that("Constructors - valid", {
    expect_true(validObject(ct2))
    expect_true(validObject(chiptsne2:::.ChIPtsne2())) # internal
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
    expect_identical(rowToRowMat(ct2), prof_mat[unique(prof_dt$region),])
    expect_identical(rownames(rowToRowMat(ct2)), rownames(ct2))
    expect_identical(names(colToRowMatCols(ct2)), colnames(ct2))
})

bam_cfg_f = system.file("extdata/bam_config.csv", package = "chiptsne2", mustWork = TRUE)
fetch_config = FetchConfig.load_config(bam_cfg_f)
fetch_config@meta_data = fetch_config@meta_data[1:2,]
query_gr = exampleQueryGR()[1:10]

suppressWarnings({
    ct2.cfg = ChIPtsne2.from_FetchConfig(fetch_config, query_gr)
})

test_that("Constructor FetchConfig", {
    expect_setequal(rownames(rowToRowMat(ct2.cfg)), names(query_gr))
    expect_equal(ncol(rowToRowMat(ct2.cfg)), 400)
})

