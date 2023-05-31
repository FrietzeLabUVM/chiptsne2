testthat::context("FetchConfig_class")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(magrittr)

bam_cfg_f = system.file("extdata/bam_config.csv", package = "chiptsne2", mustWork = TRUE)
bam_cfg_f.args = system.file("extdata/bam_config.args.csv", package = "chiptsne2", mustWork = TRUE)
bam_cfg_f.no_fragLens = system.file("extdata/bam_config.no_CFG_fragLens.csv", package = "chiptsne2", mustWork = TRUE)
bam_cfg_f.bam_fragLens = system.file("extdata/bam_config.bam_fragLens.csv", package = "chiptsne2", mustWork = TRUE)
bw_cfg_f = system.file("extdata/bigwig_config.csv", package = "chiptsne2", mustWork = TRUE)
bam_files = system.file("extdata", package = "chiptsne2", mustWork = TRUE) %>% dir(pattern = "^MCF10A_.+bam$", full.names = TRUE)
bw_files = system.file("extdata", package = "chiptsne2", mustWork = TRUE) %>% dir(pattern = "^M.+bw$", full.names = TRUE)

df = data.frame(file = bam_files)
df.no_file = data.frame(bam_files)
df.bad = data.frame(file = "asdf")


test_that("FetchConfig", {
    cfg = FetchConfig(df)
    cfg$meta_data
    expect_setequal(colnames(cfg$meta_data), c("file", "name"))
    #respect name_VAR
    cfg.sample = FetchConfig(df, name_VAR = "sample")
    cfg.sample$meta_data
    expect_setequal(colnames(cfg.sample$meta_data), c("file", "sample"))

    #works with not colnames
    cfg.no_file = FetchConfig(df.no_file)
    cfg.no_file$meta_data
    expect_setequal(colnames(cfg.no_file$meta_data), c("file", "name"))
    expect_error(FetchConfig(df.bad), "Files specified in config do not exist:")
})

test_that("FetchConfig.null", {
    cfg.null = FetchConfig.null()
    expect_true(isFetchConfigNull(cfg.null))
})

test_that("FetchConfig.files", {
    cfg = FetchConfig.files(bam_files)
    expect_setequal(colnames(cfg$meta_data), c("file", "name"))
    expect_equal(cfg$read_mode, "bam_SE")
    expect_equal(cfg$view_size, 3000)
    expect_equal(cfg$window_size, 200)
    expect_equal(cfg$fetch_options, list())

    cfg.args = FetchConfig.files(bam_files, read_mode = "bam_PE", view_size = 500, window_size = 20, fetch_options = list(win_method = "summary"), name_VAR = "sample")
    expect_setequal(colnames(cfg.args$meta_data), c("file", "sample"))
    expect_equal(cfg.args$read_mode, "bam_PE")
    expect_equal(cfg.args$view_size, 500)
    expect_equal(cfg.args$window_size, 20)
    expect_equal(cfg.args$fetch_options, list(win_method = "summary"))

    bw_files.named = bw_files
    names(bw_files.named) = sub("_.+", "", basename(bw_files))
    cfg.bw = FetchConfig.files(bw_files)
    expect_equal(cfg.bw$read_mode, "bigwig")
    expect_true(all(grepl("random100.bw", cfg.bw$meta_data$name)))

    cfg.bw_named = FetchConfig.files(bw_files.named)
    expect_equal(cfg.bw_named$read_mode, "bigwig")
    expect_true(!any(grepl("random100.bw", cfg.bw_named$meta_data$name)))
})

test_that("FetchConfig.load_config", {
    cfg.bam = FetchConfig.load_config(bam_cfg_f)
    expect_setequal(colnames(cfg.bam$meta_data), c("file", "cell", "mark", "rep", "name"))
    expect_equal(names(cfg.bam$fetch_options), c("win_method", "summary_FUN", "fragLens"))

    cfg.args = FetchConfig.load_config(bam_cfg_f.args)
    cfg.args$meta_data
    expect_setequal(colnames(cfg.args$meta_data), c("file", "cell", "mark", "rep", "sample"))
    expect_equal(cfg.args$view_size, 2000)
    expect_equal(cfg.args$window_size, 25)
    expect_equal(cfg.args$read_mode, "bam_PE")
    expect_equal(names(cfg.args$fetch_options), c("win_method", "summary_FUN"))

    cfg.bw = FetchConfig.load_config(bw_cfg_f)
    expect_setequal(colnames(cfg.bw$meta_data), c("file", "cell", "mark", "rep", "name"))
    expect_equal(levels(cfg.bw$meta_data$name), c("MCF10A_CTCF_FE_random100.bw", "MCF10AT1_CTCF_FE_random100.bw", "MCF10CA1_CTCF_FE_random100.bw"))
})

test_that("FetchConfig.save_config", {
    cfg.bam = FetchConfig.load_config(bam_cfg_f)
    t_file = tempfile()
    FetchConfig.save_config(cfg.bam, file = t_file)
    cfg.loaded = FetchConfig.load_config(t_file)

    expect_equal(cfg.bam$meta_data,
                 cfg.loaded$meta_data)
    expect_equal(cfg.bam@is_null,
                 cfg.loaded@is_null)
    expect_equal(cfg.bam@view_size,
                 cfg.loaded@view_size)
    expect_equal(cfg.bam@window_size,
                 cfg.loaded@window_size)
    expect_equal(cfg.bam@read_mode,
                 cfg.loaded@read_mode)
    expect_equal(cfg.bam@name_VAR,
                 cfg.loaded@name_VAR)
    fop.o = cfg.bam@fetch_options
    fop.o = fop.o[setdiff(names(fop.o), "summary_FUN")]
    expect_equal(cfg.loaded@fetch_options, fop.o)
})

test_that("fragLens", {
    query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
    cfg.bam = FetchConfig.load_config(bam_cfg_f)
    cfg.bam@meta_data = cfg.bam$meta_data[1,]

    cfg.bam_fragLens = FetchConfig.load_config(bam_cfg_f.bam_fragLens)
    cfg.bam_fragLens@meta_data = cfg.bam_fragLens$meta_data[1,]

    cfg.bam_no_fragLens = FetchConfig.load_config(bam_cfg_f.no_fragLens)
    cfg.bam_no_fragLens@meta_data = cfg.bam_no_fragLens$meta_data[1,]

    res.bam = fetch_signal_at_features(cfg.bam, query_gr)
    res.bam_fragLens = fetch_signal_at_features(cfg.bam_fragLens, query_gr)
    res.bam_no_fragLens = fetch_signal_at_features(cfg.bam_no_fragLens, query_gr)

    # plot(
    #     res.bam$prof_dt$y,
    #     res.bam_fragLens$prof_dt$y, pch = 16, cex = .2
    # )
    # plot(
    #     res.bam$prof_dt$y,
    #     res.bam_no_fragLens$prof_dt$y, pch = 16, cex = .2
    # )
    expect_failure(expect_equal(res.bam$prof_dt$y, res.bam_fragLens$prof_dt$y))
    expect_failure(expect_equal(res.bam$prof_dt$y, res.bam_no_fragLens$prof_dt$y))
    expect_failure(expect_equal(res.bam_no_fragLens$prof_dt$y, res.bam_fragLens$prof_dt$y))
})
