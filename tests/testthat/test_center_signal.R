testthat::context("center signal")
# flipping viewGranges
library(chiptsne2)
library(testthat)
library(ggplot2)

cfg.summary = FetchConfig.load_config(exampleBamConfigFile())
cfg.sample = cfg.summary
cfg.sample$fetch_options$win_method = "sample"
cfg.sample$window_size = 50
cfg.sample$view_size = 600
qgr = exampleQueryGR()

ct2.summary = ChIPtsne2.from_FetchConfig(cfg.summary, qgr)
p.summary = plotSignalHeatmap(ct2.summary, heatmap_format_FUN = function(p) p + labs(title = "uncentered") + geom_vline(xintercept = 0), sort_strategy = "left")

ct2.summary.c = centerProfilesAndRefetch(ct2.summary, view_size = .02, use_cache = FALSE)
p.summary.c = plotSignalHeatmap(ct2.summary.c, heatmap_format_FUN = function(p) p + labs(title = "centered .1") + geom_vline(xintercept = 0), sort_strategy = "left")

ct2.summary.c2 = centerProfilesAndRefetch(ct2.summary, use_cache = FALSE)
p.summary.c2 = plotSignalHeatmap(ct2.summary.c2, heatmap_format_FUN = function(p) p + labs(title = "centered full") + geom_vline(xintercept = 0), sort_strategy = "left")

all(rowRanges(ct2.summary.c) == rowRanges(ct2.summary))
all(rowRanges(ct2.summary.c) == rowRanges(ct2.summary.c2))

rowRanges(ct2.summary)
rowRanges(ct2.summary.c)

dim(ct2.summary)
dim(ct2.summary.c)

ct2.sample = ChIPtsne2.from_FetchConfig(cfg.sample, qgr)
p.sample = plotSignalHeatmap(ct2.sample, heatmap_format_FUN = function(p) p + labs(title = "uncentered") + geom_vline(xintercept = 0), sort_strategy = "left")

ct2.sample.c = centerProfilesAndRefetch(ct2.sample, view_size = 50, use_cache = FALSE)
p.sample.c = plotSignalHeatmap(ct2.sample.c, heatmap_format_FUN = function(p) p + labs(title = "centered 50") + geom_vline(xintercept = 0), sort_strategy = "left")

ct2.sample.c2 = centerProfilesAndRefetch(ct2.sample, view_size = 500, use_cache = FALSE)
p.sample.c2 = plotSignalHeatmap(ct2.sample.c2, heatmap_format_FUN = function(p) p + labs(title = "centered 500") + geom_vline(xintercept = 0), sort_strategy = "left")

all(rowRanges(ct2.sample.c) == rowRanges(ct2.sample))
all(rowRanges(ct2.sample.c) == rowRanges(ct2.sample.c2))

rowRanges(ct2.summary)
#flipping profiles adds strand information
#first test summary
ct2.summary = sortRegions(ct2.summary, sort_strategy = "left")
ct2.summary.flip = flipProfilesToMatch(ct2.summary)
rowRanges(ct2.summary.flip)
p.summary.flip = plotSignalHeatmap(
    ct2.summary.flip,
    heatmap_format_FUN = function(p) p + labs(title = "stranded, uncentered") + geom_vline(xintercept = 0),
    sort_strategy = "none")
rowRanges(ct2.summary.flip)

rowRanges(ct2.summary.flip)
ct2.summary.flip.c = centerProfilesAndRefetch(ct2.summary.flip, view_size = 500, use_cache = FALSE)

p.summary.flip.c = plotSignalHeatmap(
    ct2.summary.flip.c,
    heatmap_format_FUN = function(p) p + labs(title = "stranded, centered 500") + geom_vline(xintercept = 0),
    sort_strategy = "none")
#next test sample
ct2.sample = sortRegions(ct2.sample, sort_strategy = "left")
ct2.sample.flip = flipProfilesToMatch(ct2.sample)
p.sample.flip = plotSignalHeatmap(
    ct2.sample.flip,
    heatmap_format_FUN = function(p) p + labs(title = "stranded, uncentered") + geom_vline(xintercept = 0),
    sort_strategy = "none")

rowRanges(ct2.sample.flip)
ct2.sample.flip.c = centerProfilesAndRefetch(ct2.sample.flip, view_size = 500, use_cache = FALSE)

p.sample.flip.c = plotSignalHeatmap(
    ct2.sample.flip.c,
    heatmap_format_FUN = function(p) p + labs(title = "stranded, centered 500") + geom_vline(xintercept = 0),
    sort_strategy = "none")

cowplot::plot_grid(p.summary, p.summary.c, p.summary.c2, nrow = 1)

cowplot::plot_grid(p.sample, p.sample.c, p.sample.c2, nrow = 1)

cowplot::plot_grid(p.summary, p.summary.flip, p.summary.flip.c, nrow = 1)

cowplot::plot_grid(p.sample, p.sample.flip, p.sample.flip.c, nrow = 1)

plotSignalHeatmap(
    ct2.sample[1:5,],
    heatmap_format_FUN = function(p) p + labs(title = "stranded, uncentered") + geom_vline(xintercept = 0) + theme(axis.text.y = element_text()),
    sort_strategy = "none")


plotSignalHeatmap(
    ct2.sample.flip[1:5,],
    heatmap_format_FUN = function(p) p + labs(title = "stranded, uncentered") + geom_vline(xintercept = 0) + theme(axis.text.y = element_text()),
    sort_strategy = "none")

rowRanges(ct2.sample.flip[1:5,])

flipProfilesToMatch(ct2.sample[1:5,])
