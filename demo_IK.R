library(IKdata)
library(chiptsne2)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
library(tidyverse)
options(mc.cores = 20)
cache_dir = "demo_IK_cache.3x_bg"
cache_file = function(f){
    out_f = file.path(cache_dir, f)
    dir.create(dirname(out_f), showWarnings = FALSE, recursive = TRUE)
    out_f
}
out_dir = "demo_IK_output.3x_bg"
out_file = function(f){
    out_f = file.path(out_dir, f)
    dir.create(dirname(out_f), showWarnings = FALSE, recursive = TRUE)
    out_f
}
treatment_colors = c("IK" = "purple", GFP = "green")

#### setup bams ####
bam_cfg_dt = cutnrun.setup_bam_files()
bam_cfg_dt = bam_cfg_dt[mark == "ERG"]
bam_cfg_dt = bam_cfg_dt[order(treatment)]
bam_cfg_dt$name = factor(bam_cfg_dt$name, levels = unique(bam_cfg_dt$name))
bam_cfg_dt$name = droplevels(bam_cfg_dt$name)

#### setup peaks ####
np_cfg_dt = cutnrun.setup_peak_files()
np_cfg_dt = np_cfg_dt[mark == "ERG"]

#### setup csaw ####
csaw_file = "~/R_workspace.combined/preB_IK_induction/output_csaw.091923/ERG_cutnrun/IK_vs_GFP_in_all.csaw_regions.csv"
csaw_dt = fread(csaw_file)
csaw_dt = filter(csaw_dt, PValue < .05)
csaw_dt = filter(csaw_dt, direction != "mixed")
csaw_dt = csaw_dt %>% mutate(up_treatment = ifelse(direction == "up", "IK", "GFP"))
csaw_gr = GRanges(csaw_dt)
csaw_gr$original_width = width(csaw_gr)

#### fetch config ####
cfg = FetchConfig(bam_cfg_dt, read_mode = "bam_PE")

#### load and overlap peaks ####
peak_grs = easyLoad_narrowPeak(np_cfg_dt$file, file_names = np_cfg_dt$name)
olap_grs = ssvConsensusIntervalSets(peak_grs, min_number = 2, min_fraction = 0)
olap_grs = prepare_fetch_GRanges_names(olap_grs)

#### create query GRanges ####
olap_grs.csaw_ik = subsetByOverlaps(olap_grs, subset(csaw_gr, up_treatment == "IK"))
olap_grs.csaw_gfp = subsetByOverlaps(olap_grs, subset(csaw_gr, up_treatment == "GFP"))

olap_grs.csaw_ik$csaw_group = "IK"
olap_grs.csaw_gfp$csaw_group = "GFP"

diff_gr = c(
    olap_grs.csaw_ik,
    olap_grs.csaw_gfp
)

#invert selects regions NOT overlapping
notDiff_gr = subsetByOverlaps(olap_grs, diff_gr, invert = TRUE)
set.seed(0)
notDiff_gr = sample(notDiff_gr, length(diff_gr)*3)
notDiff_gr$csaw_group = "not_diff"

query_gr = c(diff_gr, notDiff_gr)

#### run fetch ####
raw_cache = cache_file("ct2.raw.Rds")
if(file.exists(raw_cache)){
    ct2.raw = readRDS(raw_cache)
}else{
    ct2.raw = ChIPtsne2.from_FetchConfig(cfg, query_gr)
    ct2.raw = setPositionVariable(ct2.raw, "bp")
    ct2.raw = setRegionVariable(ct2.raw, "ERG peak")
    saveRDS(ct2.raw, raw_cache)
}

tsne_cache = cache_file("ct2.Rds")
if(file.exists(tsne_cache)){
    ct2 = readRDS(tsne_cache)
}else{
    ct2 = ct2.raw %>%
        setSeed(0) %>%
        normalizeSignalRPM() %>%
        calculateSignalCapValue() %>%
        normalizeSignalCapValue() %>%
        centerProfilesAndRefetch() %>%
        normalizeSignalRPM() %>%
        calculateSignalCapValue() %>%
        normalizeSignalCapValue() %>%
        flipProfilesToMatch() %>%
        groupRegionsBySignalCluster() %>%
        dimReduceTSNE() %>%
        groupRegionsByDimReduceCluster()
    ct2 = setValueVariable(ct2, "normalized reads")
    saveRDS(ct2, tsne_cache)
}

prof_dt = getTidyProfile(ct2)
range(prof_dt$`normalized reads`)

colData(ct2)
rowRanges(ct2)
# undebug(plotDimReduceSummaryProfiles)

theme_update(
    panel.grid = element_blank(),
    panel.background = element_blank()
)

p_summary = plotDimReduceSummaryProfiles(
    ct2,
    color_VAR = "treatment",
    extra_vars = c("cell", "csaw_group"),
    x_bins = 8,
    value_limits = c(0, 1)
)
p_summary = p_summary +
    facet_wrap(~cell) +
    scale_color_manual(values = treatment_colors) +
    labs(title = "ERG Cut&Run")

ggsave(out_file("fig1_summary.png"), p_summary, width = 9.5, height = 4.5)


plotDimReducePoints(ct2, extra_VARS = c("cell", "mark", "treatment")) +
    facet_grid(cell~treatment)

ct2.by_csaw = split(ct2, "csaw_group")
p_points.plots = lapply(ct2.by_csaw, function(ct2.sel){
    plotDimReducePoints(ct2.sel, extra_VARS = c("cell", "mark", "treatment", "csaw_group")) +
        facet_grid(cell~treatment)
})
p_points.leg = cowplot::get_legend(p_points.plots$GFP)
p_points.plots2 = lapply(names(p_points.plots), function(nam){
    p = p_points.plots[[nam]]
    p + guides(color = "none") + labs(title = nam)
})
pg_points = cowplot::plot_grid(plotlist = c(p_points.plots2, list(p_points.leg)))
ggsave(out_file("fig2_points.png"), pg_points, width = 7.5, height = 7.5)

ct2.by_name = split(ct2, "name")
ct2.gfp = (ct2.by_name$MXP5_ERG_GFP_rep1 + ct2.by_name$PDX2_ERG_GFP_rep1)/2
ct2.ik = (ct2.by_name$MXP5_ERG_IK_rep1 + ct2.by_name$PDX2_ERG_IK_rep1)/2
ct2.diff = ct2.ik - ct2.gfp
ct2.diff$new_name = "IK - GFP"
ct2.diff = setNameVariable(ct2.diff, "new_name")

plotDimReducePoints(ct2.diff, point_size = 1)
plotDimReducePoints(ct2.diff, point_size = 1)
p_summary
plotDimReduceSummaryProfiles(ct2.diff, value_limits = c(-.2, .2)) +
    guides(color = "none")
ct2.manual = run_ChIPtsne2_shiny(ct2.diff)

centroids = calculateGroupCentroid(ct2.manual, "manual_diff")
ct2.cent = groupRegionsByCentroidDistance(ct2.manual, centroid = centroids, "centroid_diff")

diff_colors = c("up" = "purple", down = "orange", bg = "gray")

plotDimReducePoints(ct2.cent, color_VAR = "manual_diff", point_colors = diff_colors)
plotDimReducePoints(ct2.cent, color_VAR = "centroid_diff", point_colors = diff_colors)

head(getRegionMetaData(ct2.manual))

olaps = findOverlaps(
    query = rowRanges(ct2.manual),
    subject = csaw_gr
)
ct2.manual = ct2.manual %>%
    chiptsne2::addRegionAnnotation(csaw_gr, c("PValue", "num.tests", "best.logFC", "original_width"))
meta_df.manual = getRegionMetaData(ct2.manual)
meta_df.manual$manual_diff = factor(meta_df.manual$manual_diff)
meta_df.manual$original_width = as.numeric(meta_df.manual$original_width)
meta_df.manual$best.logFC = as.numeric(meta_df.manual$best.logFC)
meta_df.manual$PValue = as.numeric(meta_df.manual$PValue)
meta_df.manual$num.tests = as.numeric(meta_df.manual$num.tests)

ggplot(meta_df.manual, aes(x = csaw_group, fill = manual_diff)) +
    geom_bar(position = position_dodge(preserve = "single"))

meta_df.manual = meta_df.manual %>% mutate(manual_matches = (csaw_group == "GFP" & manual_diff == "down") | (csaw_group == "IK" & manual_diff == "up"))

ggplot(subset(meta_df.manual, csaw_group != "not_diff"), aes(x = csaw_group, fill = manual_matches, y = original_width, group = paste(csaw_group, manual_matches))) +
    geom_boxplot()

ggplot(subset(meta_df.manual, csaw_group != "not_diff"), aes(x = csaw_group, fill = manual_matches, y = num.tests, group = paste(csaw_group, manual_matches))) +
    geom_boxplot()

ggplot(subset(meta_df.manual, csaw_group != "not_diff"), aes(x = csaw_group, fill = manual_matches, y = best.logFC, group = paste(csaw_group, manual_matches))) +
    geom_boxplot()

ggplot(subset(meta_df.manual, csaw_group != "not_diff"), aes(x = csaw_group, fill = manual_matches, y = PValue, group = paste(csaw_group, manual_matches))) +
    geom_boxplot()


meta_df.manual = meta_df.manual %>% mutate(csaw_group_no_narrow = ifelse(original_width < 1000, "not_diff", csaw_group))
meta_df.manual = meta_df.manual %>% mutate(csaw_group_no_narrow = ifelse(is.na(csaw_group_no_narrow), "not_diff", csaw_group))
ct2.re_csaw = setRegionMetaData(ct2.manual, meta_df.manual)

ct2.by_csaw = split(ct2.re_csaw, "csaw_group_no_narrow")

p_points.plots = lapply(ct2.by_csaw, function(ct2.sel){
    plotDimReducePoints(ct2.sel, extra_VARS = c("cell", "mark", "treatment", "csaw_group")) +
        facet_grid(cell~treatment)
})
p_points.leg = cowplot::get_legend(p_points.plots$GFP)
p_points.plots2 = lapply(names(p_points.plots), function(nam){
    p = p_points.plots[[nam]]
    p + guides(color = "none") + labs(title = nam) +
        coord_cartesian(xlim = c(-.5, .5), ylim = c(-.5, .5))
})
pg_points = cowplot::plot_grid(plotlist = c(p_points.plots2, list(p_points.leg)))
pg_points

rowRanges(ct2.manual)
centroids = calculateGroupCentroid(ct2.re_csaw, "manual_diff")
ct2.cent = groupRegionsByCentroidDistance(ct2.re_csaw, centroid = centroids, "centroid_diff")
ct2.cent

ct2.by_csaw = split(ct2.cent, "centroid_diff")

p_points.plots = lapply(ct2.by_csaw, function(ct2.sel){
    plotDimReducePoints(ct2.sel, extra_VARS = c("cell", "mark", "treatment", "csaw_group"), point_color_limits = c(-.5, .5), has_symmetrical_limits = TRUE) +
        facet_grid(cell~treatment)
})
p_points.leg = cowplot::get_legend(p_points.plots$bg)
p_points.plots2 = lapply(names(p_points.plots), function(nam){
    p = p_points.plots[[nam]]
    p + guides(color = "none") + labs(title = nam) +
        coord_cartesian(xlim = c(-.5, .5), ylim = c(-.5, .5))
})
pg_points = cowplot::plot_grid(plotlist = c(p_points.plots2, list(p_points.leg)))
pg_points

rowRanges(ct2.cent)

plotSignalHeatmap(subsetRow(ct2.cent, centroid_diff != "bg"),
                  group_VARS = c("csaw_group",
                                 "manual_diff",
                                 "csaw_group_no_narrow",
                                 "centroid_diff"),
                  annotation_colors = list(c("GFP" = "orange", "IK" = "purple", not_diff = "gray"),
                                        c("up" = "purple", down = "orange", bg = "gray"),
                                        c("GFP" = "orange", "IK" = "purple", not_diff = "gray"),
                                        c(up = "purple", down = "orange")),
                  color_key_strategy = "all",
                  name_FUN = function(x){"IK - GFP"})

plotSignalHeatmap(subsetRow(ct2.cent, csaw_group_no_narrow != "not_diff"),
                  group_VARS = c("csaw_group",
                                 "manual_diff",
                                 "centroid_diff",
                                 "csaw_group_no_narrow"),
                  annotation_colors = list(c("GFP" = "orange", "IK" = "purple", not_diff = "gray"),
                                           c(up = "purple", down = "orange", bg = "gray"),
                                           c("up" = "purple", down = "orange", bg = "gray"),
                                           c("GFP" = "orange", "IK" = "purple")),
                  color_key_strategy = "all",
                  name_FUN = function(x){"IK - GFP"})
