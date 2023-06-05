library(chiptsne2)
library(tidyverse)

#setup bam files, simple file paths are OK too
bam_file_dt = data.BRCAprogression::chip.setup_bam_files(pooled = TRUE)
bam_file_dt = bam_file_dt %>% filter(cell %in% c("MCF10A", "MCF10AT1") & mark %in% c("H3K4me3", "H3K27ac", "H3K4me1", "input"))
bam_file_dt = bam_file_dt[!grepl("rep", file)]
bam_file_dt$fragLens = 200

# setup query_gr
np_file_dt = data.BRCAprogression::chip.setup_narrowPeak_files()
np_file_dt = np_file_dt %>% filter(cell %in% c("MCF10A", "MCF10AT1") & mark %in% c("H3K4me3", "H3K27ac", "H3K4me1"))
np_grs = seqsetvis::easyLoad_narrowPeak(np_file_dt$file, np_file_dt$group)
cols = seqsetvis::safeBrew(as.character(np_file_dt$mark))[as.character(np_file_dt$mark)]
seqsetvis::ssvFeatureBars(np_grs, bar_colors = cols)
olap_gr = seqsetvis::ssvConsensusIntervalSets(np_grs, min_number = 2, min_fraction = 0)
seqsetvis::ssvFeatureBinaryHeatmap(olap_gr)

query_gr = sampleCap(olap_gr, 500)
query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)

#fetch config is not necessarily required, only for certain operations
fetch_config = FetchConfig(bam_file_dt, name_VAR = "group_split", read_mode = "bam_SE")
options(mc.cores = 20)
ct2 = ChIPtsne2.from_FetchConfig(fetch_config, query_gr)

#extra ways to assign groups
manual_df = data.frame(id = rownames(ct2))
manual_df$group_id = sample(c("A", "B", "C"), nrow(manual_df), replace = TRUE)
manual_df$random_group_id = sample(c("A", "B", "C"), nrow(manual_df), replace = TRUE)

olap_gr.k4me3 = query_gr[, c("MCF10A H3K4me3", "MCF10AT1 H3K4me3")]
olap_gr.10a_enh = query_gr[, c("MCF10A H3K4me1", "MCF10A H3K27ac")]

ct2.a = ct2 %>%
    setSeed(1) %>%
    normalizeSignalRPM() %>%
    calculateSignalCapValue
ct2.a = ct2.a %>% normalizeSignalCapValue
ct2.a %>% centerProfilesAndRefetch

colData(ct2.a)

# workflow in a tidyverse style
ct2.final = ct2 %>%
    setSeed(1) %>%
    normalizeSignalRPM() %>%
    calculateSignalCapValue %>%
    normalizeSignalCapValue(minimum_ceiling = 1) %>%
    centerProfilesAndRefetch() %>%
    normalizeSignalRPM() %>%
    calculateSignalCapValue %>%
    normalizeSignalCapValue(minimum_ceiling = 1) %>%
    dimReduceTSNE() %>%
    groupRegionBySignalCluster(group_VAR = "cluster_id_4", n_clusters = 4) %>%
    groupRegionBySignalCluster(group_VAR = "cluster_id_6", n_clusters = 6) %>%
    groupRegionByDimReduceCluster(group_VAR = "knn_50", nearest_neighbors = 50) %>%
    groupRegionManually(assignment = manual_df) %>%
    groupRegionManually(assignment = manual_df, group_VAR = "random_group_id") %>%
    groupRegionByMembershipTable(membership = olap_gr.k4me3, group_VAR = "k4me3_overlap") %>%
    groupRegionByMembershipTable(membership = olap_gr.10a_enh, group_VAR = "10a_enhancer")

colData(ct2.final)

debug(sortRegions, "ChIPtsne2")
sortRegions(ct2.final, "cluster_id_4")

# reuse workflow
ct2.redo = rerun_history(ct2, history_source = ct2.final)
my_hist= ChIPtsne2.history(ct2.final)
names(my_hist)
#history modification
my_hist$setSeed$ARG$seed = 42

ct2.mod = rerun_history(ct2, history_source = my_hist)

# results are the same
reg1 = getRegionMetaData(ct2.final)
reg2 = getRegionMetaData(ct2.redo)
reg3 = getRegionMetaData(ct2.mod)
stopifnot(all(reg1 == reg2))
stopifnot(!all(reg1$cluster_id_4 == reg3$cluster_id_4))
# still working on plots, but here's an example
rowRanges(ct2.final)
plotDimReducePoints(ct2.final) +
    scale_color_viridis_c(limits = c(0, 1)) +
    coord_fixed()
plotDimReducePoints(ct2.final, color_VAR = c("k4me3_overlap"), point_size = 1)
plotDimReducePoints(ct2.final, color_VAR = c("cluster_id_4", "cluster_id_6", "knn_50"))

# back to seqsetvis
rowRanges(ct2.final)
prof_dt = getTidyProfile(ct2.final, sample_meta_VARS = c("cell", "mark"))
seqsetvis::ssvSignalHeatmap(prof_dt, facet_ = "group_split")

