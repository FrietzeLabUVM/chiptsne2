library(chiptsne2)
library(tidyverse)
bam_file_dt = data.BRCAprogression::chip.setup_bam_files(pooled = TRUE)
bam_file_dt = bam_file_dt %>% filter(cell %in% c("MCF10A", "MCF10AT1") & mark %in% c("H3K4me3", "H3K27ac", "H3K4me1", "input"))
bam_file_dt = bam_file_dt[!grepl("rep", file)]

np_file_dt = data.BRCAprogression::chip.setup_narrowPeak_files()
np_file_dt = np_file_dt %>% filter(cell %in% c("MCF10A", "MCF10AT1") & mark %in% c("H3K4me3", "H3K27ac", "H3K4me1"))
np_grs = seqsetvis::easyLoad_narrowPeak(np_file_dt$file, np_file_dt$group)
cols = seqsetvis::safeBrew(as.character(np_file_dt$mark))[as.character(np_file_dt$mark)]
seqsetvis::ssvFeatureBars(np_grs, bar_colors = cols)
olap_gr = seqsetvis::ssvConsensusIntervalSets(np_grs, min_number = 2, min_fraction = 0)
seqsetvis::ssvFeatureBinaryHeatmap(olap_gr)

query_gr = sampleCap(olap_gr, 500)
query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)


bam_file_dt$fragLens = 200
# bam_file_dt$fragLens = sapply(bam_file_dt$file, function(f){
#     message(f)
#     seqsetvis::fragLen_calcStranded(f, qgr = query_gr)
# })

fetch_config = FetchConfig(bam_file_dt, name_VAR = "group_split", read_mode = "bam_SE")
options(mc.cores = 20)
# debug(ChIPtsne2.from_tidy)
ct2 = ChIPtsne2.from_FetchConfig(fetch_config, query_gr)
showMethods(normalizeSignalCapValue)

ct2 %>% dimReduceTSNE()

ct2.final = ct2 %>%
    normalizeSignalRPM() %>%
    calculateSignalCapValue %>%
    normalizeSignalCapValue %>%
    centerProfilesAndRefetch() %>%
    dimReduceTSNE()

ct2.grouped = ct2.final %>%
    groupRegionBySignalCluster(group_VAR = "cluster_id_4", n_clusters = 4) %>%
    groupRegionBySignalCluster(group_VAR = "cluster_id_6", n_clusters = 6) %>%
    groupRegionByDimReduceCluster(group_VAR = "knn_50", nearest_neighbors = 50)


man_df = data.frame(id = rownames(ct2.final))
man_df$group_id = sample(c("A", "B", "C"), nrow(man_df), replace = TRUE)
man_df$random_group_id = sample(c("A", "B", "C"), nrow(man_df), replace = TRUE)
ct2.grouped2 = ct2.final %>%
    groupRegionManually(assignment = man_df) %>%
    groupRegionManually(assignment = man_df, group_VAR = "random_group_id")

rowRanges(ct2.final)
rowRanges(ct2.grouped2)

olap_gr.k4me3 = query_gr[, c("MCF10A H3K4me3", "MCF10AT1 H3K4me3")]
olap_gr.10a_enh = query_gr[, c("MCF10A H3K4me1", "MCF10A H3K27ac")]

# debug(groupRegionByMembershipTable, "ChIPtsne2")
ct2.grouped3 = ct2.final %>%
    groupRegionByMembershipTable(membership = olap_gr.k4me3, group_VAR = "k4me3_overlap") %>%
    groupRegionByMembershipTable(membership = olap_gr.10a_enh, group_VAR = "10a_enhancer")
rowRanges(ct2.grouped3)

