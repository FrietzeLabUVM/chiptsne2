library(chiptsne2)
library(tidyverse)
bam_file_dt = data.BRCAprogression::chip.setup_bam_files()
bam_file_dt = bam_file_dt %>% filter(cell %in% c("MCF10A", "MCF10AT1") & mark %in% c("H3K4me3", "H3K27ac", "H3K4me1", "input"))

np_file_dt = data.BRCAprogression::chip.setup_narrowPeak_files()
np_file_dt = np_file_dt %>% filter(cell %in% c("MCF10A", "MCF10AT1") & mark %in% c("H3K4me3", "H3K27ac", "H3K4me1"))
np_grs = seqsetvis::easyLoad_narrowPeak(np_file_dt$file, np_file_dt$group)
cols = seqsetvis::safeBrew(as.character(np_file_dt$mark))[as.character(np_file_dt$mark)]
seqsetvis::ssvFeatureBars(np_grs, bar_colors = cols)
olap_gr = seqsetvis::ssvConsensusIntervalSets(np_grs, min_number = 2, min_fraction = 0)
seqsetvis::ssvFeatureBinaryHeatmap(olap_gr)

query_gr = sampleCap(olap_gr, 500)


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

ct2$cell
# subset(ct2, cell == "MCF10A")
# dplyr::filter(ct2, cell == "MCF10A")
