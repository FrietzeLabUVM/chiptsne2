#### ChIP_tSNE for Histone4 acetylation
suppressPackageStartupMessages({
  library(data.BRCAprogression)
  library(chiptsne2)
  library(seqsetvis)
  library(GenomicRanges)
  library(data.table)
})

options(mc.cores = 20)
treatment_colors = c("H3K27ac" = "purple", "H4K5ac" = "green", "H4K8ac" = "blue")

###  data.BRCAprogression provides convenient shortcuts to setup bam and peak files with metadata.
bam_df.chip = chip.setup_bam_files()
bam_cfg_dtt = bam_df.chip
bam_cfg_dt1 = bam_cfg_dtt[mark == "H3K27ac"]
bam_cfg_dt2 = bam_cfg_dtt[mark == "H4K5ac" | mark == "H4K8ac"]
bam_cfg_dt <- rbind(bam_cfg_dt1, bam_cfg_dt2)
bam_cfg_dt = bam_cfg_dt[order(mark)] #treatment -> mark
bam_cfg_dt$name = factor(bam_cfg_dt$name, levels = unique(bam_cfg_dt$name))
bam_cfg_dt$name = droplevels(bam_cfg_dt$name)

### Load in query gr, which is the set of TE+SE in MCF10A
## need to filter on TE_grl
print(load("/slipstream_old/home/conggao/BRCA_progression/Chromtin_States/Enhancer_grl_20240115.Rdata",))

w=TE_grl$MCF10A_TE
TE_MCF10A <- w[width(w) >180]
mcols(TE_MCF10A)=NULL
SE_MCF10A=SE_grl$MCF10A_SE
olaps_grs <- ssvOverlapIntervalSets(GRangesList(TE_MCF10A, SE_MCF10A))
colnames(mcols(olaps_grs)) <- c("TE_MCF10A", "SE_MCF10A")
ssvFeatureEuler(olaps_grs)
# save(bam_cfg_dt, olaps_grs, file="./chiptsne_test_20241124.Rdata")
###
query_gr = olaps_grs
cfg = FetchConfig(bam_cfg_dt)
ct2.raw = ChIPtsne2.from_FetchConfig(cfg, query_gr)
ct2.raw = setPositionVariable(ct2.raw, "bp")

### to get cap_dt
# 1.
src_dir = "/slipstream_old/home/conggao/BRCA_progression/PTM_ChIP"
bam_files = dir(src_dir, pattern = "pooled.+Aligned.sortedByCoord.out.bam$", full.names = TRUE)
bam_dt = data.table(file = bam_files)
bam_dt[, c("cell", "mark", "rep") := tstrsplit(basename(file), split = "[_\\.]", keep = c(1,2,3))]
bam_dt[, name := paste(cell, mark, rep, sep = "_")]
bam_dt[, name_split := paste(cell, mark, rep, sep = "\n")]
bam_dt$method <- "chipseq"
bam_dt$treatment <- "native"

# bam_dt$bw  <- sub('out.bam','out_macs2_FE.bw',bam_dt$file)

todo_mark <- c("H3K27ac", "H4K5ac", "H4K8ac")
bam_dt$peak <- bam_dt$file
todo_dt <- bam_dt[is.element(bam_dt$mark, todo_mark),]
todo_dt$peak <- sub('out.bam','out_macs2_peaks.narrowPeak',todo_dt$file)
# 2.
cons_peaks =lapply(split(todo_dt, c("H3K27ac", "H4K5ac", "H4K8ac")), function(x){
  peak_qgr.sel = easyLoad_narrowPeak(x$peak)
  olap_gr = ssvConsensusIntervalSets(peak_qgr.sel, min_number = 3, min_fraction = 0, maxgap = 100)
  olap_gr
})

cap_dt = lapply(seq_len(nrow(todo_dt)), function(i){
  peak_qgr.sel = cons_peaks[[todo_dt$mark[i]]]
  qgr = resize(ssvRecipes::sampleCap(peak_qgr.sel, 5e3), 3e3, fix = "center")
  prof_dt = ssvFetchBam(todo_dt[i,], qgr, return_data.table = TRUE,
                        fragLens = 200, n_region_splits = 50)
  calc_norm_factors(prof_dt, aggFUN2 = function(x)quantile(x, .9))
}) %>% rbindlist

cap_dt
ct_cap_dt = cap_dt

colnames(ct_cap_dt) = c(
  getNameVariable(ct2.raw),
  "cap_value"
)

ct2 = ct2.raw
signal_cap_data = ct_cap_dt

### normalize chiptsne objects
ct2 = ct2.raw %>%
  setSeed(0) %>%
  normalizeSignalCapValue(signal_cap_data = ct_cap_dt)
# %>%
  centerProfilesAndRefetch() %>%
  normalizeSignalCapValue(signal_cap_data = ct_cap_dt) %>%
  flipProfilesToMatch() %>%
  groupRegionsBySignalCluster() %>%
  dimReduceTSNE() %>%
  groupRegionsByDimReduceCluster() %>%
  setValueVariable("normalized reads")



### Clustering
ct2 = exampleChIPtsne2.with_meta()

ct2 = ct2 %>% dimReduceTSNE() %>% groupRegionsByDimReduceCluster()

plotDimReducePoints(ct2) +
  scale_color_viridis_c() +
  theme_classic()

rowData(ct2)
plotDimReducePoints(ct2, color_VAR = "knn_id")

plotDimReducePoints(ct2, extra_VARS = c("knn_id", "cell")) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_grid(cell~knn_id)
# setRegionVariable(ct2.raw, "ERG peak")
