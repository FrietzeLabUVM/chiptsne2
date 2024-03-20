# Haley Greenyer
# Hi Joe, sure, that would be great!
#     The DESeq outputs are here:
#     /slipstream_old/home/hgreenyer/PROJECTS/MCF10A_RUNX1_TIMECOURSE_2023/MCF10A_TimeCourse_RNAseq/DEA_output/
#     both of the 6D samples are in there, the new ones have 6Dnew
# The 27ac ChIP peaks are here:
#     /slipstream_old/home/hgreenyer/PROJECTS/MCF10A_RUNX1_TIMECOURSE_2023/MCF10A_TimeCourse_ChIPseq/MCF10A_TimeCourse_H3K27ac_ChIP/NARROWPEAKS/
#     We also have 4me3 ChIP:
#     /slipstream_old/home/hgreenyer/PROJECTS/MCF10A_RUNX1_TIMECOURSE_2023/MCF10A_TimeCourse_ChIPseq/MCF10A_TimeCourse_H3K4me3_ChIP/NARROWPEAKS/
library(tidyverse)
library(GenomicRanges)
library(chiptsne2)
# we need a data.frame specifying file paths and metadata

#### configure bam/bw files ####
bam_files.k27 = dir("/slipstream_old/home/hgreenyer/PROJECTS/MCF10A_RUNX1_TIMECOURSE_2023/MCF10A_TimeCourse_ChIPseq/MCF10A_TimeCourse_H3K27me3_ChIP/PIPELINE_OUT/", pattern = ".+pool.+.bam$", full.names = TRUE)
bam_dirs = dir("/slipstream_old/home/hgreenyer/PROJECTS/MCF10A_RUNX1_TIMECOURSE_2023/MCF10A_TimeCourse_ChIPseq", pattern = "ChIP", full.names = TRUE)
bam_files.main = dir(file.path(bam_dirs, "BAMS"), full.names = TRUE, pattern = "pool.+bam$")
bam_files = c(bam_files.k27, bam_files.main)
cfg_df = data.frame(file = bam_files)
cfg_df = cfg_df %>% mutate(name = sub("\\..+", "", basename(file)))

cfg_df = cfg_df %>% separate(name, into = c("V1", "V2", "V3", "V4", "V5"), sep = "_", remove = FALSE)
cfg_df = cfg_df %>% mutate(mark = sub("_ChIP", "", sub("MCF10A_TimeCourse_", "", basename(dirname(dirname(file))))))

cfg_df = cfg_df %>% mutate(treatment = ifelse(V3 %in% c("DMSO", "DTAG"), V3, V2))
cfg_df = cfg_df %>% mutate(time = ifelse(V4 %in% c("pooled"), V3, V4))
cfg_df = cfg_df %>% mutate(cell = V1)
cfg_df = subset(cfg_df, treatment != "input")
cfg_df = select(cfg_df, file, cell, treatment, time, mark)
cfg_df = cfg_df %>% mutate(cell = ifelse(cell == "C7", "MCF10A", cell))
cfg_df = cfg_df %>% mutate(time = ifelse(time == "6Days", "6D", time))


table(cfg_df$cell)
table(cfg_df$mark)
table(cfg_df$treatment)
table(cfg_df$time)

cfg_df = cfg_df %>% mutate(name = paste(cell, treatment, time, mark, sep = "_"))

cfg_df$time = factor(cfg_df$time)
cfg_df$time = relevel(cfg_df$time, "3H")
cfg_df = cfg_df[order(cfg_df$time),]
cfg_df = cfg_df[order(cfg_df$mark),]
cfg_df$name = factor(cfg_df$name, levels = cfg_df$name)


qc_df = cfg_df
qc_df$mapped_reads = seqsetvis::get_mapped_reads(qc_df$file)


ggplot(qc_df, aes(x = name, y = mapped_reads, fill = mark)) +
    geom_bar(stat = "identity") +
    theme(axis.text = element_text(angle = 45, hjust = 1, vjust = 1))

#### load DE ####
# query granges
de_files = dir("/slipstream_old/home/hgreenyer/PROJECTS/MCF10A_RUNX1_TIMECOURSE_2023/MCF10A_TimeCourse_RNAseq/DEA_output/", pattern = "DESeqRes", full.names = TRUE)
de_cfg_df = data.frame(file = de_files)
de_cfg_df = de_cfg_df %>% mutate(name = sub("_DES.+", "", basename(file)))
de_cfg_df = de_cfg_df %>% separate(name, into = c("cell", "pair", "time"), sep = "_")
de_cfg_df = subset(de_cfg_df, !time %in% c("48H", "6D"))
de_cfg_df = de_cfg_df %>% mutate(name = paste(cell, pair, time, sep = '_'))
de_files = de_cfg_df$file
names(de_files) = de_cfg_df$name

# i want stats for all genes that are sig in any comparison
# load once with sig filters then load again with gene id filter
load_de_f = function(f){
    des_dt = data.table::fread(f)
    data.table::setnames(des_dt, "name", "gene_name")
    data.table::setnames(des_dt, "Gene", "gene_id")
    subset(des_dt, baseMean > 50 & padj < .01)

}


de_dt.l = lapply(de_files, load_de_f)
de_dt = data.table::rbindlist(de_dt.l, idcol = "name")

gene_keep = de_dt$gene_id

load_de_f2 = function(f){
    des_dt = data.table::fread(f)
    data.table::setnames(des_dt, "name", "gene_name")
    data.table::setnames(des_dt, "Gene", "gene_id")
    subset(des_dt, gene_id %in% gene_keep)

}

de_dt.l = lapply(de_files, load_de_f2)
de_dt = data.table::rbindlist(de_dt.l, idcol = "name")

de_dt[is.na(padj), padj := 1]
de_dt$padj %>% range
de_dt[, is_DE := padj < .01 & abs(log2FoldChange) > 1]
de_dt[, negLogpadj := -log10(padj)]

table(de_dt$is_DE)
gene_dt = data.table::dcast(de_dt, gene_id+gene_name~name, value.var = c("log2FoldChange", "negLogpadj", "is_DE"))

ref_gr = rtracklayer::import.gff("/slipstream/home/joeboyd/lab_shared/indexes/HG38/GTF/gencode.v36.annotation.gtf", feature.type = "gene")
names(ref_gr) = sub("\\..+", "", ref_gr$gene_id)

# some gene names are missing, this will fix
gene_dt[gene_name == "None"]
stopifnot(all(gene_dt$gene_id %in% names(ref_gr)))
gene_dt$gene_name = ref_gr[gene_dt$gene_id]$gene_name

#### make query gr ####
# i'm not going to get into promoter selection for different transcript etc.
# should be considered in final anlaysis
pos_dt = data.table::as.data.table(ref_gr)[, .(gene_name, seqnames, start, end, strand)]
query_dt = merge(gene_dt, pos_dt, by = "gene_name")
query_gr = GRanges(query_dt)
#center on promoters
query_gr = promoters(query_gr, upstream = 3e3, downstream = 3e3)

#### FetchConfig ####
#just a guess at fragLens to avoid slow calculation
cfg_df$fragLens = 180
fetch_cfg = FetchConfig(cfg_df, view_size = 6e3, name_VAR = "name", read_mode = "bam_SE")

#### ChIPtsne2 ####
options(mc.cores = 20)
ct2.raw = ChIPtsne2.from_FetchConfig(fetch_config = fetch_cfg, query_gr = query_gr)
ct2.raw = ct2.raw %>%
    setSeed(0) %>%
    setPositionVariable("bp from TSS") %>%
    setValueVariable("read pileup")

plotSignalHeatmap(ct2.raw)
plotSignalLinePlot(ct2.raw)

#### normalization ####

ct2 = ct2.raw %>%
    calculateSignalCapValue(cap_quantile = .95) %>%
    normalizeSignalCapValue() %>%
    setValueVariable("normalized read pileup")

ct2 = groupRegionsBySignalCluster(ct2)

plotSignalHeatmap(ct2, group_VARS = "cluster_id", relative_heatmap_width = .8, relative_heatmap_height = .8)
plotSignalLinePlot(ct2, color_VAR = "mark", group_VAR = "cluster_id")

# I'm having trouble normalizing some samples
tmp = subsetSamples(ct2.raw, mark == "H3K27ac" & treatment == "DTAG")
tmp = groupRegionsBySignalCluster(tmp, n_clusters = 3)
plotSignalHeatmap(tmp, group_VARS = "cluster_id")

tmpn = tmp %>%
    calculateSignalCapValue(cap_quantile = .5) %>%
    normalizeSignalCapValue()

p1 = plotSignalHeatmap(tmpn, group_VARS = "cluster_id",
                       relative_heatmap_width = .7, relative_heatmap_height = .8,
                       heatmap_format_FUN = function(p)p + labs(title = "6D: this looks like overamp to me\ncan't uniformly normalize") )

p2 = plotSignalLinePlot(tmpn,  group_VAR = "cluster_id", facet_VAR = "mark", color_VAR = "time")
cowplot::plot_grid(p1, p2)

p1 = plotSignalHeatmap(tmp, group_VARS = "cluster_id",
                       relative_heatmap_width = .7, relative_heatmap_height = .8,
                       heatmap_format_FUN = function(p)p + labs(title = "6D: pre normalization\nmassive signal at high signal regions") )
p2 = plotSignalLinePlot(tmp,  group_VAR = "cluster_id", facet_VAR = "mark", color_VAR = "time")
cowplot::plot_grid(p1, p2)

#### run tsne ####
ct2 = dimReduceTSNE(ct2)


#### chipseq plots ####
plotDimReduceBins(subsetSamples(ct2, mark == "H3K27ac"), x_bins = 20, extra_VARS = c("treatment", "time", "mark")) +
    facet_grid(treatment~time)

plotDimReduceBins(subsetSamples(ct2, mark == "H3K4me3"), x_bins = 20, extra_VARS = c("treatment", "time", "mark")) +
    facet_grid(treatment~time)

plotDimReduceBins(subsetSamples(ct2, mark == "H3K27me3"), x_bins = 20, extra_VARS = c("treatment", "time", "mark")) +
    facet_grid(treatment~time)


#### rnaseq plots ####
rowData(ct2)

plotDimReducePoints(ct2, color_VAR = c("is_DE_MCF10a_DMSOvDTAG_3H", "is_DE_MCF10a_DMSOvDTAG_24H", "is_DE_MCF10a_DMSOvDTAG_6Dnew"), point_size = .6)

plotDimReducePoints(ct2, color_VAR = c("log2FoldChange_MCF10a_DMSOvDTAG_3H", "log2FoldChange_MCF10a_DMSOvDTAG_24H", "log2FoldChange_MCF10a_DMSOvDTAG_6Dnew"), point_size = .6, point_color_limits = c(-2, 2))
plotDimReducePoints(ct2, color_VAR = c("negLogpadj_MCF10a_DMSOvDTAG_3H", "negLogpadj_MCF10a_DMSOvDTAG_24H", "negLogpadj_MCF10a_DMSOvDTAG_6Dnew"), point_size = .6, point_color_limits = c(0, 10))

#### differential ####
colData(ct2)
ct2 = mutateSamples(ct2, mark_time = paste(mark, time))
colData(ct2)$mark_time = factor(colData(ct2)$mark_time, levels = unique(colData(ct2)$mark_time))

ct2.by_treatment = split(ct2, "treatment")
colData(ct2.by_treatment$DMSO)$mark_time
colData(ct2.by_treatment$DTAG)$mark_time

#this is a critical step so that colnames are equal when subtracting.
#otherwise chiptsne cannot handle the data properly
ct2.by_treatment = swapNameVariable(ct2.by_treatment, "mark_time")
colnames(ct2.by_treatment$DMSO)
colnames(ct2.by_treatment$DTAG)
ct2.diff = ct2.by_treatment$DTAG - ct2.by_treatment$DMSO
colData(ct2.diff)
plotDimReducePoints(ct2.diff, point_color_limits = c(-.2, .2))

#### tsne on diff ####
ct2.diff = dimReduceTSNE(ct2.diff)
plotDimReducePoints(ct2.diff, point_color_limits = c(-.2, .2))

plotDimReducePoints(ct2.diff, color_VAR = c("is_DE_MCF10a_DMSOvDTAG_3H", "is_DE_MCF10a_DMSOvDTAG_24H", "is_DE_MCF10a_DMSOvDTAG_6Dnew"), point_size = .6)
plotDimReducePoints(ct2.diff, color_VAR = c("log2FoldChange_MCF10a_DMSOvDTAG_3H", "log2FoldChange_MCF10a_DMSOvDTAG_24H", "log2FoldChange_MCF10a_DMSOvDTAG_6Dnew"), point_size = .6, point_color_limits = c(-2, 2))
plotDimReducePoints(ct2.diff, color_VAR = c("negLogpadj_MCF10a_DMSOvDTAG_3H", "negLogpadj_MCF10a_DMSOvDTAG_24H", "negLogpadj_MCF10a_DMSOvDTAG_6Dnew"), point_size = .6, point_color_limits = c(0, 10))

hist(ct2.diff@rowToRowMat)
plotDimReduceSummaryProfiles(ct2.diff, color_VAR = "time", extra_VARS = c("mark"), value_limits = c(-.1, .1)) + facet_wrap(~mark)
plotDimReduceSummaryProfiles(ct2, color_VAR = "time", extra_VARS = c("mark", "treatment"), value_limits = c(0, .5)) + facet_wrap(treatment~mark)

plotDimReduceSummaryProfiles(subsetSamples(ct2, time == "3H"),
                             color_VAR = "treatment", extra_VARS = c("mark", "time"),
                             value_limits = c(0, .5)) +
    facet_wrap(time~mark)

plotDimReduceSummaryProfiles(subsetSamples(ct2, time == "6D"),
                             color_VAR = "treatment", extra_VARS = c("mark", "time"),
                             value_limits = c(0, .5)) +
    facet_wrap(time~mark)

plotDimReduceSummaryProfiles(ct2.diff,
                             color_VAR = "mark", extra_VARS = c("mark", "time"),
                             value_limits = c(-.5, .5)) +
    facet_wrap(time~.)
