library(seqsetvis)
library(magrittr)
library(chiptsne2)

query_gr = CTCF_in_10a_overlaps_gr
prof_dt = CTCF_in_10a_profiles_dt
meta_dt = unique(prof_dt[, .(sample)])
meta_dt[, c("cell", "mark") := data.table::tstrsplit(sample, "_")]

map_dt = unique(prof_dt[, .(x, sample)])
map_dt[, cn := paste(sample, x, sep = "_")]
cn = map_dt$cn
setequal(cn, colnames(prof_mat))
cn == colnames(prof_mat)
map_dt$cn
map_dt[, nr := seq(.N)]
map_list = split(map_dt$n, map_dt$sample)

tmp_wide = tidyr::pivot_wider(prof_dt, names_from = c("sample", "x"), values_from = "y", id_cols = "id")
prof_mat = as.matrix(tmp_wide[, -1])
rownames(prof_mat) = tmp_wide$id


prof_max = prof_dt[, .(y = max(y)), .(id, sample)] %>%
    tidyr::pivot_wider(names_from = "sample", id_cols = "id", values_from = "y")
prof_max_mat = as.matrix(prof_max[, -1])
rownames(prof_max_mat) = prof_max$id




ct2 = ChIPtsne2.from_tidy(prof_dt, query_gr)
ct2[1:5,]

ct = ChIPtsne2(assay = list(max = prof_max_mat),
               rowRanges = query_gr,
               rowToRowMat = prof_mat,
               colToRowMatCols = map_list,
               colData = meta_dt,
               metadata = list(time = date()))
library(SummarizedExperiment)
metadata(ct)
assays(ct)[[1]]
colData(ct)
rowData(ct)
rowRanges(ct)
se = SummarizedExperiment(assay = list(max = prof_max_mat),
                          rowRanges = query_gr,
                          colData = meta_dt)
rowRanges(se)
rowData(se)
se.sub = subsetByOverlaps(se, GRanges("chr1", IRanges(0, 50e6)))
dim(se.sub)
rowRanges(se.sub)
rowRanges(se)
rowRanges(ct)
rowData(ct)

ct
ct[1:5,]
ct[, 1:2]
ct[1:2, 3]

library(GenomicRanges)
subsetByOverlaps(ct, GRanges("chr1", IRanges(0, 50e6)))
findOverlaps(ct, GRanges("chr1", IRanges(0, 50e6)))
subsetByOverlaps(query_gr, GRanges("chr1", IRanges(0, 50e6)))

tidy_profiles = function(ct){
    ct@rowToRowMat
    ct@colToRowMatCols
    i = ct@colToRowMatCols[[1]]
    nam = names(ct@colToRowMatCols)[1]
    df = do.call(rbind,
                 lapply(names(ct@colToRowMatCols), function(nam){
                     i = ct@colToRowMatCols[[nam]]
                     df = reshape2::melt(ct@rowToRowMat[, i])
                     # df$Var1 = NULL
                     df$Var2 = factor(df$Var2)
                     # df$Var2 = as.numeric(sapply(strsplit(levels(df$Var2), "_"), function(x)x[length(x)]))
                     levels(df$Var2) = as.numeric(sub(paste0(nam, "_"), "", levels(df$Var2)))
                     df$name = nam
                     df
                 })
    )
    colnames(df) = c("id", "x", "y", "name")
    df
}

ct$cell
tidy_profiles(ct %>% subset(cell == "MCF10A"))
tidy_profiles(ct %>% filter(cell == "MCF10A"))

object.size(df)
object.size(ct@rowToRowMat)

query_gr$tx = ""
query_gr[rownames(xy_dt)]$tx = xy_dt$tx
query_gr$ty = ""
query_gr[rownames(xy_dt)]$ty = xy_dt$ty

dim(prof_mat)

se = SummarizedExperiment::SummarizedExperiment(colData = meta_dt, rowRanges = query_gr, assays = list(max = prof_max_mat))
se
dim(se)
se[1:3,1]
