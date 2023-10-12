library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()
ct2 = dimReduceTSNE(ct2)
ct2.by_cell = split(ct2, "cell")
rowRanges(ct2)
ct2.by_peak = split(ct2, "peak_MCF10CA1_CTCF")
args = ct2.by_cell
#### cbind ####

# rowToRowMat = rowToRowMat,
# colToRowMatCols = colToRowMatCols,
# name_VAR = name_VAR,
# position_VAR = position_VAR,
# value_VAR = value_VAR,
# region_VAR = region_VAR,
# fetch_config = fetch_config


# undebug(cbind, signature = "ChIPtsne2")
# undebug(cbind)
# debug(rbind, signature = "ChIPtsne2")
# undebug(cbind, signature = "ChIPtsne2")
# rbind = SummarizedExperiment::rbind
# debug(rbind, signature = "ChIPtsne2")
cbind(ct2.by_cell$MCF10A, ct2.by_cell$MCF10AT1, ct2.by_cell$MCF10CA1)
if(FALSE) rbind(ct2.by_cell$MCF10A, ct2.by_cell$MCF10AT1, ct2.by_cell$MCF10CA1)


if(FALSE) cbind(ct2.by_peak$`FALSE`, ct2.by_peak$`TRUE`)
rbind(ct2.by_peak$`FALSE`, ct2.by_peak$`TRUE`)

cbind(ct2.by_cell)
if(FALSE) cbind(ct2.by_cell, ct2.by_cell)
rbind(ct2.by_peak)

ct2.by_cell$MCF10A

# debug(unsplit, "ChIPtsne2List")
ct2.un_cell = unsplit(ct2.by_cell, "cell")
ct2.un_peak = unsplit(ct2.by_peak, "peak_MCF10CA1_CTCF")
#### rbind ####
dim(ct2.un_cell)
dim(ct2.by_cell$MCF10A)
dim(ct2.un_peak)
dim(ct2.by_peak[[1]])





