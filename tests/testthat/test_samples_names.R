testthat::context("sample_names")
# address some errors related to changing sample names
library(chiptsne2)
library(testthat)
library(ggplot2)

ct2 = exampleChIPtsne2.with_meta()
getSampleMetaData(ct2)
colData(ct2)$rep = "rep1"

ct2_r1 = ct2
ct2_r2 = ct2
ct2_r2$rep = "rep2"

colData(ct2_r1)
colData(ct2_r2)

colnames(ct2_r1) = paste(ct2_r1$cell, ct2_r1$rep, sep = "_")
colnames(ct2_r2) = paste(ct2_r2$cell, ct2_r2$rep, sep = "_")

ct2 = cbind(ct2_r1, ct2_r2)

colData(ct2)

ct2.agg = aggregateSamplesByGroup(ct2, group_VAR = "cell")
colData(ct2.agg)

getNameVariable(ct2)
ct2.sp = split(ct2, "sample")

ct2.10a = (ct2.sp$MCF10A_rep1 + ct2.sp$MCF10A_rep2) / 2
ct2.at1 = (ct2.sp$MCF10AT1_rep1 + ct2.sp$MCF10AT1_rep2) / 2
ct2.avg = cbind(ct2.10a, ct2.at1)

colnames(ct2.avg) = c("MCF10A", "MCF10AT1")

ct2.avg = swapNameVariable(ct2.avg, "cell")

colData(ct2.avg)
