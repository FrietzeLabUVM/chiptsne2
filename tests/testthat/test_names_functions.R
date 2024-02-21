testthat::context("all_functions")
# flipping viewGranges
library(chiptsne2)
library(testthat)

ct2_objects = list(
    ct2_noMeta = exampleChIPtsne2(),
    ct2_wMeta = exampleChIPtsne2.with_meta()
)

#all functions should adjust to renamed variables
ct2_objects = lapply(ct2_objects, function(x){
    x = setPositionVariable(x, "new_pos")
    x = setNameVariable(x, "new_name")
    x = setRegionVariable(x, "new_reg")
    x = setValueVariable(x, "new_value")
    x
})


ct2_objects.no_rowRanges = lapply(ct2_objects, function(x){
    rowRanges(x) = NULL
    x
})
names(ct2_objects.no_rowRanges) = paste0(names(ct2_objects), ".no_rowRanges")

ct2_objects.all = c(ct2_objects, ct2_objects.no_rowRanges)

ct2 = ct2_objects.all$ct2_noMeta
rownames(ct2) = paste0("row_", rownames(ct2))
colnames(ct2) = paste0("col_", colnames(ct2))
