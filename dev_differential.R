#method to subtract
#method to filter
library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()

subset(ct2, peak_MCF10CA1_CTCF, cell == "MCF10A")


# dplyr_row_slice(data, i, ...)
#
# dplyr_col_modify(data, cols)
#
# dplyr_reconstruct(data, template)

#
# filter.ChIPtsne2 = function(.data, ...){
#     subset(.data, ...)
# }

# filter(ct2, peak_MCF10CA1_CTCF, cell == "MCF10A")

library(dplyr)
filter(ct2, peak_MCF10CA1_CTCF, cell == "MCF10A")


ct2_a = filter(ct2, peak_MCF10CA1_CTCF, cell == "MCF10A")

ct2_b = filter(ct2, peak_MCF10CA1_CTCF, cell == "MCF10AT1")

ct2_a - ct2_b


ct2_sub = apply_chiptsne2_operator(ct2_a, ct2_b, "-")
ct2_plus = apply_chiptsne2_operator(ct2_a, ct2_b, "+")
ct2_div = apply_chiptsne2_operator(ct2_a, ct2_b, "/")
ct2_mult = apply_chiptsne2_operator(ct2_a, ct2_b, "*")
apply_chiptsne2_operator(ct2_a, ct2_b, "asdf")

ct2_sub.num = apply_chiptsne2_operator.num(ct2_a, 5, "-")
ct2_plus.num = apply_chiptsne2_operator.num(ct2_a, 5, "+")
ct2_div.num = apply_chiptsne2_operator.num(ct2_a, 5, "/")
ct2_mult.num = apply_chiptsne2_operator.num(ct2_a, 5, "*")

plotSignalHeatmap(ct2_sub, heatmap_colors = scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-40, 40)))

# debug(getTidyProfile, "ChIPtsne2")
getTidyProfile(ct2_sub)
getTidyProfile(ct2_plus)
getTidyProfile(ct2_div)
getTidyProfile(ct2_mult)
# debug(plotSignalHeatmap, "ChIPtsne2")

plotSignalHeatmap(ct2_sub, name_FUN = paste)
plotSignalHeatmap(ct2_plus, name_FUN = paste)
plotSignalHeatmap(ct2_div, name_FUN = paste)
plotSignalHeatmap(ct2_mult, name_FUN = paste)

plotSignalHeatmap(ct2_sub.num, name_FUN = paste)
plotSignalHeatmap(ct2_plus.num, name_FUN = paste)
plotSignalHeatmap(ct2_div.num, name_FUN = paste)
plotSignalHeatmap(ct2_mult.num, name_FUN = paste)


plotSignalHeatmap(ct2_b - ct2_a)
plotSignalHeatmap(ct2_b + ct2_a)
plotSignalHeatmap(ct2_b * ct2_a)
plotSignalHeatmap(ct2_b / ct2_a)


plotSignalHeatmap(ct2_b - 50)
plotSignalHeatmap(ct2_b + 50)
plotSignalHeatmap(ct2_b * 5)
plotSignalHeatmap(ct2_b / 5)

plotSignalHeatmap(50 - ct2_b)
plotSignalHeatmap(50 + ct2_b)
plotSignalHeatmap(5 * ct2_b)
plotSignalHeatmap(5 / ct2_b)

plotSignalHeatmap((ct2_a + 5) / (ct2_b + 5), name_FUN = paste)
plotSignalHeatmap((ct2_a) - (ct2_b), name_FUN = paste)
