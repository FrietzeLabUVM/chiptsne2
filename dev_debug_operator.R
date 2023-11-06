ct2.by_treatment = readRDS("../chiptsne2_paper/ct2.by_treatment.Rds")
e1 = ct2.by_treatment$bza
e2 = ct2.by_treatment$ctrl
operator = "-"

ct2.diff =ct2.by_treatment$bza - ct2.by_treatment$ctrl
colData(ct2.diff)
