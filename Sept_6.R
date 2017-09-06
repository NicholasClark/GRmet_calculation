# Breast cancer density project
# Level 2: LDS-1261 or HMS dataset 20256
# Level 3: LDS-1262 or HMS dataset 20257
# Level 4: LDS-1263 or HMS dataset 20258

library(readr)
library(GRmetrics)
library(dplyr)

# will incorrectly guess "Small Mol Concentration" as integer column unless you increase "guess_max"
lvl2 = read_tsv("HMS data/Breast Cancer Density/Level 2/20256.txt", guess_max = 55000)
lvl3 = read_tsv("HMS data/Breast Cancer Density/Level 3/20257.txt")
lvl4 = read_tsv("HMS data/Breast Cancer Density/Level 4/20258.txt")


.trim_mean = function(x, percent) {
  x = x[!is.na(x)]
  n = length(x)
  k = n*(percent/100)/2
  # round down if k is half an integer
  if(round(k) != k & round(k*2) == k*2) {
    lo = floor(k) + 1
    hi = n - lo + 1
  } else {
    lo = round(k) + 1
    hi = n - lo + 1
  }
  x = sort(x)[lo:hi]
  return(mean(x))
}

lvl2_gp = lvl2 %>% group_by(`Cell HMS LINCS ID`, `Cell Name`, `Small Molecule HMS LINCS ID`, `Small Molecule Name`, `Small Mol Concentration`, `Small Mol Conc Unit`, Replicate, `Seeding Density`) %>%
  dplyr::summarise(
  `Mean Total Cell Count Before Treatment` = round(mean(`Total Cell Count Before Treatment`),2),
  `Mean Total Cell Count After Treatment` = round(mean(`Total Cell Count After Treatment`),2),
  `Mean Total Control Cell Count` = round(mean(`Total Control Cell Count`),2),
  `Mean Relative Cell Count` = round(mean(`Relative Cell Count`),4),
  `Mean Normalized Growth Rate Inhibition Value` = round(mean(`Normalized Growth Rate Inhibition Value`),4) ) %>%
  dplyr::filter(`Small Mol Concentration` != 0) %>% dplyr::arrange(`Cell Name`, `Small Molecule Name`, Replicate, `Seeding Density`, `Small Mol Concentration`) %>% ungroup %>% as.data.frame() %>% `rownames<-`(seq_len(nrow(lvl3)))

lvl3_sort = lvl3 %>% arrange(`Cell Name`, `Small Molecule Name`, Replicate, `Seeding Density`, `Small Mol Concentration`) %>% as.data.frame() %>% `rownames<-`(seq_len(nrow(lvl3)))

all.equal(lvl3_sort, lvl2_gp)
# [1] "Attributes: < Component “spec”: Component “cols”: Names: 5 string mismatches >"   
# [2] "Component “Mean Total Control Cell Count”: Mean relative difference: 9.130227e-07"
# [3] "Component “Mean Relative Cell Count”: Mean relative difference: 0.0001155902"   
# only different for huge values - rounding error
# Note: their metadata says they do a 50%-trimmed mean of cell counts, but this doesn't appear to be the case

lvl3_rename = lvl3 %>% dplyr::rename(concentration = `Small Mol Concentration`, cell_count = `Mean Total Cell Count After Treatment`, cell_count__time0 = `Mean Total Cell Count Before Treatment`, cell_count__ctrl = `Mean Total Control Cell Count`)

lvl3_MATLAB_input = lvl3_rename[,1:11]

write_tsv(lvl3_MATLAB_input, path = "20257_rename.tsv")

GRmet = GRfit(lvl3_rename, c("Cell HMS LINCS ID", "Cell Name", "Small Molecule HMS LINCS ID", "Small Molecule Name", "Replicate", "Seeding Density"))
GRmet_lvl4 = GRgetMetrics(GRmet)
rownames(GRmet_lvl4) = 1:dim(GRmet_lvl4)[1]
colnames(GRmet_lvl4)[1:6] = colnames(lvl4)[1:6]
GRmet_lvl4_new = GRmet_lvl4[,c(1:8, 12:29)]
GRmet_lvl4_rename = GRmet_lvl4_new %>% rename(`Emax (IC)` = `Emax`, `AUC (IC)` = AUC, `Curve Fit (IC)` = `fit_rel_cell`, `Curve Fit (GR)` = `fit_GR`, `r^2 (IC)` = r2_rel_cell, `EC50 (IC)` = `EC50`, `Einf (IC)` = `Einf`, `Hill Coefficient (IC)` = h, `AUC (GR)` = GR_AOC, `r^2 (GR)` = r2_GR, `EC50 (GR)` = GEC50, `Hill Coefficient (GR)` = h_GR)

drop = names(GRmet_lvl4_rename) %in% names(lvl4)
GRmet_lvl4_rename = GRmet_lvl4_rename[, drop]
reorder = match(names(lvl4), names(GRmet_lvl4_rename))
GRmet_lvl4_reorder = GRmet_lvl4_rename[,reorder]

GRmet_lvl4_compare = GRmet_lvl4_reorder %>% arrange(`Cell HMS LINCS ID`, `Cell Name`, `Small Molecule HMS LINCS ID`, `Small Molecule Name`, `Replicate`, `Seeding Density`) %>% as.data.frame() %>% `rownames<-`(seq_len(nrow(lvl4)))


lvl4_compare = lvl4 %>% arrange(`Cell HMS LINCS ID`, `Cell Name`, `Small Molecule HMS LINCS ID`, `Small Molecule Name`, `Replicate`, `Seeding Density`) %>% as.data.frame() %>% `rownames<-`(seq_len(nrow(lvl4)))

all.equal(lvl4_compare, GRmet_lvl4_compare)

MATLAB_output = read_tsv("density_GRmetrics.tsv")
lvl4_compare_matlab = lvl4[, c(1:6, 15:17, 19:22)]
MATLAB_output = MATLAB_output[,c(1:4, 6:10, 14, 11:13)]
names(MATLAB_output) = names(lvl4_compare_matlab)

lvl4_compare_matlab = lvl4_compare_matlab %>% dplyr::arrange(`Cell HMS LINCS ID`, `Cell Name`, `Small Molecule HMS LINCS ID`, `Small Molecule Name`, `Replicate`, `Seeding Density`) %>% as.data.frame() %>% `rownames<-`(seq_len(nrow(lvl4)))

MATLAB_output = MATLAB_output %>% dplyr::arrange(`Cell HMS LINCS ID`, `Cell Name`, `Small Molecule HMS LINCS ID`, `Small Molecule Name`, `Replicate`, `Seeding Density`) %>% as.data.frame() %>% `rownames<-`(seq_len(nrow(lvl4)))

all.equal(MATLAB_output, lvl4_compare_matlab)
