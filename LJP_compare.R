# This script compares R output to MATLAB output for the LJP dataset
# Created 07/20/2017 Nick Clark

library(readr)
library(dplyr)
ljp_GRval_r = read_tsv("LJP_GRval_R_0726.tsv")
ljp_GRmet_r = read_tsv("LJP_GRmet_R_0726.tsv")

#ljp_GRval_avg_r = read_tsv("LJP_GRval_avg_R.tsv")
#ljp_GRmet_avg_r = read_tsv("LJP_GRmet_avg_R.tsv")

#ljp_GRval_matlab_rep = read_tsv("level 4 GR50 (no published)/LJP_GRvalues_replicates_NC.tsv")
#ljp_GRmet_matlab_rep = read_tsv("level 4 GR50 (no published)/LJP_GRmetrics_replicates_NC.tsv")

ljp_GRval_matlab_merged = read_tsv("level 4 GR50 (no published)/LJP_GRvalues_merged_NC2.tsv")
ljp_GRmet_matlab_merged = read_tsv("level 4 GR50 (no published)/LJP_GRmetrics_merged_NC2.tsv")

# rename drugs:
rename_drugs = function(input_data) {
  input_data$drug[input_data$drug == "(R)- Roscovitine"] = "Seliciclib"
  input_data$drug[input_data$drug == "BYL719"] = "Alpelisib"
  input_data$drug[input_data$drug == "Flavopiridol"] = "Alvocidib"
  input_data$drug[input_data$drug == "BEZ235"] = "Dactolisib"
  input_data$drug[input_data$drug == "BMS-387032"] = "SNS-032"
  input_data$drug[input_data$drug == "CYT387"] = "Momelotinib"
  input_data$drug[input_data$drug == "GSK2126458"] = "Omipalisib"
  input_data$drug[input_data$drug == "GDC-0941"] = "Pictilisib"
  input_data$drug[input_data$drug == "GW843682"] = "GW843682X"
  input_data$drug[input_data$drug == "NVP-AUY922"] = "Luminespib"
  input_data$drug[input_data$drug == "Rapamycin"] = "Sirolimus"
  input_data$drug[input_data$drug == "TG 101348"] = "Fedratinib"
  return(input_data)
}

names(ljp_GRmet_matlab_merged)[2] = "drug"
ljp_GRmet_matlab_merged = rename_drugs(ljp_GRmet_matlab_merged)
ljp_GRmet_matlab_merged = ljp_GRmet_matlab_merged %>% arrange(cell_line,drug)

ljp_GRmet_r = ljp_GRmet_r[,c(1,2,8,9,10,11,12,13,14,15)]
ljp_GRmet_r = as.data.frame(ljp_GRmet_r)
ljp_GRmet_matlab_merged = as.data.frame(ljp_GRmet_matlab_merged)
all.equal(ljp_GRmet_r, ljp_GRmet_matlab_merged)
fivenum(ljp_GRmet_matlab_merged$GEC50 - ljp_GRmet_r$GEC50)

noFitR = which(ljp_GRmet_r$GEC50 == 0)
noFitmatlab = which(ljp_GRmet_matlab_merged$GEC50 == 0)
noFit = union(noFitR, noFitmatlab)
ljp_GRmet_r_fit = ljp_GRmet_r[-noFit,]
ljp_GRmet_matlab_merged_fit = ljp_GRmet_matlab_merged[-noFit,]
ljp_GRmet_r_fit$log10GEC50 = log10(ljp_GRmet_r_fit$GEC50)
ljp_GRmet_matlab_merged_fit$log10GEC50 = log10(ljp_GRmet_matlab_merged_fit$GEC50)
all.equal(ljp_GRmet_r_fit, ljp_GRmet_matlab_merged_fit)

#ljp_GRval_r_avg_after = ljp_GRval_r %>% group_by(cell_line,drug,concentration) %>% summarise(cell_count__time0 = mean(cell_count__time0), cell_count = mean(cell_count), cell_count__ctrl = mean(cell_count__ctrl), GRvalue = mean(GR))
names(ljp_GRval_matlab_merged)[2] = "drug"
ljp_GRval_matlab_merged = rename_drugs(ljp_GRval_matlab_merged) %>% arrange(GRvalue)#arrange(cell_line, drug, concentration, replicate)
ljp_GRval_r = ljp_GRval_r %>% arrange(GR)
fivenum(ljp_GRval_matlab_merged$GRvalue, ljp_GRval_r$GRvalue)

fivenum(ljp_GRval_r$GR - ljp_GRval_matlab_merged$GRvalue)
# [1] -5.107026e-15 -2.220446e-16  0.000000e+00  2.220446e-16  5.329071e-15
# Success! Averaging replicate cell counts before calculating GR values produces the same GR value output

ljp_GRmet_r_sub$log10GR50 = log10(ljp_GRmet_r_sub$GR50)
ljp_GRmet_r_sub$log10GEC50 = log10(ljp_GRmet_r_sub$GEC50)

ljp_GRmet_matlab_merged$log10GR50 = log10(ljp_GRmet_matlab_merged$GR50)
ljp_GRmet_matlab_merged$log10GEC50 = log10(ljp_GRmet_matlab_merged$GEC50)

all.equal(ljp_GRmet_r_sub, ljp_GRmet_matlab_merged)
# [1] "Attributes: < Length mismatch: comparison on first 2 components >"        
# [2] "Component “GR50”: Mean absolute difference: Inf"                          
# [3] "Component “GEC50”: Mean relative difference: 0.1519975"                   
# [4] "Component “GRinf”: Mean relative difference: 0.0610278"                   
# [5] "Component “h_GR”: Mean relative difference: 0.03675003"                   
# [6] "Component “r2_GR”: Mean relative difference: 0.3870447"                   
# [7] "Component “pval_GR”: Mean relative difference: 0.2716971"                 
# [8] "Component “log10GR50”: 'is.NA' value mismatch: 81 in current 60 in target"
# [9] "Component “log10GEC50”: Mean absolute difference: Inf"  
fivenum(ljp_GRmet_r_sub$log10GR50 - ljp_GRmet_matlab_merged$log10GR50)
#[1]          -Inf -7.165006e-06  2.895582e-06  4.146045e-05           Inf
fivenum(ljp_GRmet_r_sub$log10GEC50 - ljp_GRmet_matlab_merged$log10GEC50)
#[1]          -Inf -1.347939e-05  5.752204e-06  2.033689e-04           Inf
fivenum(ljp_GRmet_r_sub$GRinf - ljp_GRmet_matlab_merged$GRinf)
#[1] -1.386795e+00 -3.807953e-05 -1.967893e-08  7.794041e-07  1.107028e+00

nofit_r = which(ljp_GRmet_r_sub$GEC50 == 0)
nofit_matlab = which(ljp_GRmet_matlab_merged$GEC50 == 0)
nofit_r
nofit_matlab
GEC_diff = union(nofit_matlab[!nofit_matlab %in% nofit_r], nofit_r[!nofit_r %in% nofit_matlab])
GEC_diff_r = ljp_GRmet_r_sub[GEC_diff,]
GEC_diff_matlab = ljp_GRmet_matlab_merged[GEC_diff,]

ljp_GRmet_r2 = read_tsv("LJP_GRmet_R_07212017.tsv")
