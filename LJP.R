# This script converts the raw data for the LINCS Joint Project in GR metrics
# Created 07/20/2017 Nick Clark

library(readr)
library(dplyr)
library(GRmetrics)

# Convert GR value data for replicates 1, 2, and 3 into averaged GR value data
# Converts datasets 20248, 20249, and 20250 into form of dataset 20251 and compares them

# Import HMS data from the 3 replicates
# http://lincs.hms.harvard.edu/db/datasets/20248
ljpHMS_rep1 = read_csv("LJP/dataset_20248_20170719191123.csv")
# http://lincs.hms.harvard.edu/db/datasets/20249
ljpHMS_rep2 = read_csv("LJP/dataset_20249_20170719191502.csv")
# http://lincs.hms.harvard.edu/db/datasets/20250
ljpHMS_rep3 = read_csv("LJP/dataset_20250_20170719191518.csv")

# Combine the 3 replicates into one table, order them, and find mean and 
# std err of GR values and relative cell count over the three replicates
ljpHMS_rep1$Replicate = 1
ljpHMS_rep2$Replicate = 2
ljpHMS_rep3$Replicate = 3
ljpHMS_allreps = rbind(ljpHMS_rep1, ljpHMS_rep2, ljpHMS_rep3)
ljpHMS_allreps_avg = ljpHMS_allreps %>% group_by(`Small Molecule Name`,`Cell Name`, `Small Mol Concentration`, `LJP Library Plate`) %>%
  summarise(meanGR = mean(`Mean Normalized Growth Rate Inhibition Value`),
            stderrGR = sd(`Mean Normalized Growth Rate Inhibition Value`)/sqrt(3),
            mean_relcell = mean(`Mean Relative Cell Count`),
            stderr_relcell = sd(`Mean Relative Cell Count`)/sqrt(3))

# Import HMS data from the already averaged replicates and order the table
# http://lincs.hms.harvard.edu/db/datasets/20251
ljpHMS_avg = read_csv("LJP/dataset_20251_20170719183318.csv")
ljpHMS_avg_order = ljpHMS_avg %>% arrange(`Small Molecule Name`,`Cell Name`, `Small Mol Concentration`, `LJP Library Plate`)

# Compare dataset 20251 to average of the other 3 datasets
max(abs(round(ljpHMS_allreps_avg$meanGR,4) - ljpHMS_avg_order$`Mean Normalized Growth Rate Inhibition Value`)) 
#[1] 1e-04
max(abs(round(ljpHMS_allreps_avg$stderrGR,4) - ljpHMS_avg_order$`Standard Error of the Mean for Normalized Growth Rate Inhibition Value`)) 
#[1] 1e-04
max(abs(round(ljpHMS_allreps_avg$mean_relcell,4) - ljpHMS_avg_order$`Mean Relative Cell Count`)) 
#[1] 1e-04
max(abs(round(ljpHMS_allreps_avg$stderr_relcell,4) - ljpHMS_avg_order$`Standard Error of the Mean for Relative Cell Count`)) 
#[1] 1e-04

##############
# Compute GR metrics based on average of replicates 1, 2, and 3 - datasets 20248, 20249, and 20250
ljpHMS_allreps_sub = ljpHMS_allreps[,c(2,4,5,7,8,10,11,12)]
colnames(ljpHMS_allreps_sub)[8] = "cell_count__ctrl"
colnames(ljpHMS_allreps_sub)[7] = "cell_count"
colnames(ljpHMS_allreps_sub)[6] = "cell_count__time0"
colnames(ljpHMS_allreps_sub)[3] = "concentration"
ljpHMS_allreps_GRmet = GRfit(ljpHMS_allreps_sub, c("Small Molecule Name", "Cell Name"))
all_reps_GRmet = GRgetMetrics(ljpHMS_allreps_GRmet) %>% arrange(`Small.Molecule.Name`,`Cell.Name`)

##############
# Import datasets 20245, 20246, and 20247
ljp1 = read_csv("LJP/dataset_20245_20170719170028.csv")
ljp2 = read_csv("LJP/dataset_20246_20170720022312.csv")
ljp3 = read_csv("LJP/dataset_20247_20170720022325.csv")

ljp1$Replicate = 1
ljp2$Replicate = 2
ljp3$Replicate = 3

ljp1$exp = paste(ljp1$`Cell Name`,ljp1$`Small Molecule Name`, ljp1$`Small Mol Concentration`, ljp1$`LJP Library Plate`)
ljp2$exp = with(ljp2, paste(`Cell Name`,`Small Molecule Name`, `Small Mol Concentration`,`LJP Library Plate`))
ljp3$exp = with(ljp3, paste(`Cell Name`,`Small Molecule Name`, `Small Mol Concentration`,`LJP Library Plate`))
length(unique(ljp1$exp))

ljp = rbind(ljp1,ljp2,ljp3)
ljp_sub = ljp[,c(2,4,5,7,8,12,13,14,16)]
colnames(ljp_sub)[8] = "cell_count__ctrl"
colnames(ljp_sub)[7] = "cell_count"
colnames(ljp_sub)[6] = "cell_count__time0"
colnames(ljp_sub)[3] = "concentration"
ljp_sub_noNA = ljp_sub[!is.na(ljp_sub$`Small Molecule Name`),]
ljpGR = GRfit(ljp_sub_noNA, c("Small Molecule Name", "Cell Name"))
ljpGRmet = GRgetMetrics(ljpGR) %>% arrange(`Small.Molecule.Name`,`Cell.Name`)
ljpGRval = GRgetValues(ljpGR) %>% arrange(`Small Molecule Name`,`Cell Name`)

##############
# Import GR metrics from HMS database, dataset 20252
# http://lincs.hms.harvard.edu/db/datasets/20252/results
ljpHMS_GRmet = read_csv("LJP/dataset_20252_20170719215358.csv") %>%
  arrange(`Small Molecule Name`, `Cell Name`)

##############
# Compare GR metrics from dataset 20252 to those calculated from 20245, 20246, 20247

sum(ljpHMS_GRmet$`Cell Name` != ljpGRmet$Cell.Name)
# [1] 0
sum(ljpHMS_GRmet$`Small Molecule Name` != ljpGRmet$Small.Molecule.Name)
# [1] 0
# All combinations of cell line and drug are the same
fivenum(ljpGRmet$GR_AOC - ljpHMS_GRmet$GR_AOC)
#[1] -0.026082069 -0.002689459  0.001959314  0.007414373  0.037939197
# not quite equal to GRmetrics R output of GR_AOC

##############
# Import GRbrowser LJP dataset
ljpGRbrowser = read_tsv("LJP/20170224_LJP_for_GRB.tsv") %>%
  arrange(`Small_Mol_HMSLID`, `Cell_Name`)

##############
# Compare GRbrowser dataset to HMS dataset GR metrics
# Note: "Fedratinib" in data from HMS website is "TG 101348" in GRbrowser data
ljpGRbrowser$Small_Mol_Name[ljpGRbrowser$Small_Mol_Name == "TG 101348"] = "Fedratinib"
ljpGRbrowser = ljpGRbrowser %>% arrange(`Small_Mol_Name`, `Cell_Name`)
which(ljpGRmet$Small.Molecule.Name != ljpGRbrowser$Small_Mol_Name)
#integer(0)
which(ljpGRmet$Cell.Name != ljpGRbrowser$Cell_Name)
#integer(0)
# drugs and cell lines match
fivenum(ljpGRbrowser$GR_AOC - ljpHMS_GRmet$GR_AOC)
# [1] -5e-05 -2e-05  0e+00  2e-05  5e-05

##############
# Compare GRbrowser GR metrics to calculated GR metrics from 20245, 20246, 20247
fivenum(ljpGRbrowser$GR_AOC - ljpGRmet$GR_AOC)
# [1] -0.037969197 -0.007410365 -0.001946717  0.002655364  0.026082069
fivenum(ljpGRbrowser$GRmax - ljpGRmet$GRmax)
# [1] 0.02591055 0.10854670 0.15423455 0.22425560 0.90388043
# GRmax should be dead on... must be some issue here with rounding or averaging.

##############
# Compare GRbrowser GR metrics to calculated GR metrics from 20248, 20249, 20250
which(all_reps_GRmet$Small.Molecule.Name != ljpGRbrowser$Small_Mol_Name)
#integer(0)
which(all_reps_GRmet$Cell.Name != ljpGRbrowser$Cell_Name)
fivenum(ljpGRbrowser$GR_AOC - all_reps_GRmet$GR_AOC)
# [1] -0.033827488 -0.006769788 -0.001510732  0.002810604  0.026212550
hist(ljpGRbrowser$GRmax - all_reps_GRmet$GRmax)
# Compare GR metrics from 20248, 20249, 20250 to GR metrics from 20245, 20246, 20247
fivenum(ljpGRmet$GR_AOC - all_reps_GRmet$GR_AOC)
# [1] -9.375279e-03  3.710637e-05  3.548758e-04  8.576722e-04  1.903413e-02
# much closer, but should be better for AOC and GRmax

#############
# Compare old LJP to new LJP
old_ljp = read_tsv("LJP/LJP_GRmetrics_merged.tsv") %>% arrange(cellLine,smallMolecule)
new_ljp = read_tsv("LJP/20170224_LJP_for_GRB.tsv") %>% arrange(`Cell_Name`,`Small_Mol_Name`)
which(old_ljp$cellLine != new_ljp$Cell_Name)
# cell lines the same
which(old_ljp$smallMolecule != new_ljp$Small_Mol_Name)
# drugs named differently
# Rename and re-arrange
old_ljp$smallMolecule[old_ljp$smallMolecule == "(R)- Roscovitine"] = "Seliciclib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "BYL719"] = "Alpelisib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "Flavopiridol"] = "Alvocidib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "BEZ235"] = "Dactolisib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "BMS-387032"] = "SNS-032"
old_ljp$smallMolecule[old_ljp$smallMolecule == "CYT387"] = "Momelotinib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "GSK2126458"] = "Omipalisib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "GDC-0941"] = "Pictilisib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "GW843682"] = "GW843682X"
old_ljp$smallMolecule[old_ljp$smallMolecule == "NVP-AUY922"] = "Luminespib"
old_ljp$smallMolecule[old_ljp$smallMolecule == "Rapamycin"] = "Sirolimus"
old_ljp = old_ljp %>% arrange(cellLine,smallMolecule)
which(old_ljp$smallMolecule != new_ljp$Small_Mol_Name)
# drugs now the same
old_ljp2 = old_ljp[,c(1:3,6,4,7,8,5,9,10)]
new_ljp2 = new_ljp[,c(2,4:12)]
names(new_ljp2)[1:2] = c("cellLine","smallMolecule")
names(new_ljp2)[7] = "Hill"
names(new_ljp2)[9] = "r2"
old_ljp2 = as.data.frame(old_ljp2)
new_ljp2 = as.data.frame(new_ljp2)
all.equal(new_ljp2, old_ljp2)
# [1] "Component “pval”: Mean relative difference: 0.0005905279"

# Old and new GRbrowser dataset for LJP are the same except for some names of
# drugs and extremely small (5e-08) p-value differences

write_tsv(ljpHMS_allreps_sub,"ljp_484950.tsv")
