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
ljpHMS_rep1$replicate = 1
ljpHMS_rep2$replicate = 2
ljpHMS_rep3$replicate = 3
ljpHMS_allreps = rbind(ljpHMS_rep1, ljpHMS_rep2, ljpHMS_rep3) %>% arrange(`Small Molecule Name`,`Cell Name`)
ljpHMS_allreps_sub = ljpHMS_allreps[,c(2,4,5,10,11,12,15)]
names(ljpHMS_allreps_sub) = c("cell_line", "drug", "concentration", "cell_count__time0", "cell_count", "cell_count__ctrl", "replicate")
#ljpHMS_allreps_sub_avg = ljpHMS_allreps_sub %>% group_by(cell_line, drug, concentration) %>% summarise(cell_count__time0 = mean(cell_count__time0), cell_count = mean(cell_count), cell_count__ctrl = mean(cell_count__ctrl))

ljpGR = GRfit(ljpHMS_allreps_sub, c("cell_line", "drug"))
ljpGRmet = GRgetMetrics(ljpGR) %>% arrange(`cell_line`,`drug`)
ljpGRval = GRgetValues(ljpGR) %>% arrange(cell_line, drug, replicate)
write_tsv(ljpGRval, "LJP_GRval_R_0726.tsv")
write_tsv(ljpGRmet, "LJP_GRmet_R_0726.tsv")

#ljpGR = GRfit(ljpHMS_allreps_sub_avg, c("cell_line", "drug"))
#ljpGRmet = GRgetMetrics(ljpGR) %>% arrange(`cell_line`,`drug`)
#ljpGRval = GRgetValues(ljpGR) %>% arrange(`cell_line`,`drug`)
#write_tsv(ljpGRval, "LJP_GRval_avg_R.tsv")
#write_tsv(ljpGRmet, "LJP_GRmet_avg_R.tsv")

