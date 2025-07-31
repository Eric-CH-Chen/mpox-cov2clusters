#!/usr/bin/env -S Rscript 
#Importation from Nextflow JMC 2023
#
# Adapted from Jessica Caleta's covflow

# Reads in clust and clean ref

library(data.table)
library(tools)

args = commandArgs(trailingOnly=TRUE)
input_interval_file <- args[1] # for file name


# Read and combine results
lfiles <- list.files(pattern = ".*\\.csv$", full.names = TRUE)
results_dt <- rbindlist(lapply(lfiles, fread))

# Read interval file
in_dt <- read.table(input_interval_file, sep = "\t", header = TRUE, quote = "")

# Combine result and interval
#  option 'all=TRUE' should be fine if theres no weirdness
final_dt <- merge(in_dt, results_dt, by.x = "RunTag", by.y = "tag", all = TRUE)

final_file_name <- paste0(tools::file_path_sans_ext(input_interval_file), "_collection.tsv")
write.table(final_dt, file=final_file_name, quote=FALSE, row.names=FALSE, sep="\t")
