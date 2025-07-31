#!/usr/bin/env -S Rscript 

library(data.table)
library(tools)

args = commandArgs(trailingOnly=TRUE)
input_transmission_file <- args[1]  # input name
output_nx_file <- args[2]			# output file name for nextstrain
output_gephi_file <- args[3]		# output file name for gephi node-time

# Read in file
raw_trans_file <- read.table(input_transmission_file, sep="\t", quote="")

# Part 1: Generate input for nextstrain
# Transform df so it can be used as meta file for nextstrain
nx_trans_file <- raw_trans_file[,-1]
colnames(nx_trans_file) <- c("strain", "date")
nx_trans_file$time <- nx_trans_file$date
nx_trans_file$date <- as.Date(nx_trans_file$date, format= "%Y-%m-%d", origin="2023-01-01")

# rbind the Ancestral Sequences
ancestral_df = data.frame("strain" = "Ancestral_Sequence", "date" = "2022-12-31", "time" = -1)
nx_trans_file = rbind(nx_trans_file, ancestral_df)

# Write output
write.table(nx_trans_file, file=output_nx_file, sep="\t", row.names = FALSE, quote = FALSE)

# Part 2: Generate time table for gephi
gephi_node_time_file <- nx_trans_file
colnames(gephi_node_time_file) <- c("Id", "date")
write.table(gephi_node_time_file, file=output_gephi_file, sep="\t", row.names = FALSE, quote = FALSE)
