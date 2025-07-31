#!/usr/bin/env -S Rscript 
#Importation from Nextflow JMC 2023
#
# Adapted from Jessica Caleta's covflow

# Reads in clust and clean ref. Saves both modified (for cluster comparison) and unmodified cluster file
#
# Note:
# The string for outlier is hard-coded: "Outlier" for reference, and "-1" for cluster file (via cov2cluster's output)

# Library
library(data.table)

# Command line args
args = commandArgs(trailingOnly=TRUE)
input_run_clust <- args[1] #just placeholder in funct definition, but default if ordered not specified
input_ref_clust <- args[2]
input_run_flag <- args[3]

# Read in data
clust_df <- read.table(input_run_clust, sep="\t", header=TRUE, quote="")
ref_df <- fread(input_ref_clust)


# Ensures the consistency in naming of the columns
colnames(clust_df) <- c("strain", "clust")
colnames(ref_df) <- c("strain", "clust")

# ensure only those that are present in both are included
clust_df <- setDT(clust_df)[strain %chin% ref_df$strain]
ref_df <- setDT(ref_df)[strain %chin% clust_df$strain]
clust_df <- clust_df[, -3]

# sort so the 'strain' are the same
setorder(clust_df, cols="strain")
setorder(ref_df, cols="strain")

# Prepare and write filtered but unmodified cluster file
clust_df_fill_f <- paste0(input_run_flag, "_clust.tsv")
ref_df_fill_f <- paste0(input_run_flag, "_ref.tsv")
write.table(clust_df, file=clust_df_fill_f, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(ref_df, file=ref_df_fill_f, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# Modify outlier so they are different
clust_df[, clust := ifelse( clust == "-1", paste0("Outlier", .I), clust)]
ref_df[, clust := ifelse( clust == "Outlier", paste0("Outlier", .I), clust)]

# Prepare and write filtered and unmodified cluster file
clust_dt_f <- paste0(input_run_flag, "_clust_rename.tsv")
ref_dt_f <- paste0(input_run_flag, "_ref_rename.tsv")
write.table(clust_df, file=clust_dt_f, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(ref_df, file=ref_dt_f, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
