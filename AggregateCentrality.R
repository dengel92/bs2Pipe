source("~/Bigscale_pipe/CentAggr.R")
source("~/Bigscale_pipe/BS2_PCA.R")
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(Rtsne)
library(rrcov)
library(gridExtra)
library(tidyverse)
library(purrr)

# File path for recorded output paths
args <- commandArgs(trailingOnly = TRUE)

# Read the paths
stime <- args[1]

print(str(stime))

print(length(stime))

pth <- args[2]

# # Ensure all paths are identical
# if (length(unique(output_paths)) > 1) {
#   stop("Error: Output paths are not consistent across runs. Check output_paths.txt.")
# }
# 
# # Proceed with aggregation using the common output path
# common_path <- unique(output_paths)[1]
# cat("All paths are consistent. Proceeding with path:", common_path, "\n")

# Your centrality aggregation logic here

print(paste("Processing data in:", paste0("BS2_",stime)))
suffixes <- rdsfilename_to_suffix(pth, stime)
print(length(suffixes))

centralities <- c("Degree", "Betweenness", "Closeness", "PAGErank")

aggrcent <- CentAggr(pth = pth, suffixes = suffixes, centralities = centralities, stime)

# str(aggrcent[[1]])
# head(aggrcent[[1]])


qc.plot <- pca.plot(aggrcent[[1]], cent = centralities)

files <- list.files(pth, full.names = TRUE, pattern = "\\.RDS$", recursive = T) %>% 
  grep(pattern = stime, . , value = T)

subout <- unique(dirname(files))

print(subout)

## Save each graph takes a separate sheet of a PDF:
pdf(paste0(subout, "/bigScale2QC_plots_output.pdf"), width = 15, height = 8.5)

# Filter out NULL values from qc.plot
filtered_qc_plot <- qc.plot[!sapply(qc.plot, is.null)]

# Use walk to print the remaining plots
walk(filtered_qc_plot, ~{
  print(marrangeGrob(.x, nrow = 1, ncol = 1))
})

dev.off()

openxlsx::write.xlsx(aggrcent[[1]], paste0(subout,"/aggrecent","_",format(Sys.time(), "%Y_%d_%m_%H_%M_%S"),".xlsx"))

openxlsx::write.xlsx(aggrcent[[2]], paste0(subout,"/Summary_stats","_",format(Sys.time(), "%Y_%d_%m_%H_%M_%S"),".xlsx"))
