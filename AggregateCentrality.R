source("~/Bigscale_pipe/CentAggr.R")

# File path for recorded output paths
args <- commandArgs(trailingOnly = TRUE)

# Read the paths
output_paths <- readLines(args)

# Ensure all paths are identical
if (length(unique(output_paths)) > 1) {
  stop("Error: Output paths are not consistent across runs. Check output_paths.txt.")
}

# Proceed with aggregation using the common output path
common_path <- unique(output_paths)[1]
cat("All paths are consistent. Proceeding with path:", common_path, "\n")

# Your centrality aggregation logic here

print(paste("Processing data in:", common_path))
suffixes <- rdsfilename_to_suffix(common_path)
aggrcent <- CentAggr(pth = common_path, suffixes = suffixes, centralities = c("Degree", "Betweenness", "Closeness", "PAGErank"))