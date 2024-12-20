library(tidyverse)

rdsfilename_to_suffix <- function(pth, stime) {
  files <- list.files(pth, full.names = TRUE, pattern = "\\.RDS$", recursive = T) %>% 
    grep(pattern = stime, . , value = T)
  print(length(files))
  # Extract the suffixes using str_extract
  suffixes <- sapply(files, function(file) {
    # Extract the filename (without path)
    filename <- basename(file)
    # print(filename)
    # Use str_extract to capture the G190Mx and Run.X parts
    g190m <- str_extract(filename, "(G\\d+M\\d+)+")  # Extract G190M and the number
    # print(length(g190m))
    # Check if the filename contains "norm" or "raw"
    data_type <- ifelse(str_detect(filename, "norm"), "norm", 
                        ifelse(str_detect(filename, "raw"), "raw", "unknown"))
    # print(paste0("data_type length is ", length(data_type)))
    print(paste0("Data Type: ", data_type))
    
    run_num <- str_extract(filename, "Run\\.\\d+")  # Extract "Run.X" (including "Run")
    # print(length(run_num))
    print(paste0("Run #: ", run_num))
    
    if (!is.na(g190m) && !is.na(run_num)) { #
      # Combine the G190M, Run number, and suffix with an underscore
      return(paste(g190m, data_type, run_num, sep = "_"))
    }
    return(paste(g190m, data_type, sep = "_"))  # Return g190m + data_type if the match is not found
  })
  
  # Remove any NAs (if any filenames did not match)
  # print(length(suffixes))
  # print(suffixes)
  suffixes <- suffixes[!is.na(suffixes)]
  
  return(suffixes)
}

load_rds_with_suffixes <- function(pth, stime, suffixes) {
  files <- list.files(pth, full.names = TRUE, pattern = "\\.RDS$", recursive = T) %>% 
  grep(pattern = stime, . , value = T)
   print(length(files))
  # Check if the number of suffixes matches the number of files
  if (length(files) != length(suffixes)) {
    stop("The number of suffixes must match the number of .RDS files.")
  }
  
  # Loop through files and suffixes to assign unique names
  lapply(seq_along(files), function(i) {
    obj_name <- paste0(make.names(gsub("\\.RDS$", "", basename(files[i]))), "_", suffixes[i])
    assign(obj_name, readRDS(files[i]), envir = .GlobalEnv)  # Assign to the global environment
  })
}


# Function to add top 100 indicators and ensure they are the last columns
add_top100_indicators <- function(df, centralities, suffix) {
  df %>%
    # Rename centrality columns with the provided suffix
    rename_with(~ paste0(., "_", suffix), .cols = all_of(centralities)) %>%
    
    # Calculate ranks and top 100 indicators for each centrality measure
    mutate(across(starts_with(paste0(centralities, "_", suffix)), 
                  list(rank = ~ rank(-.),  # Rank in descending order
                       top100 = ~ ifelse(. %in% sort(., decreasing = TRUE)[1:100], 1, 0)),
                  .names = "{.col}_{fn}")) %>%
    
    # Rearrange to ensure top 100 indicator and rank columns are at the end
    select(everything(), ends_with("_rank"), ends_with("_top100"))
}

CentAggr <- function(pth, suffixes = NULL, centralities = NULL, runtime) {
  
  load_rds_with_suffixes(pth = pth, stime = runtime, suffixes = suffixes)
  # Get the list of objects in the global environment
  all_objects <- ls(.GlobalEnv)[sapply(mget(ls(.GlobalEnv), .GlobalEnv), is.list)]
  
  #print(all_objects)
  # Exclude the function parameters (suffixes and centralities) from the list of objects
  excluded_objects <- c("suffixes", "centralities")
  
  data_sources <- mget(all_objects, mode = "list", envir = .GlobalEnv)
  
  #str(data_sources)
  
  # Use mapply to iterate over both data sources and topgenes
  dataframes <- mapply(function(data) {
    data$centrality %>%
      rownames_to_column("ids")
  }, data_sources, SIMPLIFY = FALSE)
  
  #str(dataframes)
  
  
  #print(head(dataframes))
  # Apply the function to each data frame
  dataframes_with_indicators <- mapply(function(df, suffix) {
    add_top100_indicators(df, centralities, suffix)
  }, dataframes, suffixes, SIMPLIFY = FALSE)
  
  #str(dataframes_with_indicators)
  
  combined_df <- purrr::reduce(dataframes_with_indicators, function(df1, df2) {
    full_join(df1, df2, by = "ids")  # You can use other join types like inner_join, full_join, etc.
  }) %>%
    arrange(across(ends_with("_top100"))) %>%
    select(everything(), ends_with("_top100"))
  
  # Reorder columns to place top 100 indicators at the end
  all_colnames <- names(combined_df)
  top100_columns <- grep("_top100$", all_colnames, value = TRUE)
  other_columns <- setdiff(all_colnames, top100_columns)
  
  # Select columns with non-top100 ones first, then top100 columns
  reordered_combined_df <- combined_df %>%
    select(all_of(other_columns), all_of(top100_columns))%>%
    mutate(across(ends_with("_rank"), ~ ifelse(is.na(.) | . == "", 10000, as.numeric(.))))%>%
    mutate(across(-c(ends_with("_rank"), ids), ~ ifelse(is.na(.) | . == "", 0, as.numeric(.))))
  
  
  
  sum.stats <- imap_dfr(data_sources,function(X, ids){
    
    tmp <- tibble(NetID = str_extract(ids, "(G\\d+M\\d+)+_Run\\.(\\d+|[a-zA-Z]+)"),
                  quartile.p = str_extract(ids, "\\d+\\.\\d+"),
                  cutoff.p = X$cutoff.p,
                  Node.cnt = length(igraph::V(X$graph)),
                  Clust.cnt = max(X$clusters))
  })

  
  return(list(aggrdf = reordered_combined_df, sum.stats = sum.stats))
}

