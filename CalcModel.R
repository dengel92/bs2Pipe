library(Matrix)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(knitr)
library(foreach)
library(doParallel)
library(cluster)
library(limma)
library(enrichR)
library(Rtsne)
library(Seurat)
library(SeuratObject)
source("~/bigScale2/R/Functions_v2.R")
source("~/bigScale2/R/SingleCellMethods.R")
source("~/bigScale2/R/RcppExports.R")
library(Rcpp)
sourceCpp("~/bigScale2/src/de.cpp")
sourceCpp("~/bigScale2/src/dist.cpp")
sourceCpp("~/bigScale2/src/dist_icells.cpp")
source("~/Bigscale_pipe/CentAggr.R")
source("~/Bigscale_pipe/BS2_PCA.R")
library(SingleCellExperiment)
library(tidyverse)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: 'args' cannot be empty. Please provide valid arguments.")
}

set.seed(142)  # Set the seed 

home <- args[1]

model.out <- args[2]

datetime <- args[3]

params <- read.delim2(paste0(home,"/config.tsv"))

samp <- unique(unlist(flatten(stringr::str_split(unique(params$samp.id), ","))))

datatype <- params %>% 
  select("dattype") %>% 
  mutate(dattype = case_when(
    dattype == "norm" ~ "SCT",
    dattype == "raw" ~ "RNA"
  )) %>% t() %>% as.list()


# Load and Filter Data
dat.seurat <- readRDS(unique(params$path))

#calculate Model for each data type

models <- map(datatype, function(X){
  Idents(dat.seurat) <- dat.seurat$batch
  k2 <- WhichCells(object = dat.seurat , idents = c(samp))
  dat.seurat.filt <- dat.seurat[,k2]
  Idents(dat.seurat.filt) <- dat.seurat.filt$batch
  data <- as.SingleCellExperiment(dat.seurat.filt, assay = X)
  
  model <- compute.network.model(data@assays@data$counts)
})


saveRDS(models[[1]], paste0(model.out,"/model_raw_",datetime,".RDS"))
saveRDS(models[[2]], paste0(model.out,"/model_norm_",datetime,".RDS"))