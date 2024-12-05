# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: 'args' cannot be empty. Please provide valid arguments.")
}


set_seed <- as.logical(args[4])  # Boolean value to decide if the seed should be set

if (set_seed) {
  set.seed(142)  # Set the seed if condition is TRUE
  print("Seed has been set")
} else {
  print("No seed set")
}


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
library(SingleCellExperiment)
library(tidyverse)






# Assign the Bash parameters to quantile.p and opt
quantile.p <- as.numeric(args[1])
opt <- as.numeric(args[2])
dattype <- args[3]

samp.id <- unlist(str_split(args[5], ","))
pth <- args[6]
out_dir <- args[7]

if (length(args) < 8 || is.na(args[8]) || is.na(as.numeric(args[8]))) {
  run <- "Run.Norm"
} else {
  run <- paste0("Run.", as.numeric(args[8]))
}

if (args[length(args)] %in% c("normal", "fast")) {
  speed <- args[length(args)]
} else {
  speed <- "slow"
}


print(out_dir)
print(sprintf('quantile.p: %g',quantile.p)) ###
print(sprintf('opt: %g',opt))
print(sprintf('dattype: %s',dattype))
print(sprintf('Sample ID: %s',samp.id))
print(sprintf('Speed: %s',speed))

str(samp.id)

dat.seurat <- readRDS(pth)





## Uses dat type Parameter pulled from bash code to determine whether raw or normalized data will be used.
## Also Filtering out specific Samples, then using compute.network.model like kritika does. 



# Your code or random operations can follow here


if (dattype == "norm") {
  
  if (length(samp.id)>1) {
    Idents(dat.seurat) <- dat.seurat$batch
    k2 <- WhichCells(object = dat.seurat , idents = c(samp.id))
    dat.seurat.filt <- dat.seurat[,k2]
    Idents(dat.seurat.filt) <- dat.seurat.filt$batch
    
    rsamp <- map(samp.id, function(ids){
      filt <- length(SingleCellExperiment::counts(as.SingleCellExperiment(dat.seurat.filt[, WhichCells(object = dat.seurat , idents = ids)], 
                                                                          assay = "SCT"))@Dimnames[[2]])
      #print(filt) # For troubleshooting purposes only
      
      sce <- SingleCellExperiment::counts(as.SingleCellExperiment(dat.seurat.filt[, WhichCells(object = dat.seurat , idents = ids)], assay = "SCT"))[, sample(1:filt, 3500)]
    }) %>% reduce(cbind)
    
    ## obtain all gene names and store in vector
    genes <-  rsamp@Dimnames[[1]]
    
    ## create lncRNA specific vector
    lncs <- genes[grep("lnc", genes)]
    
    ## check that lncRNAs are in lncs vector
    print(head(lncs))
    
    
    network1 <- compute.network.sparse(expr.data = rsamp, quantile.p = quantile.p, gene.names = genes, opt = opt, lncs = lncs, speed.preset=speed)
  } else {
    Idents(dat.seurat) <- dat.seurat$batch
    k2 <- WhichCells(object = dat.seurat , idents = c(samp.id))
    dat.seurat.filt <- dat.seurat[,k2]
    Idents(dat.seurat.filt) <- dat.seurat.filt$batch
    samp <- as.SingleCellExperiment(dat.seurat.filt, assay = "SCT")
    
    ## obtain all gene names and store in vector
    genes <-  samp@rowRanges@partitioning@NAMES
    
    ## create lncRNA specific vector
    lncs <- genes[grep("lnc", genes)]
    
    ## check that lncRNAs are in lncs vector
    print(head(lncs))
    
    network1 <- compute.network.sparse(expr.data = samp@assays@data$counts, quantile.p = quantile.p, gene.names = genes, opt = opt, lncs = lncs, speed.preset=speed)
  }
  
} else {
  
  if (length(samp.id)>1) {
    Idents(dat.seurat) <- dat.seurat$batch
    k2 <- WhichCells(object = dat.seurat , idents = c(samp.id))
    dat.seurat.filt <- dat.seurat[,k2]
    Idents(dat.seurat.filt) <- dat.seurat.filt$batch
    
    rsamp <- map(samp.id, function(ids){
      filt <- length(SingleCellExperiment::counts(as.SingleCellExperiment(dat.seurat.filt[, WhichCells(object = dat.seurat , idents = ids)], 
                                                                          assay = "RNA"))@Dimnames[[2]])
      #print(filt) # For troubleshooting purposes only
      
      sce <- SingleCellExperiment::counts(as.SingleCellExperiment(dat.seurat.filt[, WhichCells(object = dat.seurat , idents = ids)], assay = "RNA"))[, sample(1:filt, 3500)]
    }) %>% reduce(cbind)
    
    ## obtain all gene names and store in vector
    genes <-  rsamp@Dimnames[[1]]
    
    ## create lncRNA specific vector
    lncs <- genes[grep("lnc", genes)]
    
    ## check that lncRNAs are in lncs vector
    print(head(lncs))
    
    
    network1 <- compute.network.sparse(expr.data = rsamp, quantile.p = quantile.p, gene.names = genes, opt = opt, lncs = lncs, speed.preset=speed)
  } else {
    Idents(dat.seurat) <- dat.seurat$batch
    k2 <- WhichCells(object = dat.seurat , idents = c(samp.id))
    dat.seurat.filt <- dat.seurat[,k2]
    Idents(dat.seurat.filt) <- dat.seurat.filt$batch
    samp <- as.SingleCellExperiment(dat.seurat.filt, assay = "RNA")
    
    ## obtain all gene names and store in vector
    genes <-  samp@rowRanges@partitioning@NAMES
    
    ## create lncRNA specific vector
    lncs <- genes[grep("lnc", genes)]
    
    ## check that lncRNAs are in lncs vector
    print(head(lncs))
  
    network1 <- compute.network.sparse(expr.data = samp@assays@data$counts, quantile.p = quantile.p, gene.names = genes, opt = opt, lncs = lncs, speed.preset=speed)

    }
}

#model1 <- compute.network.model(expr.data = cbind(ctrlmsce@assays@data$counts,ctrlfsce@assays@data$counts))

#saveRDS(model1, paste0("model_bigscale2_", Sys.Date(),".RDS"))



#network2 <- compute.network(expr.data = ctrlmsce@assays@data$counts, quantile.p = 0.90, model = model1, gene.names = ctrlm.genes)

#res <- lst("model" = model1,
           #"network" = network2)
if (grepl(",", args[5])) {
  samps <- str_remove_all(args[5], ",")
} else {
  samps <- args[5]
}


saveRDS(network1, paste0(out_dir,"/bigscale2_network_nomodel_",quantile.p,"_opt",opt,"_",dattype,"_",samps,"_",run,"_",speed,"_",format(Sys.time(), "%Y_%d_%m_%H_%M_%S"),".RDS"))

#saveRDS(res, paste0("bigscale2_res_",Sys.Date(),".RDS"))
