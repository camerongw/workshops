setwd('/data/gpfs/projects/punim1236/spleen_st/data/RCTD_10.2022')
.libPaths(
  c("/data/gpfs/projects/punim1236/R_envs/R-4.2.0_Seurat-v4",
    "/usr/local/easybuild-2019/easybuild/software/mpi/gcc/11.2.0/openmpi/4.1.1/r/4.2.0/lib64/R/library"))

## import
library(Seurat)
library(spacexr)
library(plotly)
library(htmlwidgets)
library(tidyverse)
library(circular)
library(useful)
library(dplyr)
library(dbscan)
library(reshape2)

normalize <- function(x){
  x <- (x-min(x))/(max(x)-min(x))
  return(x)
}
scale_c2l <- function(x){ ## enforce row sum value at 1
  x <- x / sum(x)
  return(x)
}

# Change the presentation of decimal numbers to 4 and avoid scientific notation
options(prompt="> ", digits=4, scipen=4)

# import RCTD data
infect_res <- readRDS("/data/gpfs/projects/punim1236/spleen_st/data/RCTD_10.2022/rctd_obj_D7.rds")
infect_res <- infect_res@RCTD.reps

# set save paths and heatmap reps
reps <- 1:4

#

#

#

#

### Prepare data
## check results for infected spleen
res <- list()

for(i in 1:4){
  
  ## get cell type weights
  map = as.data.frame(as.matrix(infect_res[[i]]@results$weights))
  
  ## format
  #map = as.data.frame(t(apply(X = map, MARGIN = 1, FUN = scale_c2l)))
  map = as.data.frame(apply(X = map, MARGIN = 2, FUN = normalize))
  
  ## get coords
  coords = as.data.frame(infect_res[[i]]@spatialRNA@coords)
  coords = coords[rownames(map),]
  
  ## centre coords for rotation
  coords$x = coords$x - mean(coords$x)
  coords$y = coords$y - mean(coords$y)
  
  ## convert pixels to microns
  coords[,c('x','y')] = coords[,c('x','y')] * 0.665
  
  ## get coord data
  map$x = coords$x
  map$y = coords$y
  map = map %>% relocate(y) %>% relocate(x)
  
  ## write
  res[[i]] <- map
}

#

#

#

#

### Functions for binning data spatially ###
## taken from Seurat
## bin data
SpCor <- function(query, var1, var2, x = 'x', y = 'y', radius = 50, nbins = 2000, verbose = F){
  
  if (verbose == T){
    message("Binning spatial data")
  }
  
  ## get pos
  pos = dbscan::frNN(query[,c(x,y)], eps = radius)
  
  ## subset to variables of interest
  query = as.data.frame(
    query[,c(var1,var2)])
  
  ## add 'mean' variable
  ## this will contain the mean for each spot's bin
  query$var1_mean = NA
  query$var2_mean = NA
  
  ## randomly select 1000 spots for binning
  bins = sample(x = 1:nrow(query), size = nbins)
  
  ## do binning
  for(i in bins){
    
    ## cell type 1
    query$var1_mean[i] = mean(
      query[,var1][pos$id[[i]]])
    
    ## cell type 1
    query$var2_mean[i] = mean(
      query[,var2][pos$id[[i]]])
    
  }
  
  ## remove NAs
  query = query %>% filter(!is.na(var1_mean))
  
  ## subset to means
  query = query[,c('var1_mean','var2_mean')]
  colnames(query) = c(var1,var2)
  
  #
  
  #
  
  #
  
  ## calculate correlation
  cor.res = cor.test(
    x=query[,var1],
    y=query[,var2],
    method = 'spearman'
  )
  
  message('Correlation: ',cor.res$estimate)
  
  ## Collect results
  
  results <- list()
  results[['cor']] = cor.res$estimate
  results[['pval']] = cor.res$p.value
  results[['misc']] = cor.res
  results[['Cell type 1']] = var1
  results[['Cell type 2']] = var2
  
  ## p value from this test may not be meaningful.
  ## random shuffling test probably a better way to estimate a p value
  
  return(results)
  
}

#

#

#

#

### Do SpCor for replicates
SpCorRep <- function(query, reps = 1:length(query), var1, var2, x = 'x', y = 'y', radius = 50, nbins = 2000, verbose = F){
  
  query_reps = query
  
  pooled_res = data.frame()
  
  for(i in reps){
    
    ## get query df
    query = query_reps[[i]]
    
    ## get pos
    pos = dbscan::frNN(query[,c(x,y)], eps = radius)
    
    ## subset to variables of interest
    query = as.data.frame(
      query[,c(var1,var2)])
    
    ## add 'mean' variable
    ## this will contain the mean for each spot's bin
    query$var1_mean = NA
    query$var2_mean = NA
    
    ## randomly select 1000 spots for binning
    bins = sample(x = 1:nrow(query), size = nbins)
    
    ## do binning
    for(i in bins){
      
      ## cell type 1
      query$var1_mean[i] = mean(
        query[,var1][pos$id[[i]]])
      
      ## cell type 1
      query$var2_mean[i] = mean(
        query[,var2][pos$id[[i]]])
      
    }
    
    ## remove NAs
    query = query %>% filter(!is.na(var1_mean))
    ## subset to means
    query = query[,c('var1_mean','var2_mean')]
    colnames(query) = c(var1,var2)
    
    pooled_res = rbind(pooled_res,query)
    
  }
  
  #
  
  #
  
  #
  
  ## calculate correlation
  query = pooled_res
  
  cor.res = cor.test(
    x=query[,var1],
    y=query[,var2],
    method = 'spearman'
  )
  
  message('Correlation: ',cor.res$estimate)
  
  ## Collect results
  
  results <- list()
  results[['cor']] = cor.res$estimate
  results[['pval']] = cor.res$p.value
  results[['misc']] = cor.res
  results[['Cell type 1']] = var1
  results[['Cell type 2']] = var2
  
  ## p value from this test may not be meaningful.
  ## random shuffling test probably a better way to estimate a p value
  
  return(results)
  
}

#

#

#

#

### Plot correlations of ALL cell types ###
### create results dataframe
celltypes <- colnames(res[[1]] %>% select(-x, -y))
n <- length(celltypes)

celltype1 <- rep(celltypes, each=n)
celltype2 <- rep(celltypes, times=n)

heatmap_data <- data.frame(celltype1,celltype2)
heatmap_data$Correlation <- 0
heatmap_data$bead_correlation <- 0
heatmap_data$bead_cor_pearson <- 0
heatmap_data$pval <- NA
heatmap_data$bead_pval <- NA

## run comparisons
for(i in 1:nrow(heatmap_data)){
  
  ## NA for same cell type
  if(heatmap_data[i,1] == heatmap_data[i,2]){
    heatmap_data$Correlation[i] = NA
    heatmap_data$bead_correlation[i] = NA
  }
  
  ## run comparisons
  ## add if/else so that if comparison is already run it isn't duplicated
  else{
    celltype_1 <- heatmap_data[i,1]
    celltype_2 <- heatmap_data[i,2]
    rownumber <- which(heatmap_data[i,1] == celltype2 & heatmap_data[i,2] == celltype1)
    
    ## copy results for existing comparisons to row of interest
    ## eg if Tfh v Th1 is done then copy it to Th1 v Tfh row
    if(heatmap_data[rownumber,'Correlation'] != 0){
      heatmap_data[i,3:7] <- heatmap_data[rownumber,3:7]
    } 
    
    ## run comparison
    else{
      
      ## Run using spatial binning function
      message(as.character(heatmap_data[i,1]),' vs ',as.character(heatmap_data[i,2]))
      cor_results = SpCorRep(
        query=res, reps = reps,
        var1 = as.character(heatmap_data[i,1]), 
        var2 = as.character(heatmap_data[i,2]))
      heatmap_data$Correlation[i] = as.numeric(cor_results$cor)
      heatmap_data$pval[i] = as.numeric(cor_results$pval)
      
    }
  }
}

#

#

#

#

### Plot correlations of ALL cell types ###
### create results dataframe
celltypes <- colnames(res[[1]] %>% select(-x, -y))
n <- length(celltypes)

celltype1 <- rep(celltypes, each=n)
celltype2 <- rep(celltypes, times=n)

heatmap_data <- data.frame(celltype1,celltype2)
heatmap_data$Correlation <- 0
heatmap_data$bead_correlation <- 0
heatmap_data$bead_cor_pearson <- 0
heatmap_data$pval <- NA
heatmap_data$bead_pval <- NA

## run comparisons
for(i in 1:nrow(heatmap_data)){
  
  ## NA for same cell type
  if(heatmap_data[i,1] == heatmap_data[i,2]){
    heatmap_data$Correlation[i] = NA
    heatmap_data$bead_correlation[i] = NA
  }
  
  ## run comparisons
  ## add if/else so that if comparison is already run it isn't duplicated
  else{
    celltype_1 <- heatmap_data[i,1]
    celltype_2 <- heatmap_data[i,2]
    rownumber <- which(heatmap_data[i,1] == celltype2 & heatmap_data[i,2] == celltype1)
    
    ## copy results for existing comparisons to row of interest
    ## eg if Tfh v Th1 is done then copy it to Th1 v Tfh row
    if(heatmap_data[rownumber,'Correlation'] != 0){
      heatmap_data[i,3:7] <- heatmap_data[rownumber,3:7]
    } 
    
    ## run comparison
    else{
      
      ## Run using spatial binning function
      message(as.character(heatmap_data[i,1]),' vs ',as.character(heatmap_data[i,2]))
      cor_results = SpCorRep(
        query = res, 
        reps = 1:4,
        var1 = as.character(heatmap_data[i,1]), 
        var2 = as.character(heatmap_data[i,2]),
        radius = 50)
      heatmap_data$Correlation[i] = as.numeric(cor_results$cor)
      heatmap_data$pval[i] = as.numeric(cor_results$pval)
      
    }
  }
}

#

#

#

#

### Save data ###
write.csv(heatmap_data,'D7_spat_cor_r50.csv')

### read ###
heatmap_data <- read.csv('D7_spat_cor_r50.csv') %>% select(-X)

colnames(heatmap_data)[3] <- 'ρ'

### order per hierarchical clustering ###

## unmelt df for hierarchical clustering
hdata <- dcast(heatmap_data, formula = celltype1~celltype2, fun.aggregate=sum, value.var='ρ')
rownames(hdata) <- hdata$celltype1

hdata <- hdata %>% 
  select(-celltype1)

hdata <- hclust(dist(hdata[,1:ncol(hdata)]))

## reorder cells per hclust
heatmap_data$celltype1 <- as.factor(heatmap_data$celltype1)
heatmap_data$celltype2 <- as.factor(heatmap_data$celltype2)

heatmap_data$celltype1 <- factor(
  x = heatmap_data$celltype1,
  ordered = T,
  levels = levels(heatmap_data$celltype1)[hdata$order]
)

heatmap_data$celltype2 <- factor(
  x = heatmap_data$celltype2,
  ordered = T,
  levels = levels(heatmap_data$celltype2)[hdata$order]
)

#

#

#

#

### Plot ###
cor_plot <- ggplot(data = heatmap_data, aes(x=celltype1,y=celltype2,fill=ρ)) +
  geom_tile() + 
  scale_fill_gradient2(
    low='royalblue3',mid='white',high='red3',midpoint=0,na.value = 'grey80', limits=c(-1,1)) +
  scale_y_discrete(limits = rev(levels(heatmap_data$celltype2))) +
  xlab(NULL) + ylab(NULL) +
  theme_minimal() +
  scale_x_discrete(position = "top") + 
  theme(axis.text.x = element_text(hjust = 0, angle = 45),
        legend.position = 'right',
        legend.direction = 'vertical',
        axis.line = element_line(colour = 'black')) +
  coord_fixed()
cor_plot
ggsave(plot=cor_plot, filename = 'D7_spat_cor_r50.pdf',
       units = 'in', width = 8, height = 8, dpi = 300, bg = 'white')

library(ggdendro)
dendro <- ggdendrogram(hdata)
ggsave(plot=dendro, filename = 'D7_spat_cor_r50_dendro.pdf',
       units = 'in', width = 6, height = 4, dpi = 300, bg = 'white')




