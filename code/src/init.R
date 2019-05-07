
# Set up environment
rm(list=ls())
gc()

# Set working directory
starting_dir <- getwd()
setwd('~/Desktop/repos/Cdiff_modeling/')

# Load dependencies
#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")
deps <- c('vegan', 'randomForest', 'gplots','viridis', 'scales','Hmisc', 'AUCRF', 'RColorBrewer')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep)

# Set seed for RNG
set.seed(8619)

#-----------------#

# Set input data location string variables
# Metadata
metadata <- 'data/metadata.tsv'
# 16S
shared <- 'data/16S/cdf_models.opti_mcc.0.03.shared'
tax <- 'data/16S/cdf_models.opti_mcc.0.03.cons.format.taxonomy'

# Set output plot location string variables
fig1Plot <- 'results/figures/figure_1.pdf'

# Set color palette
heat_palette <- viridis(n=150)
resistant_col <- '#ec4902'
susceptible_col <- '#000054ff'

strep_col <- 'mediumblue'
cef_col <- 'red2'
clinda_col <- 'limegreen'
noabx_col <- 'aquamarine2'

#-----------------#

# Define frequently used functions

# Merge 2 matrices with shared row names
mergeByRow <- function(data_1, data_2) {
  clean_merged <- merge(data_1, data_2, by='row.names', all.y=TRUE)
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  return(clean_merged)
}

# Remove columns with low variance
rmLowVar <- function(data_table) {
  vars <- apply(data_table, 2, var)
  keep <- which(vars > as.vector(quantile(vars, na.rm=TRUE))[2])
  keep_table <- data_table[, keep]
  return(keep_table)
}

# Returns medians from iterative subsampling of a vector
multiRare <- function(raw_abund, subs) {
  temp_abund <- as.vector(rrarefy(raw_abund, sample=subs))
  for (x in 1:999) {
    temp_abund <- cbind(temp_abund, as.vector(rrarefy(raw_abund, sample=subs)))
  }
  final_abund <- apply(temp_abund, 1, median) 
  rm(temp_abund)
  return(final_abund)
}

# Rarefy a shared file to the lowest sample total
rarefyOTU <- function(shared) {
  subSize <- min(rowSums(shared))
  shared <- t(shared)
  for (x in 1:ncol(shared)) {
    shared[,x] <- as.vector(rrarefy(shared[,x], sample=subSize))
  }
  shared <- as.data.frame(t(shared))
  return(shared)
}

# Filter out columns that have values in minimum number of samples (ignores first column if needed)
filtOTU <- function(data, minSamp) {
  drop <- c()
  if (class(data[,1]) != 'character') {
    if (sum(data[,1] != 0) < minSamp) {
      drop <- c(drop, colnames(data)[1])
    }
  }
  for (index in 2:ncol(data)) {
    if (sum(data[,index] != 0) < minSamp) {
      drop <- c(drop, colnames(data)[index])
    }
  }
  filtered_data <- data[,!(colnames(data) %in% drop)]
  return(filtered_data)
}

# Find centroid point of given MDS coordinates based on metadata
centroidMDS <- function(mds_coords) {
  centroids <- aggregate(cbind(mds_coords$MDS1, mds_coords$MDS2) ~ mds_coords[,1], 
                         data=mds_coords, mean)
  return(centroids)
}

# Calculate median and 95% confidence for non-normal, large datasets
# Based on: Conover, W.J. (1980) Practical Nonparametric Statistics John Wiley and Sons, New York.
confInt <- function(data) {
  data_median <- median(data)
  data <- sort(unique(data))
  n <- length(data)
  q <- 0.5
  nq <- n * q
  conf_range <- 1.96 * sqrt(n * q * (1 - q))
  j <- ceiling(nq - conf_range)
  k <- ceiling(nq + conf_range)
  lower_95 <- data[j]
  upper_95 <- data[k]
  return(c(lower_95, data_median, upper_95))
}

# Plot logarithmic tick marks on axes
logAxisTicks <- function(ax, n, t.ratio=0.5, mn, mx,...){
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  axis(ax, at=major.ticks, labels=rep('', length(major.ticks)), las=1)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  axis(ax, at=minor.ticks,tcl=par("tcl")*t.ratio, labels=FALSE)
}

# AUCRF feature selection based on infection susceptibility
aucrfSus <- function(training_data){
  
  # Format levels of susceptibility for AUCRF
  colnames(training_data) <- make.names(colnames(training_data))
  levels <- as.vector(unique(training_data$susceptibility))
  training_data$susceptibility <- as.character(training_data$susceptibility)
  training_data$susceptibility[which(training_data$susceptibility==levels[1])] <- 0
  training_data$susceptibility[which(training_data$susceptibility==levels[2])] <- 1
  training_data$susceptibility <- as.factor(as.numeric(training_data$susceptibility))
  rm(levels)
  
  # Run AUCRF with reproduceable parameters
  data_AUCRF <- AUCRF(susceptibility ~ ., data=training_data, pdel=0.05, k0=5, ranking='MDA')
  print(data_AUCRF)
  
  sig_features <- data_AUCRF$RFopt
  sig_features <- as.data.frame(sig_features$importance)
  
  return(sig_features)
}

# Standard random forest feature selection
rfSus <- function(training_data){
  
  colnames(training_data) <- make.names(colnames(training_data))
  mTries <- round(sqrt(ncol(training_data) - 1))
  nTrees <- length(as.vector(levels(training_data$susceptibility))) * (ncol(training_data) - 1)

  modelRF <- randomForest(training_data$susceptibility ~ ., data=training_data[,2:ncol(training_data)], importance=TRUE, replace=FALSE, err.rate=TRUE, ntree=nTrees, mtry=mTries)
  print(modelRF)
  
  featRF <- importance(modelRF, type=1)
  featRF <- as.data.frame(featRF)
  featRF$features <- rownames(featRF)
  rm(modelRF)
  
  # Filter to significant features (Strobl 2002) and sort
  sigFeatRF <- as.data.frame(subset(featRF, featRF$MeanDecreaseAccuracy > (abs(min(featRF$MeanDecreaseAccuracy)))))
  sigFeatRF <- sigFeatRF[order(-sigFeatRF$MeanDecreaseAccuracy),]
  sigFeatRF$features <- NULL
  rm(featRF)
  
  return(sigFeatRF)
}
rfAbx <- function(training_data){
  
  colnames(training_data) <- make.names(colnames(training_data))
  training_data$abx <- droplevels(training_data$abx)
  mTries <- round(sqrt(ncol(training_data) - 1))
  nTrees <- length(as.vector(levels(training_data$abx))) * (ncol(training_data) - 1)
  
  modelRF <- randomForest(training_data$abx ~ ., data=training_data[,2:ncol(training_data)], importance=TRUE, replace=FALSE, err.rate=TRUE, ntree=nTrees, mtry=mTries)
  print(modelRF)
  
  featRF <- importance(modelRF, type=1)
  featRF <- as.data.frame(featRF)
  featRF$features <- rownames(featRF)
  rm(modelRF)
  
  # Filter to significant features (Strobl 2002) and sort
  sigFeatRF <- as.data.frame(subset(featRF, featRF$MeanDecreaseAccuracy > (abs(min(featRF$MeanDecreaseAccuracy)))))
  rm(featRF)
  
  subset_data <- training_data[,which(colnames(training_data) %in% rownames(sigFeatRF))]
  subset_data$abx <- training_data$abx

  output_data <- list('sigFeatRF'=sigFeatRF, 'subset_data'=subset_data)
  return(output_data)
}

# Filter shred files by wilcoxon and larger median
sigDiff <- function(shared1, shared2){
  pval <- c()
  for (i in 1:ncol(shared1)) {
    pval[i] <- wilcox.test(shared1[,i], shared2[,i], exact=FALSE)$p.value
    }
  pval <- as.numeric(p.adjust(pval, method='BH'))
  shared1 <- shared1[,which(pval <= 0.05)]
  shared2 <- shared2[,which(pval <= 0.05)]
  rm(pval, i)
  keep <- c()
  diffs <- c()
  y <- 1
  for (x in 1:ncol(shared1)) {
    if (median(shared1[,x]) > median(shared2[,x])) {
      diffs[y] <- (median(shared1[,x]) - median(shared2[,x]))
      keep[y] <- x
      y <- y + 1
    }
  }
  cols <- colnames(shared1[,keep])
  rm(x, y, keep, shared1, shared2)
  
  final <- as.data.frame(cbind(cols, diffs))
  colnames(final) <- c('names','diff')
  
  return(final)
}

# Generates plot for significant differences in metabolite concentrations
multiStripchart <- function(metabolome1, metabolome2, 
                            treatment1, treatment2, 
                            treatment_col1, treatment_col2, 
                            pvalues, xLabel){
  
  formattedNames <- gsub('_', ' ', colnames(metabolome1))
  formattedNames <- gsub('\\.', '-', formattedNames)
  
  layout(matrix(c(1:(ncol(metabolome1)+2)), nrow=(ncol(metabolome1)+2), ncol=1, byrow = TRUE))
  
  par(mar=c(0.2, 0, 0, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE)
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  
  legend('bottomright', legend=c(treatment1, treatment2), bty='n',
         pt.bg=c(treatment_col1, treatment_col2), pch=21, cex=1.2, pt.cex=2, ncol=2)
  
  par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
  for(i in c(1:(ncol(metabolome1)))){
    xmax <- ceiling(max(c(max(metabolome1[,i]), max(metabolome2[,i]))))
    while(xmax %% 5 != 0 ){xmax <- xmax + 1}
    if (xmax > 1000) {while(xmax %% 100 != 0 ){xmax <- xmax + 1}} else if (xmax > 70){while(xmax %% 10 != 0 ){xmax <- xmax + 1}}
    plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,1.8))
    stripchart(at=1.2, jitter(metabolome1[,i], amount=1e-5), 
               pch=21, bg=treatment_col1, method='jitter', jitter=0.12, cex=2, add=TRUE)
    stripchart(at=0.66, jitter(metabolome2[,i], amount=1e-5), 
               pch=21, bg=treatment_col2, method='jitter', jitter=0.12, cex=2, add=TRUE)
    box()
    legend('topright', legend=formattedNames[i], pch=1, cex=1.3, pt.cex=0, bty='n')
    
    #text(x=1, y=1.5, labels=formattedNames[i])
    
    if (xmax <= 10) {
      text(x=seq(0,xmax,1), y=0.42, labels=seq(0,xmax,1), cex=1)
      axis(1, at=seq(0,xmax,1), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 1000){
      text(x=seq(0,xmax,200), y=0.42, labels=seq(0,xmax,200), cex=1)
      axis(1, at=seq(0,xmax,200), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 500){
      text(x=seq(0,xmax,100), y=0.42, labels=seq(0,xmax,100), cex=1)
      axis(1, at=seq(0,xmax,100), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 100){
      text(x=seq(0,xmax,50), y=0.42, labels=seq(0,xmax,50), cex=1)
      axis(1, at=seq(0,xmax,50), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 50){
      text(x=seq(0,xmax,10), y=0.42, labels=seq(0,xmax,10), cex=1)
      axis(1, at=seq(0,xmax,10), NA, cex.axis=0.8, tck=0.015)
    } else {
      text(x=seq(0,xmax,5), y=0.42, labels=seq(0,xmax,5), cex=1)
      axis(1, at=seq(0,xmax,5), NA, cex.axis=0.8, tck=0.015)
    }
    segments(median(metabolome1[,i]), 1.03, median(metabolome1[,i]), 1.37, lwd=2.5)
    segments(median(metabolome2[,i]), 0.49, median(metabolome2[,i]), 0.83, lwd=2.5)
    if (pvalues[i] <= 0.05){
      mtext('*', side=4, font=2, cex=1.6, padj=0.6)
    } else {
      mtext('n.s.', side=4, cex=0.9)
    }
  }
  
  par(mar=c(0, 0, 0, 0))
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  text(x=0, y=3.5, labels=xLabel, cex=1.4)
  #text(x=8, y=4, labels=paste('OOB Error = ', oob ,sep=''))
}
