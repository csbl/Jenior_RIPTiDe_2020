

# Load in data
transcription <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/transcript/clinda_k12.mapped.norm.tsv', header=TRUE, row.names=1)

# Format data
transcription <- transcription$normDepth
k12_quantiles <- as.vector(quantile(transcription, c(0.5, 0.625, 0.75, 0.875)))

library(vegan)
library(ape)

# Pick which samples will be taken
small_sample <- sample(c(1:10000), 100, replace=FALSE)
large_sample <- sample(c(1:10000), 300, replace=FALSE)

# Load in  and format data
lb_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/LB_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(lb_aerobic_samples) <- make.names(colnames(lb_aerobic_samples))
lb_aerobic_samples <- lb_aerobic_samples[small_sample,]
lb_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
lb_aerobic_samples$condition <- 'in vitro'

m9_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/M9_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(m9_aerobic_samples) <- make.names(colnames(m9_aerobic_samples))
m9_aerobic_samples <- m9_aerobic_samples[small_sample,]
m9_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
m9_aerobic_samples$condition <- 'in vitro'

m9_anaerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/M9_anaerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(m9_anaerobic_samples) <- make.names(colnames(m9_anaerobic_samples))
m9_anaerobic_samples <- m9_anaerobic_samples[small_sample,]
m9_anaerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
m9_anaerobic_samples$condition <- 'in vitro'
rm(small_sample)

# Save metadata
lb_metadata <- lb_aerobic_samples[, c('sample', 'condition')]
lb_metadata$media <- 'LB'
rownames(lb_metadata) <- lb_metadata$sample
lb_metadata$sample <- NULL
rownames(lb_aerobic_samples) <- lb_aerobic_samples$sample
lb_aerobic_samples$sample <- NULL
lb_aerobic_samples$condition <- NULL
lb_aerobic_samples <- as.data.frame(t(apply(lb_aerobic_samples, 2, as.numeric)))
m9_aerobic_metadata <- m9_aerobic_samples[, c('sample', 'condition')]
m9_aerobic_metadata$media <- 'm9_aerobic'
rownames(m9_aerobic_metadata) <- m9_aerobic_metadata$sample
m9_aerobic_metadata$sample <- NULL
rownames(m9_aerobic_samples) <- m9_aerobic_samples$sample
m9_aerobic_samples$sample <- NULL
m9_aerobic_samples$condition <- NULL
m9_aerobic_samples <- as.data.frame(t(apply(m9_aerobic_samples, 2, as.numeric)))
m9_anaerobic_metadata <- m9_anaerobic_samples[, c('sample', 'condition')]
m9_anaerobic_metadata$media <- 'm9_anaerobic'
rownames(m9_anaerobic_metadata) <- m9_anaerobic_metadata$sample
m9_anaerobic_metadata$sample <- NULL
rownames(m9_anaerobic_samples) <- m9_anaerobic_samples$sample
m9_anaerobic_samples$sample <- NULL
m9_anaerobic_samples$condition <- NULL
m9_anaerobic_samples <- as.data.frame(t(apply(m9_anaerobic_samples, 2, as.numeric)))
invitro_metadata <- as.data.frame(rbind(lb_metadata, m9_aerobic_metadata, m9_anaerobic_metadata))
rm(lb_metadata, m9_aerobic_metadata, m9_anaerobic_metadata)

# Merge in vitro samples
invitro_samples <- merge(lb_aerobic_samples, m9_aerobic_samples, by='row.names')
rownames(invitro_samples) <- invitro_samples$Row.names
invitro_samples$Row.names <- NULL
invitro_samples <- merge(invitro_samples, m9_anaerobic_samples, by='row.names')
rownames(invitro_samples) <- invitro_samples$Row.names
invitro_samples$Row.names <- NULL
invitro_samples <- as.data.frame(t(invitro_samples))
rownames(invitro_samples) <- rownames(invitro_metadata)
rm(lb_aerobic_samples, m9_aerobic_samples, m9_anaerobic_samples)

# Load and format in vivo samples
invivo_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/clinda_k12.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(invivo_samples) <- make.names(colnames(invivo_samples))
invivo_samples <- invivo_samples[large_sample,]
invivo_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
invivo_samples$condition <- 'in vivo'
invivo_metadata <- invivo_samples[, c('sample', 'condition')]
invivo_metadata$media <- 'none'
rownames(invivo_metadata) <- invivo_metadata$sample
invivo_metadata$sample <- NULL
invivo_samples$sample <- NULL
invivo_samples$condition <- NULL
rownames(invivo_samples) <- rownames(invivo_metadata)

# Combine datasets
invitro_samples <- as.data.frame(t(invitro_samples))
invivo_samples <- as.data.frame(t(invivo_samples))
all_samples <- merge(invitro_samples, invivo_samples, by='row.names')
invitro_samples <- as.data.frame(t(invitro_samples))
invivo_samples <- as.data.frame(t(invivo_samples))
rownames(all_samples) <- all_samples$Row.names
all_samples$Row.names <- NULL
all_samples <- as.data.frame(t(all_samples))
all_metadata <- as.data.frame(rbind(invivo_metadata, invitro_metadata))
rm(invivo_metadata, invitro_metadata)


mod <- rda(all_samples, scale = TRUE)
biplot(mod, scaling = 3, type = c("text", "points"))

# Transform data due to negative values and calculate axes
all_samples <- all_samples + max(all_samples)
samples_dist <- vegdist(all_samples, method='bray')
samples_pcoa_obj <- pcoa(samples_dist, correction='none', rn=NULL)
samples_pcoa_points <- as.data.frame(samples_pcoa_obj$vectors[,c(1,2)])
samples_pcoa_values <- samples_pcoa_obj$values
samples_pcoa_values <- samples_pcoa_values$Relative_eig[c(1,2)] * 100.0
colnames(samples_pcoa_points) <- round(samples_pcoa_values, digits=2)
rm(samples_pcoa_obj, samples_pcoa_values)

# Combine with metadata
samples_pcoa_points <- merge(samples_pcoa_points, all_metadata, by='row.names')
rownames(samples_pcoa_points) <- samples_pcoa_points$Row.names
samples_pcoa_points$Row.names <- NULL

# Calculate differences
all_samples <- merge(all_metadata, all_samples, by='row.names')
rownames(all_samples) <- all_samples$Row.names
all_samples$Row.names <- NULL
all_samples$media <- NULL
condition_permANOVA_pval <- adonis(samples_dist ~ condition, all_samples, perm=999)$aov.tab[[6]][1]
all_samples$condition <- NULL
all_samples <- merge(all_metadata, all_samples, by='row.names')
rownames(all_samples) <- all_samples$Row.names
all_samples$Row.names <- NULL
all_samples$condition <- NULL
media_permANOVA_pval <- adonis(samples_dist ~ media, all_samples, perm=999)$aov.tab[[6]][1]
rm(all_samples, samples_dist)

# Prep for plotting
samples_pcoa_points[,2] <- samples_pcoa_points[,2] * -1
x_axis_lab <- gsub('X', '', as.character(colnames(samples_pcoa_points)[1]))
x_axis_lab <- paste ('PC1 (', x_axis_lab, '%)', sep='')
y_axis_lab <- gsub('X', '', as.character(colnames(samples_pcoa_points)[2]))
y_axis_lab <- paste ('PC2 (', y_axis_lab, '%)', sep='')

# Subset points
vivo_points <- subset(samples_pcoa_points, condition == 'in vivo')
vitro_points <- subset(samples_pcoa_points, condition == 'in vitro')
lb_points <- subset(samples_pcoa_points, media == 'LB')
m9_aer_points <- subset(samples_pcoa_points, media == 'm9_aerobic')
m9_an_points <- subset(samples_pcoa_points, media == 'm9_anaerobic')

library(scales)
library(randomForest)

# Load and format samples
pfba_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/pFBA.flux_samples.format.tsv', sep='\t', header=TRUE, row.names=1)
colnames(pfba_samples) <- make.names(colnames(pfba_samples))
pfba_samples <- pfba_samples[large_sample,]
pfba_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL

# Calculate scaled relative flux within datasets
for (x in 1:nrow(pfba_samples)) {pfba_samples[x,] <- (pfba_samples[x,] / sum(pfba_samples[x,])) * 1000}
for (x in 1:nrow(invitro_samples)) {invitro_samples[x,] <- (invitro_samples[x,] / sum(invitro_samples[x,])) * 1000}
for (x in 1:nrow(invivo_samples)) {invivo_samples[x,] <- (invivo_samples[x,] / sum(invivo_samples[x,])) * 1000}

# Combine datasets
pfba_samples$condition <- 'pfba'
invitro_samples$condition <- 'in vitro'
invivo_samples$condition <- 'in vivo'
pfba_samples <- as.data.frame(t(pfba_samples))
invitro_samples <- as.data.frame(t(invitro_samples))
invivo_samples <- as.data.frame(t(invivo_samples))
all_samples1 <- as.data.frame(merge(pfba_samples, invivo_samples, by='row.names'))
rownames(all_samples1) <- all_samples1$Row.names
all_samples1$Row.names <- NULL
all_samples1 <- as.data.frame(t(all_samples1))
all_samples2 <- merge(invitro_samples, invivo_samples, by='row.names')
rownames(all_samples2) <- all_samples2$Row.names
all_samples2$Row.names <- NULL
all_samples2 <- as.data.frame(t(all_samples2))
pfba_samples <- as.data.frame(t(pfba_samples))
invitro_samples <- as.data.frame(t(invitro_samples))
invivo_samples <- as.data.frame(t(invivo_samples))

# pFBA vs in vivo
# Prep for machine learning
# Breiman (2001). Random Forests. Machine Learning.
overlap <- intersect(colnames(lb_aerobic_samples), colnames(invivo_samples))
lb_aerobic_samples <- lb_aerobic_samples[, overlap]
invivo_samples <- invivo_samples[, overlap]
lb_aerobic_samples$condition <- 'lb'
invivo_samples$condition <- 'in vivo'
all_samples <- rbind(lb_aerobic_samples, invivo_samples)
all_samples$condition <- as.factor(as.numeric(all_samples$condition))
sub1 <- round(length(rownames(all_samples1[which(all_samples$condition=='lb'),])) * 0.623)
sub2 <- round(length(rownames(all_samples1[which(all_samples$condition=='in vivo'),])) * 0.623)
fctr <- max(c(round(sub1 / sub2), round(sub2 / sub1))) * 3
n_trees <- round(length(colnames(all_samples)) - 1) * fctr
m_tries <- round(sqrt(length(colnames(all_samples)) - 1))
condition <- all_samples$condition
all_samples$condition <- NULL
all_samples <- as.data.frame(apply(all_samples, 2, as.numeric))
rm(sub1, sub2, fctr)

# Run Random Forest and parse results
rf_obj <- randomForest(condition ~ ., data=all_samples, ntree=n_trees, mtry=m_tries,
                       importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_obj)
rf_feat <- importance(rf_obj, type=1, scale=FALSE)
pfba_invivo_mda <- subset(rf_feat, rf_feat > (50*abs(min(rf_feat)))) # ~p-value of 0.001
pfba_invivo_mda <- as.data.frame(pfba_invivo_mda)
pfba_invivo_mda$feature <- rownames(pfba_invivo_mda)
pfba_invivo_mda <- pfba_invivo_mda[order(-pfba_invivo_mda$MeanDecreaseAccuracy),]
rm(rf_obj, rf_feat, condition, all_samples1, n_trees, m_tries)






# in vitro vs in vivo
# Prep for machine learning
# Breiman (2001). Random Forests. Machine Learning.
all_samples2$condition <- as.factor(as.character(all_samples2$condition))
sub1 <- round(length(rownames(all_samples2[which(all_samples2$condition=='in vitro'),])) * 0.623)
sub2 <- round(length(rownames(all_samples2[which(all_samples2$condition=='in vivo'),])) * 0.623)
fctr <- max(c(round(sub1 / sub2), round(sub2 / sub1))) * 3
n_trees <- round(length(colnames(all_samples2)) - 1) * fctr
m_tries <- round(sqrt(length(colnames(all_samples2)) - 1))
condition <- all_samples2$condition
all_samples2$condition <- NULL
all_samples2 <- as.data.frame(apply(all_samples2, 2, as.numeric))
rm(sub1, sub2, fctr)

# Run Random Forest and parse results
rf_obj <- randomForest(condition ~ ., data=all_samples2, ntree=n_trees, mtry=m_tries,
                       importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_obj)
rf_feat <- importance(rf_obj, type=1, scale=FALSE)
invitro_invivo_mda <- subset(rf_feat, rf_feat > (50*abs(min(rf_feat)))) # ~p-value of 0.001
invitro_invivo_mda <- as.data.frame(invitro_invivo_mda)
invitro_invivo_mda$feature <- rownames(invitro_invivo_mda)
invitro_invivo_mda <- invitro_invivo_mda[order(-invitro_invivo_mda$MeanDecreaseAccuracy),]
rm(rf_obj, rf_feat, condition, all_samples2, n_trees, m_tries)


#------------------------------------------------------------------------------------------------

# Load in data
transcription <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/transcript/clinda_k12.mapped.norm.tsv', header=TRUE, row.names=1)
transcription <- transcription$normDepth
cutoffs <- as.vector(quantile(transcription, probs=c(0.5,0.625,0.75,0.875)))

# Pick which samples will be taken
sub_sample <- sample(c(1:10000), 500, replace=FALSE)

# Focus final analyses on LB vs in vivo 
# Load in data and perform basic formatting
lb_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/LB_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE, row.names=1)
colnames(lb_aerobic_samples) <- make.names(colnames(lb_aerobic_samples))
lb_aerobic_samples <- lb_aerobic_samples[sub_sample,]
lb_aerobic_samples <- as.data.frame(apply(lb_aerobic_samples, 2, as.numeric))
lb_aerobic_biomass <- lb_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M
lb_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
lb_aerobic_samples$DM_mththf_c <- NULL # Remove intracellular demand reaction

#m9_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/M9_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE, row.names=1)
#colnames(m9_aerobic_samples) <- make.names(colnames(m9_aerobic_samples))
#m9_aerobic_samples <- m9_aerobic_samples[small_sample,]
#m9_aerobic_samples <- as.data.frame(apply(m9_aerobic_samples, 2, as.numeric))
#m9_aerobic_biomass <- m9_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M
#m9_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
#m9_aerobic_samples$DM_mththf_c <- NULL # Remove intracellular demand reaction

#m9_anaerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/M9_anaerobic.flux_samples.format.tsv', sep='\t', header=TRUE, row.names=1)
#colnames(m9_anaerobic_samples) <- make.names(colnames(m9_anaerobic_samples))
#m9_anaerobic_samples <- m9_anaerobic_samples[small_sample,]
#m9_anaerobic_samples <- as.data.frame(apply(m9_anaerobic_samples, 2, as.numeric))
#m9_anaerobic_biomass <- m9_anaerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M
#m9_anaerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
#m9_anaerobic_samples$DM_mththf_c <- NULL # Remove intracellular demand reaction

invivo_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/clinda_k12.flux_samples.format.tsv', sep='\t', header=TRUE, row.names=1)
colnames(invivo_samples) <- make.names(colnames(invivo_samples))
invivo_samples <- invivo_samples[sub_sample,]
invivo_samples <- as.data.frame(apply(invivo_samples, 2, as.numeric))
invivo_biomass <- invivo_samples$BIOMASS_Ec_iJO1366_WT_53p95M
invivo_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
invivo_samples$DM_mththf_c <- NULL # Remove intracellular demand reaction
rm(sub_sample)

# Normalize data
median_biomass <- median(c(lb_aerobic_biomass, invivo_biomass))
for (x in 1:nrow(lb_aerobic_samples)) {lb_aerobic_samples[x,] <- (lb_aerobic_samples[x,] / lb_aerobic_biomass[x]) * median_biomass}
#for (x in 1:nrow(m9_aerobic_samples)) {m9_aerobic_samples[x,] <- (m9_aerobic_samples[x,] / m9_aerobic_biomass[x]) * median_biomass}
#for (x in 1:nrow(m9_anaerobic_samples)) {m9_anaerobic_samples[x,] <- (m9_anaerobic_samples[x,] / m9_anaerobic_biomass[x]) * median_biomass}
for (x in 1:nrow(invivo_samples)) {invivo_samples[x,] <- (invivo_samples[x,] / invivo_biomass[x]) * median_biomass}
rm(median_biomass, lb_aerobic_biomass, invivo_biomass)

# Subset exchanges and metabolic reactions
lb_aerobic_exchanges <- lb_aerobic_samples[, grep('EX_', colnames(lb_aerobic_samples))]
lb_aerobic_samples <- lb_aerobic_samples[, !colnames(lb_aerobic_samples) %in% colnames(lb_aerobic_exchanges)]
transporters <- c(colnames(lb_aerobic_samples)[grep('tex', colnames(lb_aerobic_samples))], 
                  colnames(lb_aerobic_samples)[grep('pp', colnames(lb_aerobic_samples))])
lb_aerobic_samples <- lb_aerobic_samples[, !colnames(lb_aerobic_samples) %in% transporters]

#m9_aerobic_exchanges <- m9_aerobic_samples[, grep('EX_', colnames(m9_aerobic_samples))]
#m9_aerobic_samples <- m9_aerobic_samples[, !colnames(m9_aerobic_samples) %in% colnames(m9_aerobic_exchanges)]
#transporters <- c(colnames(m9_aerobic_samples)[grep('tex', colnames(m9_aerobic_samples))], 
#                  colnames(m9_aerobic_samples)[grep('pp', colnames(m9_aerobic_samples))])
#m9_aerobic_samples <- m9_aerobic_samples[, !colnames(m9_aerobic_samples) %in% transporters]

#m9_anaerobic_exchanges <- m9_anaerobic_samples[, grep('EX_', colnames(m9_anaerobic_samples))]
#m9_anaerobic_samples <- m9_anaerobic_samples[, !colnames(m9_anaerobic_samples) %in% colnames(m9_anaerobic_exchanges)]
#transporters <- c(colnames(m9_anaerobic_samples)[grep('tex', colnames(m9_anaerobic_samples))], 
#                  colnames(m9_anaerobic_samples)[grep('pp', colnames(m9_anaerobic_samples))])
#m9_anaerobic_samples <- m9_anaerobic_samples[, !colnames(m9_anaerobic_samples) %in% transporters]

invivo_exchanges <- invivo_samples[, grep('EX_', colnames(invivo_samples))]
invivo_samples <- invivo_samples[, !colnames(invivo_samples) %in% colnames(invivo_exchanges)]
transporters <- c(colnames(invivo_samples)[grep('tex', colnames(invivo_samples))], 
                  colnames(invivo_samples)[grep('pp', colnames(invivo_samples))])
invivo_samples <- invivo_samples[, !colnames(invivo_samples) %in% transporters]
rm(transporters)

# Consider only net-imported metabolites
lb_aerobic_exchanges <- lb_aerobic_exchanges[, which(as.vector(colMeans(lb_aerobic_exchanges)) < 0)]
invivo_exchanges <- invivo_exchanges[, which(as.vector(colMeans(invivo_exchanges)) < 0)]

# Collect differences in exchanges across models
lb_aerobic_only <- setdiff(colnames(lb_aerobic_exchanges), colnames(invivo_exchanges))
invivo_only <- setdiff(colnames(invivo_exchanges), colnames(lb_aerobic_exchanges))
lb_aerobic_exchanges <- lb_aerobic_exchanges[, lb_aerobic_only]
invivo_exchanges <- invivo_exchanges[, invivo_only]
rm(lb_aerobic_only, invivo_only)

# Calculate absolute flux through exchanges 
lb_aerobic_exchanges <- abs(lb_aerobic_exchanges)
lb_aerobic_exchanges <- log2(lb_aerobic_exchanges + 1)
invivo_exchanges <- abs(invivo_exchanges)
invivo_exchanges <- log2(invivo_exchanges + 1)

# Subsetting to those >1
lb_aer_exch_med <- apply(lb_aerobic_exchanges, 2, median)
lb_aerobic_exchanges <- lb_aerobic_exchanges[, which(lb_aer_exch_med > 1)]
invivo_exch_med <- apply(invivo_exchanges, 2, median)
invivo_exchanges <- invivo_exchanges[, which(invivo_exch_med > 1)]
exchange_rxns <- c(colnames(invivo_exchanges), colnames(lb_aerobic_exchanges))
exchange_cpds <- c(rev(c("L-ala-D-glut-meso\ndiaminopimelate-D-ala",
                   "Glycerophosphoserine","Glycine","L-Valine")),
                   rev(c("L-ala-D-glut-meso\ndiaminopimelate",
                   "Citrate","Fructoselysine","D-Glucose 6-P","L-Glutamate",
                   "L-Leucine","L-Serine")))

# Calculate IQRs for exchange fluxes
lb_aer_exch_q25 <- apply(lb_aerobic_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.25)))
lb_aer_exch_med <- apply(lb_aerobic_exchanges, 2, median)
lb_aer_exch_q75 <- apply(lb_aerobic_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.75)))
invivo_exch_q25 <- apply(invivo_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.25)))
invivo_exch_med <- apply(invivo_exchanges, 2, median)
invivo_exch_q75 <- apply(invivo_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.75)))
exch_q25 <- c(lb_aer_exch_q25, invivo_exch_q25)
exch_median <- c(lb_aer_exch_med, invivo_exch_med)
exch_q75 <- c(lb_aer_exch_q75, invivo_exch_q75)

# Breiman (2001). Random Forests. Machine Learning.
overlap <- intersect(colnames(lb_aerobic_samples), colnames(invivo_samples))
lb_aerobic_samples <- lb_aerobic_samples[, overlap]
invivo_samples <- invivo_samples[, overlap]
lb_aerobic_samples$condition <- 'lb'
invivo_samples$condition <- 'in vivo'
all_samples <- rbind(lb_aerobic_samples, invivo_samples)
all_samples$condition <- as.factor(all_samples$condition)
sub1 <- round(length(rownames(all_samples[which(all_samples$condition=='lb'),])) * 0.623)
sub2 <- round(length(rownames(all_samples[which(all_samples$condition=='in vivo'),])) * 0.623)
fctr <- max(c(round(sub1 / sub2), round(sub2 / sub1))) * 3
n_trees <- round(ncol(all_samples) - 1) * fctr
m_tries <- round(sqrt(ncol(all_samples) - 1))
condition <- all_samples$condition
all_samples$condition <- NULL
all_samples <- as.data.frame(apply(all_samples, 2, as.numeric))
rm(sub1, sub2, fctr)

# Run Random Forest and parse results
library(randomForest)
rf_obj <- randomForest(condition ~ ., data=all_samples, ntree=n_trees, mtry=m_tries,
                       importance=TRUE, replace=FALSE, err.rate=TRUE)
print(rf_obj)
rf_feat <- importance(rf_obj, type=1, scale=FALSE)
lb_invivo_mda <- subset(rf_feat, rf_feat > (50*abs(min(rf_feat)))) # ~p-value of 0.001
lb_invivo_mda <- as.data.frame(lb_invivo_mda)
lb_invivo_mda$reaction <- rownames(lb_invivo_mda)
lb_invivo_mda <- lb_invivo_mda[order(-lb_invivo_mda$MeanDecreaseAccuracy),]
#write.table(lb_invivo_mda, file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/rf_features.tsv', row.names=FALSE, quote=FALSE, sep='\t')
rm(rf_obj, rf_feat, condition, all_samples, n_trees, m_tries)

# Read in previous machine learning results and subset data to match
lb_invivo_mda <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/rf_features.tsv', sep='\t', header=TRUE)
lb_aerobic_samples <- lb_aerobic_samples[, as.character(lb_invivo_mda$reaction)]
lb_aerobic_samples <- round(as.data.frame(apply(lb_aerobic_samples, 2, as.numeric)),3)
invivo_samples <- invivo_samples[, as.character(lb_invivo_mda$reaction)]
invivo_samples <- round(as.data.frame(apply(invivo_samples, 2, as.numeric)),3)

# Filter informative features to the largest differences
keep <- which(abs(colMeans(lb_aerobic_samples) - colMeans(invivo_samples)) > 9.0)
lb_aerobic_samples <- lb_aerobic_samples[, keep]
invivo_samples <- invivo_samples[, keep]

# Test for significant differences in flux distributions
pvals <- c()
for (x in 1:ncol(invivo_samples)) {pvals[x] <- round(wilcox.test(lb_aerobic_samples[,x], invivo_samples[,x], exact=FALSE)$p.value)}
pvals <- p.adjust(pvals, method='BH')

# Transform values with respect to sign
for (y in 1:ncol(invivo_samples)) {
  for (x in 1:nrow(invivo_samples)) {
    if (invivo_samples[x,y] < 0.0) {
      invivo_samples[x,y] <- log2(abs(invivo_samples[x,y]) + 1) * -1
    } else {
      invivo_samples[x,y] <- log2(invivo_samples[x,y] + 1)
    }
    
    if (lb_aerobic_samples[x,y] < 0.0) {
      lb_aerobic_samples[x,y] <- log2(abs(lb_aerobic_samples[x,y]) + 1) * -1
    } else {
      lb_aerobic_samples[x,y] <- log2(lb_aerobic_samples[x,y] + 1)
    }
  }
}

# Calculate summary statistics for metabolic reactions
lb_aer_met_q25 <- apply(lb_aerobic_samples, 2, function(x) as.numeric(quantile(x, probs=0.25)))
lb_aer_met_med <- apply(lb_aerobic_samples, 2, median)
lb_aer_met_q75 <- apply(lb_aerobic_samples, 2, function(x) as.numeric(quantile(x, probs=0.75)))
invivo_met_q25 <- apply(invivo_samples, 2, function(x) as.numeric(quantile(x, probs=0.25)))
invivo_met_med <- apply(invivo_samples, 2, median)
invivo_met_q75 <- apply(invivo_samples, 2, function(x) as.numeric(quantile(x, probs=0.75)))
met_median <- rbind(invivo_met_med, lb_aer_met_med)

# Assign full reaction names
rxn_names <- c('Purine\nnucleoside\nphosphorylase','Valine\ntransaminase',
               'Phospho-\npentomutase','Flavodoxin\nreductase',
               'L-alanine\ntransaminase','2-Oxogluterate\ndehydrogenase')

# Read in and format previously computed PCoA results
samples_pcoa_points <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/invitro_invivo_pcoa.tsv', header=TRUE)
x_axis_lab <- gsub('X', '', as.character(colnames(samples_pcoa_points)[1]))
x_axis_lab <- paste ('PC1 (', x_axis_lab, '%)', sep='')
y_axis_lab <- gsub('X', '', as.character(colnames(samples_pcoa_points)[2]))
y_axis_lab <- paste ('PC2 (', y_axis_lab, '%)', sep='')
vivo_points <- subset(samples_pcoa_points, condition == 'in vivo')
vitro_points <- subset(samples_pcoa_points, condition == 'in vitro')
lb_points <- subset(samples_pcoa_points, media == 'LB')
m9_aer_points <- subset(samples_pcoa_points, media == 'm9_aerobic')
m9_an_points <- subset(samples_pcoa_points, media == 'm9_anaerobic')

# Generate figure
library(scales)
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_5.png', units='in', width=7, height=8, res=300)
#pdf(file='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figure_5.pdf', width=7, height=8)
layout(matrix(c(1,1,2,2,2,
                3,3,4,4,4), nrow=2, ncol=5, byrow=TRUE))

# A
par(mar=c(4,4,1.5,1), las=1, mgp=c(2.6,1,0), xaxs='i', yaxs='i', lwd=2)
hist(transcription, main='', xlim=c(0,400), ylim=c(0,800), breaks=200, col='firebrick',
     xlab=as.expression(bquote(paste(italic('In vivo'),' Transcript Density',sep=''))), 
     ylab='Gene Frequency', cex.lab=1.2, lwd=2)
abline(v=cutoffs, lty=3, lwd=2)
#text(x=c(44.78,97.4735,116.2366,157.8262,294.2821), 
#     y=780, labels=c('5','4','3','2','1'), cex=0.75)
mtext(text=c('5','4','3','2','1'), side=3, at=c(44.78,97.4735,116.2366,157.8262,294.2821), cex=0.8)
#arrows(x0=330, x1=375, y0=780, lwd=2, length=0.05)
legend('right', pt.cex=0, bty='n',
       legend=c(as.expression(bquote(bold('Genes (Reactions)\nin ranges:'))),
                as.expression(bquote(paste(bold('1:'),' 518 (140)', sep=''))), # add reactions totals
                as.expression(bquote(paste(bold('2:'),' 517 (185)', sep=''))),
                as.expression(bquote(paste(bold('3:'),' 518 (658)', sep=''))),
                as.expression(bquote(paste(bold('4:'),' 517 (244)', sep=''))),
                as.expression(bquote(paste(bold('5:'),' 2070 (1130)', sep='')))))
box(lwd=2)
mtext('A',side=3, padj=0.5, cex=1.2, font=2, at=-80)

# B
par(mar=c(4,4,1,1.5), las=1, mgp=c(2.2,0.75,0))
plot(x=samples_pcoa_points[,2], y=samples_pcoa_points[,1], xlim=c(-0.011,0.011), ylim=c(-0.014,0.014),
     xlab=x_axis_lab, ylab=y_axis_lab, pch=19, xaxt='n', yaxt='n', cex.lab=1.4, cex=0)
axis(1, at=c(-0.01,0,0.01), labels=c('-0.01','0.0','0.01'), cex.axis=1.1, cex.lab=1.2, lwd=2)
axis(2, at=c(-0.012,0,0.012), labels=c('-0.012','0.0','0.012'), cex.axis=1.1, cex.lab=1.2, lwd=2)
points(x=vivo_points[,2], y=vivo_points[,1], bg=alpha('firebrick',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=lb_points[,2], y=lb_points[,1], bg=alpha('#ffa05d',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=m9_aer_points[,2], y=m9_aer_points[,1], bg=alpha('#4145ba',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=m9_an_points[,2], y=m9_an_points[,1], bg=alpha('white',0.75), pch=21, cex=2.4, lwd=1.5)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('in vivo'),' vs all ', italic('in vitro'), sep=''))),
                               as.expression(bquote(paste(italic('p'),' = 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('topleft', legend=c(as.expression(bquote(italic('In vivo'))),'LB aerobic','M9 aerobic','M9 anaerobic'), 
       pt.bg=c('firebrick','#ffa05d','#4145ba','white'), pch=21, pt.cex=2, pt.lwd=1.5, cex=1.1, bty='n')
legend('topright', legend='Flux Distributions of Shared Reactions', cex=1.1, pt.cex=0, bty='n')
box(lwd=2)
mtext('B',side=3, padj=0.5, cex=1.2, font=2, at=-0.014)

# C
par(mar=c(3.5,9,1,1.5), las=1, mgp=c(1.9,0.7,0), lwd=2)
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,10), ylim=c(0,length(exchange_cpds)+0.5), 
     xlab='Inv. Norm. Exchange Flux', ylab='', cex.lab=1.1)
axis(1, at=c(0,5,10), labels=c('0','30','1000'), lwd=2)
minors <- c(1.5,2.5,3.3,3.9,4.3,4.6,4.8)
axis(1, at=minors, labels=rep('', length(minors)), lwd=1, tck=-0.03)
axis(1, at=minors+5, labels=rep('', length(minors)), lwd=1, tck=-0.03)
axis(2, at=c(0.4:length(exchange_cpds)+0.4), labels=rev(exchange_cpds), tck=0, cex.axis=0.9)
bar_cols <- c(rep('#ffa05d', ncol(lb_aerobic_exchanges)), rep('firebrick', ncol(invivo_exchanges)))
barplot(exch_median, col=bar_cols, width=0.5, space=1, horiz=TRUE, add=TRUE, xaxt='n', yaxt='n')
segments(x0=exch_q25, x1=exch_q75, y0=c(0.38:length(exchange_cpds)+0.38))
abline(h=7.25)
text(x=c(8.5,9), y=c(11.2, 6.95), labels=c(as.expression(bquote(italic('in vivo'))),'LB'), cex=1.1)
mtext('C',side=3, padj=0.5, cex=1.2, font=2, at=-7.5)
par(xpd=TRUE)
text(x=10.75, y=5.5, 'Context-specific Growth Substrates', srt=270)

# D
par(mar=c(5,4,1.5,1.5), las=1, mgp=c(1.8,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(-10,10), xlim=c(0.25,11.15), 
     ylab='Normalized Flux', xlab='', xaxt='n', yaxt='n', cex.lab=1.4)
axis(2, at=c(-10,-5,0,5,10), labels=c('-1000','-30','0','30','1000'), cex.axis=0.9, lwd=2)
minors <- c(1.5,2.5,3.3,3.9,4.3,4.6,4.8)
axis(2, at=minors, labels=rep('', length(minors)), lwd=1, tck=-0.01)
axis(2, at=minors+5, labels=rep('', length(minors)), lwd=1, tck=-0.01)
axis(2, at=-minors, labels=rep('', length(minors)), lwd=1, tck=-0.01)
axis(2, at=-minors-5, labels=rep('', length(minors)), lwd=1, tck=-0.01)
abline(h=0, lty=2, lwd=1.5, col='gray40')
abline(v=c(2.1,3.9,5.7,7.5,9.3), lwd=2)
mtext('***',side=3, padj=0.5, cex=1.2, font=2, at=c(1.2,3,4.8,6.6,8.4,10.2))
legend('topleft', legend=c(as.expression(bquote(italic('in vivo'))),'LB'), 
       pt.bg=c('firebrick','#ffa05d'), pch=22, pt.cex=2.3, bty='n', pt.lwd=1.5)
barplot(met_median, col=c('firebrick','#ffa05d'), beside=TRUE, add=TRUE, 
        width=0.6, space=c(0,1), xaxt='n', yaxt='n')
segments(x0=c(1.2,3,4.8,6.6,8.4,10.2)-0.3, 
         y0=invivo_met_q25, y1=invivo_met_q75, lwd=2) 
segments(x0=c(1.2,3,4.8,6.6,8.4,10.2)+0.3, 
         y0=lb_aer_met_q25, y1=lb_aer_met_q75, lwd=2)
box(lwd=2)
mtext('D',side=3, cex=1.2, font=2, at=-1.3)
par(xpd=TRUE)
text(x=c(1.1,3,4.8,6.5,8.4,10.2), y=c(rep(-11.5,5),-12), 
     rxn_names, srt=45, cex=0.9)
text(x=11.4, y=c(4,-4), c('Increased products','Increased reactants'), srt=270, cex=1.1)
#arrows(x0=11.4, y0=c(7,-7), y1=c(8.5,-8.5), lwd=2, length=0.05)

dev.off()
