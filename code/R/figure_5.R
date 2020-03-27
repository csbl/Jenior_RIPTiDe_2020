
# Start with clean environment
rm(list=ls())
gc()

#--------------------------------------------------------------------------#
# Panel A

# Load in data
transcription <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/transcript/clinda_k12.mapped.norm.tsv', header=TRUE, row.names=1)
weights <- read.delim('~/Desktop/repos/Jenior_RIPTiDe_2019/data/clinda_k12.weights.tsv')
colnames(weights) <- c('reaction', 'weight')

# Format data
transcription <- transcription$normDepth
pruning_weights <- transcription / max(transcription)
sampling_weights <- unique(sort(weights$weight))

# Calculate density curves
transcript_density <- density(transcription)
pruning_weight_density <- density(pruning_weights)

#--------------------------------------------------------------------------#
# Panel B

# Flux sampling files
m9_aerobic_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/M9_aerobic.flux_samples.tsv'
m9_anaerobic_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/M9_anaerobic.flux_samples.tsv'
lb_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/LB_aerobic.flux_samples.tsv'
invivo_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/invivo.flux_samples.tsv'

# Read in data
m9_aerobic_samples <- read.delim(m9_aerobic_samples, sep='\t', header=TRUE)
lb_samples <- read.delim(lb_samples, sep='\t', header=TRUE)
m9_anaerobic_samples <- read.delim(m9_anaerobic_samples, sep='\t', header=TRUE)
invivo_samples <- read.delim(invivo_samples, sep='\t', header=TRUE)

# Subsample data
sub_sample <- sample(1:500, 250, replace=FALSE)
m9_aerobic_samples <- m9_aerobic_samples[sub_sample,]
lb_samples <- lb_samples[sub_sample,]
m9_anaerobic_samples <- m9_anaerobic_samples[sub_sample,]
invivo_samples <- invivo_samples[sub_sample,]
rm(sub_sample)

# Format row names
m9_aerobic_samples$X <- NULL
m9_aerobic_names <- paste('m9_aerobic_', 1:nrow(m9_aerobic_samples), sep='')
rownames(m9_aerobic_samples) <- m9_aerobic_names
lb_samples$X <- NULL
lb_names <- paste('lb_samples_', 1:nrow(lb_samples), sep='')
rownames(lb_samples) <- lb_names
m9_anaerobic_samples$X <- NULL
m9_anaerobic_names <- paste('m9_anaerobic_', 1:nrow(m9_anaerobic_samples), sep='')
rownames(m9_anaerobic_samples) <- m9_anaerobic_names
invivo_samples$X <- NULL
invivo_names <- paste('invivo_', 1:nrow(invivo_samples), sep='')
rownames(invivo_samples) <- invivo_names

# Create metadata
m9_aerobic_metadata <- cbind(m9_aerobic_names, rep('m9_aerobic', length(m9_aerobic_names)), rep('minimal', length(m9_aerobic_names)))
lb_metadata <- cbind(lb_names, rep('lb', length(lb_names)), rep('rich', length(lb_names)))
m9_anaerobic_metadata <- cbind(m9_anaerobic_names, rep('m9_anaerobic', length(m9_anaerobic_names)), rep('minimal', length(m9_anaerobic_names)))
invivo_metadata <- cbind(invivo_names, rep('invivo', length(invivo_names)), rep('rich', length(invivo_names)))
metadata <- rbind(m9_aerobic_metadata, lb_metadata, m9_anaerobic_metadata, invivo_metadata)
colnames(metadata) <- c('label', 'group', 'media')
metadata <- as.data.frame(metadata)
metadata$group <- as.factor(group)
metadata$media <- as.factor(media)
m9_metadata <- rbind(m9_aerobic_metadata, m9_anaerobic_metadata)
colnames(m9_metadata) <- c('label', 'group', 'media')
m9_metadata <- as.data.frame(m9_metadata)
m9_metadata$media <- NULL
m9_metadata$group <- as.factor(group)
rich_metadata <- rbind(lb_metadata, invivo_metadata)
colnames(rich_metadata) <- c('label', 'group', 'media')
rich_metadata <- as.data.frame(rich_metadata)
rich_metadata$media <- NULL
rich_metadata$group <- as.factor(group)
rm(m9_aerobic_metadata, lb_metadata, m9_anaerobic_metadata, invivo_metadata)

# Find overlapping reactions
shared_rxns <- intersect(colnames(invivo_samples), intersect(colnames(m9_anaerobic_samples), intersect(colnames(lb_samples), colnames(m9_aerobic_samples))))
m9_aerobic_samples <- m9_aerobic_samples[,shared_rxns]
lb_samples <- lb_samples[,shared_rxns]
m9_anaerobic_samples <- m9_anaerobic_samples[,shared_rxns]
invivo_samples <- invivo_samples[,shared_rxns]
rm(shared_rxns)

# Merge data
flux_samples <- rbind(m9_aerobic_samples, lb_samples, m9_anaerobic_samples, invivo_samples)
m9_samples <- rbind(m9_aerobic_samples, m9_anaerobic_samples)
rich_samples <- rbind(lb_samples, invivo_samples)
rm(m9_aerobic_samples, lb_samples, m9_anaerobic_samples, invivo_samples)

# PCoA analysis
library(vegan)
library(ape)
flux_samples <- flux_samples + abs(min(flux_samples))
flux_bray_dist <- vegdist(flux_samples, method='bray') # Bray-Curtis
m9_samples <- m9_samples + abs(min(m9_samples))
m9_bray_dist <- vegdist(m9_samples, method='bray') # Bray-Curtis
rich_samples <- rich_samples + abs(min(rich_samples))
rich_bray_dist <- vegdist(rich_samples, method='bray') # Bray-Curtis
pcoa_obj <- pcoa(flux_bray_dist, correction='none', rn=NULL)
pcoa_points <- as.data.frame(pcoa_obj$vectors[,c(1,2)])
pcoa_values <- pcoa_obj$values
pcoa_values <- pcoa_values$Relative_eig[c(1,2)] * 100.0
colnames(pcoa_points) <- round(pcoa_values, digits=2)
pcoa_points <- merge(pcoa_points, metadata, by.x='row.names', by.y='label')
rownames(pcoa_points) <- pcoa_points$Row.names
pcoa_points$Row.names <- NULL

# Center points
flux_x <- (abs(max(pcoa_points[,1])) - abs(min(pcoa_points[,1]))) / 2
flux_y <- (abs(max(pcoa_points[,2])) - abs(min(pcoa_points[,2]))) / 2
pcoa_points[,1] <- pcoa_points[,1] - flux_x
pcoa_points[,2] <- pcoa_points[,2] - flux_y
x_lim <- c(min(pcoa_points[,1])-0.004, max(pcoa_points[,1])+0.004)
x_lim <- round(x_lim, digits=2)
y_lim <- c(min(pcoa_points[,2])-0.004, max(pcoa_points[,2])+0.004)
y_lim <- round(y_lim, digits=2)

# Subset axes
m9_aerobic_pcoa_points <- subset(pcoa_points, rownames(pcoa_points) %in% m9_aerobic_names)
lb_samples_pcoa_points <- subset(pcoa_points, rownames(pcoa_points) %in% lb_names)
m9_anaerobic_pcoa_points <- subset(pcoa_points, rownames(pcoa_points) %in% m9_anaerobic_names)
invivo_pcoa_points <- subset(pcoa_points, rownames(pcoa_points) %in% invivo_names)
rm(m9_aerobic_names, lb_names, m9_anaerobic_names, invivo_names)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=flux_samples, by.x='label', by.y='row.names')
test$media <- NULL
group_pval <- adonis(flux_bray_dist ~ group, data=test, perm=999, method='bray')
group_pval <- group_pval$aov.tab[[6]][1]
group_pval <- as.character(round(group_pval, 3))
test <- merge(x=metadata, y=flux_samples, by.x='label', by.y='row.names')
test$group <- NULL
media_pval <- adonis(flux_bray_dist ~ media, data=test, perm=999, method='bray')
media_pval <- media_pval$aov.tab[[6]][1]
media_pval <- as.character(round(media_pval, 3))
test <- merge(x=m9_metadata, y=m9_samples, by.x='label', by.y='row.names')
m9_pval <- adonis(m9_bray_dist ~ group, data=test, perm=999, method='bray')
m9_pval <- m9_pval$aov.tab[[6]][1]
m9_pval <- as.character(round(m9_pval, 3))
test <- merge(x=rich_metadata, y=rich_samples, by.x='label', by.y='row.names')
rich_pval <- adonis(rich_bray_dist ~ group, data=test, perm=999, method='bray')
rich_pval <- rich_pval$aov.tab[[6]][1]
rich_pval <- as.character(round(rich_pval, 3))
rm(test)

# Define axis labels
x_axis_lab <- gsub('X', '', as.character(colnames(pcoa_points)[1]))
x_axis_lab <- paste ('PC1 (', x_axis_lab, '%)', sep='')
y_axis_lab <- gsub('X', '', as.character(colnames(pcoa_points)[2]))
y_axis_lab <- paste ('PC2 (', y_axis_lab, '%)', sep='')

#--------------------------------------------------------------------------#
# Panel C

# Flux sampling files
lb_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/LB_aerobic.flux_samples.tsv'
invivo_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/invivo.flux_samples.tsv'

# Read in data
lb_samples <- read.delim(lb_samples, sep='\t', header=TRUE)
invivo_samples <- read.delim(invivo_samples, sep='\t', header=TRUE)

# Subset to exchange reactions
lb_exchanges <- lb_samples[, grep('EX_', colnames(lb_samples))]
invivo_exchanges <- invivo_samples[, grep('EX_', colnames(invivo_samples))]
rm(lb_samples, invivo_samples)

# Format row names
lb_exchanges$X <- NULL
lb_names <- paste('lb_samples_', 1:nrow(lb_exchanges), sep='')
rownames(lb_exchanges) <- lb_names
invivo_exchanges$X <- NULL
invivo_names <- paste('invivo_', 1:nrow(invivo_exchanges), sep='')
rownames(invivo_exchanges) <- invivo_names

# Consider only net-imported metabolites
lb_exchanges <- lb_exchanges[, which(as.vector(colMeans(lb_exchanges)) < 0)]
invivo_exchanges <- invivo_exchanges[, which(as.vector(colMeans(invivo_exchanges)) < 0)]

# Collect differences in exchanges across models
lb_only <- setdiff(colnames(lb_exchanges), colnames(invivo_exchanges))
invivo_only <- setdiff(colnames(invivo_exchanges), colnames(lb_exchanges))
lb_exchanges <- lb_exchanges[, lb_only]
invivo_exchanges <- invivo_exchanges[, invivo_only]
rm(lb_only, invivo_only)

# Calculate absolute flux through exchanges 
lb_exchanges <- log2(abs(lb_exchanges) + 1)
invivo_exchanges <- log2(abs(invivo_exchanges) + 1)

# Subsetting to those >1
lb_exch_med <- apply(lb_exchanges, 2, median)
lb_exchanges <- lb_exchanges[, which(lb_exch_med > 1)]
invivo_exch_med <- apply(invivo_exchanges, 2, median)
invivo_exchanges <- invivo_exchanges[, which(invivo_exch_med > 1)]

# Calculate IQRs for exchange fluxes
lb_exch_q25 <- apply(lb_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.25)))
lb_exch_med <- apply(lb_exchanges, 2, median)
lb_exch_q75 <- apply(lb_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.75)))
invivo_exch_q25 <- apply(invivo_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.25)))
invivo_exch_med <- apply(invivo_exchanges, 2, median)
invivo_exch_q75 <- apply(invivo_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.75)))
exch_q25 <- c(lb_exch_q25, invivo_exch_q25)
exch_median <- c(lb_exch_med, invivo_exch_med)
exch_q75 <- c(lb_exch_q75, invivo_exch_q75)

# Collect name variables
exchange_rxns <- c(colnames(invivo_exchanges), colnames(lb_exchanges))
exchange_cpds <- c(rev(c("L-Asparagine", "D-Glucose 6-phosphate", "L-Methionine S-oxide", "Nitrite", "Thymidine")),
                   rev(c("Deoxyuridine", "L-methionine", "Nitrate", "L-valine")))

#--------------------------------------------------------------------------#
# Generate figure
library(scales)
pdf(file='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_5.pdf', width=3.54, height=3.54)
layout(matrix(c(1,2,
                3,3), nrow=2, ncol=2, byrow=TRUE))

# A
par(mar=c(3,2.2,0.5,2.2), las=1, mgp=c(2,0.45,0), xaxs='i', yaxs='i', lwd=1.5, xpd=FALSE)
hist(transcription, main='', xlim=c(0,500), ylim=c(0,800), breaks=200, col='firebrick',
     xlab='', ylab='', cex.lab=1.3, lwd=2, xaxt='n', cex.axis=0.5, tck=-0.05)
legend('topright', legend=c('Reactions: 452','Metabolites: 453','Biomass Flux: -17.74%'), cex=0.5, pt.cex=0, bty='n')
par(mgp=c(2,0.1,0))
axis(side=1, at=c(0,250,500), labels=c(0,2500,5000), lwd=1.5, cex.axis=0.5, tck=-0.05)
par(mgp=c(2,0.45,0))
axis(side=4, at=seq(0,800,160), labels=c('0.0','0.2','0.3','0.6','0.8','1.0'), lwd=1.5, cex.axis=0.5, tck=-0.05)
par(new=TRUE)
plot(rev(sampling_weights), type='l', col='dodgerblue', lwd=3, xlim=c(-120,790), xaxt='n', yaxt='n', xlab='', ylab='')
par(new=TRUE)
plot(sampling_weights, type='l', col='darkorchid4', lwd=3, xlim=c(1,790), xaxt='n', yaxt='n', xlab='', ylab='')
par(xpd=TRUE)
text(x=380, y=-0.17, 'Transcript Density', cex=0.65, font=2)
text(x=-190, y=0.5, 'Gene Frequency', srt=90, cex=0.65, col='firebrick', font=2)
text(x=1020, y=0.5, 'Reaction Weights', srt=270, cex=0.65, font=2)
text(x=970, y=0.62, 'Pruning', srt=270, cex=0.65, col='dodgerblue', font=2)
text(x=970, y=0.46, '/', srt=270, cex=0.65, font=2)
text(x=970, y=0.28, 'Sampling', srt=270, cex=0.65, col='darkorchid4', font=2)
text(x=-220, y=1, 'A', cex=1.2, font=2)
par(xpd=FALSE)
box(lwd=2)

# B
par(mar=c(2.7,2.7,0.4,0.4), las=1, mgp=c(1.6,0.55,0), lwd=1.5)
plot(x=pcoa_points[,1], y=pcoa_points[,2], xlim=c(-0.021,0.021), ylim=c(-0.023,0.021), 
     xlab=x_axis_lab, ylab=y_axis_lab, pch=19, cex.lab=0.7, cex=0, cex.axis=0.5)
points(x=invivo_pcoa_points[,2], y=invivo_pcoa_points[,1], bg=alpha('firebrick',0.75), pch=21, cex=0.8, lwd=1)
points(x=lb_samples_pcoa_points[,2], y=lb_samples_pcoa_points[,1], bg=alpha('#ffa05d',0.75), pch=21, cex=0.8, lwd=1)
points(x=m9_aerobic_pcoa_points[,2], y=m9_aerobic_pcoa_points[,1], bg=alpha('#4145ba',0.75), pch=21, cex=0.8, lwd=1)
points(x=m9_anaerobic_pcoa_points[,2], y=m9_anaerobic_pcoa_points[,1], bg=alpha('white',0.75), pch=21, cex=0.8, lwd=1)
legend('bottomleft', legend=c(as.expression(bquote(italic('in vivo'))),'LB aerobic','M9 aerobic','M9 anaerobic'), 
       pt.bg=c('firebrick','#ffa05d','#4145ba','white'), pch=21, pt.cex=0.8, pt.lwd=1, cex=0.5, bty='n')
box(lwd=1.5)
par(xpd=TRUE)
text(x=-0.032, y=0.02, 'B', cex=1.2, font=2)
text(x=-0.0115, y=0.0175, 'Core Metabolism\nDissimilarity', cex=0.5)
par(xpd=FALSE)

# C
par(mar=c(2.1,8,1,1.5), las=1, mgp=c(1.2,0.5,0), lwd=1.5, xaxs='i', yaxs='i')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,10), ylim=c(0,length(exchange_cpds)+0.5), 
     xlab=expression(paste('Inverse Exchange Flux (log'['2'],')')), ylab='', cex=0.7, cex.lab=0.8)
axis(1, at=c(0,5,10), labels=c('0','30','1000'), lwd=1.5, cex.axis=0.7)
minors <- c(1.5,2.5,3.3,3.9,4.3,4.6,4.8)
axis(1, at=minors, labels=rep('', length(minors)), tck=-0.03)
axis(1, at=minors+5, labels=rep('', length(minors)), tck=-0.03)
axis(2, at=c(0.4:length(exchange_cpds)+0.4), labels=rev(exchange_cpds), tck=0, cex.axis=0.8)
bar_cols <- c(rep('#ffa05d', ncol(lb_exchanges)), rep('firebrick', ncol(invivo_exchanges)))
barplot(exch_median, col=bar_cols, width=0.5, space=1, horiz=TRUE, add=TRUE, xaxt='n', yaxt='n')
segments(x0=exch_q25, x1=exch_q75, y0=c(0.38:length(exchange_cpds)+0.38))
abline(h=4.25)
text(x=8.75, y=5, labels=as.expression(bquote(italic('in vivo'))))
text(x=9, y=0.75, labels='LB')
par(xpd=TRUE)
text(x=5, y=10, 'Context-specific Growth Substrates', cex=0.8)
text(x=-6.4, y=10, 'C', cex=1.2, font=2)
par(xpd=FALSE)

dev.off()


