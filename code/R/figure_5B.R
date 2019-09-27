
# Start with clean environment
rm(list=ls())
gc()

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

# NMDS analysis
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

# Generate figure
library(scales)
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_5B.png', units='in', width=6, height=6, res=300)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.15,0.5,0))
plot(x=pcoa_points[,1], y=pcoa_points[,2], xlim=x_lim, ylim=y_lim,
     xlab=x_axis_lab, ylab=y_axis_lab, pch=19, cex.lab=1.4, cex=0)
points(x=invivo_pcoa_points[,2], y=invivo_pcoa_points[,1], bg=alpha('firebrick',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=lb_samples_pcoa_points[,2], y=lb_samples_pcoa_points[,1], bg=alpha('#ffa05d',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=m9_aerobic_pcoa_points[,2], y=m9_aerobic_pcoa_points[,1], bg=alpha('#4145ba',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=m9_anaerobic_pcoa_points[,2], y=m9_anaerobic_pcoa_points[,1], bg=alpha('white',0.75), pch=21, cex=2.4, lwd=1.5)
legend('topleft', legend=c(as.expression(bquote(italic('In vivo'))),'LB aerobic','M9 aerobic','M9 anaerobic'), 
       pt.bg=c('firebrick','#ffa05d','#4145ba','white'), pch=21, pt.cex=2, pt.lwd=1.5, cex=1.1, bty='n')
legend('bottomright', legend='Flux Distributions of Shared Reactions', cex=0.9, pt.cex=0, bty='n')
box(lwd=2)
dev.off()
