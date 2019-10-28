
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
raw_m9 <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/m9_anaerobic/raw_m9n.iJO1366.exchange_fluxes.tsv'
objflux <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/m9_anaerobic/max.iJO1366.exchange_fluxes.tsv'
objflux_m9 <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/m9_anaerobic/max_m9n.iJO1366.exchange_fluxes.tsv'
riptide <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/m9_anaerobic/riptide.m9n.iJO1366.exchange_fluxes.tsv'

# Read in data
raw_m9 <- read.delim(raw_m9, sep='\t', header=TRUE)
objflux <- read.delim(objflux, sep='\t', header=TRUE)
objflux_m9 <- read.delim(objflux_m9, sep='\t', header=TRUE)
riptide <- read.delim(riptide, sep='\t', header=TRUE)

# Format data
raw_m9$X <- NULL
raw_m9_samples <- paste('raw_m9_', 1:nrow(raw_m9), sep='')
rownames(raw_m9) <- raw_m9_samples
objflux$X <- NULL
objflux_samples <- paste('objflux_', 1:nrow(objflux), sep='')
rownames(objflux) <- objflux_samples
objflux_m9$X <- NULL
objflux_m9_samples <- paste('objflux_m9_', 1:nrow(objflux_m9), sep='')
rownames(objflux_m9) <- objflux_m9_samples
riptide$X <- NULL
riptide_samples <- paste('riptide_', 1:nrow(riptide), sep='')
rownames(riptide) <- riptide_samples

# Subsample data
sub_sample <- sample(1:500, 250, replace=FALSE)
raw_m9 <- raw_m9[sub_sample,]
objflux <- objflux[sub_sample,]
objflux_m9 <- objflux_m9[sub_sample,]
riptide <- riptide[sub_sample,]
rm(sub_sample)

# Create metadata
raw_m9_metadata <- cbind(raw_m9_samples, rep('raw', length(raw_m9_samples)),rep('m9', length(raw_m9_samples)), rep('none', length(raw_m9_samples)))
objflux_metadata <- cbind(objflux_samples, rep('objflux', length(objflux_samples)),rep('none', length(objflux_samples)),rep('objflux', length(objflux_samples)))
objflux_m9_metadata <- cbind(objflux_m9_samples, rep('objflux_m9', length(objflux_m9_samples)),rep('m9', length(objflux_m9_samples)),rep('objflux', length(objflux_m9_samples)))
riptide_metadata <- cbind(riptide_samples, rep('riptide', length(riptide_samples)),rep('none', length(riptide_samples)),rep('riptide', length(riptide_samples)))
m9_metadata <- rbind(raw_m9_metadata, objflux_metadata, objflux_m9_metadata, riptide_metadata)
colnames(m9_metadata) <- c('label', 'group','media','method')
m9_metadata <- as.data.frame(m9_metadata)
rm(raw_m9_metadata, objflux_metadata, objflux_m9_metadata, riptide_metadata)

# Merge data
m9_flux_samples <- rbind(raw_m9, objflux, objflux_m9, riptide)
rm(raw_m9, objflux, objflux_m9, riptide)

# NMDS analysis
library(vegan)
m9_flux_samples <- m9_flux_samples + abs(min(m9_flux_samples))
m9_flux_bray_dist <- vegdist(m9_flux_samples, method='bray') # Bray-Curtis
m9_flux_nmds <- as.data.frame(metaMDS(m9_flux_bray_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(m9_flux_nmds$MDS1)) - abs(min(m9_flux_nmds$MDS1))) / 2
flux_y <- (abs(max(m9_flux_nmds$MDS2)) - abs(min(m9_flux_nmds$MDS2))) / 2
m9_flux_nmds$MDS1 <- m9_flux_nmds$MDS1 - flux_x
m9_flux_nmds$MDS2 <- m9_flux_nmds$MDS2 - flux_y
m9_x_lim <- c(min(m9_flux_nmds$MDS1)-0.004, max(m9_flux_nmds$MDS1)+0.004)
m9_x_lim <- round(m9_x_lim, digits=2)
m9_y_lim <- c(min(m9_flux_nmds$MDS2)-0.004, max(m9_flux_nmds$MDS2)+0.004)
m9_y_lim <- round(m9_y_lim, digits=2)

# Subset axes
rownames(m9_flux_nmds) <- rownames(m9_flux_samples)
raw_m9_nmds <- subset(m9_flux_nmds, rownames(m9_flux_nmds) %in% raw_m9_samples)
objflux_no_m9_nmds <- subset(m9_flux_nmds, rownames(m9_flux_nmds) %in% objflux_samples)
objflux_m9_nmds <- subset(m9_flux_nmds, rownames(m9_flux_nmds) %in% objflux_m9_samples)
m9_riptide_nmds <- subset(m9_flux_nmds, rownames(m9_flux_nmds) %in% riptide_samples)
rm(raw_m9_samples, objflux_samples, objflux_m9_samples, riptide_samples)

# Statistical testing (permANOVA)
test <- merge(x=m9_metadata, y=m9_flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
test$group <- NULL
test$method <- NULL
m9_media_pval <- adonis(m9_flux_bray_dist ~ media, data=test, perm=999, method='bray')
m9_media_pval <- m9_media_pval$aov.tab[[6]][1]
m9_media_pval <- as.character(round(m9_media_pval, 3))
rm(test)

#----------------#

# Flux sampling files
raw_lb <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/lb_aerobic/raw_lb.iJO1366.exchange_fluxes.tsv'
objflux <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/lb_aerobic/max.iJO1366.exchange_fluxes.tsv'
objflux_lb <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/lb_aerobic/max_lb.iJO1366.exchange_fluxes.tsv'
riptide <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/lb_aerobic/riptide.lb.iJO1366.exchange_fluxes.tsv'

# Read in data
raw_lb <- read.delim(raw_lb, sep='\t', header=TRUE)
objflux <- read.delim(objflux, sep='\t', header=TRUE)
objflux_lb <- read.delim(objflux_lb, sep='\t', header=TRUE)
riptide <- read.delim(riptide, sep='\t', header=TRUE)

# Format data
raw_lb$X <- NULL
raw_lb_samples <- paste('raw_lb_', 1:nrow(raw_lb), sep='')
rownames(raw_lb) <- raw_lb_samples
objflux$X <- NULL
objflux_samples <- paste('objflux_', 1:nrow(objflux), sep='')
rownames(objflux) <- objflux_samples
objflux_lb$X <- NULL
objflux_lb_samples <- paste('objflux_lb_', 1:nrow(objflux_lb), sep='')
rownames(objflux_lb) <- objflux_lb_samples
riptide$X <- NULL
riptide_samples <- paste('riptide_', 1:nrow(riptide), sep='')
rownames(riptide) <- riptide_samples

# Subsample data
sub_sample <- sample(1:500, 250, replace=FALSE)
raw_lb <- raw_lb[sub_sample,]
objflux <- objflux[sub_sample,]
objflux_lb <- objflux_lb[sub_sample,]
riptide <- riptide[sub_sample,]
rm(sub_sample)

# Create metadata
raw_lb_metadata <- cbind(raw_lb_samples, rep('raw', length(raw_lb_samples)),rep('lb', length(raw_lb_samples)), rep('none', length(raw_lb_samples)))
objflux_metadata <- cbind(objflux_samples, rep('objflux', length(objflux_samples)),rep('none', length(objflux_samples)),rep('objflux', length(objflux_samples)))
objflux_lb_metadata <- cbind(objflux_lb_samples, rep('objflux_lb', length(objflux_lb_samples)),rep('lb', length(objflux_lb_samples)),rep('objflux', length(objflux_lb_samples)))
riptide_metadata <- cbind(riptide_samples, rep('riptide', length(riptide_samples)),rep('none', length(riptide_samples)),rep('riptide', length(riptide_samples)))
lb_metadata <- rbind(raw_lb_metadata, objflux_metadata, objflux_lb_metadata, riptide_metadata)
colnames(lb_metadata) <- c('label', 'group','media','method')
lb_metadata <- as.data.frame(lb_metadata)
rm(raw_lb_metadata, objflux_metadata, objflux_lb_metadata, riptide_metadata)

# Merge data
lb_flux_samples <- rbind(raw_lb, objflux, objflux_lb, riptide)
rm(raw_lb, objflux, objflux_lb, riptide)

# NMDS analysis
lb_flux_samples <- lb_flux_samples + abs(min(lb_flux_samples))
lb_flux_bray_dist <- vegdist(lb_flux_samples, method='bray') # Bray-Curtis
lb_flux_nmds <- as.data.frame(metaMDS(lb_flux_bray_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(lb_flux_nmds$MDS1)) - abs(min(lb_flux_nmds$MDS1))) / 2
flux_y <- (abs(max(lb_flux_nmds$MDS2)) - abs(min(lb_flux_nmds$MDS2))) / 2
lb_flux_nmds$MDS1 <- lb_flux_nmds$MDS1 - flux_x
lb_flux_nmds$MDS2 <- lb_flux_nmds$MDS2 - flux_y
lb_x_lim <- c(min(lb_flux_nmds$MDS1)-0.004, max(lb_flux_nmds$MDS1)+0.004)
lb_x_lim <- round(lb_x_lim, digits=2)
lb_y_lim <- c(min(lb_flux_nmds$MDS2)-0.004, max(lb_flux_nmds$MDS2)+0.004)
lb_y_lim <- round(lb_y_lim, digits=2)

# Subset axes
rownames(lb_flux_nmds) <- rownames(lb_flux_samples)
raw_lb_nmds <- subset(lb_flux_nmds, rownames(lb_flux_nmds) %in% raw_lb_samples)
objflux_no_lb_nmds <- subset(lb_flux_nmds, rownames(lb_flux_nmds) %in% objflux_samples)
objflux_lb_nmds <- subset(lb_flux_nmds, rownames(lb_flux_nmds) %in% objflux_lb_samples)
lb_riptide_nmds <- subset(lb_flux_nmds, rownames(lb_flux_nmds) %in% riptide_samples)
rm(raw_lb_samples, objflux_samples, objflux_lb_samples, riptide_samples)

# Statistical testing (permANOVA)
test <- merge(x=lb_metadata, y=lb_flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
test$group <- NULL
test$method <- NULL
lb_media_pval <- adonis(lb_flux_bray_dist ~ media, data=test, perm=999, method='bray')
lb_media_pval <- lb_media_pval$aov.tab[[6]][1]
lb_media_pval <- as.character(round(lb_media_pval, 3))
rm(test)

rm(flux_x, flux_y, m9_flux_bray_dist, lb_flux_bray_dist, m9_flux_samples, lb_flux_samples,
   m9_metadata, lb_metadata)

#-------------------------------------------------------------------------------------------------------------------------#

# Define palette
objflux_color <- 'gray25'
raw_m9_color <- 'deeppink3'
objflux_m9_color <- 'chartreuse2'
riptide_color <- '#4347B9'
raw_lb_color <- 'brown3'
objflux_lb_color <- 'cyan3'

# Plot - ordination
library(scales)
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S5.png', units='in', width=9, height=4.5, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE))
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(2.15,0.75,0))

plot(x=m9_flux_nmds$MDS1, y=m9_flux_nmds$MDS2, xlim=m9_x_lim, ylim=m9_y_lim,
     xlab='NMDS axis 1', ylab='NMDS axis 2', cex=0, cex.axis=0.9)
points(x=objflux_no_m9_nmds$MDS1, y=objflux_no_m9_nmds$MDS2, bg=objflux_color, pch=21, cex=2, lwd=2)
points(x=raw_m9_nmds$MDS1, y=raw_m9_nmds$MDS2, bg=raw_m9_color, pch=21, cex=2, lwd=2)
points(x=objflux_m9_nmds$MDS1, y=objflux_m9_nmds$MDS2, bg=objflux_m9_color, pch=21, cex=2, lwd=2)
points(x=m9_riptide_nmds$MDS1, y=m9_riptide_nmds$MDS2, bg=riptide_color, pch=21, cex=2, lwd=2)
legend('bottomright', legend='p = 0.001', bty='n', pt.cex=0, cex=0.8)
legend('topleft', legend='M9 + glc anaerobic', bty='n', pt.cex=0)
legend('bottomleft', legend=c('>80% objective flux only','Media exchanges only','>80% objective flux + Media exchanges','RIPTiDe'), 
       pt.bg=c(objflux_color,raw_m9_color,objflux_m9_color,riptide_color), 
       pch=21, pt.cex=1.5, pt.lwd=1.5, cex=0.7, bty='n')
mtext('A', side=2, line=2, las=2, adj=1, padj=-11, cex=1.4, font=2)
box(lwd=3)

plot(x=lb_flux_nmds$MDS1, y=lb_flux_nmds$MDS2, xlim=lb_x_lim, ylim=lb_y_lim,
     xlab='NMDS axis 1', ylab='NMDS axis 2', cex=0, cex.axis=0.9)
points(x=objflux_no_lb_nmds$MDS1, y=objflux_no_lb_nmds$MDS2, bg=objflux_color, pch=21, cex=2, lwd=2)
points(x=raw_lb_nmds$MDS1, y=raw_lb_nmds$MDS2, bg=raw_lb_color, pch=21, cex=2, lwd=2)
points(x=objflux_lb_nmds$MDS1, y=objflux_lb_nmds$MDS2, bg=objflux_lb_color, pch=21, cex=2, lwd=2)
points(x=lb_riptide_nmds$MDS1, y=lb_riptide_nmds$MDS2, bg=riptide_color, pch=21, cex=2, lwd=2)
legend('bottomright', legend='p = 0.001', bty='n', pt.cex=0, cex=0.8)
legend('topleft', legend='LB aerobic', bty='n', pt.cex=0)
legend('topright', legend=c('>80% objective flux only','Media exchanges only','>80% objective flux + Media exchanges','RIPTiDe'), 
       pt.bg=c(objflux_color,raw_lb_color,objflux_lb_color,riptide_color), 
       pch=21, pt.cex=1.5, pt.lwd=1.5, cex=0.7, bty='n')
mtext('B', side=2, line=2, las=2, adj=1, padj=-11, cex=1.4, font=2)
box(lwd=3)

dev.off()


