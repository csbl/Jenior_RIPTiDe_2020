
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
raw_m9 <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/raw_m9.iJO1366.exchange_fluxes.tsv'
objflux <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/max.iJO1366.exchange_fluxes.tsv'
objflux_m9 <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/max_m9.iJO1366.exchange_fluxes.tsv'
riptide <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/riptide.iJO1366.exchange_fluxes.tsv'

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
metadata <- rbind(raw_m9_metadata, objflux_metadata, objflux_m9_metadata, riptide_metadata)
colnames(metadata) <- c('label', 'group','media','method')
metadata <- as.data.frame(metadata)
sub_metadata <- rbind(raw_m9_metadata, objflux_m9_metadata, riptide_metadata)
colnames(sub_metadata) <- c('label', 'group','media','method')
sub_metadata <- as.data.frame(sub_metadata)
rm(raw_m9_metadata, objflux_metadata, objflux_m9_metadata, riptide_metadata)

# Merge data
flux_samples <- rbind(raw_m9, objflux, objflux_m9, riptide)
sub_flux_samples <- rbind(raw_m9, objflux_m9, riptide)
rm(raw_m9, objflux, objflux_m9, riptide)

# NMDS analysis
library(vegan)
flux_samples <- flux_samples + abs(min(flux_samples))
flux_bray_dist <- vegdist(flux_samples, method='bray') # Bray-Curtis
sub_flux_samples <- sub_flux_samples + abs(min(sub_flux_samples))
sub_flux_bray_dist <- vegdist(sub_flux_samples, method='bray') # Bray-Curtis
#flux_thetayc_dist <- designdist(flux_samples, method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
flux_nmds <- as.data.frame(metaMDS(flux_bray_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
x_lim <- c(min(flux_nmds$MDS1)-0.004, max(flux_nmds$MDS1)+0.004)
x_lim <- round(x_lim, digits=2)
y_lim <- c(min(flux_nmds$MDS2)-0.004, max(flux_nmds$MDS2)+0.004)
y_lim <- round(y_lim, digits=2)

# Subset axes
rownames(flux_nmds) <- rownames(flux_samples)
raw_m9_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% raw_m9_samples)
objflux_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% objflux_samples)
objflux_m9_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% objflux_m9_samples)
riptide_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% riptide_samples)
rm(raw_m9_samples, objflux_samples, objflux_m9_samples, riptide_samples)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
test$group <- NULL
test$method <- NULL
media_pval <- adonis(flux_bray_dist ~ media, data=test, perm=999, method='bray')
media_pval <- media_pval$aov.tab[[6]][1]
media_pval <- as.character(round(media_pval, 3))
rm(test)
test <- merge(x=sub_metadata, y=sub_flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
test$group <- NULL
test$method <- NULL
sub_media_pval <- adonis(sub_flux_bray_dist ~ media, data=test, perm=999, method='bray')
sub_media_pval <- sub_media_pval$aov.tab[[6]][1]
sub_media_pval <- as.character(round(sub_media_pval, 3))
rm(test)


# Create supplementary table
table_S5 <- as.matrix(flux_bray_dist)
write.table(table_S5, file='~/Desktop/repos/Jenior_RIPTiDe_2019/results/tables/table_S5.tsv', 
            row.names=TRUE, col.names=TRUE, quote=FALSE, sep='\t')


#-------------------------------------------------------------------------------------------------------------------------#

# Define palette
objflux_color <- 'gray25'
raw_m9_color <- 'deeppink3'
objflux_m9_color <- 'chartreuse2'
riptide_color <- '#4347B9'

# Plot - ordination
library(scales)
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_3A.png', units='in', width=7, height=3.5, res=300)
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(2.15,0.75,0))
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=x_lim, ylim=y_lim,
     xlab='NMDS axis 1', ylab='NMDS axis 2', cex=0, cex.axis=0.9)
points(x=objflux_nmds$MDS1, y=objflux_nmds$MDS2, bg=objflux_color, pch=21, cex=2, lwd=2)
points(x=raw_m9_nmds$MDS1, y=raw_m9_nmds$MDS2, bg=raw_m9_color, pch=21, cex=2, lwd=2)
points(x=objflux_m9_nmds$MDS1, y=objflux_m9_nmds$MDS2, bg=objflux_m9_color, pch=21, cex=2, lwd=2)
points(x=riptide_nmds$MDS1, y=riptide_nmds$MDS2, bg=riptide_color, pch=21, cex=2, lwd=2)
legend('bottomright', legend='Shared Exchange Reaction Flux Dissimilarity', bty='n', pt.cex=0, cex=0.8)
box(lwd=3)
dev.off()


