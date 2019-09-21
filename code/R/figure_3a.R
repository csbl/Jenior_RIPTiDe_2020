
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
#raw <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/raw.iJO1366.exchange_fluxes.tsv'
raw_m9 <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/raw_m9.iJO1366.exchange_fluxes.tsv'
objflux <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/max.iJO1366.exchange_fluxes.tsv'
objflux_m9 <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/max_m9.iJO1366.exchange_fluxes.tsv'
riptide <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/riptide.iJO1366.exchange_fluxes.tsv'
#riptide_m9 <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/exchanges/riptide_m9.iJO1366.exchange_fluxes.tsv'

# Read in data
#raw <- read.delim(raw, sep='\t', header=TRUE)
raw_m9 <- read.delim(raw_m9, sep='\t', header=TRUE)
objflux <- read.delim(objflux, sep='\t', header=TRUE)
objflux_m9 <- read.delim(objflux_m9, sep='\t', header=TRUE)
riptide <- read.delim(riptide, sep='\t', header=TRUE)
#riptide_m9 <- read.delim(riptide_m9, sep='\t', header=TRUE)

# Format data
#raw$X <- NULL
#raw_samples <- paste('raw_', 1:nrow(raw), sep='')
#rownames(raw) <- raw_samples
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
#riptide_m9$X <- NULL
#riptide_m9_samples <- paste('riptide_m9_', 1:nrow(riptide_m9), sep='')
#rownames(riptide_m9) <- riptide_m9_samples

# Subsample data
sub_sample <- sample(1:500, 500, replace=FALSE)
#raw <- raw[sub_sample,]
raw_m9 <- raw_m9[sub_sample,]
objflux <- objflux[sub_sample,]
objflux_m9 <- objflux_m9[sub_sample,]
riptide <- riptide[sub_sample,]
#riptide_m9 <- riptide_m9[sub_sample,]
rm(sub_sample)

# Create metadata
#raw_metadata <- cbind(raw_samples, rep('raw', length(raw_samples)),rep('none', length(raw_samples)), rep('none', length(raw_samples)))
raw_m9_metadata <- cbind(raw_m9_samples, rep('raw', length(raw_m9_samples)),rep('m9', length(raw_m9_samples)), rep('none', length(raw_m9_samples)))
objflux_metadata <- cbind(objflux_samples, rep('objflux', length(objflux_samples)),rep('none', length(objflux_samples)),rep('objflux', length(objflux_samples)))
objflux_m9_metadata <- cbind(objflux_m9_samples, rep('objflux_m9', length(objflux_m9_samples)),rep('m9', length(objflux_m9_samples)),rep('objflux', length(objflux_m9_samples)))
riptide_metadata <- cbind(riptide_samples, rep('riptide', length(riptide_samples)),rep('none', length(riptide_samples)),rep('riptide', length(riptide_samples)))
#riptide_m9_metadata <- cbind(riptide_m9_samples, rep('riptide_m9', length(riptide_m9_samples)),rep('m9', length(raw_samples)),rep('riptide', length(raw_samples)))
metadata <- rbind(raw_m9_metadata, objflux_metadata, objflux_m9_metadata, riptide_metadata)
colnames(metadata) <- c('label', 'group','media','method')
metadata <- as.data.frame(metadata)
rm(raw_m9_metadata, objflux_metadata, objflux_m9_metadata, riptide_metadata)

# Merge data
flux_samples <- rbind(raw_m9, objflux, objflux_m9, riptide)
rm(raw_m9, objflux, objflux_m9, riptide)

# NMDS analysis
library(vegan)
flux_samples <- flux_samples + abs(min(flux_samples))
flux_bray_dist <- vegdist(flux_samples, method='bray') # Bray-Curtis
#flux_thetayc_dist <- designdist(flux_samples, method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
flux_nmds <- as.data.frame(metaMDS(flux_bray_dist, k=2, trymax=50)$points)

# Collect distances between points








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
#raw_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% raw_samples)
raw_m9_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% raw_m9_samples)
objflux_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% objflux_samples)
objflux_m9_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% objflux_m9_samples)
riptide_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% riptide_samples)
#riptide_m9_nmds <- subset(flux_nmds, rownames(flux_nmds) %in% riptide_m9_samples)
rm(raw_m9_samples, objflux_samples, objflux_m9_samples, riptide_samples)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
test$media <- NULL
test$method <- NULL
group_pval <- adonis(flux_bray_dist ~ group, data=test, perm=999, method='bray')
group_pval <- group_pval$aov.tab[[6]][1]
group_pval <- as.character(round(group_pval, 3))
test <- merge(x=metadata, y=flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
test$group <- NULL
test$method <- NULL
media_pval <- adonis(flux_bray_dist ~ media, data=test, perm=999, method='bray')
media_pval <- media_pval$aov.tab[[6]][1]
media_pval <- as.character(round(media_pval, 3))
test <- merge(x=metadata, y=flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
test$group <- NULL
test$media <- NULL
method_pval <- adonis(flux_bray_dist ~ method, data=test, perm=999, method='bray')
method_pval <- method_pval$aov.tab[[6]][1]
method_pval <- as.character(round(method_pval, 3))
rm(test)



#-------------------------------------------------------------------------------------------------------------------------#

# Define palette
objflux_color <- '#B1B1B0'
raw_m9_color <- 'paleturquoise1'
objflux_m9_color <- 'orangered'
riptide_color <- '#4347B9'

# Plot - ordination
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_3A.png', units='in', width=7, height=4, res=300)
par(mar=c(4,4,1,1), las=1, mgp=c(2.8,0.75,0))
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=x_lim, ylim=y_lim,
     xlab='NMDS axis 1', ylab='NMDS axis 2', cex=0, cex.axis=1.2, cex.lab=1.2)
points(x=objflux_nmds$MDS1, y=objflux_nmds$MDS2, bg=objflux_color, pch=21, cex=2, lwd=1.2)
points(x=raw_m9_nmds$MDS1, y=raw_m9_nmds$MDS2, bg=raw_m9_color, pch=21, cex=2, lwd=1.2)
points(x=objflux_m9_nmds$MDS1, y=objflux_m9_nmds$MDS2, bg=objflux_m9_color, pch=21, cex=2, lwd=1.2)
points(x=riptide_nmds$MDS1, y=riptide_nmds$MDS2, bg=riptide_color, pch=21, cex=2, lwd=1.2)
legend('topright', legend=c('Hi flux','Media only','Hi flux + media', 'RIPTiDe only'), box.lwd=2,
       pt.bg=c(objflux_color,raw_m9_color,objflux_m9_color,riptide_color), pch=21, cex=1.1, pt.cex=2)
legend('bottomright', legend='Shared Exchange Flux Distributions', bty='n', cex=1.2, pt.cex=0)
box(lwd=2)
dev.off()


