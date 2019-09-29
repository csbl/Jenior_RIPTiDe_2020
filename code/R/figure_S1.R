
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
m9_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/M9_aerobic.flux_samples.tsv'
lb_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/LB_aerobic.flux_samples.tsv'

# Read in data
m9_samples <- read.delim(file=m9_samples, sep='\t', header=TRUE)
lb_samples <- read.delim(file=lb_samples, sep='\t', header=TRUE)

# Biomass flux
m9_growth <- as.numeric(m9_samples$BIOMASS_Ec_iJO1366_WT_53p95M)
lb_growth <- as.numeric(lb_samples$BIOMASS_Ec_iJO1366_WT_53p95M)

# G3PD2 flux samples
m9_samples <- as.numeric(m9_samples$G3PD2)
lb_samples <- as.numeric(lb_samples$G3PD2)

# Convert to growth rates
m9_growth <- (1/m9_growth) * 3600
lb_growth <- (1/lb_growth) * 3600

# Test differences
wilcox.test(m9_samples, lb_samples, exact=F)
wilcox.test(m9_growth, lb_growth, exact=F)

# Find maxima for plotting
flux_max <- max(c(max(m9_samples), max(lb_samples)))
growth_max <- max(c(max(m9_growth), max(lb_growth)))

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S1.png', 
    units='in', width=4, height=5, res=300)
par(mar=c(3,4,0.5,0.5), las=1, mgp=c(2.5,0.7,0))
plot(0, type='n', ylim=c(-100,1000), xlim=c(0,2), 
     ylab='Reaction Flux Samples', xlab='', xaxt='n', yaxt='n', cex.lab=1.4)
axis(2, at=seq(-100,1000,100), lwd=2, cex.axis=1.2)
stripchart(m9_samples, at=0.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col='#4145ba', method='jitter', jitter=0.2)
segments(0.25, median(m9_samples), 0.75, lwd=4)
stripchart(lb_samples, at=1.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col='#ffa05d', method='jitter', jitter=0.2)
segments(1.25, median(lb_samples), 1.75, lwd=4)
mtext(c('M9 (aerobic)', 'LB (aerobic)'), side=1, at=c(0.5,1.5), padj=1, cex=1.2)
legend('topleft', legend='Glycerol-3-P Dehydrogenase', cex=1.2, pt.cex=0, bty='n')
text(x=0.5, y=100, labels='***', cex=2, font=2)
abline(h=0, lty=5, lwd=1.5, col='gray40')
box(lwd=2)
dev.off()

