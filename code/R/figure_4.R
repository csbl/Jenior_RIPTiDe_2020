
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
m9_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/M9_aerobic.flux_samples.tsv'
lb_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/LB_aerobic.flux_samples.tsv'

# Read in data
m9_samples <- read.delim(file=m9_samples, sep='\t', header=TRUE)
lb_samples <- read.delim(file=lb_samples, sep='\t', header=TRUE)

# G3PD2 fluxes (Glycerol-3-phosphate dehydrogenase)
m9_g3pd2 <- as.numeric(m9_samples$G3PD2)
lb_g3pd2 <- as.numeric(lb_samples$G3PD2)

# THD2pp fluxes (NADP transhydrogenase)
m9_thd2pp <- as.numeric(m9_samples$THD2pp)
# pruned from LB-specific model
rm(m9_samples, lb_samples)

# Test differences
g3pd2_pval <- round(wilcox.test(m9_g3pd2, lb_g3pd2, exact=F)[[3]], 3)

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/sub_panels/figure_4A.png', 
    units='in', width=4, height=5, res=300)
par(mar=c(3,4,0.5,2), las=1, mgp=c(2.5,0.7,0), xpd=FALSE)
plot(0, type='n', ylim=c(0,1100), xlim=c(0,2), 
     ylab='Reaction Flux', xlab='', xaxt='n', yaxt='n', cex.lab=1.4)
axis(2, at=seq(0,1000,100), lwd=2, cex.axis=1.2)
boxplot(m9_thd2pp, cex=0, lwd=4, at=0.5, col='#4145ba',ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=4, xaxt='n', yaxt='n', add=TRUE)
segments(0.25, median(m9_thd2pp), 0.75, lwd=4)
mtext(c('M9+glc (aerobic)', 'LB (aerobic)'), side=1, at=c(0.5,1.5), padj=1, cex=1.2)
legend('topright', legend='NADP Transhydrogenase', cex=1.2, pt.cex=0, bty='n')
text(x=1.5, y=50, labels='pruned', cex=1.2, font=2)
abline(h=0, lty=5, lwd=1.5, col='gray40')
box(lwd=2)
par(xpd=TRUE)
text(x=2.2, y=550, 'Primary NADPH Source', srt=270, cex=1.5)
dev.off()

png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/sub_panels/figure_4B.png', 
    units='in', width=4, height=5, res=300)
par(mar=c(3,4,0.5,2), las=1, mgp=c(2.5,0.7,0), xpd=FALSE)
plot(0, type='n', ylim=c(-100,1100), xlim=c(0,2), 
     ylab='Reaction Flux', xlab='', xaxt='n', yaxt='n', cex.lab=1.4)
axis(2, at=seq(-100,1000,100), lwd=2, cex.axis=1.2)
boxplot(m9_g3pd2, cex=0, lwd=4, at=0.5, col='#4145ba',ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=4, xaxt='n', yaxt='n', add=TRUE)
boxplot(lb_g3pd2, cex=0, lwd=4, at=1.5, col='#ffa05d',ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=4, xaxt='n', yaxt='n', add=TRUE)
mtext(c('M9+glc (aerobic)', 'LB (aerobic)'), side=1, at=c(0.5,1.5), padj=1, cex=1.2)
legend('topright', legend='Glycerol-3-P Dehydrogenase', cex=1.2, pt.cex=0, bty='n')
text(x=0.5, y=100, labels='***', cex=2, font=2)
abline(h=0, lty=5, lwd=1.5, col='gray40')
box(lwd=2)
par(xpd=TRUE)
text(x=2.2, y=550, 'Primary NADPH Source', srt=270, cex=1.5)
dev.off()

