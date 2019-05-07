#------------Read in data--------------

library(vioplot)
library(plotrix)
sub_sample <- sample(c(1:10000), 500, replace=FALSE)

# Growth rates
invivo_rates <- as.data.frame(read.delim(file='~/Desktop/repos/Jenior_RIPTiDe_2019/data/invivo_growth_rates.tsv', 
                                         sep='\t', header=FALSE))[,1][sub_sample]
pfba_rates <- as.data.frame(read.delim(file='~/Desktop/repos/Jenior_RIPTiDe_2019/data/pfba_growth_rates.tsv', 
                                         sep='\t', header=FALSE))[,1][sub_sample]

# G3PD2 flux samples
m9_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/M9_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
m9_aerobic_samples <- m9_aerobic_samples[, c('G3PD2','G6PDH2r','GND')]
m9_aerobic_samples <- m9_aerobic_samples[sub_sample,]
m9_aerobic_samples <- as.data.frame(apply(m9_aerobic_samples, 2, as.numeric))
lb_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/LB_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
lb_aerobic_samples <- lb_aerobic_samples[, c('G3PD2','ICDHyr')]
lb_aerobic_samples <- lb_aerobic_samples[sub_sample,]
lb_aerobic_samples <- as.data.frame(apply(lb_aerobic_samples, 2, as.numeric))
rm(sub_sample)

# Test differences
wilcox.test(m9_aerobic_samples$G3PD2, lb_aerobic_samples$G3PD2, exact=F)
wilcox.test(invivo_rates, pfba_rates, exact=F)

#------------Generate figure--------------

png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S1.png', 
    units='in', width=8, height=5, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

par(mar=c(3,4,0.5,0.5), las=1, mgp=c(2.5,0.7,0))
plot(0, type='n', ylim=c(-100,700), xlim=c(0,2), 
     ylab='Reaction Flux Samples', xlab='', xaxt='n', yaxt='n', cex.lab=1.4)
axis(2, at=seq(-100,700,100), lwd=2, cex.axis=1.2)
stripchart(m9_aerobic_samples$G3PD2, at=0.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col='#4145ba', method='jitter', jitter=0.2)
segments(0.25, median(m9_aerobic_samples$G3PD2), 0.75, lwd=4)
stripchart(lb_aerobic_samples$G3PD2, at=1.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col='#ffa05d', method='jitter', jitter=0.2)
segments(1.25, median(lb_aerobic_samples$G3PD2), 1.75, lwd=4)
mtext(c('M9 (aerobic)', 'LB (aerobic)'), side=1, at=c(0.5,1.5), padj=1, cex=1.2)
legend('topright', legend='Glycerol-3-P Dehydrogenase', cex=1.2, pt.cex=0, bty='n')
segments(x0=0.5, y0=550, x1=1.5, lwd=4)
text(x=1, y=580, labels='***', cex=2, font=2)
legend('bottomleft', legend='p-value << 0.001', pt.cex=0, bty='n')
abline(h=0, lty=5, lwd=1.5, col='gray40')
box(lwd=2)
mtext('A', side=3, padj=1, cex=1.5, font=2, at=-0.5)

par(mar=c(2,3.3,1,1), xpd=FALSE, las=1, mgp=c(2,0.75,0), lwd=2)
plot(0,0,type="n",xlim=c(0,2), ylim=c(0,1.5),  xaxt = 'n', xlab='', xaxt='n', 
     yaxt='n', ylab=expression(paste('Calculated Growth Rates (hr'^'-1',')')), cex.lab=1.4)
axis(2, at=c(0,0.5,1,1.5), labels=c('0.0','0.5','1.0','1.5'), lwd=2) 
vioplot(pfba_rates, at=0.5, col='gray50', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(invivo_rates, at=1.5, col='firebrick', lwd=2, drawRect=FALSE, add=TRUE)
legend('topright', legend='Experimentally measured', 
       pch=1, cex=0.75, pt.cex=0, box.lwd=2, lwd=2, lty=3, col='firebrick2')
text(x=1, y=1.2, '***', font=2, cex=2)
segments(x0=0.5, x1=1.5, y0=1.1, lwd=3)
legend('bottomleft', legend='p-value << 0.001', pt.cex=0, bty='n')
mtext('pFBA', side=1, at=0.5, padj=1, cex=1.2)
mtext('in vivo', side=1, at=1.5, padj=1, cex=1.2, font=3)
segments(x0=1.65, x1=1.75, y0=0.67, lwd=3, col='firebrick2')
segments(x0=1.7, y0=0.67, y1=1.33, lwd=3, col='firebrick2', lty=2)
segments(x0=1.65, x1=1.74, y0=1.33, lwd=3, col='firebrick2')
mtext('B',side=3, padj=0.5, cex=1.5, font=2, at=-0.4)

dev.off()

#------------Clean up----------------

# Clean up
#rm(list=ls())
#gc()
