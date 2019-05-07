
library(scales)

# Read in and format data
sub_sample <- sample(c(1:10000), 500, replace=FALSE)
m9_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/M9_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
m9_aerobic_samples <- m9_aerobic_samples[, c('G3PD2','G6PDH2r','GND')]
m9_aerobic_samples <- m9_aerobic_samples[sub_sample,]
m9_aerobic_samples <- as.data.frame(apply(m9_aerobic_samples, 2, as.numeric))
lb_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/LB_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
lb_aerobic_samples <- lb_aerobic_samples[, c('G3PD2','ICDHyr')]
lb_aerobic_samples <- lb_aerobic_samples[sub_sample,]
lb_aerobic_samples <- as.data.frame(apply(lb_aerobic_samples, 2, as.numeric))
rm(sub_sample)


# Generate figures

# M9
png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_4A.png', units='in', width=4, height=5, res=300)
par(mar=c(3,4,0.5,1.5), las=1, mgp=c(2.5,0.7,0), xpd=FALSE)
plot(0, type='n', ylim=c(0,0.1), xlim=c(0,2), 
     ylab='Reaction Flux Samples', xlab='', xaxt='n', yaxt='n', cex.lab=1.2)
axis(2, at=seq(0,0.1,0.01), lwd=2)
abline(h=0, lty=5, lwd=1.5, col='gray40')
stripchart(m9_aerobic_samples$G6PDH2r, at=0.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col=alpha('#4145ba',0.75), method='jitter', jitter=0.2)
segments(0.25, median(m9_aerobic_samples$G6PDH2r), 0.75, lwd=4)
stripchart(m9_aerobic_samples$GND, at=1.5, vertical=TRUE, add=TRUE,  cex=1.3,
           pch=16, col=alpha('#4145ba',0.75), method='jitter', jitter=0.2)
segments(1.25, median(m9_aerobic_samples$GND), 1.75, lwd=4)
mtext(c('Glucose-6-P\nDehydrogenase','Phosphogluconate\nDehydrogenase'), 
      side=1, at=c(0.5,1.5), padj=1, cex=1)
par(xpd=TRUE)
text(x=2.2, y=0.05, 'Pentose Phosphate Pathway (NADPH Source)', srt=270, cex=1.1)
box(lwd=2)
dev.off()

# LB
png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_4B.png', units='in', width=4, height=5, res=300)
par(mar=c(3,4,0.5,1.5), las=1, mgp=c(2.5,0.7,0), xpd=FALSE)
plot(0, type='n', ylim=c(-100,500), xlim=c(0,1), 
     ylab='Reaction Flux Samples', xlab='', xaxt='n', yaxt='n', cex.lab=1.2)
axis(2, at=seq(-100,500,100), lwd=2)
abline(h=0, lty=5, lwd=1.5, col='gray40')
stripchart(lb_aerobic_samples$ICDHyr, at=0.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col=alpha('#ffa05d',0.75), method='jitter', jitter=0.2)
segments(0.25, median(lb_aerobic_samples$ICDHyr), 0.75, lwd=4)
mtext('Isocitrate\nDehydrogenase', side=1, at=0.5, padj=1, cex=1.1)
par(xpd=TRUE)
text(x=1.1, y=200, 'Citric Acid Cycle (NADPH Source)', srt=270, cex=1.1)
box(lwd=2)
dev.off()

# Combined glycerol-3-phosphate dehydrogenase
pdf(file='~/Desktop/repos/Cdiff_modeling/results/figure_S__.pdf', width=4, height=5)
par(mar=c(3,4,0.5,0.5), las=1, mgp=c(2.5,0.7,0))
plot(0, type='n', ylim=c(-100,600), xlim=c(0,2), 
     ylab='Reaction Flux Samples', xlab='', xaxt='n', yaxt='n', cex.lab=1.4)
axis(2, at=seq(-100,600,100), lwd=2, cex.axis=1.2)
abline(h=0, lty=5, lwd=1.5, col='gray40')
stripchart(m9_aerobic_samples$G3PD2, at=0.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col=alpha('#4145ba',0.75), method='jitter', jitter=0.2)
segments(0.25, median(m9_aerobic_samples$G3PD2), 0.75, lwd=4)
stripchart(lb_aerobic_samples$G3PD2, at=1.5, vertical=TRUE, add=TRUE, cex=1.3,
           pch=16, col=alpha('#ffa05d',0.75), method='jitter', jitter=0.2)
segments(1.25, median(lb_aerobic_samples$G3PD2), 1.75, lwd=4)
mtext(c('M9 (aerobic)', 'LB (aerobic)'), side=1, at=c(0.5,1.5), padj=1, cex=1.2)
legend('topright', legend='Glycerol-3-P Dehydrogenase', cex=1.2, pt.cex=0, bty='n')
box(lwd=2)
dev.off()

