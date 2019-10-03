
# Start with clean environment
rm(list=ls())
gc()

# Assign objective values
unconstrained = 105.77
pfba = 86.89
m9_a = 85.9
m9_n = 85.64
lb = 86.58

# Calculate doubling time
unconstrained = (1 / unconstrained) * 3600
pfba = (1 / pfba) * 3600
m9_a = (1 / m9_a) * 3600
m9_n = (1 / m9_n) * 3600
lb = (1 / lb) * 3600

# Generate figure
library(plotrix)
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S1.png', units='in', 
    width=5, height=4, res=300)
par(mar=c(3,3.4,1.5,1), xpd=FALSE, las=1, mgp=c(2.5,0.75,0), lwd=2)
plot(0,0,type="n",xlim=c(0.5,4.5), ylim=c(41,42.2),  xaxt = 'n', xlab='', yaxt='n', 
     ylab='Calculated Doubling Time (minutes)')
axis(2, at=c(41.0,41.2,41.4,41.6,41.8,42.0,42.2), labels=c(0,41.2,41.4,41.6,41.8,42.0,42.2), lwd=2)
mtext('* Derived from biomass reaction flux levels', side=3, cex=0.8, adj=1, padj=-0.5)
mtext(c('Maximum\nPasimony','LB\nAerobic','M9 + glc\nAerobic','M9 + glc\nAnaerobic'), 
      side=1, at=c(1,2,3,4), padj=1)
legend('topleft', legend='Unconstrained iJO1366: ~34.04 minutes', pt.cex=0, cex=0.8, bty='n')
rect(xleft=0.75, ybottom=0, xright=1.25, ytop=pfba, col='#b2b2b1', lwd=3)
rect(xleft=1.75, ybottom=0, xright=2.25, ytop=lb, col='chocolate2', lwd=3)
rect(xleft=2.75, ybottom=0, xright=3.25, ytop=m9_a, col='blue3', lwd=3)
rect(xleft=3.75, ybottom=0, xright=4.25, ytop=m9_n, col='white', lwd=3)
box(lwd=2)
axis.break(2, 41.1, style='slash')
dev.off()
