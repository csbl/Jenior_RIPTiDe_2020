
# Campos M et al. (2014) A constant size extension drives bacterial cell size homeostasis. 
# Cell 159:1433–46. doi: 10.1016/j.cell.2014.11.022. p.1439 right column top paragraph
# bionumbers.hms.harvard.edu/bionumber.aspx?id=111767

library(vioplot)
library(scales)
library(vegan)
library(plotrix)

#------------------------------------------------------------------------------

# Load in data and format
growth_rates <- as.data.frame(t(read.delim(file='~/Desktop/repos/Jenior_RIPTiDe_2019/data/new_growth_rates2.tsv', sep='\t', header=FALSE, row.names=1)))
base_pfba <- as.numeric(growth_rates$base_pfba)
base_pfba  <- base_pfba[!is.na(base_pfba)]
lb_aerobic <- as.numeric(growth_rates$lb_aerobic)
lb_aerobic  <- lb_aerobic[!is.na(lb_aerobic)]
m9_gluc_aerobic <- as.numeric(growth_rates$m9_gluc_aerobic)
m9_gluc_aerobic  <- m9_gluc_aerobic[!is.na(m9_gluc_aerobic)]
m9_gluc_anaerobic <- as.numeric(growth_rates$m9_gluc_anaerobic)
m9_gluc_anaerobic  <- m9_gluc_anaerobic[!is.na(m9_gluc_anaerobic)]
rm(growth_rates)
unconstrained <- read.delim(file='~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/growth_rates/unconstrained_growth_rates.tsv', header=FALSE)[,1]
sub_sample <- sample(c(1:length(unconstrained)), length(m9_gluc_anaerobic), replace=FALSE)
unconstrained <- unconstrained[sub_sample]
unconstrained <- unconstrained[!is.na(unconstrained)]

# Calculate significant differences
pvals <- p.adjust(c(round(wilcox.test(base_pfba, m9_gluc_aerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(base_pfba, m9_gluc_anaerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(base_pfba, lb_aerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(m9_gluc_aerobic, m9_gluc_anaerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(m9_gluc_aerobic, lb_aerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(m9_gluc_anaerobic, lb_aerobic, exact=FALSE)$p.value,5)), method='BH')

# Log Transform to plot
base_pfba  <- log10(base_pfba + 1)
lb_aerobic  <- log10(lb_aerobic + 1)
m9_gluc_aerobic  <- log10(m9_gluc_aerobic + 1)
m9_gluc_anaerobic  <- log10(m9_gluc_anaerobic + 1)
unconstrained <- log10(unconstrained + 1)

# Find ceiling for data
max(c(max(unconstrained),max(base_pfba),max(lb_aerobic),max(m9_gluc_aerobic),max(m9_gluc_anaerobic)))

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_3c.png', units='in', 
    width=5.5, height=3.5, res=300)
par(mar=c(3,3.3,1.5,1), xpd=FALSE, las=1, mgp=c(2,0.75,0), lwd=2)
plot(0,0,type="n",xlim=c(0.5,6), ylim=c(0,1),  xaxt = 'n', xlab='', yaxt='n', 
     ylab=expression(paste('Calculated Doubling Time (hr'^'-1',')')))
axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c('0.0','1.2','2.4','3.6','4.8','6.0'), lwd=2) 
abline(v=2.75, lwd=1.5, lty=2)
abline(v=1.3, lwd=2)
text(x=4, y=0.08, labels='RIPTiDe-contextualized\ntranscriptomes', cex=0.9)
legend('topright', legend='Experimentally measured', 
       pch=1, cex=0.75, pt.cex=0, box.lwd=2, lwd=2, lty=3, col='firebrick2')
mtext('* Derived from biomass reaction flux levels', side=3, cex=0.8, adj=1, padj=-0.5)

# Add data
vioplot(unconstrained, at=0.8, col='darkorchid4', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(base_pfba, at=2, col='#b2b2b1', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(lb_aerobic, at=3.5, col='chocolate2', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(m9_gluc_aerobic, at=4.5, col='blue3', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(m9_gluc_anaerobic, at=5.5, col='white', lwd=2, drawRect=FALSE, add=TRUE)
mtext(c('Unconstrained\niJO1366','Maximum\nPasimony','LB\nAerobic','M9 + glc\nAerobic','M9 + glc\nAnaerobic'), 
      side=1, at=c(0.8,2,3.5,4.5,5.5), padj=1, cex=c(0.8,1,1,1,1))

# LB
segments(x0=3.4, x1=3.6, y0=0.2787536, lwd=3, col='firebrick2')
segments(x0=3.5, y0=0.2787536, y1=0.39794, lty=3, lwd=3, col='firebrick2')
segments(x0=3.4, x1=3.6, y0=0.39794, lwd=3, col='firebrick2')
# M9 - aerobic
segments(x0=4.5, y0=0.2121876, y1=0.2380461, lty=3, lwd=3, col='firebrick2')
segments(x0=4.4, x1=4.6, y0=0.2121876, lwd=3, col='firebrick2')
segments(x0=4.4, x1=4.6, y0=0.2380461, lwd=3, col='firebrick2')
# M9 - anaerobic
segments(x0=5.5, y0=0.07918125, y1=0.20412, lty=3, lwd=3, col='firebrick2')
segments(x0=5.4, x1=5.6, y0=0.07918125, lwd=3, col='firebrick2')
segments(x0=5.4, x1=5.6, y0=0.20412, lwd=3, col='firebrick2')

# Add stats
segments(x0=3.6, x1=4.4, y0=0.8, lwd=2)
text(x=4, y=0.85, '***', font=2, cex=1.5)
segments(x0=3.6, x1=5.4, y0=0.65, lwd=2)
text(x=4.5, y=0.7, '***', font=2, cex=1.5)
segments(x0=4.6, x1=5.4, y0=0.5, lwd=2)
text(x=5, y=0.55, '***', font=2, cex=1.5)

dev.off()

#------------Clean up----------------

# Clean up
rm(list=ls())
gc()

