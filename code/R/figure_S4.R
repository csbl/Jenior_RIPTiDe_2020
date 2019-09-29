
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
lb_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/LB_aerobic.flux_samples.tsv'
m9_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/M9_aerobic.flux_samples.tsv'

# Read in data
lb_samples <- read.delim(lb_samples, sep='\t', header=TRUE)
m9_samples <- read.delim(m9_samples, sep='\t', header=TRUE)

m9_samples <- read.delim(m9_samples, sep='\t', header=TRUE)

# Format data
overlap <- intersect(colnames(lb_samples), colnames(m9_samples))
lb_samples <- lb_samples[, overlap]
m9_samples <- m9_samples[, overlap]
rm(overlap)

# Subset data
sub_sample <- sample(1:500, 250, replace=FALSE)
lb_samples <- lb_samples[sub_sample,]
m9_samples <- m9_samples[sub_sample,]
rm(sub_sample)

# Reduce data by sizes of median difference
#lb_med <- apply(lb_samples, 2, median)
#m9_med <- apply(m9_samples, 2, median)
#diffs <- abs(lb_med - m9_med)
#diffs <- names(subset(diffs, diffs > 0.0))
#lb_samples <- lb_samples[, diffs]
#m9_samples <- m9_samples[, diffs]
#rm(diffs)

# Merge data
lb_samples$condition <- 1
m9_samples$condition <- 0
all_samples <- rbind(lb_samples, m9_samples)
all_samples$condition <- as.factor(all_samples$condition)

# Run AUCRF and obtain feature lists
library(AUCRF)
set.seed(906801)
top <- 10
all_aucrf <- AUCRF(condition ~ ., data=all_samples, pdel=0, k0=top)
print(all_aucrf)
top_rxns_importance <- all_aucrf$ranking[1:all_aucrf$Kopt]
top_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
colnames(top_rxns) <- c('reaction','MDA')
rm(all_aucrf)

# Get reaction names
rxn_ids <- top_rxns$reaction
rxn_mdas <- round(as.numeric(as.character(top_rxns$MDA)), 1)
rxn_names <- c("Putrescine\ntransport", "Ribose-5-phosphate\nisomerase", 
               "Adenosine\nexchange", "Phosphopentomutase", 
               "Adenosine\ntransport\n(extracellular)","Adenosine\ntransport\n(periplasm)", 
               "Purine-nucleoside\nphosphorylase", "Dihydrofolate\nreductase", 
               "Thioredoxin\nreductase", "	Glycine\nhydroxymethyltransferase")

# Subset flux samples to top reactions
lb_samples <- lb_samples[, rxn_ids]
m9_samples <- m9_samples[, rxn_ids]

# Transform values with respect to sign
for (y in 1:ncol(m9_samples)) {
  for (x in 1:nrow(m9_samples)) {
    if (m9_samples[x,y] < 0.0) {
      m9_samples[x,y] <- log2(abs(m9_samples[x,y]) + 1) * -1
    } else {
      m9_samples[x,y] <- log2(m9_samples[x,y] + 1)
    }
    
    if (lb_samples[x,y] < 0.0) {
      lb_samples[x,y] <- log2(abs(lb_samples[x,y]) + 1) * -1
    } else {
      lb_samples[x,y] <- log2(lb_samples[x,y] + 1)
    }
  }
}

# Calculate summary statistics for metabolic reactions
lb_q25 <- apply(lb_samples, 2, function(x) as.numeric(quantile(x, probs=0.25)))
lb_med <- apply(lb_samples, 2, median)
lb_q75 <- apply(lb_samples, 2, function(x) as.numeric(quantile(x, probs=0.75)))
m9_q25 <- apply(m9_samples, 2, function(x) as.numeric(quantile(x, probs=0.25)))
m9_med <- apply(m9_samples, 2, median)
m9_q75 <- apply(m9_samples, 2, function(x) as.numeric(quantile(x, probs=0.75)))
rxn_medians <- rbind(m9_med, lb_med)
rm(m9_med, lb_med)

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S4.png', 
    units='in', width=12, height=6, res=300)
par(mar=c(6,4,1.5,1.5), las=1, mgp=c(1.8,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(-10,10), xlim=c(1,17.6), 
     ylab='Reaction Flux', xlab='', xaxt='n', cex.lab=1.4)
abline(h=0, lty=2, lwd=1.5, col='gray40')
abline(v=seq(2.1,16.5,1.8), lwd=2)
legend('topleft', legend=c(as.expression(bquote(italic('in vivo'))),'LB'), 
       pt.bg=c('firebrick','#ffa05d'), pch=22, pt.cex=2.3, bty='n', pt.lwd=1.5)
barplot(rxn_medians, col=c('firebrick','#ffa05d'), beside=TRUE, add=TRUE, 
        width=0.6, space=c(0,1), xaxt='n', yaxt='n')
segments(x0=seq(1.2,17.4,1.8)-0.3, y0=m9_q25, y1=m9_q75, lwd=2) 
segments(x0=seq(1.2,17.4,1.8)+0.3, y0=lb_q25, y1=lb_q75, lwd=2)
box(lwd=2)
par(xpd=TRUE)
text(x=seq(1.1,17.3,1.8), y=rep(-40,length(rxn_names)), rxn_names, srt=45, cex=0.8)
text(x=18.5, y=c(15,-15), c('Increased products','Increased reactants'), srt=270, cex=1.1)
dev.off()

