
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
lb_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/LB_aerobic.flux_samples.tsv'
invivo_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/invivo.flux_samples.tsv'

# Read in data
lb_samples <- read.delim(lb_samples, sep='\t', header=TRUE)
invivo_samples <- read.delim(invivo_samples, sep='\t', header=TRUE)

# Format data
overlap <- intersect(colnames(lb_samples), colnames(invivo_samples))
lb_samples <- lb_samples[, overlap]
invivo_samples <- invivo_samples[, overlap]
rm(overlap)

# Subset data
sub_sample <- sample(1:500, 250, replace=FALSE)
lb_samples <- lb_samples[sub_sample,]
invivo_samples <- invivo_samples[sub_sample,]
rm(sub_sample)

# Merge data
lb_samples$condition <- 1
invivo_samples$condition <- 0
all_samples <- rbind(lb_samples, invivo_samples)
all_samples$condition <- as.factor(all_samples$condition)
rm(lb_samples, invivo_samples)

# Run AUCRF and obtain feature lists
#library(AUCRF)
#set.seed(906801)
#all_aucrf <- AUCRF(condition ~ ., data=all_samples, pdel=0, k0=20)
#print(all_aucrf)
rm(all_samples)

# Assemble feature table
#top_rxns_importance <- all_aucrf$ranking[1:all_aucrf$Kopt]
#rf_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
#colnames(rf_rxns) <- c('id','mda')
#rm(all_aucrf, top_rxns_importance)

# Read in data
rf_rxns <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/MDA.tsv'
rf_rxns <- read.delim(rf_rxns, sep='\t', header=TRUE)
rf_rxns$name <- gsub('_', ' ', rf_rxns$name )

# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_S4.png', 
    units='in', width=6, height=7, res=300)
par(mar=c(3,3,1,1), xpd=FALSE, mgp=c(2,1,0), xaxt='n')
dotchart(rev(rf_rxns$mda), labels=rev(rf_rxns$name), pt.cex=1.5, bg='dodgerblue3', 
         xlab='Mean Decrease Accuracy', xlim=c(0,16), pch=21)
par(xaxt='s')
axis(1, at=seq(0,16,4), labels=seq(0,16,4), cex=1.2)
text(x=7, y=20.5, labels='in vivo vs LB aerobic')
dev.off()

