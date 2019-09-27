
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
lb_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/LB_aerobic.flux_samples.tsv'
invivo_samples <- '~/Desktop/repos/Jenior_RIPTiDe_2019/data/flux_samples/media_conditions/invivo.flux_samples.tsv'

# Read in data
lb_samples <- read.delim(lb_samples, sep='\t', header=TRUE)
invivo_samples <- read.delim(invivo_samples, sep='\t', header=TRUE)

# Subset to exchange reactions
lb_exchanges <- lb_samples[, grep('EX_', colnames(lb_samples))]
invivo_exchanges <- invivo_samples[, grep('EX_', colnames(invivo_samples))]
rm(lb_samples, invivo_samples)

# Format row names
lb_exchanges$X <- NULL
lb_names <- paste('lb_samples_', 1:nrow(lb_exchanges), sep='')
rownames(lb_exchanges) <- lb_names
invivo_exchanges$X <- NULL
invivo_names <- paste('invivo_', 1:nrow(invivo_exchanges), sep='')
rownames(invivo_exchanges) <- invivo_names

# Consider only net-imported metabolites
lb_exchanges <- lb_exchanges[, which(as.vector(colMeans(lb_exchanges)) < 0)]
invivo_exchanges <- invivo_exchanges[, which(as.vector(colMeans(invivo_exchanges)) < 0)]

# Collect differences in exchanges across models
lb_only <- setdiff(colnames(lb_exchanges), colnames(invivo_exchanges))
invivo_only <- setdiff(colnames(invivo_exchanges), colnames(lb_exchanges))
lb_exchanges <- lb_exchanges[, lb_only]
invivo_exchanges <- invivo_exchanges[, invivo_only]
rm(lb_only, invivo_only)

# Calculate absolute flux through exchanges 
lb_exchanges <- log2(abs(lb_exchanges) + 1)
invivo_exchanges <- log2(abs(invivo_exchanges) + 1)

# Subsetting to those >1
lb_exch_med <- apply(lb_exchanges, 2, median)
lb_exchanges <- lb_exchanges[, which(lb_exch_med > 1)]
invivo_exch_med <- apply(invivo_exchanges, 2, median)
invivo_exchanges <- invivo_exchanges[, which(invivo_exch_med > 1)]

# Calculate IQRs for exchange fluxes
lb_exch_q25 <- apply(lb_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.25)))
lb_exch_med <- apply(lb_exchanges, 2, median)
lb_exch_q75 <- apply(lb_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.75)))
invivo_exch_q25 <- apply(invivo_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.25)))
invivo_exch_med <- apply(invivo_exchanges, 2, median)
invivo_exch_q75 <- apply(invivo_exchanges, 2, function(x) as.numeric(quantile(x, probs=0.75)))
exch_q25 <- c(lb_exch_q25, invivo_exch_q25)
exch_median <- c(lb_exch_med, invivo_exch_med)
exch_q75 <- c(lb_exch_q75, invivo_exch_q75)

# Collect name variables
exchange_rxns <- c(colnames(invivo_exchanges), colnames(lb_exchanges))
exchange_cpds <- c(rev(c("L-Asparagine", "D-Glucose 6-phosphate", "L-Methionine S-oxide", "Nitrite", "Thymidine")),
                   rev(c("Deoxyuridine", "L-methionine", "Nitrate", "L-valine")))


# Generate figure
png(filename='~/Desktop/repos/Jenior_RIPTiDe_2019/results/figures/figure_5C.png', 
    units='in', width=5, height=6, res=300)
par(mar=c(3.5,9,1,1.5), las=1, mgp=c(1.9,0.7,0), lwd=2, xaxs='i', yaxs='i')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,10), ylim=c(0,length(exchange_cpds)+0.5), 
     xlab='Inv. Norm. Exchange Flux', ylab='', cex.lab=1.1)
axis(1, at=c(0,5,10), labels=c('0','30','1000'), lwd=2)
minors <- c(1.5,2.5,3.3,3.9,4.3,4.6,4.8)
axis(1, at=minors, labels=rep('', length(minors)), lwd=1, tck=-0.03)
axis(1, at=minors+5, labels=rep('', length(minors)), lwd=1, tck=-0.03)
axis(2, at=c(0.4:length(exchange_cpds)+0.4), labels=rev(exchange_cpds), tck=0, cex.axis=0.9)
bar_cols <- c(rep('#ffa05d', ncol(lb_exchanges)), rep('firebrick', ncol(invivo_exchanges)))
barplot(exch_median, col=bar_cols, width=0.5, space=1, horiz=TRUE, add=TRUE, xaxt='n', yaxt='n')
segments(x0=exch_q25, x1=exch_q75, y0=c(0.38:length(exchange_cpds)+0.38))
abline(h=4.25)
text(x=c(8.5,9), y=c(9.2, 3.95), labels=c(as.expression(bquote(italic('in vivo'))),'LB'), cex=1.1)
par(xpd=TRUE)
text(x=10.75, y=5.5, 'Context-specific Growth Substrates', srt=270)
dev.off()

