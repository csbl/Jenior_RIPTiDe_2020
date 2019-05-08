# Parsimonious transcript data integration improves context-specific predictions of bacterial metabolism in complex environments

### Abstract

The metabolic responses of bacteria to dynamic extracellular conditions drives not only the behavior of single species, but also entire communities of microbes. Over the last decade, genome-scale metabolic network reconstructions have assisted in our appreciation of important metabolic determinants of bacterial physiology. These network models have been a powerful force in understanding the metabolic capacity that species may utilize in order to succeed in an environment. However, previous approaches to use these platforms in concert with omics data to better characterize experimental systems have met challenges due to assumptions necessary by the various integration platforms or due to large input data requirements. With these challenges in mind we developed RIPTiDe (Reaction Inclusion by Parsimony and Transcript Distribution) which uses both parsimony of overall flux and transcriptomic abundances to identify the most cost-effective usage of metabolism that also best reflects the cell’s investments into transcription. This approach is a fundamental shift from previous omic integration algorithms which sought to strictly maximize concordance in reaction usage with the given transcriptome. Utilizing a metabolic network reconstruction for the model organism *Escherichia coli* str. K-12 substr. MG1655, we found that RIPTiDe correctly identifies context-specific metabolic pathway activity without supervision. We also assessed the application of RIPTiDe to in vivo metatranscriptomic data where *E. coli* was present at high abundances, and found that our approach also correctly predicts metabolic behaviors of host-associated bacteria with high accuracy. In the context of human health, metabolic changes in bacteria have downstream impacts at multiple levels including the expression of virulence factors from pathogens or inducing disease-associated dysbiosis across the microbiota. Toward this point, RIPTiDe has far reaching potential for understanding context-specific bacterial metabolism within complex communities.

### Overview

	project
	|- README          		# the top level description of content
	|- LICENSE         		# the license for this project
	|
	|- doc/					# additional documents associated with the study
	|
	|- data/          		# raw and primary data
	| |- /reconstructions	# genome-scale metabolic network reconstructions
	| |- /flux_samples		# results from contextualized flux sampling of GENREs
	| |- /transcript		# mapped transcript abundances
	| +- /metabolome		# results from untargeted LCMS analysis
	|
	|- code/				# all programmatic code (python & R)
	|
	|- results/				# all output from workflows and analyses
	| |- figures/			# manuscript figures
	| +- tables/			# supplementary tables
	|
	|- notebooks/			# jupyter notebooks for the analyses performed during this study

