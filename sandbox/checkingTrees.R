#
library(NELSI)
setwd('~/Dropbox/projects_WORKING/speciesJumpSimulations/HostJump_Sims/sandbox/')
dir()

trees <- read.nexus('SARS_lowMinkSample.nexus.tree')
tree_heights <- sapply(trees, function(x) max(allnode.times(x, keeproot = T, reverse = T)))
tree_sizes <- sapply(trees, function(x) length(x$tip.label))

hist(tree_heights)
hist(tree_sizes)

