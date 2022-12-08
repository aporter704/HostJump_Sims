setwd('/home/sebastiand/Dropbox/projects_WORKING/speciesJumpSimulations/HostJump_Sims/sandbox')
library(NELSI)


tr <- read.tree('set_test_template_r1_nMig19_migAge0.62_splitAge1.29.newick.tree')
tr

tip_ages <- data.frame(tips = tr$tip.label,
                       ages = as.numeric(gsub('.+_|\'', '', tr$tip.label)))
head(tip_ages)

# FIrst human case:
range( tip_ages$ages[grep('type1', tip_ages$tips)] )

#human is type1

# In complete tree



