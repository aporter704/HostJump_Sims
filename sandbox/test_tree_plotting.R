library(NELSI)
library(RColorBrewer)


tr <- read.annotated.nexus('set1_MPXV_common_r0_nMig16_migAge0.19_splitAge1.91.nexus.tree')
tr
colours <- rep(brewer.pal(8, 'Set1'), length.out = length(tr$edge.length))


plot(tr, edge.col = colours, edge.width = 4)

str(tr)

