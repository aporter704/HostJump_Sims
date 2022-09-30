#install.packages("rjson")
library(rjson)

trajectory <- fromJSON(file= "~/Downloads/Hendra_lowBatSampling.json")

plot(trajectory$t, trajectory$I[[1]], type = 'l', col = 'red') #infections for human
lines(trajectory$t, trajectory$I[[2]], type = 'l', col='blue')#infections for bat
lines(trajectory$t, trajectory$I[[3]], type = 'l', col='green')#infections for horse

lines(trajectory$t, trajectory$I_sampled[[1]], type = 'l',   col='red', lty=2)#sampling for human
lines(trajectory$t, trajectory$I_sampled[[2]], type = 'l', col='blue', lty=2)#sampling for bat
lines(trajectory$t, trajectory$I_sampled[[3]], type = 'l', col='green', lty=2)#sampling for horse


#install.packages('phangorn')
library(phangorn)

#generating tree
tr <- rtree(10)
plot(tr)

#generating alignments
aln <- as.DNAbin(simSeq(tr, rate=1e-4))
aln

q = matrix(NA, 4, 4)
#transversions
q[1,2] <- 1/8 #a-c 
q[1,4] <- 1/8 #a-t
q[2,1] <- 1/8 #c-a
q[2,3] <- 1/8 #c-g
q[3,2] <- 1/8 
q[3,4] <- 1/8 
q[4,1] <- 1/8 
q[4,3] <- 1/8 
#transition
q[1,3] <- 1/4 #a-g
q[2,4] <- 1/4 #c-t
q[3,1] <- 1/4 
q[4,2] <- 1/4 

q[2, 2] <- -sum(q[2,], na.rm = T)
q[3, 3] <- -sum(q[3,], na.rm = T)
q[4, 4] <- -sum(q[4,], na.rm = T)

aln <- as.DNAbin(simSeq(tr, rate=1e-4, l = 150000))




