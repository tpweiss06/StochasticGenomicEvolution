# This script explores the probability of a single allele being sampled by at 
#    least k landscapes out of 22 (the number of structured populations) if it is
#    present in the source population at p frequency.

library(plot3D)

# First make a function to perform the calculation for a given p and k
ProbUnSampled <- function(p, k){
     # First calculate the frequency of individuals with at least one copy of
     #    the allele, assuming Hardy-Weinberg equilibrium
     f <- 2*p*(1-p) + p^2
     # Now calculate the probability of a single landscape not sampling the 
     #    allele
     m <- (1 - f) ^ 20
     # Now use that to calculate the probability of k  structured populations
     #    not sampling the allele and return that value
     return(m ^ (k))
}

ParamSpace <- matrix(NA, nrow = 21, ncol = 1000)

kSeq <- 1:21
pSeq <- seq(0, 0.05, length.out = 1000)

for(k in 1:length(kSeq)){
     for(p in 1:length(pSeq)){
          ParamSpace[k,p] <- ProbUnSampled(p = pSeq[p], k = kSeq[k])
     }
}

pdf(file = "SamplingAlleles1.pdf", width = 7, height = 4, onefile = FALSE, paper = "special")
     image2D(z = ParamSpace, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             main = "")
     axis(1, at = seq(0, 1, length.out = 21), labels = 1:21, cex.axis = 1)
     axis(2, at = seq(0, 1, by = 0.2), labels = seq(0, 0.05, by = 0.01), las = 1, cex.axis = 1)
     axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
     mtext("Landscapes not sampling allele", side = 1, line = 3, cex = 1.25)
     mtext("Allele Frequency", side = 2, line = 3, cex = 1.25)
dev.off()
