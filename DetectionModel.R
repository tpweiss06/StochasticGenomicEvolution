# This script will perform a mechanistic simulation of allele sampling during
#    the processes of landscape founding, pooling, and sequencing to determine
#    the probability of identifying a rare allele as a function of its initial
#    frequency in the source population.
#    Adjusted from Silas Tittes' script

library(tidyverse)

# Read in the observed sampling depths
dep <- read_delim(
     "pop_depths_freq-table.txt", 
     col_names = F, 
     delim = " ") %>% 
     gather(pop, count, -X1)
Z <- read.csv("Z.csv") %>% 
     gather(binary, land, -prefix) %>%
     filter(land == 1) %>%
     arrange(prefix) %>%
     mutate(pop = sort(unique(dep$pop)))

dep_join <- full_join(dep, Z, "pop")

mean_founder <- dep_join %>% 
     filter(binary == "founder") %>%
     group_by(X1) %>%
     summarise(mean_dep = mean(count)) %>%
     filter(X1 <= 22 & X1 >= 4)

n_founders <- sum(Z$binary == "founder")
n_pools <- n_founders
#n_pools <- 4
n_beetles <- 20
chromosomes <- 2*n_beetles
freq_rare <- 0.01

#generate founder alleles according to defined frequency
pop_sample <- function(n_pools, chromosomes, freq_rare){
     PoolAlleles <- vector(mode = "list", length = n_pools)
     for(i in 1:n_pools){
          PoolAlleles[[i]] <- sample(c("A", "a"), size = chromosomes,
                                     prob = c(1 - freq_rare, freq_rare),
                                     replace = TRUE)
     }
     return(PoolAlleles)
}

#For each pool, draw sample depths from empirical founder frequencies
depth_sample <- function(mean_founder, n_pools){
     SampleDepths <- sample(mean_founder$X1, size = n_pools, prob = mean_founder$mean_dep,
                    replace = TRUE)
     return(SampleDepths)
}

#sample founder alleles to simulate pooling
pool_sample <- function(PoolAlleles, SampleDepths){
     DomAlleleFreqs <- rep(NA, length(PoolAlleles))
     for(i in 1:length(PoolAlleles)){
          SampledAlleles <- sample(PoolAlleles[[i]], size = SampleDepths[i],
                                   replace = TRUE)
          DomAlleleFreqs[i] <- mean(SampledAlleles == "A")
     }
     return(DomAlleleFreqs)
}

# Now perform the simulation
Nsim <- 10000
Nfreqs <- 50
InitFreqs <- seq(5e-5, 0.01, length.out = Nfreqs)
est <- rep(NA, Nfreqs)
lwr <- rep(NA, Nfreqs)
upr <- rep(NA, Nfreqs)
DetectionProbability <- data.frame(InitFreq = InitFreqs, est = est, lwr = lwr, upr = upr)

for(i in 1:Nfreqs){
     SimPops <- rep(NA, Nsim)
     for(j in 1:Nsim){
          PoolAlleles <- pop_sample(n_pools = n_pools, chromosomes = chromosomes, 
                                    freq_rare = DetectionProbability$InitFreq[i])
          DepthSample <- depth_sample(mean_founder = mean_founder, n_pools = n_pools)
          FreqPools <- pool_sample(PoolAlleles = PoolAlleles, SampleDepths = DepthSample)
          SimPops[j] <- sum(FreqPools < 1)
     }
     MeanProb <- mean(SimPops != 0)
     IntVal <- 1.96 * sqrt((MeanProb * (1 - MeanProb)) / Nsim)
     DetectionProbability$est[i] <- MeanProb
     DetectionProbability$lwr[i] <- MeanProb - IntVal
     DetectionProbability$upr[i] <- MeanProb + IntVal
     print( i / Nfreqs)
}

pdf(file = "SamplingProbability.pdf", width = 7, height = 5, onefile = FALSE, paper = "special")
     plot(DetectionProbability$InitFreq, DetectionProbability$est, las = 1,
          xlab = "Allele frequency", ylab = "Detection probability", 
          main = "", cex.axis = 1.25, cex.lab = 1.25)
     segments(x0 = DetectionProbability$InitFreq, y0 = DetectionProbability$lwr,
              x1 = DetectionProbability$InitFreq, y1 = DetectionProbability$upr)
dev.off()

# Figure out at which allele frequency value does the detection probability cross
#    the 95% probability threshold
AboveThreshold <- which(DetectionProbability$est >= 0.95)
DetectionProbability[(AboveThreshold[1] - 1):AboveThreshold[2],]
