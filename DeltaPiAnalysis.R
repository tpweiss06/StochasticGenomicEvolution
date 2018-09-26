# Fit linear mixed models with autocorrelated random effects to the average, 
#    standard deviation, and mu^2/sigma (SelIndex in the data frame) of the delta
#    pi values for all windows with complete observations across all landscapes.
#    The scipt calculates p values for the fixed effects and generates confidence
#    intervals for each chromosome at each level of spatial structure.

setwd("/home/topher/RangeExpansionGenetics/FinalAnalyses/DeltaPi/")
library(nlme)
load("DeltaData.rdata")

# Create objects to hold the output values
Pvals <- expand.grid(statistic = c("mu", "sigma"), WinSize = c(5000, 7500, 10000, 12500, 15000))
pIntercept <- rep(NA, nrow(Pvals))
pLoc <- rep(NA, nrow(Pvals))
pChrom <- rep(NA, nrow(Pvals))
pInteraction <- rep(NA, nrow(Pvals))
Pvals <- cbind(Pvals, pIntercept, pLoc, pChrom, pInteraction)
MuInts <- vector(length = 5, mode = "list")
SigmaInts <- vector(length = 5, mode = "list")
MuPercChange <- vector(length = 5, mode = "list")
SigmaPercChange <- vector(length = 5, mode = "list")

# Run through each window size, fit the models, record the p values, and
#    calculate the relevant confidence intervals
for(i in 1:5){
     # Prepare the data for subsequent analyses by converting the chromosome
     #    column to a factor and adding a useful interaction term
     CurData <- DeltaData[[i]]
     CurData$Chrom <- as.factor(CurData$Chrom)
     CurData$ChromLoc <- interaction(CurData$Chrom, CurData$location)
     
     # Fit the standard models for the mean and standard deviation
     MuMod <- lme(fixed = mu ~ location + Chrom + location*Chrom, data = CurData,
                  random = ~ 1|window, correlation = corExp(), method = "ML")
     SigmaMod <- lme(fixed = sigma ~ location + Chrom + location*Chrom, data = CurData,
                     random = ~ 1|window, correlation = corExp(), method = "ML")
     
     # Calculate and extract p values for the model terms
     MuAnova <- anova(MuMod, test = "LRT")
     SigmaAnova <- anova(SigmaMod, test = "LRT")
     Pvals[1 + 2*(i-1), 3:6] <- MuAnova$'p-value'
     Pvals[2 + 2*(i-1), 3:6] <- SigmaAnova$'p-value'
     
     # Now refit the models with the interaction term to calculate confidence
     #    intervals
     MuMod <- lme(fixed = mu ~ -1 + ChromLoc, data = CurData, random = ~ 1|window, 
                  correlation = corExp(), method = "ML")
     SigmaMod <- lme(fixed = sigma ~ -1 + ChromLoc, data = CurData, random = ~ 1|window, 
                  correlation = corExp(), method = "ML")
     MuInts[[i]] <- intervals(MuMod, which = "fixed")
     SigmaInts[[i]] <- intervals(SigmaMod, which = "fixed")
     
     # Now calculate the percent change in edge values from core and shuffled
     #    values
     CoreVals <- 1:10
     EdgeVals <- 11:20
     ShufVals <- 21:30
     MuEdgeCore <- (MuInts[[i]]$fixed[CoreVals,2] - MuInts[[i]]$fixed[EdgeVals,2]) / 
          MuInts[[i]]$fixed[CoreVals,2]
     MuEdgeShuf <- (MuInts[[i]]$fixed[ShufVals,2] - MuInts[[i]]$fixed[EdgeVals,2]) / 
          MuInts[[i]]$fixed[ShufVals,2]
     SigmaEdgeCore <- (SigmaInts[[i]]$fixed[CoreVals,2] - SigmaInts[[i]]$fixed[EdgeVals,2]) / 
          SigmaInts[[i]]$fixed[CoreVals,2]
     SigmaEdgeShuf <- (SigmaInts[[i]]$fixed[ShufVals,2] - SigmaInts[[i]]$fixed[EdgeVals,2]) / 
          SigmaInts[[i]]$fixed[ShufVals,2]
     MuPercChange[[i]] <- list(EdgeCore = MuEdgeCore, EdgeShuf = MuEdgeShuf)
     SigmaPercChange[[i]] <- list(EdgeCore = SigmaEdgeCore, EdgeShuf = SigmaEdgeShuf)
}

save(Pvals, MuInts, SigmaInts, MuPercChange, SigmaPercChange,
     file = "DeltaPiResults.rdata")
