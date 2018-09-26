# Bootstrap the p value for the effect of location on the correlation of eighth
#    generation pi values among replicates. 

setwd("/home/topher/RangeExpansionGenetics/FinalAnalyses/DeltaPi/")
load("PiWinVals.rdata")
reps <- 10000

# Create a data frame to hold the p values and confidence intervals for
#    each window size
PermResults <- expand.grid(location = c("core", "edge", "shuf"), 
                           WinSize = c(5000, 7500, 10000, 12500, 15000))
est <- rep(NA, 15)
lwr <- rep(NA, 15)
upr <- rep(NA, 15)
pVal <- rep(NA, 15)
PermResults <- cbind(PermResults, est, lwr, upr, Pval)

# Create a list to hold all the rho vals for later graphing
RhoVals <- vector(mode = "list", length = 5)

# Now loop through each window size and fill in the results data frame
for(i in 1:5){
     CurPiMat <- PiValsMats[[i]]
     # Calculate the correlation matrix for each sub group and the full set
     FullRhoMat <- cor(CurPiMat, use = "complete.obs")
     FoundRhoMat <- cor(CurPiMat[,c(ShuffFound, StructFound)], use = "complete.obs")
     CoreRhoMat <- cor(CurPiMat[,CorePops], use = "complete.obs")
     EdgeRhoMat <- cor(CurPiMat[,EdgePops], use = "complete.obs")
     ShufRhoMat <- cor(CurPiMat[,ShufPops], use = "complete.obs")
     # Extract the rho values as vectors and then remove the matrices in the
     #    interest of space
     FoundRho <- FoundRhoMat[lower.tri(FoundRhoMat)]
     CoreRho <- CoreRhoMat[lower.tri(CoreRhoMat)]
     EdgeRho <- EdgeRhoMat[lower.tri(EdgeRhoMat)]
     ShufRho <- ShufRhoMat[lower.tri(ShufRhoMat)]
     rm(CurPiMat, FoundRhoMat, CoreRhoMat, EdgeRhoMat, ShufRhoMat)
     
     # Save the rho values in the list
     RhoVals[[i]] <- list(Found = FoundRho, Core = CoreRho, Edge = EdgeRho, Shuf = ShufRho)
     
     # Save the mean values for each group
     ObsTrtMeans <- c(mean(CoreRhoVals), mean(EdgeRhoVals), mean(ShufRhoVals))

     # Calculate the mean absolute difference between groups for the test statistic
     obsMAD <- sum(abs(c(diff(ObsTrtMeans), diff(ObsTrtMeans, lag = 2)))) / 3
     
     # Perform the permutation test accounting for the dependence among values
     MAD <- rep(NA, reps)
     for(j in 1:reps){
          # randomize the treatments
          Rand <- sample(1:nrow(FullRhoMat))
          RandMat <- FullRhoMat[Rand,]
          RandMat <- RandMat[,Rand]
          CoreRandMat <- randRhoMat[1:22, 1:22]
          EdgeRandMat <- randRhoMat[23:44, 23:44]
          ShufRandMat <- randRhoMat[45:59, 45:59]
          CoreRandVals <- CoreRandMat[lower.tri(CoreRandMat)]
          EdgeRandVals <- EdgeRandMat[lower.tri(EdgeRandMat)]
          ShufRandVals <- ShufRandMat[lower.tri(ShufRandMat)]
          
          # Calculate the test statistic for the randomized values
          RandTrtMeans <- c(mean(CoreRandVals), mean(EdgeRandVals), mean(ShufRandVals))
          MAD[i] <- sum(abs(c(diff(RandTrtMeans), diff(RandTrtMeans, lag=2)))) / 3
     }
     # Calculate how many simulations have a value greater than the observed
     PermResults$pVal[1:3 + 3*(i-1)] <- sum(MAD >= obsMAD) / reps
     
     # Now calculate the percent reduction in correlations from founders to
     #    each treatment group, along with confidence intervals
     FoundMean <- mean(FoundRho)
     PermResults$est[1 + 3*(i-1)] <- (FoundMean - ObsTrtMeans[1]) / FoundMean
     PermResults$est[2 + 3*(i-1)] <- (FoundMean - ObsTrtMeans[2]) / FoundMean
     PermResults$est[3 + 3*(i-1)] <- (FoundMean - ObsTrtMeans[3]) / FoundMean
     
     # Now, permute the data again to calculate bootstrapped confidence intervals
     SimCoreReduct <- rep(NA, reps)
     SimEdgeReduct <- rep(NA, reps)
     SimShufReduct <- rep(NA, reps)
     for(j in 1:reps){
          SimEdge <- sample(EdgeRho, size = length(EdgeRho), replace = TRUE)
          SimCore <- sample(CoreRho, size = length(CoreRho), replace = TRUE)
          SimShuf <- sample(ShufRho, size = length(ShufRho), replace = TRUE)
          SimFound <- sample(FounderRho, size = length(FounderRho), replace = TRUE)
          
          SimCoreReduct[j] <- (mean(SimFound) - mean(SimCore)) / mean(SimFound)
          SimEdgeReduct[j] <- (mean(SimFound) - mean(SimEdge)) / mean(SimFound)
          SimShufReduct[j] <- (mean(SimFound) - mean(SimShuf)) / mean(SimFound)
     }
     # Calculate the delta statistic for each (bootstrap - empirical)
     CoreDelta <- SimCoreReduct - PermResults$est[1 + 3*(i-1)]
     EdgeDelta <- SimEdgeReduct - PermResults$est[2 + 3*(i-1)]
     ShufDelta <- SimShufReduct - PermResults$est[3 + 3*(i-1)]
     
     # Finally, use that to calculate the 95% confidence intervals
     CoreDeltaQuant <- quantile(CoreDelta, c(0.025, 0.975))
     EdgeDeltaQuant <- quantile(EdgeDelta, c(0.025, 0.975))
     ShufDeltaQuant <- quantile(ShufDelta, c(0.025, 0.975))
     PermResults$lwr[1 + 3*(i-1)] <- PermResults$est[1 + 3*(i-1)] + CoreDeltaQuant[1]
     PermResults$upr[1 + 3*(i-1)] <- PermResults$est[1 + 3*(i-1)] + CoreDeltaQuant[2]
     PermResults$lwr[2 + 3*(i-1)] <- PermResults$est[2 + 3*(i-1)] + EdgeDeltaQuant[1]
     PermResults$upr[2 + 3*(i-1)] <- PermResults$est[2 + 3*(i-1)] + EdgeDeltaQuant[2]
     PermResults$lwr[3 + 3*(i-1)] <- PermResults$est[3 + 3*(i-1)] + ShufDeltaQuant[1]
     PermResults$upr[3 + 3*(i-1)] <- PermResults$est[3 + 3*(i-1)] + ShufDeltaQuant[2]
}

save(PermResults, RhoVals, file = "PermTestResults.rdata")