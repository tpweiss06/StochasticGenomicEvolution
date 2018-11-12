# Test the effect of spatial structure on the average reduction in pi values
#    for both autosomes and the X chromosome and generate the relevant confidence
#    intervals.

setwd("/home/topher/RangeExpansionGenetics/FinalAnalyses/MeanPi")
library(lme4)
library(boot)
library(pbkrtest)
PopPi <- read.csv("PiMeans.csv")
PopPi$block <- factor(PopPi$block)

# Subset and format the data
FounderData <- subset(PopPi, gen == 0)
Gen8Data <- subset(PopPi, gen == 8)
Gen8Data$location <- ifelse(is.na(Gen8Data$location), "S", Gen8Data$location)
Gen8Data$location <- as.factor(Gen8Data$location)
levels(Gen8Data$location) <- c("C", "E", "S")

# Check how many landscapes in each block and treatment we have
NumLandscapes <- matrix(NA, nrow = 3, ncol = 2)
colnames(NumLandscapes) <- c("struct", "shuf")
rownames(NumLandscapes) <- c("B2", "B3", "B4")
NumLandscapes[1,1] <- nrow(subset(FounderData, (block == 2) & (treat == "struct")))
NumLandscapes[1,2] <- nrow(subset(FounderData, (block == 2) & (treat == "shuff")))
NumLandscapes[2,1] <- nrow(subset(FounderData, (block == 3) & (treat == "struct")))
NumLandscapes[2,2] <- nrow(subset(FounderData, (block == 3) & (treat == "shuff")))
NumLandscapes[3,1] <- nrow(subset(FounderData, (block == 4) & (treat == "struct")))
NumLandscapes[3,2] <- nrow(subset(FounderData, (block == 4) & (treat == "shuff")))
NumLandscapes

# Create the data frames to hold all the results from the statistical analyses
WinSizes <- c(5000, 7500, 10000, 12500, 15000)

MeanBlockAuto <- data.frame(window = WinSizes, P = rep(NA, 5), FoundP = rep(NA, 5), 
                            EC = rep(NA,5), EClwr = rep(NA,5), ECupr = rep(NA,5),
                            ES = rep(NA, 5), ESlwr = rep(NA,5), ESupr = rep(NA,5),
                            E = rep(NA, 5), Elwr = rep(NA, 5), Eupr = rep(NA, 5), 
                            C = rep(NA, 5), Clwr = rep(NA, 5), Cupr = rep(NA, 5), 
                            S = rep(NA, 5), Slwr = rep(NA, 5), Supr = rep(NA, 5))
MeanBlockSex <- data.frame(window = WinSizes, P = rep(NA, 5), FoundP = rep(NA, 5), 
                           EC = rep(NA,5), EClwr = rep(NA,5), ECupr = rep(NA,5),
                           ES = rep(NA, 5), ESlwr = rep(NA,5), ESupr = rep(NA,5),
                           E = rep(NA, 5), Elwr = rep(NA, 5), Eupr = rep(NA, 5), 
                           C = rep(NA, 5), Clwr = rep(NA, 5), Cupr = rep(NA, 5), 
                           S = rep(NA, 5), Slwr = rep(NA, 5), Supr = rep(NA, 5))
MeanNoBlockAuto <- data.frame(window = WinSizes, P = rep(NA, 5), FoundP = rep(NA, 5), 
                              EC = rep(NA,5), EClwr = rep(NA,5), ECupr = rep(NA,5),
                              ES = rep(NA, 5), ESlwr = rep(NA,5), ESupr = rep(NA,5),
                              E = rep(NA, 5), Elwr = rep(NA, 5), Eupr = rep(NA, 5), 
                              C = rep(NA, 5), Clwr = rep(NA, 5), Cupr = rep(NA, 5), 
                              S = rep(NA, 5), Slwr = rep(NA, 5), Supr = rep(NA, 5))
MeanNoBlockSex <- data.frame(window = WinSizes, P = rep(NA, 5), FoundP = rep(NA, 5), 
                             EC = rep(NA,5), EClwr = rep(NA,5), ECupr = rep(NA,5),
                             ES = rep(NA, 5), ESlwr = rep(NA,5), ESupr = rep(NA,5),
                             E = rep(NA, 5), Elwr = rep(NA, 5), Eupr = rep(NA, 5), 
                             C = rep(NA, 5), Clwr = rep(NA, 5), Cupr = rep(NA, 5), 
                             S = rep(NA, 5), Slwr = rep(NA, 5), Supr = rep(NA, 5))

# Create a function to use for bootstrapping in the models with block effects
ModEffects <- function(fitobj){
     # Get the fixed effect coefficients from the fit object
     FixedCoefs <- summary(fitobj)$coefficients[,1]
     # Calculate the core, edge, and shuffled values
     Core <- FixedCoefs[1]
     Edge <- FixedCoefs[1] + FixedCoefs[2]
     Shuf <- FixedCoefs[1] + FixedCoefs[3]
     # Calculate the percent reduction in edge populations from core and shuffled
     EC <- (Core-Edge)/Core
     ES <- (Shuf-Edge)/Shuf
     # Return all the relevant quantities
     return(c(Core, Edge, Shuf, EC, ES))
}


# Now loop through the different window sizes to calculate all the necessary
#    quantitites
for(i in 1:5){
     # Calculate the difference in mean pi values between generation 8 and the founders
     Gen8Data$AutoDiff <- rep(NA, dim(Gen8Data)[1])
     Gen8Data$SexDiff <- rep(NA, dim(Gen8Data)[1])
     
     for(j in 1:dim(Gen8Data)[1]){
          cur_landscape <- Gen8Data$landscape[j]
          founder <- subset(FounderData, landscape == cur_landscape)
          Gen8Data$AutoDiff[j] <- Gen8Data[j,6+4*(i-1)] - founder[1,6+4*(i-1)]
          Gen8Data$SexDiff[j] <- Gen8Data[j,7+4*(i-1)] - founder[1,7+4*(i-1)]
     }
     
     ################# Now calculate the statistics without using a block effect
     ########## Autosomal results
     FullMod <- lm(AutoDiff ~ location, data = Gen8Data)
     NullMod <- lm(AutoDiff ~ 1, data = Gen8Data)
     MeanNoBlockAuto$P[i] <- anova(FullMod, NullMod)$"Pr(>F)"[2]
     MeanNoBlockAuto$EC[i] <- (-1*FullMod$coefficients[2]) / FullMod$coefficients[1]
     MeanNoBlockAuto$ES[i] <- (FullMod$coefficients[3] - FullMod$coefficients[2]) / 
          (FullMod$coefficients[1] + FullMod$coefficients[3])
     ##### CIs
     CIDat <- expand.grid(location=factor(c("C", "E", "S")))
     temp <- predict(FullMod, CIDat, interval = "confidence")
     MeanNoBlockAuto$C[i] <- temp[1,1]
     MeanNoBlockAuto$Clwr[i] <- temp[1,2]
     MeanNoBlockAuto$Cupr[i] <- temp[1,3]
     MeanNoBlockAuto$E[i] <- temp[2,1]
     MeanNoBlockAuto$Elwr[i] <- temp[2,2]
     MeanNoBlockAuto$Eupr[i] <- temp[2,3]
     MeanNoBlockAuto$S[i] <- temp[3,1]
     MeanNoBlockAuto$Slwr[i] <- temp[3,2]
     MeanNoBlockAuto$Supr[i] <- temp[3,3]
     MeanNoBlockAuto$EClwr <- (MeanNoBlockAuto$Cupr[i] - MeanNoBlockAuto$Elwr[i]) / MeanNoBlockAuto$Cupr[i]
     MeanNoBlockAuto$ECupr <- (MeanNoBlockAuto$Clwr[i] - MeanNoBlockAuto$Eupr[i]) / MeanNoBlockAuto$Clwr[i]
     MeanNoBlockAuto$ESlwr <- (MeanNoBlockAuto$Supr[i] - MeanNoBlockAuto$Elwr[i]) / MeanNoBlockAuto$Supr[i]
     MeanNoBlockAuto$ESupr <- (MeanNoBlockAuto$Slwr[i] - MeanNoBlockAuto$Eupr[i]) / MeanNoBlockAuto$Slwr[i]
     ##### Founders
     FullMod <- lm(FounderData[,6+4*(i-1)] ~ FounderData$treat)
     NullMod <- lm(FounderData[,6+4*(i-1)] ~ 1)
     MeanNoBlockAuto$FoundP[i] <- anova(FullMod, NullMod)$"Pr(>F)"[2]
     ########## Sex chromosome results
     FullMod <- lm(SexDiff ~ location, data = Gen8Data)
     NullMod <- lm(SexDiff ~ 1, data = Gen8Data)
     MeanNoBlockSex$P[i] <- anova(FullMod, NullMod)$"Pr(>F)"[2]
     MeanNoBlockSex$EC[i] <- (-1*FullMod$coefficients[2]) / FullMod$coefficients[1]
     MeanNoBlockSex$ES[i] <- (FullMod$coefficients[3] - FullMod$coefficients[2]) / 
          (FullMod$coefficients[1] + FullMod$coefficients[3])
     ##### CIs
     CIDat <- expand.grid(location=factor(c("C", "E", "S")))
     temp <- predict(FullMod, CIDat, interval = "confidence")
     MeanNoBlockSex$C[i] <- temp[1,1]
     MeanNoBlockSex$Clwr[i] <- temp[1,2]
     MeanNoBlockSex$Cupr[i] <- temp[1,3]
     MeanNoBlockSex$E[i] <- temp[2,1]
     MeanNoBlockSex$Elwr[i] <- temp[2,2]
     MeanNoBlockSex$Eupr[i] <- temp[2,3]
     MeanNoBlockSex$S[i] <- temp[3,1]
     MeanNoBlockSex$Slwr[i] <- temp[3,2]
     MeanNoBlockSex$Supr[i] <- temp[3,3]
     MeanNoBlockSex$EClwr <- (MeanNoBlockSex$Cupr[i] - MeanNoBlockSex$Elwr[i]) / MeanNoBlockSex$Cupr[i]
     MeanNoBlockSex$ECupr <- (MeanNoBlockSex$Clwr[i] - MeanNoBlockSex$Eupr[i]) / MeanNoBlockSex$Clwr[i]
     MeanNoBlockSex$ESlwr <- (MeanNoBlockSex$Supr[i] - MeanNoBlockSex$Elwr[i]) / MeanNoBlockSex$Supr[i]
     MeanNoBlockSex$ESupr <- (MeanNoBlockSex$Slwr[i] - MeanNoBlockSex$Eupr[i]) / MeanNoBlockSex$Slwr[i]
     ##### Founders
     FullMod <- lm(FounderData[,7+4*(i-1)] ~ FounderData$treat)
     NullMod <- lm(FounderData[,7+4*(i-1)] ~ 1)
     MeanNoBlockSex$FoundP[i] <- anova(FullMod, NullMod)$"Pr(>F)"[2]
     
     ########################## Now calculate the statistics with a block effect
     ########## Autosomal results
     FullMod <- lmer(AutoDiff ~ location + (1|block), data = Gen8Data)
     NullMod <- lmer(AutoDiff ~ 1 + (1|block), data = Gen8Data)
     PermTest <- PBmodcomp(FullMod, NullMod, nsim = 10000)
     MeanBlockAuto$P[i] <- PermTest$test$p.value[2]
     PointEsts <- ModEffects(FullMod)
     MeanBlockAuto$C[i] <- PointEsts[1]
     MeanBlockAuto$E[i] <- PointEsts[2]
     MeanBlockAuto$S[i] <- PointEsts[3]
     MeanBlockAuto$EC[i] <- PointEsts[4]
     MeanBlockAuto$ES[i] <- PointEsts[5]
     ##### CIs
     # Perform the bootstrap simulations
     bootpreds <- bootMer(FullMod, ModEffects, nsim = 10000)
     # Calculate the resultant intervals
     ParamCIs <- vector("list", 5)
     for (j in 1:5){
          ParamCIs[[j]] <- boot.ci(bootpreds, type = "perc", index = j)
     }
     MeanBlockAuto$Clwr <- ParamCIs[[1]]$percent[1,4]
     MeanBlockAuto$Cupr <- ParamCIs[[1]]$percent[1,5]
     MeanBlockAuto$Elwr <- ParamCIs[[2]]$percent[1,4]
     MeanBlockAuto$Eupr <- ParamCIs[[2]]$percent[1,5]
     MeanBlockAuto$Slwr <- ParamCIs[[3]]$percent[1,4]
     MeanBlockAuto$Supr <- ParamCIs[[3]]$percent[1,5]
     MeanBlockAuto$EClwr <- ParamCIs[[4]]$percent[1,4]
     MeanBlockAuto$ECupr <- ParamCIs[[4]]$percent[1,5]
     MeanBlockAuto$ESlwr <- ParamCIs[[5]]$percent[1,4]
     MeanBlockAuto$ESupr <- ParamCIs[[5]]$percent[1,5]
     ##### Founders
     FullMod <- lmer(FounderData[,6+4*(i-1)] ~ FounderData$treat + (1|FounderData$block))
     NullMod <- lmer(FounderData[,6+4*(i-1)] ~ 1+ (1|FounderData$block))
     PermTest <- PBmodcomp(FullMod, NullMod, nsim = 10000)
     MeanBlockAuto$FoundP[i] <- PermTest$test$p.value[2]
     
     ########## Sex chromosome results
     FullMod <- lmer(SexDiff ~ location + (1|block), data = Gen8Data)
     NullMod <- lmer(SexDiff ~ 1 + (1|block), data = Gen8Data)
     PermTest <- PBmodcomp(FullMod, NullMod, nsim = 10000)
     MeanBlockSex$P[i] <- PermTest$test$p.value[2]
     PointEsts <- ModEffects(FullMod)
     MeanBlockSex$C[i] <- PointEsts[1]
     MeanBlockSex$E[i] <- PointEsts[2]
     MeanBlockSex$S[i] <- PointEsts[3]
     MeanBlockSex$EC[i] <- PointEsts[4]
     MeanBlockSex$ES[i] <- PointEsts[5]
     ##### CIs
     # Perform the bootstrap simulations
     bootpreds <- bootMer(FullMod, ModEffects, nsim = 10000)
     # Calculate the resultant intervals
     ParamCIs <- vector("list", 5)
     for (j in 1:5){
          ParamCIs[[j]] <- boot.ci(bootpreds, type = "perc", index = j)
     }
     MeanBlockSex$Clwr <- ParamCIs[[1]]$percent[1,4]
     MeanBlockSex$Cupr <- ParamCIs[[1]]$percent[1,5]
     MeanBlockSex$Elwr <- ParamCIs[[2]]$percent[1,4]
     MeanBlockSex$Eupr <- ParamCIs[[2]]$percent[1,5]
     MeanBlockSex$Slwr <- ParamCIs[[3]]$percent[1,4]
     MeanBlockSex$Supr <- ParamCIs[[3]]$percent[1,5]
     MeanBlockSex$EClwr <- ParamCIs[[4]]$percent[1,4]
     MeanBlockSex$ECupr <- ParamCIs[[4]]$percent[1,5]
     MeanBlockSex$ESlwr <- ParamCIs[[5]]$percent[1,4]
     MeanBlockSex$ESupr <- ParamCIs[[5]]$percent[1,5]
     ##### Founders
     FullMod <- lmer(FounderData[,7+4*(i-1)] ~ FounderData$treat + (1|FounderData$block))
     NullMod <- lmer(FounderData[,7+4*(i-1)] ~ 1+ (1|FounderData$block))
     PermTest <- PBmodcomp(FullMod, NullMod, nsim = 10000)
     MeanBlockSex$FoundP[i] <- PermTest$test$p.value[2]
}

save(MeanNoBlockAuto, MeanNoBlockSex, MeanBlockAuto, MeanBlockSex, 
     file = "MeanPiResults.rdata")

