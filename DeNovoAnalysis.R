# This script will analyze the delta pi values associated with de novo mutations

setwd("/home/topher/RangeExpansionGenetics/FinalAnalyses/DeNovo/")
DeNovoResults <- read.csv("DeNovoResults.csv")

# First run the models for the effect of location
FullModel <- lm(DeltaPi ~ location, data = DeNovoResults)
NullModel <- lm(DeltaPi ~ 1, data = DeNovoResults)
anova(FullModel, NullModel) # p = 0.21

# Check out the model diagnostics
plot(FullModel)

# Generate 95% confidence intervals
ConfDat <- expand.grid(location=factor(c("C", "E", "S")))
temp <- predict(FullModel, ConfDat, interval='confidence') 
DeNovoInts <- cbind(ConfDat, temp)

# Save the results for plotting
save(DeNovoInts, file = "DeNovoReslts.rdata")

