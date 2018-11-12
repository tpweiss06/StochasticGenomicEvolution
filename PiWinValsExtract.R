# Extract the pi values at individual windows (for each window size) and save
#    them for downstream analysis

setwd("/mnt/md0raid1_2tb/SpatialStructSeqData/IndPileUps")
FileNames <- read.table("file_names.txt")
NumPops <- dim(FileNames)[1]
struct <- c(31:40, 51:60, 71:80)
SamplePrefix <- strsplit(as.character(FileNames$V1[1]), ".pi")[[1]]

# Create lists to hold all the values for each window value
PiValsMats <- vector(mode = "list", length = 5)
EdgePops <- NULL
CorePops <- NULL
ShufPops <- NULL
StructFound <- NULL
ShuffFound <- NULL

# Get chromosome information in separate vector as populating matrix

# Loop through all the input files and sort them into the relevant vectors above
for(i in 1:NumPops){
     CurData <- as.character(FileNames$V1[i])
     MetaData <- strsplit(CurData, split = "_")[[1]]
     if(MetaData[2] == "0"){
          if(as.numeric(MetaData[1]) %in% struct){
               StructFound <- c(StructFound, i)
          } else{
               ShuffFound <- c(ShuffFound, i)
          }
     } else{
          if(MetaData[3] == "E"){
               EdgePops <- c(EdgePops, i)
          } else if(MetaData[3] == "C"){
               CorePops <- c(CorePops, i)
          } else if(MetaData[3] == "NA"){
               ShufPops <- c(ShufPops, i)
          }
     }
}

# Create a vector with all the necessary prefixes for the pi calculations.
# NOTE: The 10,000 window value has no prefix as that is the default.
WinVals <- c("_5000", "_7500", "", "_12500", "_15000")
ChromVals <- vector(mode = "list", length = 5)
for(i in 1:5){
     # First, Load in an example data file to extract the number of windows
     #    and use that to set up the PiValsMat for the current window size
     SampleFile <- paste(SamplePrefix, WinVals[i], ".pi", sep = "")
     SampleData <- read.table(SampleFile, na.strings = "na")
     NumWins <- dim(SampleData)[1]
     PiValsMats[[i]] <- matrix(NA, nrow = NumWins, ncol = NumPops)
     
     # Get the chromosome information for each window
     ChromVals[[i]] <- SampleData$V1
     
     # Now loop through each file and read in the pi values
     for(j in 1:NumPops){
          CurData <- as.character(FileNames$V1[j])
          Prefix <- strsplit(CurData, split = ".pi")[[1]]
          InFile <- paste(Prefix, WinVals[i], ".pi", sep = "")
          PiVals <- read.table(InFile, na.strings = "na")
          PiValsMats[[i]][,j] <- PiVals$V5
     }
}

save(PiValsMats, StructFound, ShuffFound, EdgePops, CorePops, ShufPops,
     file = "/home/topher/RangeExpansionGenetics/FinalAnalyses/DeltaPi/PiWinVals.rdata")

# Now calculate the mean and standard deviation of delta pi values at each window
#    for each window size
DeltaData <- vector(mode = "list", length = 5)
for(i in 1:5){
     # First subset for only the data without any missing values
     GoodWins <- which( rowSums(is.na(PiValsMats[[i]])) == 0 )
     WinChroms <- ChromVals[[i]][GoodWins]
     CompletePiVals <- PiValsMats[[i]][GoodWins,]
     
     # Next create matrices to hold all the delta pi values for each category
     #    (core, edge, and shuffled)
     DeltaPiCore <- matrix(NA, nrow = nrow(CompletePiVals), ncol = length(CorePops))
     DeltaPiEdge <- matrix(NA, nrow = nrow(CompletePiVals), ncol = length(EdgePops))
     DeltaPiShuf <- matrix(NA, nrow = nrow(CompletePiVals), ncol = length(ShufPops))
     
     # Now populate those matrices
     for(j in 1:length(CorePops)){
          for(l in 1:nrow(CompletePiVals)){
               DeltaPiCore[l,j] <- CompletePiVals[l, CorePops[j]] - CompletePiVals[l, StructFound[j]]
               DeltaPiEdge[l,j] <- CompletePiVals[l, EdgePops[j]] - CompletePiVals[l, StructFound[j]]
               if(j <= length(ShufPops)){
                    DeltaPiShuf[l,j] <- CompletePiVals[l, ShufPops[j]] - CompletePiVals[l, ShuffFound[j]]
               }
          }
     }
     
     # Finally, create a data frame to hold the mean and standard deviation of
     #    the delta pi values at each window. 
     # NOTE: Since this data represents aggregates over all the individual
     #    landscapes, it is not possible to incorporate a block effect
     location <- c("core", "edge", "shuf")
     CurDeltaData <- expand.grid(window = GoodWins, location = location)
     sigma <- rep(NA, nrow(CurDeltaData))
     mu <- rep(NA, nrow(CurDeltaData))
     Chrom <- rep(WinChroms, 3)
     CurDeltaData <- cbind(CurDeltaData, mu, sigma, Chrom)
     # Create an index for each location to keep track
     CoreIndex <- 1
     EdgeIndex <- 1
     ShufIndex <- 1
     for(j in 1:nrow(CurDeltaData)){
          if(CurDeltaData$location[j] == "core"){
               CurDeltaData$mu[j] <- mean(DeltaPiCore[CoreIndex,])
               CurDeltaData$sigma[j] <- sd(DeltaPiCore[CoreIndex,])
               CoreIndex <- CoreIndex + 1
          } else if(CurDeltaData$location[j] == "edge"){
               CurDeltaData$mu[j] <- mean(DeltaPiEdge[EdgeIndex,])
               CurDeltaData$sigma[j] <- sd(DeltaPiEdge[EdgeIndex,])
               EdgeIndex <- EdgeIndex + 1
          } else{
               CurDeltaData$mu[j] <- mean(DeltaPiShuf[ShufIndex,])
               CurDeltaData$sigma[j] <- sd(DeltaPiShuf[ShufIndex,])
               ShufIndex <- ShufIndex + 1
          }
     }
     DeltaData[[i]] <- CurDeltaData
}

save(DeltaData, file = "/home/topher/RangeExpansionGenetics/FinalAnalyses/DeltaPi/DeltaData.rdata")
