# This script will calculate the mean pi value across all autosomes and the sex
#    chromosome for each window size and save the output in a single csv file.

setwd("/mnt/md0raid1_2tb/SpatialStructSeqData/IndPileUps")

# Set up the data frame to hold all the results
FileNames <- read.table("file_names.txt")
NumFiles <- dim(FileNames)[1]
PiMeans <- data.frame(landscape = rep(NA, NumFiles), block = rep(NA, NumFiles),
                      gen = rep(NA, NumFiles), treat = rep(NA, NumFiles), 
                      location = rep(NA, NumFiles),
                      AutoPi_5 = rep(NA, NumFiles), SexPi_5 = rep(NA, NumFiles),
                      AutoSE_5 = rep(NA, NumFiles), SexSE_5 = rep(NA, NumFiles),
                      AutoPi_7.5 = rep(NA, NumFiles), SexPi_7.5 = rep(NA, NumFiles),
                      AutoSE_7.5 = rep(NA, NumFiles), SexSE_7.5 = rep(NA, NumFiles),
                      AutoPi_10 = rep(NA, NumFiles), SexPi_10 = rep(NA, NumFiles),
                      AutoSE_10 = rep(NA, NumFiles), SexSE_10 = rep(NA, NumFiles),
                      AutoPi_12.5 = rep(NA, NumFiles), SexPi_12.5 = rep(NA, NumFiles),
                      AutoSE_12.5 = rep(NA, NumFiles), SexSE_12.5 = rep(NA, NumFiles),
                      AutoPi_15 = rep(NA, NumFiles), SexPi_15 = rep(NA, NumFiles),
                      AutoSE_15 = rep(NA, NumFiles), SexSE_15 = rep(NA, NumFiles))

# Create a vector of all the replicates in each block and treatment to allow
#    for easy processing
B2 <- 21:40
B3 <- 41:60
B4 <- 61:80
Shuff <- c(21:30, 41:50, 61:70)

# Create a vector with all the necessary prefixes for the pi calculations.
# NOTE: The 10,000 window value has no prefix as that is the default.
WinVals <- c("_5000", "_7500", "", "_12500", "_15000")

# Got through each file and fill in the appropriate values in the data frame
for(i in 1:NumFiles){
     # Extract the prefix and relevant metadata from the current population
     CurData <- as.character(FileNames$V1[i])
     Prefix <- strsplit(CurData, split = ".pi")[[1]]
     MetaData <- strsplit(CurData, split = "_")[[1]]
     PiMeans$landscape[i] <- as.numeric(MetaData[1])
     if(PiMeans$landscape[i] %in% B2){
          CurBlock <- 2
     } else if (PiMeans$landscape[i] %in% B3){
          CurBlock <- 3
     } else{
          CurBlock <- 4
     }
     PiMeans$block[i] <- CurBlock
     PiMeans$gen[i] <- as.numeric(MetaData[2])
     PiMeans$treat[i] <- ifelse(PiMeans$landscape[i] %in% Shuff, "shuff", "struct")
     if((PiMeans$gen[i] == 0) | PiMeans$treat[i] == "shuff"){
          PiMeans$location[i] <- NA
     } else{
          PiMeans$location[i] <- MetaData[3]
     }
     # Now loop through the pi values for each window size and populate the data frame
     for(j in 1:5){
          if(j == 3){
               InFile <- CurData
          } else{
               InFile <- paste(Prefix, WinVals[j], ".pi", sep = "")
          }
          PiVals <- read.table(InFile, na.strings = "na")
          Autosomes <- subset(PiVals, V1 != "CM000276.3")
          SexChromosomes <- subset(PiVals, V1 == "CM000276.3")
          PiMeans[i,6+4*(j-1)] <- mean(Autosomes$V5, na.rm = TRUE)
          PiMeans[i,7+4*(j-1)] <- mean(SexChromosomes$V5, na.rm = TRUE)
          PiMeans[i,8+4*(j-1)] <- sd(Autosomes$V5, na.rm = TRUE) / sqrt(sum(!is.na(Autosomes$V5)))
          PiMeans[i,9+4*(j-1)] <- sd(SexChromosomes$V5, na.rm = TRUE) / sqrt(sum(!is.na(SexChromosomes$V5)))
     }
}

write.csv(PiMeans, file = "/home/topher/RangeExpansionGenetics/FinalAnalyses/MeanPi/PiMeans.csv",
          row.names = FALSE)
