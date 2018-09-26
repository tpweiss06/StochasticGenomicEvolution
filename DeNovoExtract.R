# This script will extract the values for pi, delta_pi, and allele frequencies for
#    the de novo mutations from the data as Silas extracted it on the server.
# Topher Weiss-Lehman

setwd("/home/topher/RangeExpansionGenetics/FinalAnalyses/DeNovo/")

# Load in the denovo_freqs data and the de_novos object
load("denovo_data.rdata")

# Explore the structure of these data objects
class(denovo_freqs)
str(denovo_freqs)
head(denovo_freqs)
dim(denovo_freqs)
# denovo_freqs is a data frame with a number of rows presumably equal to the
#    number of denovo mutations we found. The number of columns corresponds
#    to the number of core, edge, and shuffled populations (founders are not
#    included)

class(de_novos)
length(de_novos)
str(de_novos)
# de_novos is a list of data frames with a list entry for each experimental
#    group (founders, core, edge, shuffled). Each individual data frame is
#    composed of 66 rows for the 66 de novo mutations found and a number of columns
#    equal to the number of landscapes in each experimental group plus two 
#    (chromosome and position)

# ------------------------------------------------------------------------------
# What I need from this is a single data frame with columns for the mutation index,
#    landscape, location, pi, delta pi, and the frequency of the novel allele.

MutIndex <- NULL
landscape <- NULL
location <- NULL
pi <- NULL
DeltaPi <- NULL
Gen8Freq <- NULL

for(m in 1:nrow(denovo_freqs)){
     cur_MutIndex <- m
     for(p in 1:ncol(denovo_freqs)){
          col_info <- unlist(strsplit(names(denovo_freqs)[p], split = "_"))
          cur_location <- col_info[1]
          cur_landscape <- as.numeric(col_info[2])
          if(denovo_freqs[m,p] < 1){
               MutIndex <- c(MutIndex, cur_MutIndex)
               landscape <- c(landscape, cur_landscape)
               location <- c(location, cur_location)
               Gen8Freq <- c(Gen8Freq, denovo_freqs[m,p])
               if(cur_location == "C"){
                    cur_pi <- de_novos$core[[names(denovo_freqs)[p]]][m]
               } else if(cur_location == "E"){
                    cur_pi <- de_novos$edge[[names(denovo_freqs)[p]]][m]
               } else if(cur_location == "S"){
                    cur_pi <- de_novos$shuffled[[names(denovo_freqs)[p]]][m]
               }
               FounderPi <- de_novos$founder[[paste("F", cur_landscape, sep = "_")]][m]
               pi <- c(pi, cur_pi)
               DeltaPi <- c(DeltaPi, cur_pi - FounderPi)
          }
     }
}

DeNovoResults <- data.frame(MutIndex, landscape, location, pi, DeltaPi, Gen8Freq)

write.csv(DeNovoResults, file = "DeNovoResults.csv", row.names = FALSE)
