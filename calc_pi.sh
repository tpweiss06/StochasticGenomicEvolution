#!/bin/bash
# This script will run the sliding window variance Popoolation script for all 
#	the individual pileup files

# First set the windo and step sizes
WinSizes=(5000 7500 12500 15000)

for win in ${WinSizes[@]}; do
	for file in *.gz; do

		# First unzip the file
		#gunzip $file

		# Get the unzipped file name
		NewFile=`echo $file | awk -F "." '{print $1"."$2}'`
		OutFile=`echo $file | awk -F "." '{print $1}'`"_"$win".pi"

		# Now run the sliding window pi calculation
		#perl /home/topher/RangeExpansionGenetics/popoolation-code/Variance-sliding.pl --input $NewFile --output $OutFile --measure pi --window-size $win --step-size $win --min-count 2 —-max-count 2 --min-coverage 4 --max-coverage 22 --min-qual 20 --pool-size 40 --fastq-type sanger

		# Finally, re-zip the mpileup file to keep the size down
		#gzip $NewFile
        	echo $NewFile
        	echo $OutFile
	done
done
