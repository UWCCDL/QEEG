#!/bin/bash
## Convert_to_emotiv.sh
## ---------------------------------------------------------------- ##
## (c) 2016, Andrea Stocco
##     University of Washington, 
##     Seattle, WA 98195
##     Email: stocco@uw.edu
## ---------------------------------------------------------------- ##
## Converts EDF-style data to Emotiv-like datasets for EEG. The 
## Emotiv format used at CCDL is already ready for analysis in R.
## ---------------------------------------------------------------- ## 


if [ -e header.txt ]; then
    # If we have a header file, then we write it into
    # the temp file
    cat header.txt >> file.tmp
fi

# First we count the lines in the file 
n=`grep -v "%" $1 | wc | awk '{print $1}'`
echo $n
# Remove the first 5 seconds (@ 256Hz)
l=$((n - 1280))
echo $l
# Write the last (N-5) seconds of recording on temp file
grep -v "%" ${1} | tail -${l} | tr ',' '\t' >> file.tmp

# Overwrite the original file
mv file.tmp $1
