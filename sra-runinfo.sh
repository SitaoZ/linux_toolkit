#!/usr/bin/env bash
#
# Downloads data from SRA for a given month and year.
#
# Example usage: bash sra-runinfo.sh 2010 5
#
# The script can be automated with GNU parallel like so:
#
# seq 2007 2019 > years.txt
# seq 1 12 > months.txt
# wget http://data.biostarhandbook.com/sra/sra-runinfo.sh
# chmod +x sra-runinfo.sh
# parallel -j 1 --eta --progress --verbose ./sra-runinfo.sh {1} {2} :::: years.txt :::: months.txt
#

# Bail out on errors.
set -ue

# Shortcut to parameters.
YEAR=$1

# The month needs zero padding.
MONTH=$(printf %02d $2)

# This will be the output filename.
FNAME=$YEAR-$MONTH.csv

echo "*** Downloading SRA run info for Year=$YEAR Month=$MONTH."

# Build the command esearch command.
esearch -db sra -query "\"$1/$MONTH\"[Publication Date]" | efetch -format runinfo > $FNAME

# Print the finished file name and size.
SIZE=`cat $FNAME | wc -l | tr -d ' '`
echo "*** Saved results to $FNAME, $SIZE lines."
