#!/bin/bash

# Use this to collect mutation location data for all genes in directory of choice
# Take in a list of desired gene names
# Designed for use on small gene subsets (pseudo genes, etc)
# Gives out a normalised position (codon) for each gene, and when I write it, also generates a plot.

GENE_LIST=$1
echo $GENE_LIST

TIME_STAMP=$(date +"%y%m%d")
FULL_PATH=$(realpath $0)
SCRIPT_PATH=$(dirname $FULL_PATH)

if [ -z $GENE_LIST ]; then
    echo "Usage: $0 <list of genes.txt to study>"
    echo "Code assumes execution in analysis directory"
    echo "This means that each directory it accesses should be an ST"
    exit
fi

# Starting from analysis folder
MAIN_PATH=$(pwd)

for SUBDIR in `ls -d */`
do
    SUBDIR=${SUBDIR%/}
    for GENE_NAME in `cat $GENE_LIST`
#    input="$GENE_LIST"
#    while IFS= read -r GENE_NAME
    do
        GZ_FILE="$SUBDIR-$GENE_NAME-removed.macse-codon-dist.tsv.gz"
        IN_FILE="$SUBDIR-$GENE_NAME-removed.macse-codon-dist.tsv"
        GZ_SWITCH="F"

        if [ -f $MAIN_PATH/$SUBDIR/$GENE_NAME/$GZ_FILE ]; then
            gunzip $GZ_FILE
            GZ_SWITCH="T"
            echo "Will not work on gzipped file. Fix first."
        fi
        
        perl $SCRIPT_PATH/normalise-mut-pos.pl "$SUBDIR/$GENE_NAME/$IN_FILE"

        if [ $GZ_SWITCH == "T" ]; then 
            gzip $IN_FILE
        fi
    done
    mv "normalised-pos-aamuts.txt" "$SUBDIR-norm-pos-aamuts.txt"
done
    


