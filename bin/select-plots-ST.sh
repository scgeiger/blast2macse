#!/bin/bash

INPUT_FILE=$1
TIME_STAMP=$(date +"%C%m%d")
ERROR_FILE="$TIME_STAMP-visual_error.log"

if [ -z $INPUT_FILE ]; then
    echo "Usage: $0 <File with gene_name\n>"
    echo "Run this code in the ST directory\n";
    echo "As in, you should ls and see geneIDs\n";
    exit
fi

while read line
do
    GENE=$line
    ST=$(basename "$PWD")
    base_path="$PWD"
    ID="${ST}"-"${GENE}"
    script_path=$0
    OUTDIR="/mnt/blast2macse/all-cds/analysis/topFuFs/."
    file_path="${base_path}/${GENE}/${ID}-removed.macse-aa-uniq.tsv"
    echo "$file_path"

    # Visualize uniqalleles
    if [ ! -f "$ID-removed-macse-aa-uniq.png" ]; then
            Rscript /mnt/projects/devspace/blast2macse/bin/plot-macse-aa-uniq.R "$file_path"
            outfile="${base_path}/${GENE}/${ID}-removed-macse-aa-uniq.png"
            mv "$outfile" "$OUTDIR"
            if [ ! -f "$ID-removed-macse-aa-uniq.png" ]; then
            echo "$ID   aa uniq" >> "$ERROR_FILE"
            fi
        else
            echo "$ID: uniq figure here already"
    fi
    
    file_path="${base_path}/${GENE}/${ID}.macse-aa-uniq.tsv"
    # Visualize uniqalleles (not removed)
    if [ ! -f "$ID-macse-aa-uniq.png" ]; then
            Rscript /mnt/projects/devspace/blast2macse/bin/plot-macse-aa-uniq.R "$file_path"
            outfile="${base_path}/${GENE}/${ID}-macse-aa-uniq.png"
            mv "$outfile" "$OUTDIR"
            if [ ! -f "$ID-macse-aa-uniq.png" ]; then
            echo "$ID   aa uniq" >> "$ERROR_FILE"
            fi
        else
            echo "$ID: uniq figure here already"
     fi
    
    # Visualize codon dist(aa)
    file_path="${base_path}/${GENE}/${ID}-removed.macse-codon-dist.tsv"
    if [ ! -f "$ID-removed-macse-dec-aa-dist.png" ]; then
            Rscript /mnt/projects/devspace/blast2macse/bin/plot-codon-distribution.R "$file_path"
            temp="${base_path}/${GENE}/${ID}-removed-macse-aa-dist.png"
            mv "$temp" "$OUTDIR"
            temp="${base_path}/${GENE}/${ID}-removed-macse-nt-dist.png"
            mv "$temp" "$OUTDIR"
            temp="${base_path}/${GENE}/${ID}-removed-macse-dec-aa-dist.png"
            mv "$temp" "$OUTDIR"
            temp="${base_path}/${GENE}/${ID}-removed-macse-dec-nt-dist.png"
            mv "$temp" "$OUTDIR"
            if [ ! -f "$ID-removed-macse-aa-dist.png" ]; then
                echo "$ID   codon dist" >> "$ERROR_FILE"
            fi
        else
            echo "$ID: uniq figure here already"
     fi   
 
done < $INPUT_FILE
