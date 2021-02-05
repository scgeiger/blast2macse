#!/bin/bash

INPUT_FILE=$1
TIME_STAMP=$(date +"%C%m%d")
ERROR_FILE="$TIME_STAMP-visual_error.log"

if [ -z $INPUT_FILE ]; then
    echo "Usage: $0 <File with ST\tgene_name>"
    echo "Run this code in the dir where you want the plots"
    exit
fi

while read line
do
    GENE=$(echo $line | cut -f2 -d' ')
    ST=$(echo $line | cut -f1 -d' ')
    base_path="/mnt/projects/EC_ST131/201001/analysis/$ST/$GENE/"   
    ID="${ST}"-"${GENE}"


    file_path="${base_path}${ID}-removed.macse-aa-uniq.tsv"
    # Visualize uniqalleles
    if [ ! -f "$ID-removed-macse-aa-uniq.png" ]; then
            Rscript /mnt/projects/devspace/blast2macse/bin/plot-macse-aa-uniq.R "$file_path"
            outfile="${base_path}${ID}-removed-macse-aa-uniq.png"
            mv "$outfile" "/mnt/projects/EC_ST131/201001/summary/top_FuFs/"
            if [ ! -f "$ID-removed-macse-aa-uniq.png" ]; then
            echo "$ID   aa uniq" >> "$ERROR_FILE"
            fi
        else
            echo "$ID: uniq figure here already"
     fi
    
    # Visualize codon dist(aa)
    file_path="${base_path}${ID}-removed.macse-codon-dist.tsv"
    if [ ! -f "$ID-removed-macse-dec-aa-dist.png" ]; then
            Rscript /mnt/projects/devspace/blast2macse/bin/plot-codon-distribution.R "$file_path"
            temp="${base_path}${ID}-removed-macse-aa-dist.png"
            mv "$temp" "/mnt/projects/EC_ST131/201001/summary/top_FuFs/"
            temp="${base_path}${ID}-removed-macse-nt-dist.png"
            mv "$temp" "/mnt/projects/EC_ST131/201001/summary/top_FuFs/"
            temp="${base_path}${ID}-removed-macse-dec-aa-dist.png"
            mv "$temp" "/mnt/projects/EC_ST131/201001/summary/top_FuFs/"
            temp="${base_path}${ID}-removed-macse-dec-nt-dist.png"
            mv "$temp" "/mnt/projects/EC_ST131/201001/summary/top_FuFs/"
            if [ ! -f "$ID-removed-macse-aa-dist.png" ]; then
                echo "$ID   codon dist" >> "$ERROR_FILE"
            fi
        else
            echo "$ID: uniq figure here already"
     fi   
 
done < $INPUT_FILE
