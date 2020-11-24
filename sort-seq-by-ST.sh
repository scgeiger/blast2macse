#!/bin/bash

echo "loaded"
TIME_STAMP=$(date +"%C%m%d")
ERROR_FILE="$TIME_STAMP-sort-seq-by-ST.log"

input="../ST-seq-list.txt"
while read -r line
do
    SEQID="$(cut -f1 <<< $line)"
    ST="$(cut -f2 <<< $line)"
    if [ "$ST" == "NF" ] || [ "$ST" == "NF*" ]; then
        echo "llama"
    else
        if [[ "$ST" =~ .*"*" ]]; then
            ST="${ST//"*"}"
        fi
        if [ -d "/mnt/projects/EC_ST131/200923/analysis/$ST" ]; then
            cp "/mnt/seqs/Ecoli/20023-assemblies/fna/$SEQID.fna.gz" "/mnt/projects/EC_ST131/200923/analysis/$ST"
        else
            mkdir "/mnt/projects/EC_ST131/200923/analysis/$ST"
            cp "/mnt/seqs/Ecoli/20023-assemblies/fna/$SEQID.fna.gz" "/mnt/projects/EC_ST131/200923/analysis/$ST"
        fi
    fi        
done < $input
