#!/bin/bash

TIME_STAMP=$(date +"%C%m%d")
MAIN_PATH=$("/mnt/projects/EC_ST131/200923/analysis")
ERROR_FILE="$MAIN_PATH/$TIME_STAMP-make-blast-db-error.log"

# Execute this within /analysis
for SUBDIR in `ls -d */`
do
    SUBDIR=${SUBDIR%/}
    cd "/mnt/projects/EC_ST131/200923/analysis/$SUBDIR"
    ls -lR *fna* | wc -l >> "$SUBDIR-count.txt" 
    gunzip *.gz
    cat *.fna >> all.fasta
    makeblastdb -in all.fasta -parse_seqids -dbtype nucl -out "$SUBDIR.db"
    RESULT=$?
    if [ $RESULT -ne 0 ]; then
        echo "$SUBDIR: error making blastdb" >> "$ERROR_FILE"
        continue
    fi 
    rm -f all.fasta
    mv *.db* "/mnt/projects/EC_ST131/200923/db"
    cd "/mnt/projects/EC_ST131/200923/analysis"

done
    
