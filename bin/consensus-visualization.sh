#!/bin/bash

TIME_STAMP=$(date +"%C%m%d")
ERROR_FILE="$TIME_STAMP-consensusanalysis_error.log"

MAIN_PATH=$(pwd) #this leads to starting directory
ERROR_FILE="$MAIN_PATH/$TIME_STAMP-consensusanalysis_error.log"
MAIN_DIR=${PWD##*/}
# For each subdirectory in this directory
for SUBDIR in `ls -d */`
do
    ID=${SUBDIR%/}
    echo $ID
    cd $ID
    java -jar ~/macse_v2.03.jar -prog alignSequences -seq "all-$ID-consensus.fasta"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running MACSE" >> "$ERROR_FILE"
                continue
            fi
        mv "all-$ID-consensus_NT.fasta" "$ID.consensus_NT.aln"
        mv "all-$ID-consensus_AA.fasta" "$ID.consensus_AA.aln"
           
        perl /mnt/projects/devspace/consensus/consensus.pl "$ID.consensus_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running consensus.pl" >> "$ERROR_FILE"
                continue
            fi
        mv consensus.fasta "$ID.consensus.consensus"

        # MACSE visualization
        perl /mnt/projects/EC_ST131/200923/scripts/codon-distribution.pl -ref "$ID.consensus.consensus" -query "$ID.consensus_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running codon-distribution.pl" >> "$ERROR_FILE"
                continue
            fi
        mv macse-codon-distribution.tsv "$ID.consensus.macse-codon-dist.tsv"
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-codon-distribution.R "$ID.consensus.macse-codon-dist.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running Rscript plot-codon-distribution.R" >> "$ERROR_FILE"
            fi
        perl /mnt/projects/EC_ST131/200923/scripts/nu-macse-uniq.pl -query "$ID.consensus_NT.aln" -ref "$ID.consensus.consensus" -greedy TRUE
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running nu-macse-uniq.pl" >> "$ERROR_FILE"
                continue
            fi
        mv nt-uniq-macse-output.tsv "$ID.consensus.macse-nt-uniq.tsv"
        mv aa-uniq-macse-output.tsv "$ID.consensus.macse-aa-uniq.tsv"
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-macse-nt-uniq.R "$ID.consensus.macse-nt-uniq.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running Rscript plot-macse-nt-uniq.R" >> "$ERROR_FILE"
            fi
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-macse-aa-uniq.R "$ID.consensus.macse-aa-uniq.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running Rscript plot-macse-aa-uniq.R" >> "$ERROR_FILE"
            fi
        rm -f Rplots.pdf
        cd $MAIN_PATH
done



        
        

        
