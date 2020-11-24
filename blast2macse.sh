#!/bin/bash

REF_PATH=$1
TIME_STAMP=$(date +"%C%m%d")
ERROR_FILE="$TIME_STAMP-blast2macse_error.log"

if [ -z $REF_PATH ]; then
    echo "Usage: $0 <directory with reference IDs>"
    echo "Code will assume execution in the parent directory for a datasset"
    echo "This means that it expects each directory here to be of subclades"
    exit
fi

# start this from analysis folder
MAIN_PATH=$(pwd) #this leads to starting directory

ERROR_FILE="$MAIN_PATH/$TIME_STAMP-blast2macse_error.log"
# For each subdirectory in this directory
for SUBDIR in `ls -d */`
do
    SUBDIR=${SUBDIR%/}
    BLAST_DB="$SUBDIR.db"

    # For each gene in ref gene directory 
    for GENE_NAME in `ls "$REF_PATH"` 
    do
        GENE_NAME=${GENE_NAME%.*} 
        mkdir "$MAIN_PATH/$SUBDIR/$GENE_NAME"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error making directory" >> "$ERROR_FILE"
                continue
            fi
        cd "$MAIN_PATH/$SUBDIR/$GENE_NAME"
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error moving to new directory" >> "$ERROR_FILE"
                continue
            fi 
        ID="$SUBDIR-$GENE_NAME"
        echo $ID

        # Run BLAST analysis
        blastn -db "/mnt/projects/EC_ST131/200923/db/$BLAST_DB" -query "$REF_PATH/$GENE_NAME.nt" -out "$MAIN_PATH/$SUBDIR/$GENE_NAME/$ID.raw-blast" -outfmt "6 sseqid sseq" -perc_identity 90 -qcov_hsp_perc 90 -num_alignments 100000
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running BLASTcmd" >> "$ERROR_FILE"
                continue
            fi

        perl /mnt/projects/devspace/blast2clustal/blast2clustal.pl "$ID.raw-blast"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running blast2clustal.pl" >> "$ERROR_FILE"
                continue                
            fi

        mv "formatted-$ID.blast" "$ID.blast"
   
        # Get unique sequences ready and formatted
        perl /mnt/projects/devspace/uniqallele/seq_uniq/uniqseqs.pl "$ID.blast"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running uniqseqs.pl" >> "$ERROR_FILE"
                continue
            fi

        mv uniqseqs.fasta "$ID.uniqseqs"
        perl /mnt/projects/EC_ST131/200923/scripts/dealign.pl "$ID.uniqseqs"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running dealign.pl" >> "$ERROR_FILE"
                continue
            fi

        # Start allelic analysis
        perl /mnt/projects/EC_ST131/200923/scripts/get-major-allele.pl "$ID.uniqseqs.dealign"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running get-major-allele.pl" >> "$ERROR_FILE"
                continue
            fi
        java -jar ~/macse_v2.03.jar -prog alignSequences -seq "$ID.uniqseqs.dealign"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running MACSE" >> "$ERROR_FILE"
                continue
            fi
        mv "$ID.uniqseqs_NT.dealign" "$ID.uniqseqs_NT.aln"
        mv "$ID.uniqseqs_AA.dealign" "$ID.uniqseqs_AA.aln"
        perl /mnt/projects/devspace/expand-uniqseqs/expand-uniqseqs.pl "$ID.uniqseqs_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running expand-uniqseqs.pl" >> "$ERROR_FILE"
                continue
            fi
        perl /mnt/projects/devspace/consensus/consensus.pl "expanded-$ID.uniqseqs_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running consensus.pl" >> "$ERROR_FILE"
                continue
            fi
        mv consensus.fasta "$ID.consensus"

        # MACSE visualization
        perl /mnt/projects/EC_ST131/200923/scripts/codon-distribution.pl -ref "$ID.consensus" -query "expanded-$ID.uniqseqs_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running codon-distribution.pl" >> "$ERROR_FILE"
                continue
            fi
        mv macse-codon-distribution.tsv "$ID.macse-codon-dist.tsv"
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-codon-distribution.R "$ID.macse-codon-dist.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running Rscript plot-codon-distribution.R" >> "$ERROR_FILE"
            fi
        perl /mnt/projects/EC_ST131/200923/scripts/nu-macse-uniq.pl -query "expanded-$ID.uniqseqs_NT.aln" -ref "$ID.consensus" -greedy TRUE
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running nu-macse-uniq.pl" >> "$ERROR_FILE"
                continue
            fi
        mv nt-uniq-macse-output.tsv "$ID.macse-nt-uniq.tsv"
        mv aa-uniq-macse-output.tsv "$ID.macse-aa-uniq.tsv"
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-macse-nt-uniq.R "$ID.macse-nt-uniq.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running Rscript plot-macse-nt-uniq.R" >> "$ERROR_FILE"
            fi
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-macse-aa-uniq.R "$ID.macse-aa-uniq.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running Rscript plot-macse-aa-uniq.R" >> "$ERROR_FILE"
            fi

        rm -f Rplots.pdf
    done # looping through genes
done #looping through subdirs
        

        # SNAP Calculations
#        perl /mnt/projects/devspace/snapconvert/snapconvert.pl "expanded-$ID.uniqseqs_NT.aln"
#            RESULT=$?
#            if [ $RESULT -ne 0 ]; then
#                echo "$ID: Error running snapconvert.pl" >> "$ERROR_FILE"
#                continue
#            fi
#        perl /mnt/projects/EC_ST131/200923/scripts/SNAP.pl "snap-fmt-expanded-$ID.fasta"
#            RESULT=$?
#            if [ $RESULT -ne 0 ]; then
#                echo "$ID: Error running SNAP.pl" >> "$ERROR_FILE"
#                continue
#            fi



        
        

        
