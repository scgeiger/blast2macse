#!/usr/bin/bash
XREF_PATH=$1
XTIME_STAMP=$(date +"%C%m%d")

if [ -z $XREF_PATH ]; then
    echo "Usage: $0 <directory with reference IDs>"
    echo "Code will assume execution in the parent directory for a datasset"
    echo "This means that it expects each directory here to be of subclades"
    exit
fi

XMAIN_PATH=$(pwd) #this leads to starting directory
XERROR_FILE="$XMAIN_PATH/$XTIME_STAMP-popgenome_error.log"
XMAIN_DIR=${PWD##*/}

#for XSUBDIR in `ls -d */`
#do
    XSUBDIR="10"
    echo "subdir is $XSUBDIR"
    for XGENE_NAME in `ls "$XREF_PATH"`
    do
        XGENE_NAME=${XGENE_NAME%.*}
        echo "gene $XGENE_NAME"
        XID="$XSUBDIR-$XGENE_NAME"
        echo "id $XID"
        cd "$XMAIN_PATH/$XSUBDIR/$XGENE_NAME"
        echo $PWD
        perl /mnt/projects/EC_ST131/200923/scripts/rm-delN.pl -aln "expanded-$XID.uniqseqs_NT.aln" -dist "$XID.macse-codon-dist.tsv"
        RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$XID: Error running rm-delN.pl" >> "$XERROR_FILE"
                continue
            fi
        perl /mnt/projects/EC_ST131/200923/scripts/codon-distribution.pl -ref "$XID.consensus" -query "$XID-removed-macse.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$XID: Error running codon-distribution.pl" >> "$ERROR_FILE"
                continue
            fi
        mv macse-codon-distribution.tsv "$XID-removed.macse-codon-dist.tsv"
         Rscript /mnt/projects/EC_ST131/200923/scripts/plot-codon-distribution.R "$XID-removed.macse-codon-dist.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$XID: Error running Rscript plot-codon-distribution.R" >> "$ERROR_FILE"
            fi
        perl /mnt/projects/EC_ST131/200923/scripts/nu-macse-uniq.pl -query "$XID-removed-macse.aln" -ref "$XID.consensus" -greedy TRUE
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$XID: Error running nu-macse-uniq.pl" >> "$ERROR_FILE"
                continue
            fi
        mv nt-uniq-macse-output.tsv "$XID-removed.macse-nt-uniq.tsv"
        mv aa-uniq-macse-output.tsv "$XID-removed.macse-aa-uniq.tsv"
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-macse-nt-uniq.R "$XID-removed.macse-nt-uniq.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$XID: Error running Rscript plot-macse-nt-uniq.R" >> "$ERROR_FILE"
            fi
        Rscript /mnt/projects/EC_ST131/200923/scripts/plot-macse-aa-uniq.R "$XID-removed.macse-aa-uniq.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$XID: Error running Rscript plot-macse-aa-uniq.R" >> "$ERROR_FILE"
            fi

        rm -f Rplots.pdf

        cp "$XID-removed-macse.aln" "pg-fmt-$XID-removed-macse.aln"
        sed -i 's/\!/-/g' "pg-fmt-$XID-removed-macse.aln"
        sed -i '/^#/d' "pg-fmt-$XID-removed-macse.aln"
        Rscript /mnt/projects/EC_ST131/200923/scripts/mini-PopGenome.R "pg-fmt-$XID-removed-macse.aln"
#        echo ("Error message: In mean (as.numeric(x)): NAs introduced by coercion is expected. Just ignore")
    done
#done

