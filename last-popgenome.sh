#!/usr/bin/bash
REF_PATH=$1
TIME_STAMP=$(date +"%C%m%d")

if [ -z $REF_PATH ]; then
    echo "Usage: $0 <directory with reference IDs>"
    echo "Code will assume execution in the parent directory for a datasset"
    echo "This means that it expects each directory here to be of subclades"
    exit
fi

MAIN_PATH=$(pwd) #this leads to starting directory
ERROR_FILE="$MAIN_PATH/$TIME_STAMP-popgenome_error.log"
MAIN_DIR=${PWD##*/}

#for SUBDIR in `ls -d */`
#do
    SUBDIR="10"
    echo "subdir is $SUBDIR"
    for GENE_NAME in `ls "$REF_PATH"`
    do
        GENE_NAME=${GENE_NAME%.*}
        echo "gene $GENE_NAME"
        ID="$SUBDIR-$GENE_NAME"
        echo "id $ID"
        cd "$MAIN_PATH/$SUBDIR/$GENE_NAME"
        cp "expanded-$ID.uniqseqs_NT.aln" "pg-fmt-expanded-$ID.uniqseqs_NT.aln"
        sed -i 's/\!/-/g' "pg-fmt-expanded-$ID.uniqseqs_NT.aln"
        sed -i '/^#/d' "pg-fmt-expanded-$ID.uniqseqs_NT.aln"
        Rscript /mnt/projects/EC_ST131/200923/scripts/mini-PopGenome.R "pg-fmt-expanded-$ID.uniqseqs_NT.aln"
#        echo ("Error message: In mean (as.numeric(x)): NAs introduced by coercion is expected. Just ignore")
    done
#done

