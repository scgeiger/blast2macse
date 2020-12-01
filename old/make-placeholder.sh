#!/bin/bash

REF_PATH=$1
TIME_STAMP=$(date +"%C%m%d")

if [ -z $REF_PATH ]; then
    echo "Usage: $0 <directory with reference IDs>"
    echo "Code will assume execution in the parent directory for a datasset"
    echo "This means that it expects each directory here to be of subclades"
    exit
fi

MAIN_PATH=$(pwd) #this leads to starting directory

ERROR_FILE="$MAIN_PATH/$TIME_STAMP-make-placeholder_error.log"
for SUBDIR in `ls -d */`
do
    SUBDIR=${SUBDIR%/}
    for GENE_NAME in `ls "$REF_PATH"`
    do
        GENE_NAME=${GENE_NAME%.*}
        cd "$MAIN_PATH/$SUBDIR/$GENE_NAME"
        if [ -f "$SUBDIR-$GENE_NAME-macse-aa-dist.png" ]; then
            echo "$SUBDIR-$GENE_NAME-macse-aa-dist.png found"
        else
            convert -size 11338x3023 xc:\#9197AE "$SUBDIR-$GENE_NAME-macse-aa-dist.png"
            echo "triggered 1"
        fi
        if [ -f "$SUBDIR-$GENE_NAME-macse-nt-dist.png" ]; then
            echo "$SUBDIR-$GENE_NAME-macse-nt-dist.png found"
        else
            convert -size 11338x3023 xc:\#9197AE "$SUBDIR-$GENE_NAME-macse-nt-dist.png"
            echo "triggered 2"
        fi
        if [ -f "$SUBDIR-$GENE_NAME-macse-dec-aa-dist.png" ]; then
            echo "$SUBDIR-$GENE_NAME-macse-dec-aa-dist.png found"
        else
            convert -size 6299x6299 xc:\#9197AE "$SUBDIR-$GENE_NAME-macse-dec-aa-dist.png"
            echo "triggered 3"
        fi
        if [ -f "$SUBDIR-$GENE_NAME-macse-dec-nt-dist.png" ]; then
            echo "$SUBDIR-$GENE_NAME-macse-dec-nt-dist.png found"
        else
            convert -size 6299x6299 xc:\#9197AE "$SUBDIR-$GENE_NAME-macse-dec-nt-dist.png"
            echo "triggered 4"
        fi
    done
done




