#!/bin/bash

# 
REF_PATH=$1

if [ -z $REF_PATH ]; then
    echo "Usage: $0 <directory with reference IDs>"
    echo "Code will assume execution in the parent directory for a datasset"
    echo "This means that it expects each directory here to be of subclades"
    exit
fi

# For each subdirectory in this directory
for GENE_NAME in `ls "$REF_PATH"`
do
    GENE_NAME=${GENE_NAME%.*}
    echo "working on $GENE_NAME"
    montage -mode concatenate -tile 5x22 *"$GENE_NAME-macse-nt-dist.png" "$GENE_NAME-fig-macse-nt-dist.png"
    echo "done with macse-nt-dist"
    montage -mode concatenate -tile 5x22 *"$GENE_NAME-macse-aa-dist.png" "$GENE_NAME-fig-macse-aa-dist.png"
    echo "done with macse-aa-dist"
    montage -mode concatenate -tile 10x11 *"$GENE_NAME-macse-dec-nt-dist.png" "$GENE_NAME-fig-dec-nt-dist.png"
    montage -mode concatenate -tile 10x11 *"$GENE_NAME-macse-dec-aa-dist.png" "$GENE_NAME-fig-dec-aa-dist.png"
    montage -mode concatenate -tile 10x11 *"$GENE_NAME-macse-nt-uniq.png" "$GENE_NAME-fig-macse-nt-uniq.png" 
    montage -mode concatenate -tile 10x11 *"$GENE_NAME-macse-aa-uniq.png" "$GENE_NAME-fig-macse-aa-uniq.png"

    mv "$GENE_NAME"* ../final-panels
done
