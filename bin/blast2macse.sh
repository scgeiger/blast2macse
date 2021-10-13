#!/bin/bash

# This script is for a gene-based analysis across sequence types, in which a list of ref genes has been compiled
# These genes are located in REF_PATH
# The assumed structure is main/path/ST/gene_names
# I still need to figure out a better way to handle error files

REF_PATH=$1
TIME_STAMP=$(date +"%y%m%d")
ERROR_FILE="$TIME_STAMP-blast2macse_error.log"
FULL_PATH=$(realpath $0)
SCRIPT_PATH=$(dirname $FULL_PATH)


if [ -z $REF_PATH ]; then
    echo "Usage: $0 <directory with reference IDs>"
    echo "Code will assume execution in the parent directory for a dataset"
    echo "This means that it expects each directory here to be of subpopulations (ST, etc)"
    exit
fi

# start this from analysis folder
MAIN_PATH=$(pwd) #this leads to starting directory. This is assumed to be a directory containing subdirs for each sequence type
ERROR_FILE="$MAIN_PATH/$TIME_STAMP-blast2macse_error.log"
echo "#ID   Category    Error" > $ERROR_FILE

# For each subdirectory in this directory
for SUBDIR in `ls -d */`
do
    SUBDIR=${SUBDIR%/} 
    BLAST_DB="$SUBDIR.db"  # Blast dbs were named with this convention when generated using assemblies

    # For each gene in ref gene directory 
    for GENE_NAME in `ls "$REF_PATH"` 
    do
        GENE_NAME=${GENE_NAME%.*}
        ID="$SUBDIR-$GENE_NAME" #this will give a labeling tag unique to each ST and gene

        echo ""
        echo "looking at $GENE_NAME" 
       
        # Set up Directory
        if [ -d "$MAIN_PATH/$SUBDIR/$GENE_NAME" ]; then
            echo "$ID: Dir exists. Moving to dir."
        else
            mkdir "$MAIN_PATH/$SUBDIR/$GENE_NAME"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   CMD Error making directory" >> "$ERROR_FILE"
            fi 
            echo "$ID: Made dir."
        fi

        # Move to the directory
        cd "$MAIN_PATH/$SUBDIR/$GENE_NAME"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   CMD Error moving to new directory" >> "$ERROR_FILE"
            fi 
        echo "$ID: Moved directories."

        #######################
        #     BLAST BLOCK     # 
        #######################
    
        # Run BLAST analysis
        # Inputs: ref gene path, db location, output location, blast specs
        # Output: $ID.raw-blast
        if [ ! -f "$ID.raw-blast" ] && [ ! -f "$ID.blast" ]; then
        blastn -db "/mnt/seqs/db/$BLAST_DB" -query "$REF_PATH/$GENE_NAME.nt" -out "$MAIN_PATH/$SUBDIR/$GENE_NAME/$ID.raw-blast" -outfmt "6 sseqid sseq" -perc_identity 90 -qcov_hsp_perc 90 -num_alignments 100000
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   BLAST   Error running cmd" >> "$ERROR_FILE"
                echo "$ID: Error running blast command"
                echo "$ID: Make sure the db is correct"
                continue
            fi
            if [ ! -s "$ID.raw-blast" ]; then
                echo "$ID   BLAST   No matches found" >> "$ERROR_FILE"
                echo "$ID: No BLAST hits found"
                continue
            fi
        else
            if [ ! -s "$ID.raw-blast" ]; then
                echo "$ID   BLAST   No matches found" >> "$ERROR_FILE"
                echo "$ID: No BLAST hits found"
                continue
            fi
            echo "$ID: BLAST file already exists. Checking format."
        fi

        # Format the blast script to clustal format
        # Make sure blast output exists and has contents. If no matches were found, empty file will be made
        if [ ! -f "$ID.blast" ]; then
            if [ -f "$ID.raw-blast" ] && [ -s "$ID.raw-blast" ]; then
                perl $SCRIPT_PATH/blast2clustal.pl "$ID.raw-blast"
                RESULT=$?
                if [ $RESULT -ne 0 ]; then
                    echo "$ID   BLAST2CLUSTAL unknown" >> "$ERROR_FILE"
                    echo "$ID: Error running blast2clustal.pl" 
                    continue                
                fi
                mv "formatted-$ID.blast" "$ID.blast"
                echo "$ID: formatted BLAST file."
                rm -f "$ID.raw-blast"
                echo "removed $ID.raw-blast"
            elif [ ! -f $ID.raw-blast ]; then
                echo "$ID   BLAST   file not generated" >> "$ERROR_FILE"
                echo "$ID: BLAST file was not generated. Exiting loop."
                continue
            elif [ ! -s $ID.raw-blast ]; then
                echo "$ID   BLAST   output file empty" >> "$ERROR_FILE"
                echo "$ID: BLAST file was empty. No hits presumed. Exiting loop";
                continue
            fi
        else
            echo "$ID: File $ID.blast identified. Skipping format conversion." 
        fi

        ##########################
        #   UNIQSEQS & MACSE     #
        ##########################

        # get unique sequences ready and formatted
        if [ ! -f "$ID.uniqseqs" ]; then
            perl $SCRIPT_PATH/uniqseqs.pl "$ID.blast"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   UNIQSEQ Error running uniqseqs.pl" >> "$ERROR_FILE"
                echo "$ID: Error running uniqseqs.pl"
                continue
            fi
            mv uniqseqs.fasta "$ID.uniqseqs"
        else
            echo "$ID: File $ID.uniqseqs identified. Skipping uniqalleles."
        fi

        # Clean up uniqseqs to make sure all gaps are removed prior to alignment
        if [ ! -f "$ID.uniqseqs.dealign" ]; then
            perl $SCRIPT_PATH/dealign.pl "$ID.uniqseqs"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running dealign.pl"    
                echo "$ID   DEALIGN Error running dealign.pl" >> "$ERROR_FILE" 
                continue
            fi
        else
            echo "$ID: Dealigned file here. Skipping dealignment."
        fi
 
        # Running MACSE
        if [ ! -f "$ID.uniqseqs_NT.aln" ] || [ ! -f "$ID.uniqseqs_AA.aln" ]; then
            java -jar ~/macse_v2.03.jar -prog alignSequences -seq "$ID.uniqseqs.dealign"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   MACSE   Error running alignment" >> "$ERROR_FILE"
                echo "$ID: Error running MACSE"
                continue
            fi
            if [ -f "$ID.uniqseqs_NT.dealign" ] && [ -f "$ID.uniqseqs_AA.dealign" ]; then
                mv "$ID.uniqseqs_NT.dealign" "$ID.uniqseqs_NT.aln"
                mv "$ID.uniqseqs_AA.dealign" "$ID.uniqseqs_AA.aln"        
            else
                echo "$ID   MACSE   Error running alignment" >> "$ERROR_FILE"
                echo "$ID: Error running MACSE"
                continue
            fi 
        else
            echo "$ID: uniqseqs_NT/AA.aln files here. Skipping MACSE alignment."
        fi

        # Expanding unique sequences into whole dataset
        if [ ! -f "expanded-$ID.uniqseqs_NT.aln" ]; then
            perl $SCRIPT_PATH/expand-uniqseqs.pl "$ID.uniqseqs_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running expand-uniqseqs.pl"
                echo "$ID   UNIQSEQ Error expanding alignment" >> "$ERROR_FILE"
                continue
            fi
        else 
            echo "$ID: expanded alignments are present. Skipping expand-uniqseqs.pl"
        fi

        # Making the consensus
        if [ ! -f "$ID.consensus" ]; then
            perl $SCRIPT_PATH/consensus.pl "expanded-$ID.uniqseqs_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running consensus.pl"
                echo "$ID   CONSENSUS   Error running consensus.pl" >> "$ERROR_FILE"
                continue
            fi
            mv consensus.fasta "$ID.consensus"
        else
            echo "$ID: consensus file already found. Skipping consensus generation"
        fi

        ############################
        #  FIRST ALIGN ANALYSIS    #
        ############################

        # Getting codon use and distribution
        if [ ! -f "$ID.macse-codon-dist.tsv" ]; then
            perl $SCRIPT_PATH/codon-distribution.pl -ref "$ID.consensus" -query "expanded-$ID.uniqseqs_NT.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running codon-distribution.pl" 
                echo "$ID   CODONDIST   Error running codon-distribution.pl" >> "$ERROR_FILE"
                continue
            fi
            mv macse-codon-distribution.tsv "$ID.macse-codon-dist.tsv"
        else
            echo "$ID: codon dist file found. Skipping codon-distribution.pl"
        fi

        # Getting uniqallele data for latter plotting
        if [ ! -f "$ID.macse-nt-uniq.tsv" ] || [ ! -f "$ID.macse-aa-uniq.tsv" ]; then
            perl $SCRIPT_PATH/macse-uniq.pl -query "expanded-$ID.uniqseqs_NT.aln" -ref "$ID.consensus" -greedy TRUE
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running macse-uniq.pl" 
                echo "$ID   MACSE   Error running macse-uniq.pl" >> "$ERROR_FILE"
                continue
            fi
            mv nt-uniq-macse-output.tsv "$ID.macse-nt-uniq.tsv"
            mv aa-uniq-macse-output.tsv "$ID.macse-aa-uniq.tsv"
        else
            echo "$ID: macse uniqallele files found. Skipping macse-uniq.pl"
        fi 

        # Running PopGenome analysis
        # Need to format file first
        cp "expanded-$ID.uniqseqs_NT.aln" "pg-fmt-expanded-$ID.uniqseqs_NT.aln"
        sed -i 's/\!/-/g' "pg-fmt-expanded-$ID.uniqseqs_NT.aln"
        sed -i '/^#/d' "pg-fmt-expanded-$ID.uniqseqs_NT.aln"

        if [ ! -f "$ID-PopGenome.tsv" ]; then
            Rscript $SCRIPT_PATH/mini-PopGenome.R "pg-fmt-expanded-$ID.uniqseqs_NT.aln"
            # echo ("Error message: In mean (as.numeric(x)): NAs introduced by coercion is expected. Just ignore")
            # That just means that we're likely seeing only one allele present, or not enough info to calc some stats
            # Regardless, this error will trigger the usual error message, so it's been changed here. 
            mv "PopGenome.tsv" "$ID-PopGenome.tsv"
            rm -f "pg-fmt-expanded-$ID.uniqseqs_NT.aln"

            if [ ! -f "$ID-PopGenome.tsv" ]; then
                echo "$ID: Error running mini-PopGenome.R"
                echo "$ID   POPGENOME   Error running mini-PopGenome.R" >> "$ERROR_FILE"
                continue
            fi
        else
            echo "$ID: PopGenome file was detected. Skipping mini-PopGenome.R"
        fi 
        
        if [ ! -f "$ID.majority-allele.fasta" ]; then
            perl $SCRIPT_PATH/get-major-allele.pl "$ID.uniqseqs.dealign"
            if [ ! -f "$ID.majority-allele.fasta" ]; then
                echo "$ID: Could not get majority allele file"
                echo "$ID   M_ALLELE   Error getting major allele" >> "$ERROR_FILE"
            fi
        else
            echo "$ID: Majority allele file was detected. Skipping get."
        fi

    ######################
    ######################
    ## SEQ REM ANALYSIS ##
    ######################
    ######################

        echo "Now looking at what happens when problem sequences are removed."
        
        # Removing seqs with deletions (1 codon at least) and Ns
        if [ ! -f "$ID-removed-macse.aln" ]; then
            perl $SCRIPT_PATH/rm-delN.pl -aln "expanded-$ID.uniqseqs_NT.aln" -dist "$ID.macse-codon-dist.tsv"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID: Error running rm-delN.pl" >> "$ERROR_FILE"
                continue
            fi
        else
            echo "$ID: removed alignment file found. Skipping rm-delN.pl"
        fi

        # Getting codon use and distribution
        if [ ! -f "$ID-removed.macse-codon-dist.tsv" ]; then
            perl $SCRIPT_PATH/codon-distribution.pl -ref "$ID.consensus" -query "$ID-removed-macse.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID-removed: Error running codon-distribution.pl"
                echo "$ID   R_CODONDIST   Error running codon-distribution.pl" >> "$ERROR_FILE"
                continue
            fi
            mv macse-codon-distribution.tsv "$ID-removed.macse-codon-dist.tsv"
        else
            echo "$ID-removed: codon dist file found. Skipping codon-distribution.pl"
        fi

        # Getting uniqallele data for latter plotting
        if [ ! -f "$ID-removed.macse-nt-uniq.tsv" ] || [ ! -f "$ID-removed.macse-aa-uniq.tsv" ]; then
            perl $SCRIPT_PATH/macse-uniq.pl -query "$ID-removed-macse.aln" -ref "$ID.consensus" -greedy TRUE
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID-removed: Error running macse-uniq.pl"
                echo "$ID   R_MACSE   Error running macse-uniq.pl" >> "$ERROR_FILE"
                continue
            fi
            mv nt-uniq-macse-output.tsv "$ID-removed.macse-nt-uniq.tsv"
            mv aa-uniq-macse-output.tsv "$ID-removed.macse-aa-uniq.tsv"
        else
            echo "$ID-removed: macse uniqallele files found. Skipping macse-uniq.pl"
        fi

        # Running PopGenome analysis
        # Format file first
        cp "$ID-removed-macse.aln" "pg-fmt-$ID-removed-macse.aln"
        sed -i 's/\!/-/g' "pg-fmt-$ID-removed-macse.aln"
        sed -i '/^#/d' "pg-fmt-$ID-removed-macse.aln"
        
        if [ ! -f "$ID-removed-macse-PopGenome.tsv" ]; then
            Rscript $SCRIPT_PATH/mini-PopGenome.R "pg-fmt-$ID-removed-macse.aln"
            # echo ("Error message: In mean (as.numeric(x)): NAs introduced by coercion is expected. Just ignore")
            # That just means that we're likely seeing only one allele present, or not enough info to calc some stats
            # Regardless, this error will trigger the usual error message, so it's been changed here.
            mv "PopGenome.tsv" "$ID-removed-macse-PopGenome.tsv" 
            rm -f "pg-fmt-$ID-removed-macse.aln"

            if [ ! -f "$ID-removed-macse-PopGenome.tsv" ]; then
                echo "$ID-removed: Error running mini-PopGenome.R"
                echo "$ID   R_POPGENOME   Error running mini-PopGenome.R" >> "$ERROR_FILE"
                continue
            fi
        else
            echo "$ID-removed: PopGenome file was detected. Skipping mini-PopGenome.R"
        fi
#### INsert running dnds here
# get uniqseqs of removed- macse aln
        if [ ! -f "$ID-removed-macse.uniqseqs" ]; then
            perl $SCRIPT_PATH/uniqseqs.pl "$ID-removed-macse.aln"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   UNIQSEQ Error running uniqseqs.pl" >> "$ERROR_FILE"
                echo "$ID: Error running uniqseqs.pl"
                continue
            fi
            mv uniqseqs.fasta "$ID-removed-macse.uniqseqs"
        else
            echo "$ID: File $ID-removed-macse.uniqseqs identified. Skipping uniqalleles."
        fi
# run allele mapping
        if [ ! -f "am-$ID-removed-macse.fasta" ]; then
            perl $SCRIPT_PATH/allele-map.pl "$ID-removed-macse.uniqseqs"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   ALLELEMAP Error running allele-map.pl" >> "$ERROR_FILE"
                echo "$ID: Error running allele-map.pl"
                continue
            fi
        else
            echo "$ID: File am-$ID-removed-macse.fasta identified. Skipping allelemapping."
        fi 
# run yn00 wrapper
        if [ ! -f "am-$ID-removed-macse-dnds-matrix.tsv" ]; then
            perl $SCRIPT_PATH/wrap-run-yn.pl "am-$ID-removed-macse.fasta"
            RESULT=$?
            if [ $RESULT -ne 0 ]; then
                echo "$ID   YNWRAP Error running wrap-run-yn.pl" >> "$ERROR_FILE"
                echo "$ID: Error running wrap-run-yn.pl"
                continue
            fi
        else
            echo "$ID: File am-$ID-removed-macse-dnds-matrix.tsv. Skipping dnds calcs."
        fi       
 
        rm -f "$ID-removed-macse.aln"
        rm -f "expanded-$ID.uniqseqs_NT.aln"

        echo ""
        echo "done looping through $GENE_NAME"    
    done # looping through genes
    echo ""
    echo "Done looping through genes for $SUBDIR"

done #looping through subdirs
echo ""
echo "Done looping through all subdirs"
        
