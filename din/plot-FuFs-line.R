pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, patchwork, GGally, ggstance)
# take in gcf cds from genomic.tsv, make into dataframe
GCF.FILE <- "GCF_000285655.3_EC958.v1_cds_from_genomic.tsv"
GCF.df <- read.table(GCF.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#", quote = "")

gene.loc <- GCF.df[, c("gene", "location")]

#if $location containts join, separate by comma and make new row with other location

gene.sep <- separate(data = gene.loc, col = location, into = c("start", "stop"), sep = "\\..")

# now have df with gene / start / stop. Need to blend with Fu's Fs score. 

SUMM.FILE <- "/mnt/projects/EC_ST131/210220/summary/210309-all-ST131-summary.tsv"

SUMM.df <- read.table(SUMM.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#", quote = "")
