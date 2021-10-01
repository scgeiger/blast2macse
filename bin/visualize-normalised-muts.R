# This script is for plotting normalised mutations (codon) with aa
# position on the x-axis and frequency on the y-axis. Colors are coded
# According to mut type.
# Updated 210907
# Not ready to be automated yet!

pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, patchwork, GGally, ggstance)

IN.FILE      <- "normalised-pos-aamuts.txt"
PROT.ID.FILE <- "2ln-GCF_000285655-protID.txt"
SUMMARY.FILE <- "131-all-summary.tsv"

d <- read.table(IN.FILE, 
                header = FALSE,  
                sep = "\t", 
                as.is = TRUE, 
                comment.char = "#", 
                quote = "", 
                col.names = c("geneID", "mutID", "normPos")
)

b <- read.table(PROT.ID.FILE, 
                header = FALSE, 
                sep = "\t",
                strip.white = TRUE,
                as.is = TRUE, 
                comment.char = "#", 
                check.names = TRUE,
                quote = "",
                col.names = c("geneID", "protID")
)

s <- read.table(SUMMARY.FILE,
                header = TRUE,
                sep = "\t",
                strip.white = TRUE,
                as.is = TRUE,
                comment.char = "#",
                check.names = TRUE,
                quote = "",
)

b$geneID <- str_trim(b$geneID)
d        <- left_join(d, b, by = "geneID")
s        <- subset(s, s$gene %in% d$geneID)

mut.cols <- c("Missense"   = "#90BD6D",
              "Nonsense"   = "#2A9D8F",   
              "Synonymous" = "#7C2C8C",
              "Frameshift" = "#F3722C",
              "Insertion"  = "#F9C74F"
             )
mut.types <- unique(d$mutID)

# Setting basic plot params for mutation plot
norm.mutplot.params = list(
     scale_color_manual ("muttype", values = mut.cols), 
     scale_fill_manual ("muttype", values = mut.cols), 
     geom_histogram(aes(color = mutID, fill = mutID), binwidth = .025),
     theme_light(),
     ggtitle(plot.title),
     labs(x = "Normalised position",
          y = "Frequency of mutation",
          fill = "Mutation Type")
)

### Get summary of overall data
protID.count <- count(d, protID)
plot.title <- "Pseudogenes with most mutations"
# Which genes have the most mutations? (all)
    ggplot(protID.count, aes(x = protID, y = n)) +
    geom_bar(stat = "identity") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust=1)) +
    ggtitle(plot.title)
# What happens if we remove the transposases?
d.noT      <- d[!grepl("transposase", d$protID),]
d.noT.ct   <- count(d.noT, protID)
d.onlyT    <- d[grepl("transposase", d$protID),]
d.onlyT.ct <- count(d.onlyT, protID)

noT <- ggplot(d.noT.ct, aes(x = protID, y = n)) +
       geom_bar(stat = "identity") +
       theme_light() +
       theme(axis.text.x = element_text(angle = 90,
                                        vjust = 0.5,
                                        hjust=1)) +
       ggtitle("Pseudogene mutations by gene, no transposase")

onlyT <- ggplot(d.onlyT.ct, aes(x = protID, y = n)) +
         geom_bar(stat = "identity") +
         theme_light() +
         theme(axis.text.x = element_text(angle = 90,
                                          vjust = 0.5,
                                          hjust=1)) +
         ggtitle("Pseudogene mutations by gene, only transposase")

# outliers appear to be gamma-glutamyltransferase EC958_RS27210,
# a hypothetical protein category, and IS66 element that wasn't tagged
# transferase. Of the transferases, IS1 family and IS66 had the most mutations.

# Supplementary plot- what does the gamma-glutamyltransferase look like?
plot.title = "Supplementary: gamma-glutamyltransferase mutations"
s.gg <- ggplot((subset(d, protID == "gamma-glutamyltransferase")), 
                aes(normPos, fill = mutID), col = mut.cols) +
        norm.mutplot.params

# What do the transposase only mutations look like?
plot.title = "Supplementary: transposase only mutations"
s.to <- ggplot(d.onlyT, aes(normPos, fill = mutID), col = mut.cols) +
        norm.mutplot.params

# For our control, what does acrR look like?
plot.title = "Main: ST131-acrR normalised mutations"
acrR.nm <- ggplot((subset(d, geneID == "acrR")), 
           aes(normPos, fill = mutID), col = mut.cols) +
           norm.mutplot.params

# We still need to limit the dataset; bias toward different frequencies of muts
nmut.gene <- count(d, geneID)
# ggplot(nmut.gene, aes(x = geneID, y = n)) + geom_bar(stat="identity")+
#    scale_y_log10() # duplicate of previous plot but in log scale

# I'll limit each of these plots by number of mutations
plot.title <- "Supplementary: only genes > 10,000 mutations"
nmut.10k <- subset(nmut.gene, nmut.gene$n >= 10000)
d.10k    <- subset(d, geneID %in% nmut.10k$geneID)
s.10km   <- ggplot(d.10k, aes(normPos, fill = mutID), col = mut.cols) +
            norm.mutplot.params

plot.title <- "Supplementary: only genes 1000 =< x <10,000 mutations"
nmut.1k <- subset(nmut.gene, n < 10000 & n >= 1000) 
d.1k    <- subset(d, geneID %in% nmut.1k$geneID)
s.1km   <- ggplot(d.1k, aes(normPos, fill = mutID), col = mut.cols) +
           norm.mutplot.params

plot.title <- "Supplementary: only genes < 1,000 mutations"
nmut.k <- subset(nmut.gene, nmut.gene$n < 1000)
d.k    <- subset(d, geneID %in% nmut.k$geneID)
s.km   <- ggplot(d.k, aes(normPos, fill = mutID), col = mut.cols) +
          norm.mutplot.params

# We're seeing more insertions, frameshifts at this level. 
# In particular, there are four regions with frameshifts accumulating.
# Are there any specific genes accounting for these frameshifts?

s.km   <- ggplot(d.k, aes(normPos, fill = mutID), col = geneID) +
          norm.mutplot.params

fs.only <- subset(d.k, mutID == "Frameshift")
fs.only <- count(fs.only, geneID)

fs.only <- ggplot(fs.only, aes(x = geneID, y = n)) + 
        <- geom_bar(stat = "identity") + 
        <- theme(axis.text.x = element_text(angle = 90,
                                            vjust = 0.5,
                                            hjust=1))

high.fs <- subset(fs.only, n > 50)
high.fs.s <- subset(s, gene %in% high.fs$geneID)
high.fs.d <- subset(d, geneID %in% high.fs$geneID)

# list of high.fs genes:
#       geneID   n
#          acrR  95
# EC958_RS06930 413
# EC958_RS11765 448
# EC958_RS11780 429
# EC958_RS16775 370

# make allele plots 
# /mnt/projects/devspace/blast2macse/bin/plot-macse-aa-uniq.R *-aa-uniq.tsv

# make mut distribution plots
# /mnt/projects/devspace/blast2macse/bin/plot-codon-distribution.R *-codon-dist.tsv

# all plots are in pseudo-figs



######################Below current workspace


p <- ggplot(d.rem, aes(V3, fill = V1), col = mut.cols) +
     geom_histogram(aes(color = V1, fill = V1), binwidth = .025)





# All mut types
p <- ggplot(d, aes(V3, fill = V2.x), col = mut.cols) +
     scale_color_manual ("muttype", values = mut.cols) +
     scale_fill_manual ("muttype", values = mut.cols) +
     geom_histogram(aes(color = V2.x, fill = V2.x), binwidth = .025)

NSFS.only <- subset(d, V2.x == "Frameshift" | V2.x == "Nonsense")
p <- ggplot(NSFS.only, aes(V3, fill = V2.x), col = mut.cols) +
     scale_color_manual ("muttype", values = mut.cols) +
     scale_fill_manual ("muttype", values = mut.cols) +
     geom_histogram(aes(color = V2, fill = V2), binwidth = .025)


s <- subset(d, V2 == "Missense")
p <- ggplot(s, aes(V3, fill = V2), col = mut.cols) +
     scale_color_manual ("muttype", values = mut.cols) +
     scale_fill_manual ("muttype", values = mut.cols) +
     geom_histogram(aes(color = V2, fill = V2), binwidth = .025)


s <- subset(d, V2 == "Nonsense")
p <- ggplot(s, aes(V3, fill = V2), col = mut.cols) +
     scale_color_manual ("muttype", values = mut.cols) +
     scale_fill_manual ("muttype", values = mut.cols) +
     geom_histogram(aes(color = V2, fill = V2), binwidth = .025)


s <- subset(d, V2 == "Synonymous")
p <- ggplot(s, aes(V3, fill = V2), col = mut.cols) +
     scale_color_manual ("muttype", values = mut.cols) +
     scale_fill_manual ("muttype", values = mut.cols) +
     geom_histogram(aes(color = V2, fill = V2), binwidth = .025)


s <- subset(d, V2 == "Frameshift")
p <- ggplot(s, aes(V3, fill = V2), col = mut.cols) +
     scale_color_manual ("muttype", values = mut.cols) +
     scale_fill_manual ("muttype", values = mut.cols) +
     geom_histogram(aes(color = V2, fill = V2), binwidth = .025)


s <- subset(d, V2 == "Insertion")
p <- ggplot(s, aes(V3, fill = V2), col = mut.cols) +
     scale_color_manual ("muttype", values = mut.cols) +
     scale_fill_manual ("muttype", values = mut.cols) +
     geom_histogram(aes(color = V2, fill = V2), binwidth = .025)













