# source("/mnt/projects/devspace/macse-toolkit/plot-codon-distribution.R")
pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra)
#summary(warnings())
args <-commandArgs(TRUE)

#Make sure input file is specified
if (length(args)==0) {
    stop ("Need to specify macse codon distribution file for analysis",
        call.=FALSE)
}

CODON.FILE <- args[1]
INDEX.NO <- gsub("\\..*","",CODON.FILE)
print(INDEX.NO)

d <- read.table(CODON.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#")
if (dim(d)[1] == 0) {
    stop ("No comparison data. Likely single allele is present.",
        call.=FALSE)
} 
nSeq <- max(d$nseqs)
nCodon <- ((max(d$length)+1) / 3)
seqLength <-(max(d$length)+1)
n.codon <- data.frame("codonpos" = (0:2), "freq" = (0:2)) 
#Number of codons in length:
if ((max(d$length)+1) %% 3 == 0) {
    n.codon[1,"freq"] = (max(d$length)+1) / 3
    n.codon[2,"freq"] = (max(d$length)+1) / 3
    n.codon[3,"freq"] = (max(d$length)+1) / 3 
} else if ((max(d$length)+1) %% 3 == 1) {
    n.codon[1,"freq"] = ((max(d$length)+1) / 3) + 1
    n.codon[2,"freq"] =  (max(d$length)+1) / 3
    n.codon[3,"freq"] =  (max(d$length)+1) / 3
} else if ((max(d$length)+1) %% 3 == 2) {
    n.codon[1,"freq"] = ((max(d$length)+1) / 3) + 1
    n.codon[2,"freq"] = ((max(d$length)+1) / 3) + 1
    n.codon[3,"freq"] =  (max(d$length)+1) / 3
}

### Plot initiation

main.title <- INDEX.NO

### Nucleotide analysis
# frequency of different types of mutations
title <- paste("Distribution of nt mut types: ", main.title, sep = " ")
g <- d %>% tabyl(pos, nt_type)
col.names <- as.list(colnames(g))
nt.muts <- col.names[-1]
data <- reshape2::melt(g, id.vars = "pos", 
        measure.vars = nt.muts ) 
data$prop <- with(data, value / nSeq)

nta <- ggplot(data, aes(x = pos, y = prop, fill = variable)) +
     geom_bar (position = "stack", stat="identity") +     
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Position (nt)") +
     xlim(-1,seqLength) +
#     ylim(0,0.01) +
     ylab("Proportion of sequences with muttype (nt)") +
     scale_fill_manual("Mutation Types", values = c(
                                         "Deletion"     = "#CF113A",
                                         "Frameshift"   = "#F3722C",
                                         "Insertion"    = "#F9C74F",
                                         "Transition"   = "#90BD6D",
                                         "Transversion" = "#2A9D8F",
                                         "Ambiguous"    = "#9197AE"
                                        ))

# frequency of position in codon mutations
# Doesn't really work for indels in this context
title <- paste("Distribution of codon pos muts:", main.title, sep = " ")
g <- d %>% tabyl(pos, codonpos)
col.names <- as.list(colnames(g))
categories <- col.names[-1]

data <- reshape2::melt(g, id.vars = "pos",
        measure.vars = categories )
data$prop <- with(data, value / nSeq)

ntb <- ggplot(data, aes(x = pos, y = prop, fill = variable)) +
     geom_bar (position = "stack", stat="identity") +
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Position (nt)") +
     xlim(-1,seqLength+1) +
     #ylim(0,0.01) +
     ylab("Proportion of mutation in codon position") +
     scale_fill_manual("Codon Positions", values = c(
                                         "0"  = "#CF113A",
                                         "1"  = "#F9C74F",
                                         "2"  = "#2A9D8F"
                                        ))


# position in codon (0, 1, 2)
    # frequency of mutation type
title <- paste("proportion of nt mut type at each codon pos:", main.title, sep = " ")
g <- d %>% tabyl(codonpos, nt_type)
col.names <- as.list(colnames(g))
categories <- col.names[-1]

data <- reshape2::melt(g, id.vars = "codonpos",
        measure.vars = categories )
data2 <- merge(data, n.codon, by="codonpos")

data2$prop <- with(data2, (value / (nSeq * freq)))

ntc <- ggplot(data2, aes(x = codonpos, y = prop, fill = variable)) +
     geom_bar (position = "stack", stat="identity") +
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Codon Position") +
     ylab("Proportion of Mutation Type (nt)") +
     scale_fill_manual("Mutation Types", values = c(
                                         "Deletion"     = "#CF113A",
                                         "Frameshift"   = "#F3722C",
                                         "Insertion"    = "#F9C74F",
                                         "Transition"   = "#90BD6D",
                                         "Transversion" = "#2A9D8F",
                                         "Ambiguous"    = "#9197AE"
                                        ))
    
# position in amino acid translation
    # frequency of different types of mutations
title <- paste("Distribution of aa mut types: ", main.title, sep = " ")

g <- d %>% tabyl(codon_index, aa_type)
col.names <- as.list(colnames(g))
categories <- col.names[-1]

data <- reshape2::melt(g, id.vars = "codon_index",
        measure.vars = categories ) 
data$prop <- with(data, value / nSeq)

aaa <- ggplot(data, aes(x = codon_index, y = prop, fill = variable)) +
     geom_bar (position = "stack", stat="identity") +
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Position (aa)") +
     xlim(-1,nCodon) +
     #ylim(0,0.01) +
     ylab("Proportion of sequences with muttype (aa)") +
     scale_fill_manual("Mutation Types", values = c(
                                         "Deletion"     = "#CF113A",
                                         "Frameshift"   = "#F3722C",
                                         "Insertion"    = "#F9C74F",
                                         "Missense"     = "#90BD6D",
                                         "Nonsense"     = "#2A9D8F", 
                                         "Synonymous"   = "#7C2C8C",
                                         "Ambiguous"    = "#9197AE"
                                        ))

# frequency of position in codon mutations
g <- d %>% tabyl(codon_index, codonpos)
col.names <- as.list(colnames(g))
categories <- col.names[-1]

title <- paste("Distribution of codon pos muts:", main.title, sep = " ")

data <- reshape2::melt(g, id.vars = "codon_index",
        measure.vars = categories)
data$prop <- with(data, value / nSeq)

aab <- ggplot(data, aes(x = codon_index, y = prop, fill = variable)) +
     geom_bar (position = "stack", stat="identity") +
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Position (aa)") +
     xlim(-1,nCodon) +
     #ylim(0,0.01) +
     ylab("Proportion of mutation in codon position") +
     scale_fill_manual("Codon Positions", values = c(
                                         "0"  = "#CF113A",
                                         "1"  = "#F9C74F",
                                         "2"  = "#2A9D8F"
                                        ))


# codon index
    # Frequency of types of mutation    
g <- d %>% tabyl(codonpos, aa_type)
col.names <- as.list(colnames(g))
categories <- col.names[-1]

title <- paste("proportion of aa mut type at each codon pos:", main.title, " n=", nSeq, sep = " ")

data <- reshape2::melt(g, id.vars = "codonpos",
        measure.vars = categories )
data2 <- merge(data, n.codon, by="codonpos")

data2$prop <- with(data2, (value / (nSeq * freq)))
aac <- ggplot(data2, aes(x = codonpos, y = prop, fill = variable)) +
     geom_bar (position = "stack", stat="identity") +
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Codon Position") +
     ylab("Proportion of Mutation Type (aa)") +
     scale_fill_manual("Mutation Types", values = c(
                                         "Deletion"     = "#CF113A",
                                         "Frameshift"   = "#F3722C",
                                         "Insertion"    = "#F9C74F",
                                         "Missense"     = "#90BD6D",
                                         "Nonsense"     = "#2A9D8F",
                                         "Synonymous"   = "#7C2C8C",
                                         "Ambiguous"    = "#9197AE"
                                        )) 

# Plotting nt mutations on own axes
title <- paste("Deconstructed distribution of nt mut types: ", main.title, sep = " ")
g <- d %>% tabyl(pos, nt_type)
col.names <- as.list(colnames(g))
nt.muts <- col.names[-1]
data <- reshape2::melt(g, id.vars = "pos",
        measure.vars = nt.muts )
data$prop <- with(data, value / nSeq)

ntd <- ggplot(data, aes(x = pos, y = prop, fill = variable)) +
     geom_bar (stat="identity") +
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Position (nt)") +
     xlim(-1,seqLength+1) +
#     ylim(0,0.01) +
     ylab("Proportion of sequences with muttype (nt)") +
     scale_fill_manual("Mutation Types", values = c(
                                         "Deletion"     = "#CF113A",
                                         "Frameshift"   = "#F3722C",
                                         "Insertion"    = "#F9C74F",
                                         "Transition"   = "#90BD6D",
                                         "Transversion" = "#2A9D8F",
                                         "Ambiguous"    = "#9197AE"
                                        )) +
    facet_wrap(~ variable, ncol = 1, strip.position = "right", scales = "free")

# Plotting aa mutations on own axes
title <- paste("Deconstructed distribution of aa mut types: ", main.title, sep = " ")
g <- d %>% tabyl(codon_index, aa_type)
col.names <- as.list(colnames(g))
categories <- col.names[-1]

data <- reshape2::melt(g, id.vars = "codon_index",
        measure.vars = categories )
data$prop <- with(data, value / nSeq)

aad <- ggplot(data, aes(x = codon_index, y = prop, fill = variable)) +
     geom_bar (stat="identity") +
     theme_classic() +
     ggtitle(title) +
     theme(plot.title = element_text(size = 7, face = "bold")) +
     xlab("Position (aa)") +
     xlim(-1,nCodon) +
     #ylim(0,0.01) +
     ylab("Proportion of sequences with muttype (aa)") +
     scale_fill_manual("Mutation Types", values = c(
                                         "Deletion"     = "#CF113A",
                                         "Frameshift"   = "#F3722C",
                                         "Insertion"    = "#F9C74F",
                                         "Missense"     = "#90BD6D",
                                         "Nonsense"     = "#2A9D8F",
                                         "Synonymous"   = "#7C2C8C",
                                         "Ambiguous"    = "#9197AE"
                                        )) +
    facet_wrap(~ variable, ncol = 1, strip.position = "right", scales = "free")


nt.figure <- grid.arrange(nta, ntb, ntc, nrow = 1)
OUTFILE <- paste(INDEX.NO, "-macse-nt-dist.png", sep = "") 
ggsave(OUTFILE, plot = nt.figure, device = "png", width = 90, height = 24, units = "cm", dpi = "retina")

aa.figure <- grid.arrange(aaa, aab, aac, nrow = 1 )
OUTFILE <- paste(INDEX.NO, "-macse-aa-dist.png", sep = "")
ggsave(OUTFILE, plot = aa.figure, device = "png", width = 90, height = 24, units = "cm", dpi = "retina")

OUTFILE <- paste(INDEX.NO, "-macse-dec-nt-dist.png", sep = "")
ggsave(OUTFILE, plot = ntd, device = "png", width = 50, height = 50, units = "cm", dpi = "retina")

OUTFILE <- paste(INDEX.NO, "-macse-dec-aa-dist.png", sep = "")
ggsave(OUTFILE, plot = aad, device = "png", width = 50, height = 50, units = "cm", dpi = "retina")

print("All plots saved to *.png")
