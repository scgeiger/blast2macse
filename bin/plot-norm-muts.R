# This script is for plotting normalised mutations (codon) with aa
# position on the x-axis and frequency on the y-axis. Colors are coded
# According to mut type.
# Updated 210907
# Not ready to be automated yet!

pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, patchwork, GGally, ggstance)

args <-commandArgs(TRUE)
if (length(args)==0) {
    stop ("specify a file in dataset-gene-removed-nromalised-pos-muts.txt fmt",
        call.=FALSE)
}
IN.FILE      <- args[1]
ID           <- str_remove(IN.FILE, "-removed-normalised-pos-muts.txt")

d <- read.table(IN.FILE, 
                header = FALSE,  
                sep = "\t", 
                as.is = TRUE, 
                comment.char = "#", 
                quote = "", 
                col.names = c("geneID", "ntPosN", "ntMutID", "aaPosN", "aaMutID", "codonIndex")
)

# Looking at amino acid mutations
aa.mut.cols <- c("Missense"   = "#90BD6D",
                "Nonsense"   = "#2A9D8F",   
                "Synonymous" = "#7C2C8C",
                "Frameshift" = "#F3722C",
                "Insertion"  = "#F9C74F"
               )
aa.mut.types <- unique(d$aaMutID)

# Setting basic plot params for aa mutation plot
norm.mutplot.params = list(
     scale_color_manual ("muttype", values = aa.mut.cols), 
     scale_fill_manual  ("muttype", values = aa.mut.cols), 
     geom_histogram(aes(color = aaMutID, fill = aaMutID), binwidth = .025),
     theme_light(),
     xlim(-.05,1.05),
     labs(x = "Normalised position",
          y = "Frequency of mutation",
          fill = "Amino Acid Mutation Type")
)

plot.title <- paste0(ID, " aa mutation normalised positions")
p          <- ggplot(d, aes(aaPosN, fill = aaMutID), col = geneID) +
              ggtitle(plot.title) + 
              norm.mutplot.params
OUTFILE <- paste(ID, "-aa-norm-pos.png", sep = "")
ggsave(OUTFILE, plot = p, device = "png")


# Looking at nucleotide acid mutations
nt.mut.cols <- c("Deletion"     = "#CF113A",
                 "Frameshift"   = "#F3722C",
                 "Insertion"    = "#F9C74F",
                 "Transition"   = "#90BD6D",
                 "Transversion" = "#2A9D8F",
                 "Ambiguous"    = "#9197AE"
               )
nt.mut.types <- unique(d$ntMutID)

# Setting basic plot params for nt mutation plot
norm.mutplot.params = list(
     scale_color_manual ("muttype", values = nt.mut.cols),
     scale_fill_manual  ("muttype", values = nt.mut.cols),
     geom_histogram(aes(color = ntMutID, fill = ntMutID), binwidth = .025),
     theme_light(),
     xlim(-0.05,1.05),
     labs(x = "Normalised position",
          y = "Frequency of mutation",
          fill = "Nucleotide Mutation Type")
)

plot.title <- paste0(ID, " nt mutation normalised positions")
p          <- ggplot(d, aes(ntPosN, fill = ntMutID), col = geneID) +
              ggtitle(plot.title) +
              norm.mutplot.params

OUTFILE <- paste(ID, "-nt-norm-pos.png", sep = "")
ggsave(OUTFILE, plot = p, device = "png")

print(paste0(ID, ": normalised plots finished"))
