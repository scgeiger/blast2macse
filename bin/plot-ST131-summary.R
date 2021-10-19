# This script is explicitly to plot summary figures for score distribution across sequence types, and recommend genes for further use
# Output figures will be saved as png files in analysis/ directory, along with file names with convention ST-look-further.tsv
# look-further.tsv will be organized as follows:
# ST    gene    plot to generate
# A script will be made to process these plots 
# These plots will be saved in the /analysis dir under dir /plots, with potential subdirs for the plot type or by gene
# plot names will be ST-gene-plottype.png
# Updated 211001

pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, patchwork, GGally, ggstance, ggrepel)

#args <-commandArgs(TRUE)

#if (length(args)==0) {
#    stop ("Need to specify summary table for analysis",
#        call.=FALSE)
#}

#SUMMARY.FILE <- args[1] 
CO <- 0.01 #Cut off for distribution grab
SUMMARY.FILE   <- "all-summary.tsv"
SUMMARY.PREFIX <- "all-summary"
PROT.ID.FILE   <- "EC958-refseq-210825-protID.txt"
PSEUDO.ID.FILE <- ""
GENE.POS.FILE  <- "2ln-EC958-refseq-210825-genepos.txt"
EC958_LEN      <- 5109767 #ref seq NZ_HG941718.1

d <- read.table(SUMMARY.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#", quote = "")

b <- read.table(PROT.ID.FILE,
                header = FALSE,
                sep = "\t",
                strip.white = TRUE,
                as.is = TRUE,
                comment.char = "#",
                check.names = TRUE,
                quote = "",
                col.names = c("gene", "protID")
)

b$gene <- str_trim(b$gene)
b$gene <- gsub(" ", "", b$gene)
d      <- left_join(d, b, by = "gene")

g <- read.table(GENE.POS.FILE,
                header = FALSE,
                sep = "\t",
                strip.white = TRUE,
                as.is = TRUE,
                comment.char = "#",
                check.names = TRUE,
                quote = "",
                col.names = c("gene", "start", "stop", "ori")
)

g$gene <- str_trim(g$gene)
g$gene <- gsub(" ", "", g$gene)
d      <- left_join(d, g, by = "gene")
d$start <- gsub('join\\(', "", d$start)
d$start <- gsub("<", "", d$start)
d$stop <- gsub(">", "", d$stop)

d$start<- as.numeric(d$start)
d$stop <- as.numeric(d$stop)

ds.cols <- c("ST10" = "#90bd6d", # should make this more flexible based on unique values in ds; assign by order?
             "ST11" = "#2a9d8f",
             "ST131"= "#f94144"
            )

OUTFILE          <- paste(SUMMARY.PREFIX, "-figlist.tsv", sep = "") # A list of figures to be made for top hits
QC.OUT           <- paste(SUMMARY.PREFIX, "-qc-rem.tsv", sep = "")  # A list of gene IDs that were removed by QC and why
data.order       <- c("ST10", "ST11", "ST131") # Should make this more flexible
d$dataset        <- as.factor(d$dataset)
gene.of.interest <- "acrR" # build into args
ST.of.interest   <- "131"    # build into args
GOI.data         <- subset(d, dataset == ST.of.interest & gene == gene.of.interest)

ds.list 	<-	unique(d$dataset)
gene.list 	<- 	unique(d$gene)
prop_R.list <-  c("prop_MS_R", "prop_SYN_R", "prop_NSFS_R")
stat_R.list <-  c("PG_TajimaD_R","PG_RozasR2_R","PG_FuliF_R","PG_FuliD_R","PG_FuFs_R","PG_pi_R")
stat.list   <-  c("PG_TajimaD","PG_RozasR2","PG_FuliF","PG_FuliD","PG_FuFs","PG_pi")

# Calculate Proportion columns
d <- d %>% mutate( prop_MS_R   = ((aaMS_R) / (aaNS_R + aaFS_R + aaMS_R + aaSYN_R)))
d <- d %>% mutate( prop_SYN_R  = ((aaSYN_R) / (aaNS_R + aaFS_R + aaMS_R + aaSYN_R)))
d <- d %>% mutate( prop_NSFS_R = ((aaNS_R + aaFS_R) / (aaNS_R + aaFS_R + aaMS_R + aaSYN_R)))
d <- d %>% mutate( rat_hitdb   = (hit_nseqs/db_nseqs))
d <- d %>% mutate( rat_rem     = (nRem/hit_nseqs))

# Finding genes common between all three datasets
# If no sequence is found, nSegSite = 0 (or NA? Verify)
c <- subset(d, nSegSite != "NA")
gene.freq <- count(c, gene)
gene.freq <- subset(gene.freq, n == 3) # make this flexible to fit multiple n datasets
c <- c[c$gene %in% gene.freq$gene,]
#m <- as_tibble(c)
#fufs.only <- select(m, gene, dataset, PG_FuFs_R)
#fufs.only <- fufs.only %>% pivot_wider( names_from = c(dataset), values_from = c(PG_FuFs_R) )


############## Distribution and spread for individual measurements
# PGS for all removed metrics before cleaning
stat.plots <- list()
for(stat_ in stat.list) {
    title <- paste0(stat_, "-all:")
    p     <- ggplot(d, aes_string(x = "dataset", y = stat_, col = "dataset")) +
                 geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
                 geom_violin(alpha = .5) +
                 theme_light() +
                 theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
                 ggtitle(title) +
                 scale_color_manual ("dataset", values = ds.cols) +
                 scale_x_discrete(limits = data.order) +
                 geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

    if (stat_ == "PG_RozasR2" | stat_ == "PG_pi") { # Need to change y-axis to log10 scale
        p <- p + scale_y_continuous(trans = 'log10') +
                 annotation_logticks(sides="l")
    }
    stat.plots[[stat_]] = p
    ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "-distribution.png"))
}

stat.plots <- list()
for(stat_ in stat_R.list) {
    title <- paste0(stat_, "-all:")
    p     <- ggplot(d, aes_string(x = "dataset", y = stat_, col = "dataset")) +
                 geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
                 geom_violin(alpha = .5) +
                 theme_light() +
                 theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
                 ggtitle(title) +
                 scale_color_manual ("dataset", values = ds.cols) +
                 scale_x_discrete(limits = data.order) +
                 geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

    if (stat_ == "PG_RozasR2_R" | stat_ == "PG_pi_R") { # Need to change y-axis to log10 scale
        p <- p + scale_y_continuous(trans = 'log10') +
                 annotation_logticks(sides="l")
    }
    stat.plots[[stat_]] = p
    ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "-distribution.png"))
}

###
# Tables of top hits
# Fu's Fs Table (looping through datasets and statistics)
temp <- d[(d$hit_nseqs != 0),]
temp <- temp[order(temp$PG_FuFs_R),]
only.10 <- subset(temp, dataset == "ST10")
row.cutoff <- ceiling(CO * nrow(only.10))
only.10 <- only.10[1:row.cutoff, ]

only.11 <- subset(temp, dataset == "ST11")
row.cutoff <- ceiling(CO * nrow(only.11))
only.11 <- only.11[1:row.cutoff, ]

only.131 <- subset(temp, dataset == "ST131")
row.cutoff <- ceiling(CO * nrow(only.131))
only.131 <- only.131[1:row.cutoff, ]

# Quality control (overall)
p <- ggplot(d, aes(x = (hit_nseqs/db_nseqs), y = (nRem/hit_nseqs), col = dataset)) +
     geom_point() +
#    geom_label(d=d %>% filter((nRem/hit_nseqs > 0.5) & (hit_nseqs/db_nseqs < 0.5)), mapping = aes(label=gene), nudge_x = .25) +
    scale_color_manual ("dataset", values = ds.cols) +
    theme_light() +
    ggtitle("Overlap of both quality measures") +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(p, file = paste0("qc-11.10.131-distribution.png"))

# How many coding sequences have hits in each ST?
nrow( subset( d, dataset == "ST131" & hit_nseqs > 0) ) #4828
nrow( subset( d, dataset == "ST11"  & hit_nseqs > 0) ) #3961
nrow( subset( d, dataset == "ST10"  & hit_nseqs > 0) ) #4421

# How many genes are found in 90% of the sequences?
nrow( subset( d, dataset == "ST131" & hit_nseqs > 0 & rat_hitdb > 0.9 ) ) #4211
nrow( subset( d, dataset == "ST11"  & hit_nseqs > 0 & rat_hitdb> 0.9 ) ) #3613
nrow( subset( d, dataset == "ST10"  & hit_nseqs > 0 & rat_hitdb > 0.9 ) ) #3547

# What sorts of genes are present at > 1 copy?
nrow( subset(d, dataset == "ST131" & rat_hitdb > 1)) # 146
nrow( subset(d, dataset == "ST11" & rat_hitdb > 1)) # 88 
nrow( subset(d, dataset == "ST10" & rat_hitdb > 1)) # 39

many       <- subset(d, rat_hitdb > 1 )
many       <- subset(many, dataset == "ST131")
many       <- many %>% group_by(protID) %>% summarise(counts = n())

plot.title <- "ST131 Protein IDs with > 1 copy" 
p    <- ggplot(many,
        aes(x = protID, y = counts)) +
        geom_bar(stat = "identity") +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust=1)) +
        ggtitle(plot.title)

ggsave(p, file = paste0("ST131-genes_greater_than_1.png"))

many       <- subset(d, rat_hitdb > 1 )
many       <- subset(many, dataset == "ST11")
many       <- many %>% group_by(protID) %>% summarise(counts = n())

plot.title <- "ST 11 Protein IDs with > 1 copy"
p    <- ggplot(many,
        aes(x = protID, y = counts)) +
        geom_bar(stat = "identity") +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust=1)) +
        ggtitle(plot.title)
ggsave(p, file = paste0("ST11-genes_greater_than_1.png"))

##
many       <- subset(d, rat_hitdb > 1 )
many       <- subset(many, dataset == "ST10")
many       <- many %>% group_by(protID) %>% summarise(counts = n())

plot.title <- "ST10 Protein IDs with > 1 copy"
p    <- ggplot(many,
        aes(x = protID, y = counts)) +
        geom_bar(stat = "identity") +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust=1)) +
        ggtitle(plot.title)

ggsave(p, file = paste0("ST10-genes_greater_than_1.png"))

many       <- subset(d, rat_hitdb > 1)
many       <- many %>% group_by(protID) %>% summarise(counts = n())

plot.title <- "Protein IDs with > 1 copy"
p    <- ggplot(many,
        aes(x = protID, y = counts)) +
        geom_bar(stat = "identity") +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust=1)) +
        ggtitle(plot.title)

ggsave(p, file = paste0("allST-genes_greater_than_1.png"))

# Next, we're interested solely in genes that are well-represented in the pop
# I'll use 0.9 as my cutoff for hits / db
# I'll use 0.1 as my cutoff for prop of sequences removed due to gaps
# Remove sequences with > 1 hit

# nrow d = 14490
f <- subset(d, rat_hitdb >= 0.9) # nrow 11373
f <- subset(f, rat_rem <= .1)    # nrow 11184
f <- subset(f, rat_hitdb <= 1)   # nrow 11010

p <- ggplot(f, aes(x = (hit_nseqs/db_nseqs), y = (nRem/hit_nseqs), col = dataset)) +
     geom_point() +
#    geom_label(d=d %>% filter((nRem/hit_nseqs > 0.5) & (hit_nseqs/db_nseqs < 0.5)), mapping = aes(label=gene), nudge_x = .25) +
    scale_color_manual ("dataset", values = ds.cols) +
    theme_light() +
    ggtitle("Overlap of both quality measures") +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(p, file = paste0("trimmed-qc-11.10.131-distribution.png"))

# I'll use this dataset to move forward in my overall analysis
d <- f

# Making between-dataset stat distribution comparison plots for all calculations
stat.plots <- list()
for(stat_ in stat.list) {
    title <- paste0(stat_, "-all:")
    p     <- ggplot(d, aes_string(x = "dataset", y = stat_, col = "dataset")) +
                 geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
                 geom_violin(alpha = .5) +
                 theme_light() +
                 theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
                 ggtitle(title) +
                 scale_color_manual ("dataset", values = ds.cols) +
                 scale_x_discrete(limits = data.order) +
                 geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

    if (stat_ == "PG_RozasR2" | stat_ == "PG_pi") { # Need to change y-axis to log10 scale
        p <- p + scale_y_continuous(trans = 'log10') +
                 annotation_logticks(sides="l")
    }
    stat.plots[[stat_]] = p
    ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "-cleaned-distribution.png"))
}

# Making between-dataset stat distribution comparison plots for removed calculations
stat.plots <- list()
for(stat_ in stat_R.list) {
	title <- paste0(stat_, "-all:")
	p 	  <- ggplot(d, aes_string(x = "dataset", y = stat_, col = "dataset")) +
				 geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           		 geom_violin(alpha = .5) +
           		 theme_light() +
           		 theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
          		 ggtitle(title) +
            	 scale_color_manual ("dataset", values = ds.cols) +
           		 scale_x_discrete(limits = data.order) +
           		 geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

	if (stat_ == "PG_RozasR2_R" | stat_ == "PG_pi_R") { # Need to change y-axis to log10 scale
		p <- p + scale_y_continuous(trans = 'log10') +
           		 annotation_logticks(sides="l") 
	}
	stat.plots[[stat_]] = p
	ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "cleaned-distribution.png"))
}

# Getting top hits after cleaning
temp <- d[(d$hit_nseqs != 0),]
temp <- temp[order(temp$PG_FuFs_R),]
only.10 <- subset(temp, dataset == "ST10")
row.cutoff <- ceiling(CO * nrow(only.10))
only.10 <- only.10[1:row.cutoff, ]

only.11 <- subset(temp, dataset == "ST11")
row.cutoff <- ceiling(CO * nrow(only.11))
only.11 <- only.11[1:row.cutoff, ]

only.131 <- subset(temp, dataset == "ST131")
row.cutoff <- ceiling(CO * nrow(only.131))
only.131 <- only.131[1:row.cutoff, ]
short.131 <- only.131 %>% select(gene, ref_length, PG_FuFs_R, protID)

write.table(only.131, file = "Top_ST131_FuFs_R.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(short.131, file = "Top_ST131_FuFs_less_R.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Investigating top ST131 hits dot plots
# use only.131, only.10, only.11
p <- ggplot(d, aes(x = ref_length, y = PG_FuFs_R, col = dataset)) +
     geom_point() +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("Gene length vs Fu's Fs")

p <- ggplot(only.131, aes(x = PG_FuFs_R, y = ((aaNS_R + aaFS_R) / (aaNS_R + aaFS_R + aaMS_R + aaSYN_R + aaIN_R)), col = dataset, label = gene)) +
     geom_point() +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("Nonsyn & FS ratio") +
     geom_label_repel(aes(label=gene))
p
ggsave(p, file = paste0("allST-genes_greater_than_1.png"))

# There's more after this. Skipp comments
####################### Common genes only
#stat.plots <- list()
#for(stat_ in stat_R.list) {
#    title <- paste0(stat_, "-common:")
#   p     <- ggplot(c, aes_string(x = "dataset", y = stat_, col = "dataset")) +
#                geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
#                 geom_violin(alpha = .5) +
#                theme_light() +
#                theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
#                ggtitle(title) +
#                scale_color_manual ("dataset", values = ds.cols) +
#                scale_x_discrete(limits = data.order) +
#                geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

#   if (stat_ == "PG_RozasR2_R" | stat_ == "PG_pi_R") { # Need to change y-axis to log10 scale
#       p <- p + scale_y_continuous(trans = 'log10') +
#                annotation_logticks(sides="l")
#   }
#   stat.plots[[stat_]] = p
    #ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "-common-distribution.png"))
#}

#test.plot <- arrangeGrob(grobs = c(stat.plots), nrow = 1)

#pg.figure <- ( FuFs | TajiD | FuliF | FuliD | RozasR | Pi ) & theme(legend.position = "right")
#outfig    <- pg.figure + plot_layout(guides = "collect")
#OUTFILE <- paste(SUMMARY.PREFIX, "-common-popgen.png", sep = "")
#ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

##############Correlation plot
# This section answers the question - did removing sequences with gaps >3 nt in length affect scores?

#FuFs
title <- "PG_FuFs vs PG_FuFs_R"
FuFs  <- ggplot(d, aes(x = PG_FuFs, y = PG_FuFs_R, col = dataset)) +
           geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_abline(intercept = 0, slope = 1, color = "#f89613") +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

#TajiD
title <- "PG_TajimaD vs PG_TajimaD_R"
TajiD <- ggplot(d, aes(x = PG_TajimaD, y = PG_TajimaD_R, col = dataset)) +
           geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_abline(intercept = 0, slope = 1, color = "#f89613")  +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

#FuliF
title <- "PG_FuliF vs PG_FuliF_R"
FuliF <- ggplot(d, aes(x = PG_FuliF, y = PG_FuliF_R, col = dataset)) +
           geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_abline(intercept = 0, slope = 1, color = "#f89613")  +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

#FuliD
title <- "PG_FuliD vs PG_FuliD_R"
FuliD <- ggplot(d, aes(x = PG_FuliD, y = PG_FuliD_R, col = dataset)) +
           geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_abline(intercept = 0, slope = 1, color = "#f89613")  +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

#RozasR
title <- "PG_RozasR2 vs PG_RozasR2_R"
RozasR <- ggplot(d, aes(x = PG_RozasR2, y = PG_RozasR2_R, col = dataset)) +
           geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_abline(intercept = 0, slope = 1, color = "#f89613") +
           scale_y_continuous(trans = 'log10') +
           scale_x_continuous(trans = 'log10')  +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene)) +
           annotation_logticks(sides="lb")

title <- "PG_pi vs PG_pi_R"
Pi <- ggplot(d, aes(x = PG_pi, y = PG_pi_R, col = dataset)) +
           geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+ 
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_abline(intercept = 0, slope = 1, color = "#f89613") +
           scale_y_continuous(trans = 'log10') +
           scale_x_continuous(trans = 'log10')  +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene)) +
           annotation_logticks(sides="lb")

pg.figure <- ( FuFs | TajiD | FuliF | FuliD | RozasR | Pi ) & theme(legend.position = "right")
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "-PGvsR-correlation.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

#######################Mutation properties
# What is the relationship between Fu's Fs and Pi?
#PiVsFuFs
title <- "Pi's relationship to Fu's Fs"
p     <- ggplot(d, aes(y = PG_FuFs_R, x = PG_pi_R, col = dataset)) +
            geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_continuous(trans = 'log10') +
           annotation_logticks(sides="b") +#trbl
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

ggsave(p, file = paste0("allST-pivsfu_R.png"))

########### Investigating the top hits for ST131

# Investigating top ST131 hits dot plots
# use only.131, only.10, only.11
p <- ggplot(only.131, aes(x = prop_NSFS_R, y = PG_FuFs_R, col = dataset, label = gene)) +
     geom_point(alpha = .7) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("Nonsyn & FS aa mutation proportion vs Fu's Fs") +
     geom_label_repel(aes(label=gene))

ggsave(p, file = paste0("ST131-propNSFS_RvsFusFs_R.png"))

p <- ggplot(only.131, aes(x = ref_length, y = PG_FuFs_R, col = dataset, label = gene)) +
     geom_point(alpha = .7) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("ref length vs Fu's Fs") +
     geom_label_repel(aes(label=gene))

ggsave(p, file = paste0("ST131-reflengthvsFusFs_R.png"))

p <- ggplot(only.131, aes(x = PG_pi_R, y = PG_FuFs_R, col = dataset, label = gene)) +
     geom_point(alpha = .7) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("Pi vs Fu's Fs") +
     geom_label_repel(aes(label=gene)) +
     scale_x_continuous(trans = 'log10') +
     annotation_logticks(sides="b") 

ggsave(p, file = paste0("ST131-PGPi_RvsFusFs_R.png"))

p <- ggplot(only.131, aes(x = (aaFS_R/(ref_length/3)), y = PG_FuFs_R, col = dataset, label = gene)) +
     geom_point(alpha = .7) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("aaFS_R normalised to length vs Fu's Fs") +
     geom_label_repel(aes(label=gene))

ggsave(p, file = paste0("ST131-nFS_RRatvsFusFs_R.png"))

p <- ggplot(only.131, aes(x = ((aaFS_R+aaNS_R)/(ref_length/3)), y = PG_FuFs_R, col = dataset, label = gene)) +
     geom_point(alpha = .7) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("aaFS_R + aaNS_R normalised to length vs Fu's Fs") +
     geom_label_repel(aes(label=gene))

ggsave(p, file = paste0("ST131-aaFSNS_RRatvsFusFs_R.png"))

p <- ggplot(only.131, aes(x = (aaSYN_R)/(ref_length/3), y = PG_FuFs_R, col = dataset, label = gene)) +
     geom_point(alpha = .7) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("aaSYN_R normalised to length vs Fu's Fs") +
     geom_label_repel(aes(label=gene))

ggsave(p, file = paste0("ST131-aaSYN_RatvsFusFs_R.png"))

p <- ggplot(only.131, aes(x = ((aaFS_R+aaNS_R+aaMS_R)/(ref_length/3)), y = PG_FuFs_R, col = dataset, label = gene)) +
     geom_point(alpha = .7) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("aaFS_R + aaNS_R normalised to length vs Fu's Fs") +
     geom_label_repel(aes(label=gene))

ggsave(p, file = paste0("ST131-aaFSNSMS_RRatvsFusFs_R.png"))

# Top hits bar chart
many <- rbind(only.131, only.10, only.11)
many <- many %>% group_by(protID) %>% summarise(counts = n())

plot.title <- "All Top 50 Protein IDs"
p    <- ggplot(many,
        aes(x = counts, y = protID)) +
        geom_bar(stat = "identity") +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust=1)) +
        ggtitle(plot.title)

ggsave(p, file = paste0("all-top-FusFs-protIDs.png"))


many <- rbind(only.131, only.10, only.11)
many <- many %>% group_by(gene) %>% summarise(counts = n())

plot.title <- "All Top 50 Protein IDs"
p    <- ggplot(many,
        aes(x = counts, y = gene)) +
        geom_bar(stat = "identity") +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust=1)) +
        ggtitle(plot.title)

ggsave(p, file = paste0("all-top-FusFs-gene.png"))



### plotting Fu's Fs vs location

temp <- f[(f$hit_nseqs != 0),]
temp <- temp[order(temp$PG_FuFs_R),]
only.131 <- subset(temp, dataset == "ST131")

p    <- ggplot(only.131, aes(EC958_LEN, PG_FuFs_R)) +
        geom_point(aes(x = ((stop + start) / 2), y = PG_FuFs_R)) +
        xlim(0, EC958_LEN) +
        ylim(min(only.131$PG_FuFs_R), (max(only.131$PG_FuFs_R) - .1 *(max(only.131$PG_FuFs_R)))) +
#        geom_segment(aes(x = start, y = PG_FuFs_R, xend = stop, yend = PG_FuFs_R, size = 5, show.legend=F))

ggsave(p, file = paste0("test-fusfs-pos.png"))

p    <- ggplot(only.131, aes(EC958_LEN, PG_FuFs_R)) +
        geom_point(aes(x = ((stop + start) / 2), y = PG_FuFs_R)) +
        xlim(0, EC958_LEN) +
        ylim(min(only.131$PG_FuFs_R), (max(only.131$PG_FuFs_R) - .1 *(max(only.131$PG_FuFs_R)))) 
#        geom_segment(aes(x = start, y = PG_FuFs_R, xend = stop, yend = PG_FuFs_R, size = 5, show.legend=F))

ggsave(p, file = paste0("test-fusfs-pos.png"))
