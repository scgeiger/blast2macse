# This script is explicitly to plot summary figures for score distribution across sequence types, and recommend genes for further use
# Output figures will be saved as png files in analysis/ directory, along with file names with convention ST-look-further.tsv
# look-further.tsv will be organized as follows:
# ST    gene    plot to generate
# A script will be made to process these plots 
# These plots will be saved in the /analysis dir under dir /plots, with potential subdirs for the plot type or by gene
# plot names will be ST-gene-plottype.png
# Updated 210113

pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, patchwork, GGally, ggstance)

#args <-commandArgs(TRUE)

#if (length(args)==0) {
#    stop ("Need to specify summary table for analysis",
#        call.=FALSE)
#}

#SUMMARY.FILE <- args[1] 
CO <- 0.01 #Cut off for distribution grab
SUMMARY.FILE <- "all-ST-summary.tsv"
SUMMARY.PREFIX <- "all-ST-summary"
d <- read.table(SUMMARY.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#", quote = "")

ds.cols <- c("10" = "#90bd6d", # should make this more flexible based on unique values in ds; assign by order?
             "11" = "#2a9d8f",
             "131"= "#f94144"
            )

OUTFILE          <- paste(SUMMARY.PREFIX, "-figlist.tsv", sep = "") # A list of figures to be made for top hits
QC.OUT           <- paste(SUMMARY.PREFIX, "-qc-rem.tsv", sep = "")  # A list of gene IDs that were removed by QC and why
data.order       <- c("10", "11", "131") # Should make this more flexible
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
d <- d %>% mutate( prop_MS_R = ((aaMS_R) / (aaNS_R + aaFS_R + aaMS_R + aaSYN_R)))
d <- d %>% mutate( prop_SYN_R = ((aaSYN_R) / (aaNS_R + aaFS_R + aaMS_R + aaSYN_R)))
d <- d %>% mutate( prop_NSFS_R = ((aaNS_R + aaFS_R) / (aaNS_R + aaFS_R + aaMS_R + aaSYN_R)))

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

# Tables of top hits
# Fu's Fs Table (looping through datasets and statistics)
temp <- d[(d$hit_nseqs != 0),]
temp <- temp[order(temp$PG_FuFs_R),]
only.10 <- subset(temp, dataset == "10")
row.cutoff <- ceiling(CO * nrow(only.10))
only.10 <- only.10[1:row.cutoff, ]

only.11 <- subset(temp, dataset == "11")
row.cutoff <- ceiling(CO * nrow(only.11))
only.11 <- only.11[1:row.cutoff, ]

only.131 <- subset(temp, dataset == "131")
row.cutoff <- ceiling(CO * nrow(only.131))
only.131 <- only.131[1:row.cutoff, ]

# Quality control (overall)
p <- ggplot(d, aes(x = dataset, y = (hit_nseqs/db_nseqs), col = dataset)) +
     geom_boxplot(alpha = .5, width = .1 ) +
     geom_violin(alpha = .5) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("Ratio of hits over total db seqs") +
     scale_x_discrete(limits = data.order) +
     theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))

p <- ggplot(d, aes(x = dataset, y = (nRem/hit_nseqs), col = dataset)) +
     geom_boxplot(alpha = .5, width = .1 ) +
     geom_violin( alpha = .5) +
     scale_color_manual ("dataset", values = ds.cols) +
     theme_light() +
     ggtitle("Ratio of removed over total hits") +
     scale_x_discrete(limits = data.order) +
     theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))     

p <- ggplot(d, aes(x = (hit_nseqs/db_nseqs), y = (nRem/hit_nseqs), col = dataset)) +
     geom_point() +
#    geom_label(d=d %>% filter((nRem/hit_nseqs > 0.5) & (hit_nseqs/db_nseqs < 0.5)), mapping = aes(label=gene), nudge_x = .25) +
    scale_color_manual ("dataset", values = ds.cols) +
    theme_light() +
    ggtitle("Overlap of both quality measures") +
    theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))


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
    ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "-distribution.png"))
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
	ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "-distribution.png"))
}

pg.figure <- ( FuFs | TajiD | FuliF | FuliD | RozasR | Pi ) & theme(legend.position = "right")
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "-popgen.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

############################## Removed

pg.figure <- ( FuFs | TajiD | FuliF | FuliD | RozasR | Pi ) & theme(legend.position = "right")
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "_R-popgen.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

####################### Common genes only
stat.plots <- list()
for(stat_ in stat_R.list) {
    title <- paste0(stat_, "-common:")
    p     <- ggplot(c, aes_string(x = "dataset", y = stat_, col = "dataset")) +
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
    #ggsave(stat.plots[[stat_]], file = paste0("comp-", stat_, "-common-distribution.png"))
}

test.plot <- arrangeGrob(grobs = c(stat.plots), nrow = 1)

pg.figure <- ( FuFs | TajiD | FuliF | FuliD | RozasR | Pi ) & theme(legend.position = "right")
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "-common-popgen.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")


##############Correlation plot

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

#Fs/Length
title <- "Frameshift to total length ratio"
FSVsLength <- ggplot(d, aes(y = PG_FuFs_R, x = ntFS_R/aln_length, col = dataset)) +
            geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene)) +
           scale_x_continuous(trans = 'log10') +
            annotation_logticks(sides="b")

title <- "Nonsense to total length ratio"
NSVsLength <- ggplot(d, aes(y = PG_FuFs_R, x = aaNS_R/(aln_length/3), col = dataset)) +
            geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene)) +
           scale_x_continuous(trans = 'log10') +
            annotation_logticks(sides="b")

title <- "Nonsense and Frameshift to total length ratio"
NSFSVsLength <- ggplot(d, aes(y = PG_FuFs_R, x = (aaNS_R + aaFS_R)/(aln_length/3), col = dataset)) +
            geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene)) +
           scale_x_continuous(trans = 'log10') +
            annotation_logticks(sides="b")

#PiVsFuFs
title <- "Overall Diversity"
PiVsFuFs <- ggplot(d, aes(y = PG_FuFs_R, x = PG_pi_R, col = dataset)) +
            geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_continuous(trans = 'log10') +
           annotation_logticks(sides="b") +#trbl
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

#SingHapsVsFuFs
title <- "Ratio of NS/MS to all muts +  Coeff."
MutIdentityVsFuFs <- ggplot(data=subset(d, !is.na(ntHap_R)), aes(y = PG_FuFs_R, x = ((as.numeric(ntHap_R)/(hit_nseqs-nRem))*((aaMS_R+aaNS_R)/(aaMS_R+aaNS_R+aaSYN_R))), col = dataset)) +
            geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           #scale_x_continuous(trans = 'log10') +
           #annotation_logticks(sides="b") +#trbl
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

title <- "Ratio of NS/MS to all muts +  Coeff. 2"
MutIdentityVsFuFs2 <- ggplot(data=subset(d, !is.na(ntHap_R)), aes(y = PG_FuFs_R, x = ((as.numeric(ntHap_R)-1/as.numeric(ntHap_R))*((aaMS_R+aaNS_R)/(aaMS_R+aaNS_R+aaSYN_R))), col = dataset)) +
            geom_point(alpha = .7) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           #scale_x_continuous(trans = 'log10') +
           #annotation_logticks(sides="b") +#trbl
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene))

####################Similar to ST131
#same pi
same.pi <- subset(d, PG_pi_R < GOI.data$PG_pi_R + .001 & PG_pi_R > GOI.data$PG_pi - .001) 
same.FuFs <- subset(d, PG_FuFs_R < GOI.data$PG_FuFs_R + .1 &  PG_FuFs_R > GOI.data$PG_FuFs_R - .1)
same.SingHaps <- subset(d, (ntSingHap_R / (hit_nseqs-nRem)) == (GOI.data$ntSingHap_R / (GOI.data$hit_nseqs-GOI.data$nRem)))
