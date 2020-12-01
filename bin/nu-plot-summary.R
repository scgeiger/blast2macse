# This script is explicitly to plot summary figures for score distribution across sequence types, and recommend genes for further use
# Output figures will be saved as png files in analysis/ directory, along with file names with convention ST-look-further.tsv
# look-further.tsv will be organized as follows:
# ST    gene    plot to generate
# A script will be made to process these plots 
# These plots will be saved in the /analysis dir under dir /plots, with potential subdirs for the plot type or by gene
# plot names will be ST-gene-plottype.png

pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, patchwork)
#args <-commandArgs(TRUE)

#if (length(args)==0) {
#    stop ("Need to specify summary table for analysis",
#        call.=FALSE)
#}

#SUMMARY.FILE <- args[1] 
SUMMARY.FILE <- "all-ST-summary.tsv"
SUMMARY.PREFIX <- "all-ST-summary"
d <- read.table(SUMMARY.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#", quote = "")

ds.cols <- c("10" = "#90bd6d",
             "11" = "#2a9d8f",
             "131"= "#f94144"
            )

data.order <- c("10", "11", "131")

d$dataset <- as.factor(d$dataset)

gene.of.interest <- "acrR" # build into args
ST.of.interest <- "131"    # build into args
GOI.data <- subset(d, dataset == ST.of.interest & gene == gene.of.interest)

ST131.df <- subset(d, dataset == "131" & nSegSite != "NA")
ST11.df  <- subset(d, dataset == "11" & nSegSite != "NA")
ST10.df  <- subset(d, dataset == "10" & nSegSite != "NA")

s <- subset(d, nSegSite != "NA")
gene.freq <- count(s, gene)
gene.freq <- subset(gene.freq, n == 3)

s <- s[s$gene %in% gene.freq$gene,]

############## Plot summary table
# If no sequence is found, nSegSite = 0
ds.summary <- data.frame(data.order)
ds.summary <- ds.summary %>% rename(ST = data.order)
#ds.summary <- mutate(ds.summary, nGenes = nrow(subset(d, dataset == ds.summary$ST)))

############## Distribution and spread for individual measurements
#Fu's Fs
title <- paste("PG_FuFs-all:")
FuFs <- ggplot(d, aes( x = dataset, y = PG_FuFs, col = dataset)) + 
           geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light()+
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) + 
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

#TajiD
title <- paste("PG_TajimaD-all:")
TajiD <- ggplot(d, aes( x = dataset, y = PG_TajimaD, col = dataset)) +
           geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

#FuliF
title <- paste("PG_FuliF-all:")
FuliF <- ggplot(d, aes( x = dataset, y = PG_FuliF, col = dataset)) +
           geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

#FuliD
title <- paste("PG_FuliD-all:")
FuliD <- ggplot(d, aes( x = dataset, y = PG_FuliD, col = dataset)) +
           geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

#RozasR2
title <- paste("PG_RozasR2-all:")
RozasR <- ggplot(d, aes( x = dataset, y = PG_RozasR2, col = dataset)) +
           geom_boxplot(alpha = .5, notch = TRUE, width = .1) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           scale_y_continuous(trans = 'log10') +
           annotation_logticks(sides="l") + 
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)#trbl
#Pi
title <- paste("PG_Pi-all:")
Pi   <- ggplot(d, aes( x = dataset, y = PG_pi, col = dataset)) +
           theme_light() +
           geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5)) +
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           scale_y_continuous(trans = 'log10') +
           annotation_logticks(sides="l") + #trbl
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)

pg.figure <- ( FuFs | TajiD | FuliF | FuliD | RozasR | Pi ) & theme(legend.position = "right")
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "-popgen.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

############################## Removed

#FuFs
title <- paste("PG_FuFs_R-all")
FuFs <- ggplot(d, aes( x = dataset, y = PG_FuFs_R, col = dataset)) +
        geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#TajiD
title <- paste("PG_TajimaD_R-all:")
TajiD <- ggplot(d, aes( x = dataset, y = PG_TajimaD_R, col = dataset)) +
           geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#FuliF
title <- paste("PG_FuliF_R-all:")
FuliF <- ggplot(d, aes( x = dataset, y = PG_FuliF_R, col = dataset)) +
        geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#FuliD
title <- paste("PG_FuliD_R-all:")
FuliD <- ggplot(d, aes( x = dataset, y = PG_FuliD_R, col = dataset)) +
          geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#RozasR2
title <- paste("PG_RozasR2_R-all:")
RozasR <- ggplot(d, aes( x = dataset, y = PG_RozasR2_R, col = dataset)) +
            geom_boxplot(alpha = .5, notch = TRUE, width = .1) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           scale_y_continuous(trans = 'log10') +
           annotation_logticks(sides="l") +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)#trbl
#Pi
title <- paste("PG_Pi_R-all:")
Pi   <- ggplot(d, aes( x = dataset, y = PG_pi_R, col = dataset)) +
            geom_boxplot(alpha = .5, notch = TRUE, width = .1) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           scale_y_continuous(trans = 'log10') +
           annotation_logticks(sides="l") +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)#trbl

pg.figure <- ( FuFs | TajiD | FuliF | FuliD | RozasR | Pi ) & theme(legend.position = "right")
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "_R-popgen.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

####################### Common genes only

title <- paste("PG_FuFs-common")
FuFs <- ggplot(s, aes( x = dataset, y = PG_FuFs_R, col = dataset)) +
        geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#TajiD
title <- paste("PG_TajimaD-common:")
TajiD <- ggplot(s, aes( x = dataset, y = PG_TajimaD_R, col = dataset)) +
           geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#FuliF
title <- paste("PG_FuliF-common:")
FuliF <- ggplot(s, aes( x = dataset, y = PG_FuliF_R, col = dataset)) +
        geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#FuliD
title <- paste("PG_FuliD-common:")
FuliD <- ggplot(s, aes( x = dataset, y = PG_FuliD_R, col = dataset)) +
          geom_boxplot(alpha = .5, notch = TRUE, width = .1 ) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)
#RozasR2
title <- paste("PG_RozasR2-common")
RozasR <- ggplot(s, aes( x = dataset, y = PG_RozasR2_R, col = dataset)) +
            geom_boxplot(alpha = .5, notch = TRUE, width = .1) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           scale_y_continuous(trans = 'log10') +
           annotation_logticks(sides="l") +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)#trbl
#Pi
title <- paste("PG_Pi-common")
Pi   <- ggplot(s, aes( x = dataset, y = PG_pi_R, col = dataset)) +
            geom_boxplot(alpha = .5, notch = TRUE, width = .1) +
           geom_violin(alpha = .5) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45), plot.title = element_text(face = "bold", hjust = 0.5))+
           ggtitle(title) +
           scale_color_manual ("dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) +
           scale_y_continuous(trans = 'log10') +
           annotation_logticks(sides="l") +
           geom_label(d=d %>% filter(gene == gene.of.interest), aes(label=gene), nudge_x = .25)#trbl

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
