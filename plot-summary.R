pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, plotly)

#args <-commandArgs(TRUE)

#Make sure input file is specified
#if (length(args)==0) {
#    stop ("Need to specify summary file for analysis",
#        call.=FALSE)
#}

#SUMMARY.FILE <- args[1]

SUMMARY.FILE <- "all-summary.tsv"
SUMMARY.PREFIX <- gsub("\\..*","",SUMMARY.FILE)
print(SUMMARY.PREFIX)

d <- read.table(SUMMARY.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#", quote = "")

# Setting up key aes vars
gene.shapes <- c("acrR" = 0,
                 "acrA" = 1, 
                 "acrB" = 2, 
                 "gyrA" = 3,
                 "parC" = 4,
                 "marA" = 5,
                 "soxS" = 6,
                 "rob"  = 7,
                 "mprA" = 8,
                 "envR" = 9,
                 "betI" = 10,
                 "yjdC" = 11,
                 "uidR" = 12
                )

gene.order <- c("acrR", "acrA", "acrB", "gyrA", "parC", 
                "marA", "soxS", "rob", 
                "mprA", "envR", 
                "betI", "yjdC", "uidR")

gene.colors <- c("acrR" = "#9E0142",
                 "acrA" = "#D53E4F",
                 "acrB" = "#F46D43",
                 "gyrA" = "#FDAE61",
                 "parC" = "#FEE08B",
                 "marA" = "#FFFFBF",
                 "soxS" = "#E6F598",
                 "rob"  = "#ABDDA4",
                 "mprA" = "#66C2A5",
                 "envR" = "#3288BD",
                 "betI" = "#5E4FA2",
                 "yjdC" = "#9197AE",
                 "uidR" = "#2A2D34"
                )

ds.cols <- c("#9E0142",
                 "#D53E4F",
                 "#F46D43",
                 "#FDAE61",
                 "#FEE08B",
                 "#FFFFBF",
                 "#E6F598",
                 "#ABDDA4",
                 "#66C2A5",
                 "#3288BD",
                 "#5E4FA2",
                 "#9197AE",
                 "#2A2D34"
                )

#################
#   Old vs PG   #
#################
title <- paste("PG vs PG_R: SingHap: ", SUMMARY.PREFIX, sep = "")
CompSingHap <- ggplot(d, aes(x = PG_SingHaps_R, y = PG_SingHaps))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=gene), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144") 
x11()

title <- paste("my vs PG_R: SingHap: ", SUMMARY.PREFIX, sep = "")
CompSingHap <- ggplot(d, aes(x = ntSingHap, y = PG_SingHaps_R))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

title <- paste("PG vs PG_R: FuFs: ", SUMMARY.PREFIX, sep = "")
CompSingHap <- ggplot(d, aes(x = PG_FuFs_R, y = PG_FuFs))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=gene), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

#################
# Summary Stats #
#################
title <- paste("Seqs per ST: ", SUMMARY.PREFIX, sep = "")
seqinds <- ggplot(d, aes(x = reorder(dataset, -db_nseqs), y = db_nseqs))+
            theme_light() +
            geom_point() +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
            ggtitle(title)

title <- paste("Singleton Haps / hitNseqs: ", SUMMARY.PREFIX, sep = "")
singhap.dbnseq <- ggplot(d), aes(x = reorder(dataset, -db_nseqs), y = (ntSingHap / hit_nseqs)))+
            theme_light() +
            geom_point(size = 2, alpha = 0.6) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            ggtitle(title)

title <- paste("Haps / hitNseqs: ", SUMMARY.PREFIX, sep = "")
hap.dbnseq <- ggplot(subset(d, gene %in% c("acrR")), aes(x = reorder(dataset, -db_nseqs), y = (ntHap / hit_nseqs)))+
            theme_light() +
            geom_point(size = 2, alpha = 0.6) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            ggtitle(title)

title <- paste("PG_Haplo / PG_nSeqs: ", SUMMARY.PREFIX, sep = "")
pg.hap.dbnseq <- ggplot(subset(d, gene %in% c("acrR")), aes(x = reorder(dataset, -PG_nSeqs), y = (PG_Haplo / PG_nSeqs)))+
            theme_light() +
            geom_point(size = 2, alpha = 0.6) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            ggtitle(title)

title <- paste("Mine vs PG: SegSite: ", SUMMARY.PREFIX, sep = "")
CompSegSite <- ggplot(subset(d, gene %in% c("acrR")), aes(x = nSegSite, y = PG_SegSite))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

title <- paste("Mine vs PG: Haplo: ", SUMMARY.PREFIX, sep = "")
CompHaplo <- ggplot(subset(d, gene %in% c("acrR")), aes(x = ntHap, y = PG_Haplo))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

title <- paste("Mine vs PG: SingHap: ", SUMMARY.PREFIX, sep = "")
CompSingHap <- ggplot(subset(d, gene %in% c("acrR")), aes(x = ntSingHap, y = PG_SingHaps))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

title <- paste("Mine vs PG: nSeq: ", SUMMARY.PREFIX, sep = "")
CompnSeq <- ggplot(subset(d, gene %in% c("acrR")), aes(x = hit_nseqs, y = PG_nSeqs))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

title <- paste("Mine vs PG: All Gene SegSite: ", SUMMARY.PREFIX, sep = "")
allCompSegSite <- ggplot(d, aes(x = nSegSite, y = PG_SegSite))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144") +
            scale_shape_manual ("Gene", values = gene.shapes)

title <- paste("Mine vs PG: All Gene Haplo: ", SUMMARY.PREFIX, sep = "")
allCompHaplo <- ggplot(d, aes(x = ntHap, y = PG_Haplo))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")+
            scale_shape_manual ("Gene", values = gene.shapes)

title <- paste("Mine vs PG: All Gene SingHap: ", SUMMARY.PREFIX, sep = "")
allCompSingHap <- ggplot(d, aes(x = ntSingHap, y = PG_SingHaps))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")+
            scale_shape_manual ("Gene", values = gene.shapes)

title <- paste("Mine vs PG: All Gene nSeq: ", SUMMARY.PREFIX, sep = "")
allCompnSeq <- ggplot(d, aes(x = hit_nseqs, y = PG_nSeqs))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")+
            scale_shape_manual ("Gene", values = gene.shapes)
##### PG_R vs mine 
title <- paste("Mine vs PG_R: All Gene nSeq: ", SUMMARY.PREFIX, sep = "")
allCompnSeq <- ggplot(d, aes(x = hit_nseqs, y = PG_nSeqs_R))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")+
            scale_shape_manual ("Gene", values = gene.shapes)

title <- paste("Mine vs PG: SegSite: ", SUMMARY.PREFIX, sep = "")
CompSegSite <- ggplot(subset(d, gene %in% c("acrR")), aes(x = nSegSite, y = PG_SegSite))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

title <- paste("Mine vs PG: Haplo: ", SUMMARY.PREFIX, sep = "")
CompHaplo <- ggplot(subset(d, gene %in% c("acrR")), aes(x = ntHap, y = PG_Haplo))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")

title <- paste("Mine vs PG: SingHap: ", SUMMARY.PREFIX, sep = "")
CompSingHap <- ggplot(subset(d, gene %in% c("acrR")), aes(x = ntSingHap, y = PG_SingHaps_R))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) +
            geom_abline(intercept = 0, slope = 1, color = "#F94144")






#######################################################################
title <- paste("PG_FuFs: ", SUMMARY.PREFIX, sep = "")
FuFs <- ggplot(d, aes( x = hit_nseqs, y = PG_FuFs)) +
           theme_light()+
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           geom_text(aes(label=dataset), hjust = 0, vjust = 0)

######################################################################

alleles.over.nseq.ds <- ggplot(subset(d, gene %in% c("acrR")), aes(x = as.character(dataset), y = (ntSingHap / hit_nseqs)))+
            theme_light() +
            geom_point(size = 2, alpha = 0.6) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45))

title <- paste("singHaps / nseq by nseqs: ", SUMMARY.PREFIX, sep = "")
combo <- ggplot(subset(d, gene %in% c("acrR")), aes(x = hit_nseqs, y = (ntSingHap / hit_nseqs)))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title) 
            
all.combo <- ggplot(d, aes(x = hit_nseqs, y = (ntSingHap / hit_nseqs)))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title)+
            scale_shape_manual ("Gene", values = gene.shapes) +
            ylim(0, 0.1)

all.combo <- ggplot(d, aes(x = hit_nseqs, y = PG_FuFs))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title)+
            scale_shape_manual ("Gene", values = gene.shapes) 

nsyn.syn <- ggplot(d, aes(x = hit_nseqs, y = (ntSingHap/ntHap)*((aaMS+aaNS)/(aaMS+aaNS+aaSYN))))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene)) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title)+
            scale_shape_manual ("Gene", values = gene.shapes)

nsyn.syn <- ggplot(d, aes(x = ntHap/hit_nseqs, y = ((aaMS+aaNS)/(aaMS+aaNS+aaSYN)), color = gene))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene )) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title)+
            scale_shape_manual ("Gene", values = gene.shapes) +
            scale_fill_manual ("Gene", values = gene.colors)

title <- paste("Del per Seq: ", SUMMARY.PREFIX, sep = "")
nDel <- ggplot(d, aes(x = reorder(dataset, -hit_nseqs), y = aaDEL/hit_nseqs))+
            theme_light() +
            geom_point(size = 2, alpha = 0.5, aes(shape = gene )) +
            theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
            geom_text(aes(label=dataset), hjust = 0, vjust = 0) +
            ggtitle(title)+
            scale_shape_manual ("Gene", values = gene.shapes) 
###################
# PopGenome Stats #
###################

title <- paste("PG_FuFs: ", SUMMARY.PREFIX, sep = "")
FuFs <- ggplot(d, aes( x = dataset, y = PG_FuFs)) +
           theme_light()+
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(col = dataset)) +
#           scale_shape_manual ("Gene", values = gene.shapes) +
#           scale_color_manual ("Dataset", values = ds.cols) #+
#           scale_x_discrete(limits = gene.order)

title <- paste("PG_RozasR2: ", SUMMARY.PREFIX, sep = "")
RozasR2 <- ggplot(d, aes( x = gene, y = PG_RozasR2)) +
           theme_light()+
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           ylim(0, .2) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = gene.order)

title <- paste("PG_TajimaD_all: ", SUMMARY.PREFIX, sep = "")
TajiD <- ggplot(d, aes( x = dataset, y = PG_TajimaD)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = gene.order)

title <- paste("PG_FuliD: ", SUMMARY.PREFIX, sep = "")
FuliD <- ggplot(d, aes( x = dataset, y = PG_FuliD)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = gene.order)

title <- paste("PG_FuliF: ", SUMMARY.PREFIX, sep = "")
FuliF <- ggplot(d, aes( x = dataset, y = PG_FuliF)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = gene.order)

title <- paste("PG_Pi: ", SUMMARY.PREFIX, sep = "")
pi_PG <- ggplot(d, aes( x = dataset, y = PG_pi)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = gene.order)

pg.figure <- (FuFs | RozasR2 | TajiD | FuliD | FuliF | pi_PG) & theme(legend.position = "right")
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "-popgen.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

#######################
# Singletons vs major #
#######################

title <- "SinVsMaj_all"
SinVsMaj_all <- ggplot(d, aes( x = n_ntMajorHap/hit_nseqs, y = n_ntSingHap/hit_nseqs)) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

#######################
# Mut nature vs FuFS  #
#######################

title <- "HapsperSeq"
nthap_nseq <- ggplot(d, aes( x = gene, y = ntHap/hit_nseqs)) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "aaMSVsFuFs_ST131"
aaMSVsFuFs_ST131 <- ggplot(subset(d, dataset %in% c("tree-A-ST131", "tree-B-ST131", "tree-C-ST131", "tree-all-ST131")), aes( x = PG_FuFs, y = aaMS/(aaMS+aaNS+aaSYN))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "aaNSVsFuFs_all"
aaNSVsFuFs_all <- ggplot(d, aes( x = PG_FuFs, y = aaNS/(aaMS+aaNS+aaSYN))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "aaNSVsFuFs_ST131"
aaNSVsFuFs_ST131 <- ggplot(subset(d, dataset %in% c("tree-A-ST131", "tree-B-ST131", "tree-C-ST131", "tree-all-ST131")), aes( x = PG_FuFs, y = aaNS/(aaMS+aaNS+aaSYN))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "aaSYNVsFuFs_all"
aaSYNVsFuFs_all <- ggplot(d, aes( x = PG_FuFs, y = aaSYN/(aaMS+aaNS+aaSYN))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "aaSYNVsFuFs_ST131"
aaSYNVsFuFs_ST131 <- ggplot(subset(d, dataset %in% c("tree-A-ST131", "tree-B-ST131", "tree-C-ST131", "tree-all-ST131")), aes( x = PG_FuFs, y = aaSYN/(aaMS+aaNS+aaSYN))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "aaNSMSVsFuFs_all"
aaNSMSVsFuFs_all <- ggplot(d, aes( x = PG_FuFs, y = (aaNS+aaMS)/(aaSYN))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "aaNSMSVsFuFs_ST131"
aaNSMSVsFuFs_ST131 <- ggplot(subset(d, dataset %in% c("tree-A-ST131", "tree-B-ST131", "tree-C-ST131", "tree-all-ST131")), aes( x = PG_FuFs, y = (aaNS + aaMS) / (aaSYN))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)


## Ds/Dn
title <- "aaNSMSVsFuFs_ST131"
aaNSMSVsFuFs_ST131 <- ggplot(subset(d, dataset %in% c("tree-A-ST131", "tree-B-ST131", "tree-C-ST131", "tree-all-ST131")), aes( x = PG_FuFs, y = ((aaSYN) / (aaNS + aaMS)))) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

#### Experimenting with coefficients
title <- "h.nseqco_all"
h.nseqco_all <- ggplot(d, aes( x = PG_FuFs, y = (ntHap/hit_nseqs) * ((aaMS+aaNS)/(aaMS+aaNS+aaSYN)) )) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "h.nseqco_ST131"
h.nseqco_ST131 <- ggplot(subset(d, dataset %in% c("tree-A-ST131", "tree-B-ST131", "tree-C-ST131", "tree-all-ST131")), aes( x = PG_FuFs, y = (ntHap/hit_nseqs) * ((aaMS+aaNS)/(aaMS+aaNS+aaSYN)) )) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

#### Experimenting with coefficients
title <- "h.hco_all"
h.hco_all <- ggplot(d, aes( x = PG_FuFs, y = (ntHap-1/ntHap) * ((aaMS+aaNS)/(aaMS+aaNS+aaSYN)) )) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "h.hco_ST131"
h.hco_ST131 <- ggplot(subset(d, dataset %in% c("tree-A-ST131", "tree-B-ST131", "tree-C-ST131", "tree-all-ST131")), aes( x = PG_FuFs, y = (ntHap-1/ntHap) * ((aaMS+aaNS)/(aaMS+aaNS+aaSYN)) )) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "TetR_HapsperSeq"
tetR_nthap_nseq <- ggplot(subset(d, gene %in% c("acrR", "envR", "betI", "yjdC", "uidR")), aes( x = gene, y = ntHap/hit_nseqs)) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "Upreg_HapsperSeq"
upreg_nthap_nseq <- ggplot(subset(d, gene %in% c("marA", "soxS", "rob")), aes( x = gene, y = ntHap/hit_nseqs)) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "Downreg_HapsperSeq"
Downreg_nthap_nseq <- ggplot(subset(d, gene %in% c("acrR", "mprA", "envR")), aes( x = gene, y = ntHap/hit_nseqs)) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "FQR_HapsperSeq"
FQR_nthap_nseq <- ggplot(subset(d, gene %in% c("acrR", "acrA", "acrB", "gyrA", "parC")), aes( x = gene, y = ntHap/hit_nseqs)) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- "Haplotypes over nSequences"
hap_nseq <- ggplot(d, aes( x = gene, y = ntHap/hit_nseqs)) +
           theme_light() +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

OUTFILE <- paste(SUMMARY.PREFIX, "-hapnseq.png", sep = "")
ggsave(hap_nseq, plot = outfig, device = "png", width = 25, height = 12.5, units = "cm", dpi = "retina")
