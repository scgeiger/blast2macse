pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, patchwork)

args <-commandArgs(TRUE)

#Make sure input file is specified
if (length(args)==0) {
    stop ("Need to specify summary file for analysis",
        call.=FALSE)
}

SUMMARY.FILE <- args[1]

SUMMARY.PREFIX <- gsub("\\..*","",SUMMARY.FILE)
print(SUMMARY.PREFIX)

d <- read.table(SUMMARY.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#")

ds.cols <- c("tree-A-ST131"   = "#CF113A",
             "tree-B-ST131"   = "#F3722C",
             "tree-C-ST131"   = "#F9C74F",
             "tree-all-ST131" = "#90BD6D",
             "tree-non-ST131" = "#2A9D8F"
            )
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

data.order <- c("tree-A-ST131","tree-B-ST131","tree-C-ST131","tree-all-ST131","tree-non-ST131")

title <- paste("PG_FuFs: ", SUMMARY.PREFIX, sep = "")
FuFs <- ggplot(d, aes( x = dataset, y = PG_FuFs)) +
           theme_light()+
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) 

title <- paste("PG_TajimaD_all: ", SUMMARY.PREFIX, sep = "")
TajiD <- ggplot(d, aes( x = dataset, y = PG_TajimaD)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) 

title <- paste("PG_FuliF: ", SUMMARY.PREFIX, sep = "")
FuliF <- ggplot(d, aes( x = dataset, y = PG_FuliF)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order)

title <- paste("PG_FuliD: ", SUMMARY.PREFIX, sep = "")
FuliD <- ggplot(d, aes( x = dataset, y = PG_FuliD)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order)

title <- paste("PG_Pi: ", SUMMARY.PREFIX, sep = "")
pi_PG <- ggplot(d, aes( x = dataset, y = PG_pi)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols) +
           scale_x_discrete(limits = data.order) 

pg.figure <- (FuFs | TajiD | FuliF | FuliD | pi_PG) & theme(legend.position = "right") 
outfig    <- pg.figure + plot_layout(guides = "collect")
OUTFILE <- paste(SUMMARY.PREFIX, "-popgen.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")


########### Allelic breakdown

title <- paste("Alleles over Sequence Ratio: ", SUMMARY.PREFIX, sep = "")
hapseq <- ggplot(d, aes( x = dataset, y = ntHap/hit_nseqs)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45)) +
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)+
           scale_x_discrete(limits = data.order)

title <- paste("nSeqs in Major allele vs nSingletons: ", SUMMARY.PREFIX, sep = "")
MajSin <- ggplot(d, aes( x = n_ntMajorHap/hit_nseqs, y = n_ntSingHap/hit_nseqs)) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- paste("h/nseqs * nonsyn ratio: ", SUMMARY.PREFIX, sep = "")
hapnseq <- ggplot(d, aes( x = PG_FuFs, y = (ntHap/hit_nseqs) * ((aaMS+aaNS)/(aaMS+aaNS+aaSYN)) )) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)

title <- paste("h-1/h * nonsyn ratio: ", SUMMARY.PREFIX, sep = "")
hless1h <- ggplot(d, aes( x = PG_FuFs, y = (ntHap-1/ntHap) * ((aaMS+aaNS)/(aaMS+aaNS+aaSYN)) )) +
           theme_light() +
           theme(text=element_text(size=8), axis.text.x=element_text(angle=45))+
           ggtitle(title) +
           geom_point(aes(shape = gene, col = dataset)) +
           scale_shape_manual ("Gene", values = gene.shapes) +
           scale_color_manual ("Dataset", values = ds.cols)
hap.figure <- (hapseq | MajSin | hapnseq | hless1h) & theme(legend.position = "right")
outfig     <- hap.figure + plot_layout(guides = "collect")
OUTFILE    <- paste(SUMMARY.PREFIX, "-haplo.png", sep = "")
ggsave(OUTFILE, plot = outfig, device = "png", width = 45, height = 12, units = "cm", dpi = "retina")

