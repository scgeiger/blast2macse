pacman::p_load(hfufs, PopGenome, plyr, dplyr, ggplot2)

args <- commandArgs(TRUE)

if (length(args)==0) {
    stop ("Need to specify alignment file for analysis",
        call.=FALSE)
}

INPUT.FILE <- ("pg-fmt-ST131-C-acrR-removed-macse.aln") #args[1]  #expanded alignmenta
INDEX.NO <- gsub("\\..*","",INPUT.FILE)
INDEX.NO <- gsub("pg-fmt-expanded-","",INDEX.NO)
print(INDEX.NO)

# Run your calculations on your gene
d <- hf.readData(INPUT.FILE)
d <- diversity.stats(d, pi=T)
d <- neutrality.stats(d, detail=T, do.R2=T)
n <- data.frame(get.neutrality(d)[[1]])
n$pi <- get.diversity(d)[[1]][,3]
n$x <- xaxis <- sapply(strsplit(sub(" :", "", d@region.names),split=" - "), function(x){return(mean(as.numeric(x)))})
h <- d@region.stats@haplotype.counts
n$n.sequences <- as.numeric(lapply(h, sum))
templist <- lapply(h, ncol)
templist[sapply(templist, is.null)] <- NA
n$n.haplotypes <- as.numeric(templist)
n$n.singleton.haplotypes <- as.numeric(lapply(h, function(x) length(which(x==1))))
for(i in which(is.nan(n$Fu.F_S))) {
  n$Fu.F_S[i] <- hfufs(n$n.sequences[i], n$n.haplotypes[i], n$pi[i])
}

####
#params <- new("test.params")
#params@theta <- rep(n$pi, d@genelength)
# Check, but n$pi is the same as d@theta_Tajima

ms.class <- MS(d, niter = 100, thetaID = "Tajima", neutrality=TRUE)
#####

# Run calculations on your simulation
ms.sim <- readMS("ms.out")
ms.sim <- diversity.stats(ms.sim, pi=T)
ms.sim <- neutrality.stats(ms.sim, detail=T, do.R2=T)
m <- data.frame(get.neutrality(ms.sim)[[1]])
m$pi <- get.diversity(ms.sim)[[1]][,3]
m$x <- xaxis <- sapply(strsplit(sub(" :", "", ms.sim@region.names),split=" - "), function(x){return(mean(as.numeric(x)))})
h <- ms.sim@region.stats@haplotype.counts
m$n.sequences <- as.numeric(lapply(h, sum))
templist <- lapply(h, ncol)
templist[sapply(templist, is.null)] <- NA
m$n.haplotypes <- as.numeric(templist)
m$n.singleton.haplotypes <- as.numeric(lapply(h, function(x) length(which(x==1))))
for(i in which(is.nan(m$Fu.F_S))) {
  m$Fu.F_S[i] <- hfufs(m$n.sequences[i], m$n.haplotypes[i], m$pi[i])
}
dens.FuFs <- density(m$Fu.F_S)
plot(dens.FuFs)
x11()
ggplot(m, aes(x=Fu.F_S)) +
    geom_density()

out.df <- select(n, -c(Fay.Wu.H, Zeng.E, Strobeck.S, x))
OUTFILE <- "PopGenome.tsv"
write.table(out.df, file=OUTFILE, sep="\t", row.names = FALSE, quote = FALSE)
