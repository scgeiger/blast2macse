pacman::p_load(hfufs, PopGenome, plyr, dplyr)

args <- commandArgs(TRUE)

if (length(args)==0) {
    stop ("Need to specify alignment file for analysis",
        call.=FALSE)
}

INPUT.FILE <- args[1]  #expanded alignmenta
INDEX.NO <- gsub("\\..*","",INPUT.FILE)
INDEX.NO <- gsub("pg-fmt-expanded-","",INDEX.NO)
print(INDEX.NO)

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

out.df <- select(n, -c(Fay.Wu.H, Zeng.E, Strobeck.S, x))
OUTFILE <- paste(INDEX.NO, "-PopGenome.tsv", sep = "")
write.table(out.df, file=OUTFILE, sep="\t", row.names = FALSE, quote = FALSE)
