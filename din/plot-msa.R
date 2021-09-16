cat *.dealign >> all-yjdC-consensus.fasta
#change seq names

library(msa)

file <- "all-yjdC-consensus.fasta"

seq <- readAAStringSet(file)
aln <- msaClustalOmega(seq)

sink("yjdC.aln")
print(aln, show="complete")
sink()
q()
n

