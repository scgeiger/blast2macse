pacman::p_load(tidyverse, janitor)
args <-commandArgs(TRUE)

if (length(args)==0) {
    stop ("Need to specify dnds matrix file for analysis",
        call.=FALSE)
}

INPUT.FILE <- "am-ST131-acrR-removed-macse-uniqseqs-dnds-matrix.tsv" #args[1]
INDEX.NO <- gsub("\\..*","",INPUT.FILE)
print(INDEX.NO)

input.data <- read.table(INPUT.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#", stringsAsFactors = FALSE)

input.data <- input.data %>% remove_rownames %>% column_to_rownames(var="X")
freq.dnds <- as.data.frame(table(unlist(input.data)))
freq.dnds <- freq.dnds[!(freq.dnds$Var1=="-"),]
total.comps <- sum(freq.dnds$Freq)

dsiszero <- subset(freq.dnds, Var1 == "-1")
dnanddsiszero <- subset(freq.dnds, Var1 == "-2")
errs <- subset(freq.dnds, Var1 == "ERR")

freq.scores <- freq.dnds[!(freq.dnds$Var1=="-"  |
                         freq.dnds$Var1=="-1" |
                         freq.dnds$Var1=="-2" |
                         freq.dnds$Var1=="ERR"),]

freq.scores$Var1 <- as.numeric(paste(freq.scores$Var1))

dndsless1 <- subset(freq.scores, Var1 <= 1)
dndsless1 <- sum(dndsless1$Freq)
dndsmore1 <- subset(freq.scores, Var1 > 1 )
dndsmore1 <- sum(dndsmore1$Freq)

p <- ggplot(freq.scores, aes(x = Var1)) +
     geom_histogram(binwidth = .1)

# Proportion of hits:
# categories: ERR, ds is zero, dnandds is zero, scores <= 1, scores > 1

divide <- function(x, factor) {
    x / factor
}

ST <- c('ST131', 'ST131', 'ST131', 'ST131', 'ST131')
cat <- c('errs','dsiszero','dnanddsiszero','dndsless1','dndsmore1')
frequency <- c(errs$Freq, dsiszero$Freq, dnanddsiszero$Freq, dndsless1, dndsmore1)
proportion <- c(lapply(frequency, divide, factor = total.comps))
proportion <- as.character(proportion)
summary.values <- data.frame(ST, cat, frequency, proportion)

p <- ggplot(summary.values, aes(x = ST, y = frequency, fill = cat)) +
     geom_col(position = "fill") 


