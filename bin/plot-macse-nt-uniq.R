pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra)
args <-commandArgs(TRUE)

if (length(args)==0) {
    stop ("Need to specify macse codon distribution file for analysis",
        call.=FALSE)
}

INPUT.FILE <- args[1]
INDEX.NO <- gsub("\\..*","",INPUT.FILE)
print(INDEX.NO)

# Defining input variables
input.data <- read.table(INPUT.FILE, header = TRUE, sep = "\t", as.is = TRUE, comment.char = "#")

# Prepping the data
# This will give the x,y axes for the gene bars
#if (input.data$bstart[1] == input.data$start[1] - 1 && input.data$bstop[1] == input.data$bstop[1] - 1) {
#    bufferZone = "OFF"
#    print ("bufferZone is off")
#} else { 
#    bufferZone = "ON"
#    print ("bufferZone is on")
#}

input.data <- transform(input.data, x1 = input.data$bstart)
input.data <- transform(input.data, x2 = input.data$bstop)
geneLength <- input.data$bstop[1] - input.data$bstart[1]
spacingParam <- geneLength / 400
input.data <- transform(input.data, y1 = 1.9 + (input.data$allele - 1))
input.data <- transform(input.data, y2 = 2.1 + (input.data$allele - 1))

# This gives the length of the gene bar
x <- c(input.data$start[1]:input.data$stop[1])

print("Done calculating coordinates")

# Determining figure format	
plot.title <- paste("Unique Alleles (nt):", INDEX.NO, sep = " ")
png.title <- paste(INDEX.NO, "-macse-nt-uniq.png", sep = "")
png(png.title, width = 25, height = 25, units = "cm", res = 300) 

# Setting up plot skeleton
plot(
    1,
    type = "n",
    main = plot.title,
    xlab = "Genome Coordinate",
    ylab = "Unique Alleles",
    ylim = range(0 : max(input.data$allele) + 1),
    xlim = range(input.data$bstart[1] : input.data$bstop[1]),
    yaxt = 'n',		#will remove the y-axis
    frame = FALSE
);

# Determining if there's a buffer
#if (bufferZone == "ON") {
#rect(input.data$start[1], #xleft (gene start)
#    min(input.data$y1), #ybottom (y min)
#      input.data$stop[1], #xright (gene stop)
#    max(input.data$y2), #ytop (y max)
#      col="lightcyan", border=FALSE)
#}

counter = 0

# Prints the SNP frequency at the bottom of the plot
#for (i in 1 : nrow(snp.df)) {
#    text ( x = snp.df$Var1[i],
#           y = 1,
#           cex = 0.4,
#           paste(snp.df$Freq[i])
#    )
#}

for (i in 1 : nrow(input.data)) {
    if (counter != input.data$allele[i]) {
        counter = counter + 1

        # Setting the gene lines and mutations
        segments(input.data$x1[i],
                 input.data$y1[i] + .1,
                 input.data$x2[i],
                 input.data$y1[i] + .1, 
                 col = "snow3"
        )
    }
    if (isTRUE(input.data$type[i] == "TS")) {
        rect(input.data$mstart[i]-.05,
            input.data$y1[i],
            input.data$mstop[i]+.05,
            input.data$y2[i],
            col="#90BD6D", border="#90BD6D"
        )
    }
    else if (isTRUE(input.data$type[i] == "TV")) {
        rect(input.data$mstart[i]-.05,
            input.data$y1[i],
            input.data$mstop[i]+.05,
            input.data$y2[i],
            col="#2A9D8F", border="#2A9D8F"
        )
    }
    else if (isTRUE(input.data$type[i] == "AM")) {
        rect(input.data$mstart[i]-.05,
            input.data$y1[i],
            input.data$mstop[i]+.05,
            input.data$y2[i],
            col="#9197AE", border="#9197AE"
        )
    }
    else if (isTRUE(input.data$type[i] == "IN")) {
        rect(input.data$mstart[i]-.05,
            input.data$y1[i],
            input.data$mstop[i]+.05,
            input.data$y2[i],
            col = "#F9C74F", border="#F9C74F"
        )
    }
    else if (isTRUE(input.data$type[i] == "DEL")) {
        rect(input.data$mstart[i]-.05,
            input.data$y1[i],
            input.data$mstop[i]+.05,
            input.data$y2[i],
            col = "#CF113A", border="#CF113A"
        )
    }
    else if (isTRUE(input.data$type[i] == "FS")) {
        rect(input.data$mstart[i]-.05,
            input.data$y1[i],
            input.data$mstop[i]+.05,
            input.data$y2[i],
            col = "#F3722C", border="#F3722C"
        )
    }
       
    # Putting in allele frequency 
    text(y = input.data$allele[i] + 1,
        x = input.data$bstop[i] + spacingParam, #need a calc value for this
        labels = input.data$freq[i],
        cex = 0.5
    )
        # Putting in the total # mutations
    text(y = input.data$allele[i] + 1,
         x = input.data$bstop[i] + (spacingParam * 10),#need a calc value for this
         labels = input.data$tmut[i],
         cex = 0.5
    )
}

text(y = max(input.data$allele) + 1.5,
     x = input.data$bstop[1] + spacingParam,
     labels = "freq",
     cex = 0.5
)

text(y = max(input.data$allele) + 1.5,
     x = input.data$bstop[1] + (spacingParam * 10),
     labels = "#mut",
     cex = 0.5
)

dev.off()
print("Figure has been printed to file")

#End script
