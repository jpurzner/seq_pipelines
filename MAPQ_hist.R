# generates a histogram from MAPQ file  

require(ggplot2)
require(stringr)

args <- commandArgs(trailingOnly = TRUE)

print(args)

qual_vals <- c(0:45)
qual_df <- as.data.frame(qual_vals)

file_df <- read.table(args[1], sep = "", header = F)
qual_df <- merge(x = qual_df, y = file_df, by.x = 1, by.y = 2, all = T)

colnames(qual_df) <- c("MAPQ", "count")
name = gsub("_mapq.txt", "", args[1])

pdf(str_c(name, "_mapq.pdf"))

ggplot(qual_df, aes(x = MAPQ, y = count)) + geom_bar(stat="identity") + scale_x_discrete(breaks = c(0,5,10,15,20,25,30,35,40,45), labels = c("0", "5", "10", "15", "20", "25", "30", "35", "40", "45"), name = "Mapping quality") + scale_y_continuous(name = "total reads") + ggtitle(basename(name))


dev.off()
