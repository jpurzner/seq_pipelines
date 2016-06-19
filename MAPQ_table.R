args <- commandArgs(trailingOnly = TRUE)

qual_vals <- c(0:45)

qual_df <- as.data.frame(qual_vals)

for (i in 1:length(args)) {
    file_df <- read.table(args[i], sep = "", header = F)
    colnames(file_df) <- c(args[i], "MAPQ")
    qual_df <- merge(x = qual_df, y = file_df, by.x = 1, by.y = 2, all = T) 
}
qual_df[is.na(qual_df)] <- 0
colnames(qual_df) <- gsub("_sort_map.txt", "", colnames(qual_df))
write.table(qual_df, "MAPQ_table.txt", sep ="\t", col.names = T, row.names = F)
