### This was implemented based on Ching-Heng's request

if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}
if(!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2")
  require(ggplot2, quietly = TRUE)
}

### nFeature_RNA is the number of genes detected in each cell
### nCount_RNA is the total number of molecules detected within a cell

CH_table <- data.frame(Library=unique(SeuratObj_meta$library),
                       stringsAsFactors = FALSE, check.names = FALSE)

CH_table$Cell_Num <- sapply(unique(SeuratObj_meta$library), function(x) length(which(SeuratObj_meta$library == x)))

CH_table$Mean_nCount_RNA <- sapply(unique(SeuratObj_meta$library), function(x) {
  return(mean(as.numeric(SeuratObj_meta$nCount_RNA[which(SeuratObj_meta$library == x)])))
})

CH_table$Mean_nFeature_RNA <- sapply(unique(SeuratObj_meta$library), function(x) {
  return(mean(as.numeric(SeuratObj_meta$nFeature_RNA[which(SeuratObj_meta$library == x)])))
})

write.xlsx(CH_table, file = "./CH_table.xlsx", row.names = FALSE)

### visualizing
df <- data.frame(Library=rep(CH_table$Library, 3),
                 Value=c(CH_table$Cell_Num, CH_table$Mean_nCount_RNA, CH_table$Mean_nFeature_RNA),
                 Group=c(rep("Cell_Num", nrow(CH_table)), rep("Mean_nCount_RNA", nrow(CH_table)), rep("Mean_nFeature_RNA", nrow(CH_table))),
                 stringsAsFactors = FALSE, check.names = FALSE)
p <- ggplot(data = df, aes_string(x = "Library", y = "Value", group = "Group")) +
  geom_line(aes_string(color = "Group"), size = 2) +
  labs(colour="", x= "", y="Numbers") +
  theme_classic(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
        axis.text.y = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
        axis.title = element_text(size = 35, color = "black", face = "bold"),
        legend.title = element_text(size = 30, color = "black", face = "bold"),
        legend.text = element_text(size = 25, color = "black", face = "bold"),
        legend.key.size = unit(0.7, 'cm'))
ggsave(file = "./CH_Figure.png", plot = p,
       width = 25, height = 10, dpi = 350)

