###
#   File name : DE_Results_Comparison.R
#   Author    : Hyunjin Kim
#   Date      : Jun 15, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Compare the two DE results based on logFC and on the FDR
#               1. GMP CAR+ cells last vs do not last
#               2. Responder vs Non-responder
#
#   Instruction
#               1. Source("DE_Results_Comparison.R")
#               2. Run the function "de_results_comparison" - specify the input file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DE_Results_Comparison.R/DE_Results_Comparison.R")
#               > de_results_comparison(de_result1_path="./results/PROTO/DEEP/DE_Results/ALL_Patients_DE_Genes_DESeq2.xlsx",
#                                       de_result2_path="./results/PROTO/DEEP/DE_Results/JCC212_Aggreg21Feb_GMP_CARPOS__Markers.txt",
#                                       outputDir="./results/PROTO/DEEP/DE_Results/")
###

de_results_comparison <- function(de_result1_path="./results/PROTO/DEEP/DE_Results/ALL_Patients_DE_Genes_DESeq2.xlsx",
                                  de_result2_path="./results/PROTO/DEEP/DE_Results/JCC212_Aggreg21Feb_GMP_CARPOS__Markers.txt",
                                  outputDir="./results/PROTO/DEEP/DE_Results/") {
  
  ### load libraries
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(ggrepel, quietly = TRUE)) {
    install.packages("ggrepel")
    require(ggrepel, quietly = TRUE)
  }
  
  ### load DE results
  de_result1 <- read.xlsx2(file = de_result1_path, sheetIndex = 1,
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(de_result1) <- de_result1$Gene_Symbol
  de_result2 <- read.table(file = de_result2_path, header = TRUE, sep = "\t", row.names = 1,
                           stringsAsFactors = FALSE, check.names = FALSE)
  
  ### get shared gene symbols between the two
  ### In the DE result2, there are identical gene symbols and we use any first ones
  shared_genes <- intersect(rownames(de_result1), rownames(de_result2))
  
  ### combine two DE results
  combined_result <- data.frame(Gene_Symbol=shared_genes,
                                logFC1=de_result1[shared_genes,"avg_logFC"],
                                logFC2=de_result2[shared_genes,"avg_logFC"],
                                adj_p1=de_result1[shared_genes,"FDR"],
                                adj_p2=de_result2[shared_genes,"p_val_adj"],
                                Label="",
                                stringsAsFactors = FALSE, check.names = FALSE)
  combined_result <- combined_result[which(combined_result$adj_p1 != ""),]
  combined_result[2:5] <- sapply(combined_result[2:5], as.numeric)
  
  ### mark the interesting genes
  interesting_genes <- c("RPS26", "CD40LG", "TIMP1", "PPDPF", "YBX1", "RPS2", "HLA-DRB1")
  interesting_gene_idx <- which(combined_result$Gene_Symbol %in% interesting_genes)
  combined_result$Label[interesting_gene_idx] <- combined_result$Gene_Symbol[interesting_gene_idx]
  
  ### comparison based on logFC
  ggplot(data = combined_result, aes(x=logFC1, y=logFC2)) +
    geom_point(color = "black", size = 1) +
    geom_label_repel(aes(logFC1, logFC2, label = Label), color = "red", box.padding = unit(0.45, "lines")) +
    labs(title="logFC Comparison",
         subtitle=sprintf("S.Cor = %s, p-value = %s",
                          round(cor(combined_result$logFC1, combined_result$logFC2,
                                    use = "pairwise.complete.obs", method = "spearman"), 5),
                          signif(cor.test(combined_result$logFC1, combined_result$logFC2,
                                          method = "spearman")$p.value, 5))
    ) +
    xlab("logFC of (GMP CAR+ Last vs NOT)") +
    ylab("logFC of (Responder vs Non-responder)") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  
  
  
  
}
