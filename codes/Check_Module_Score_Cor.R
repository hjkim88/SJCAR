### check correlations of effector precursor signature (module score percentage)
### to any car expansion numbers
### with varying
### 1. module score cut-off
### 2. input DE gene #
### 3. log scale
### 4. include px12
### 5. CD8 definition by EXP
### 6. TIGIT > 1
### 7. cluster3 or cluster8 independently
### 8. using DE genes from real CD8 GMP precursor vs Other CD8 GMPs
### 9 . min.pct, logFC, test method change in DE analysis
### 10. nbin, ctrl, seed change

### load libraries
if(!require(Seurat, quietly = TRUE)) {
  install.packages("Seurat")
  require(Seurat, quietly = TRUE)
}
if(!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2")
  require(ggplot2, quietly = TRUE)
}
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}
if(!require(ggrepel, quietly = TRUE)) {
  install.packages("ggrepel")
  require(ggrepel, quietly = TRUE)
}

### set variables
module_score_cut_off <- c(-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2)
input_de_gene_num <- c(3, 5, 10, 20, 30, 50, 100, 200) # TIGIT, SELL, CD27 or all significant DE genes
nbins <- 24 # default: 24
nctrl_genes <- 100 # default: 100
seed <- 1234 # default: 1

### load Jeremy's object
JCC_Seurat_Obj <- readRDS(file = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds")

### check whether the orders are the same
print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
print(identical(names(Idents(object = JCC_Seurat_Obj)), rownames(JCC_Seurat_Obj@meta.data)))

### make a correlation plot with b cell recovery data
b_cell_recovery_time <- rep(NA, length(unique(JCC_Seurat_Obj$px)))
names(b_cell_recovery_time) <- unique(JCC_Seurat_Obj$px)
b_cell_recovery_time["SJCAR19-02"] <- 73
b_cell_recovery_time["SJCAR19-03"] <- 85
b_cell_recovery_time["SJCAR19-05"] <- 174
b_cell_recovery_time["SJCAR19-06"] <- 131
b_cell_recovery_time["SJCAR19-15"] <- 44

### tumor burden
Tumor_Burden <- c(98, 10, 1, 0, 12, 1, 80, 72, 78, 84, 0, 2, 51, 43, 0)
names(Tumor_Burden) <- unique(JCC_Seurat_Obj$px)[2:16]

### correlate with peak expansion and wk1 expansion
peakcar_ug <- c(199054, 42149, 61777, 288670, 224445, 167092, 4806, 142422,
                301705, 356424, 64212, 2193043, 87863, 18532)
peakcar_ml <- c(82939, 561988, 947250, 2165026, 1122227, 3369685, 76891, 569687,
                5732386, 2024664, 2825332, 17544341, 1640116, 520669)
wk1car_ug <- c(3114, 13912, 61777, 3745, 224445, 1011, 657, 46046,
               268449, 356424, 8071, 6403, 87863, 18532)
wk2car_ug <- c(199054, 42149, 10100, 288670, 169344, 167092, 565, 142422,
               301705, 68633, 64212, 2193043, 3635, 15465)
wk3car_ug <- c(13744, 8992, 681, 38325, 7013, 47074, 4806, 55085,
               218312, 15701, 36722, 341921, 22611, 3786)
wk1car_ml <- c(9341, 296787, 947250, 24967, 1122227, 50700, 7561, 314646,
               268449, 1960334, 236747, 49091, 1640116, 395350)
wk2car_ml <- c(82939, 561988, 254181, 2165026, 254016, 3369685, 16755, 569687,
               5732386, 2024664, 2825332, 17544341, 69673, 520669)
wk3car_ml <- c(6414, 149865, 16461, 421573, 31556, 894403, 76891, 321331,
               2328662, 578331, 1621885, 4843874, 471052, 81407)
names(peakcar_ug) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                       "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")
names(peakcar_ml) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                       "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")
names(wk1car_ug) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                      "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                      "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")
names(wk2car_ug) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                      "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                      "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")
names(wk3car_ug) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                      "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                      "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")
names(wk1car_ml) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                      "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                      "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")
names(wk2car_ml) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                      "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                      "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")
names(wk3car_ml) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                      "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                      "SJCAR19-11", "SJCAR19-13", "SJCAR19-14", "SJCAR19-15")

### pairs add-ons
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x, plot = FALSE)
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "#D39F3A", ...)
  # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}
panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                      cex = 1, col.smooth = "#487A8F", span = 2/3, iter = 3, ...)  {
  reg <- function(x, y, col) abline(lm(y~x), col=col)
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  Cor_PV <- round(cor.test(x, y, method = "spearman", use = "complete.obs")$p.value, 2)
  if(!is.na(Cor_PV)) {
    if(Cor_PV < 0.05) col.smooth = "red"
    if (any(ok)) reg(x[ok], y[ok], col.smooth)
  }
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- cor(x, y, method = "spearman", use = "complete.obs")
  Cor_PV <- round(cor.test(x, y, method = "spearman", use = "complete.obs")$p.value, 2)
  if(!is.na(Cor_PV)) {
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
      cex.cor <- 0.4 / strwidth(txt)
    }
    if(Cor_PV < 0.05) {
      text(0.5, 0.5, paste0("r=", txt, "\np=", Cor_PV),
           cex = 1 + cex.cor * abs(Cor),
           col = "red")
    } else {
      text(0.5, 0.5, paste0("r=", txt, "\np=", Cor_PV),
           cex = 0.8 + cex.cor * abs(Cor))
    }
  }
}

### CD4/CD8 annotation
JCC_Seurat_Obj$CD4_CD8_by_Clusters <- "NA"
JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "2", "9", "10", "11", "14", "15", "18"))] <- "CD4"
JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20"))] <- "CD8"

### 1. the normalized number of TIGIT+ cells in each patient
print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
plot(density(JCC_Seurat_Obj@assays$RNA@counts["TIGIT",]))
gmp_tigit_pos_idx <- intersect(which(JCC_Seurat_Obj@assays$RNA@counts["TIGIT",] > 0),
                               which(JCC_Seurat_Obj$time2 == "GMP"))
TIGIT_Pos_Cell_Num <- rep(0, length(unique(JCC_Seurat_Obj$px)))
names(TIGIT_Pos_Cell_Num) <- unique(JCC_Seurat_Obj$px)
TIGIT_Pos_CD8_Cell_Num <- rep(0, length(unique(JCC_Seurat_Obj$px)))
names(TIGIT_Pos_CD8_Cell_Num) <- unique(JCC_Seurat_Obj$px)
for(px in unique(JCC_Seurat_Obj$px)) {
  TIGIT_Pos_Cell_Num[px] <- length(intersect(which(JCC_Seurat_Obj$px == px),
                                             gmp_tigit_pos_idx)) / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                    which(JCC_Seurat_Obj$time2 == "GMP")))
  TIGIT_Pos_CD8_Cell_Num[px] <- length(intersect(which(JCC_Seurat_Obj$px == px),
                                                 intersect(gmp_tigit_pos_idx,
                                                           which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8")))) / length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                                              which(JCC_Seurat_Obj$time2 == "GMP")),
                                                                                                                                    which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8")))
}

### PI subsisters in the cluster3&8
cluster3_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                 which(JCC_Seurat_Obj$AllSeuratClusters == "3")),
                                                                       which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
cluster8_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                 which(JCC_Seurat_Obj$AllSeuratClusters == "8")),
                                                                       which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
cluster38_pi_subsisters <- c(cluster3_pi_subsisters, cluster8_pi_subsisters)

### PI subsister clones in the cluster3&8
cluster3_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster3_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
cluster8_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster8_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
cluster38_pi_subsister_clones <- unique(c(cluster3_pi_subsister_clones, cluster8_pi_subsister_clones))

### remove NA clones
cluster3_pi_subsister_clones <- cluster3_pi_subsister_clones[which(!is.na(cluster3_pi_subsister_clones))]
cluster8_pi_subsister_clones <- cluster8_pi_subsister_clones[which(!is.na(cluster8_pi_subsister_clones))]
cluster38_pi_subsister_clones <- cluster38_pi_subsister_clones[which(!is.na(cluster38_pi_subsister_clones))]

### get the GMP subsister indicies that end up in PI cluster 3 & 8
GMP_Subsisters_PI_Cluster38_idx <- intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                             which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% cluster38_pi_subsister_clones))

### add a column for the info
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 <- "Others"
JCC_Seurat_Obj@meta.data[cluster38_pi_subsisters, "GMP_Subsisters_End_Up_In_Cluster38"] <- "PI_Subsisters_In_Cluster_3_And_8"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[GMP_Subsisters_PI_Cluster38_idx] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[intersect(which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "Others"),
                                                            which(JCC_Seurat_Obj$GMP_CARpos_Persister == "YES"))] <- "Other_GMP_Subsisters"

### GMP subsisters end up in cluster3 and 8 vs other CD8 GMP subsisters
JCC_Seurat_Obj@meta.data$GMP_Subsisters_End_Up_In_Cluster38_CD8 <- "Others"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_CD8[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_CD8[intersect(which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"),
                                                                which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "Other_GMP_Subsisters"))] <- "Other_CD8_GMP_Subsisters"

### GMP subsisters end up in cluster3 and 8 vs all other GMPs
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 <- "Others"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[setdiff(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                            which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8"))] <- "Other_GMPs"

### GMP subsisters end up in cluster3 and 8 vs all other CD8 GMPs
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 <- "Others"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[intersect(which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"),
                                                                  which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "Other_GMPs"))] <- "Other_CD8_GMPs"

### GMP subsisters end up in cluster3 and 8 vs other GMP subsisters
JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                           cells = rownames(JCC_Seurat_Obj@meta.data),
                           value = JCC_Seurat_Obj@meta.data$GMP_Subsisters_End_Up_In_Cluster38_2_CD8)
de_result <- FindMarkers(JCC_Seurat_Obj,
                         ident.1 = "GMP_Subsisters_End_Up_In_Cluster_3_And_8",
                         ident.2 = "Other_CD8_GMPs",
                         min.pct = 0.2,
                         logfc.threshold = 0.2,
                         test.use = "wilcox")

### with varying the DE gene #
test_set <- list()
test_cnt <- 1
for(de_gene_num in input_de_gene_num) {
  
  ### print progress
  writeLines(paste(de_gene_num))
  
  ### add module scores
  JCC_Seurat_Obj <- AddModuleScore(JCC_Seurat_Obj,
                                   features = list(rownames(de_result)[which(de_result$avg_log2FC > 0)][1:de_gene_num]),
                                   nbin = nbins,
                                   ctrl = nctrl_genes,
                                   seed = seed,
                                   name="PM_positive")
  JCC_Seurat_Obj <- AddModuleScore(JCC_Seurat_Obj,
                                   features = list(rownames(de_result)[which(de_result$avg_log2FC < 0)][1:de_gene_num]),
                                   nbin = nbins,
                                   ctrl = nctrl_genes,
                                   seed = seed,
                                   name="PM_negative")
  
  ### pairwise factor set
  factor1_list <- c("TIGIT_Cell_Num", "TIGIT_CD8_Cell_Num", "Precursor_Pcnt", "Precursor_CD8_Pcnt",
                    "Precursor_Pcnt2", "Precursor_CD8_Pcnt2", "Precursor_Pcnt3", "Precursor_CD8_Pcnt3")
  factor2_list <- c("B_Cell_Recovery_Time", "PeakCAR_ug", "Wk1CAR_ug", "Wk2CAR_ug", "Wk3CAR_ug",
                    "PeakCAR_ml", "Wk1CAR_ml", "Wk2CAR_ml", "Wk3CAR_ml", "Tumor_Burden")
  test_df <- data.frame(Variable1=rep("", length(factor1_list)*length(factor2_list)*length(module_score_cut_off)),
                        Variable2=rep("", length(factor1_list)*length(factor2_list)*length(module_score_cut_off)),
                        Threshold=rep("", length(factor1_list)*length(factor2_list)*length(module_score_cut_off)),
                        Cor=NA,
                        PVal=NA,
                        Adj.Pval=NA,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  cnt <- 1
  for(thresh in module_score_cut_off) {
    
    ### get the percentage
    Precursor_Pcnt <- rep(0, length(unique(JCC_Seurat_Obj$px)))
    names(Precursor_Pcnt) <- unique(JCC_Seurat_Obj$px)
    Precursor_CD8_Pcnt <- rep(0, length(unique(JCC_Seurat_Obj$px)))
    names(Precursor_CD8_Pcnt) <- unique(JCC_Seurat_Obj$px)
    Precursor_Pcnt2 <- rep(0, length(unique(JCC_Seurat_Obj$px)))
    names(Precursor_Pcnt2) <- unique(JCC_Seurat_Obj$px)
    Precursor_CD8_Pcnt2 <- rep(0, length(unique(JCC_Seurat_Obj$px)))
    names(Precursor_CD8_Pcnt2) <- unique(JCC_Seurat_Obj$px)
    Precursor_Pcnt3 <- rep(0, length(unique(JCC_Seurat_Obj$px)))
    names(Precursor_Pcnt3) <- unique(JCC_Seurat_Obj$px)
    Precursor_CD8_Pcnt3 <- rep(0, length(unique(JCC_Seurat_Obj$px)))
    names(Precursor_CD8_Pcnt3) <- unique(JCC_Seurat_Obj$px)
    
    for(px in unique(JCC_Seurat_Obj$px)) {
      Precursor_Pcnt[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                       which(JCC_Seurat_Obj$time2 == "GMP")),
                                             which(JCC_Seurat_Obj$PM_positive1 > thresh))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                    which(JCC_Seurat_Obj$time2 == "GMP")))
      
      Precursor_CD8_Pcnt[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                           which(JCC_Seurat_Obj$time2 == "GMP")),
                                                 intersect(which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"),
                                                           which(JCC_Seurat_Obj$PM_positive1 > thresh)))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                                   intersect(which(JCC_Seurat_Obj$time2 == "GMP"),
                                                                                                                                             which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"))))
      
      Precursor_Pcnt2[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                        which(JCC_Seurat_Obj$time2 == "GMP")),
                                              which(JCC_Seurat_Obj$PM_negative1 < thresh))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                     which(JCC_Seurat_Obj$time2 == "GMP")))
      
      Precursor_CD8_Pcnt2[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                            which(JCC_Seurat_Obj$time2 == "GMP")),
                                                  intersect(which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"),
                                                            which(JCC_Seurat_Obj$PM_negative1 < thresh)))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                                    intersect(which(JCC_Seurat_Obj$time2 == "GMP"),
                                                                                                                                              which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"))))
      
      Precursor_Pcnt3[px] <- (Precursor_Pcnt[px] + Precursor_Pcnt2[px]) / 2
      Precursor_CD8_Pcnt3[px] <- (Precursor_CD8_Pcnt[px] + Precursor_CD8_Pcnt2[px]) / 2
      
    }
    
    plot_df <- data.frame(Patient=names(peakcar_ug),
                          TIGIT_Cell_Num=as.numeric(TIGIT_Pos_Cell_Num[names(peakcar_ug)]),
                          TIGIT_CD8_Cell_Num=as.numeric(TIGIT_Pos_CD8_Cell_Num[names(peakcar_ug)]),
                          Precursor_Pcnt=as.numeric(Precursor_Pcnt[names(peakcar_ug)]),
                          Precursor_CD8_Pcnt=as.numeric(Precursor_CD8_Pcnt[names(peakcar_ug)]),
                          Precursor_Pcnt2=as.numeric(Precursor_Pcnt2[names(peakcar_ug)]),
                          Precursor_CD8_Pcnt2=as.numeric(Precursor_CD8_Pcnt2[names(peakcar_ug)]),
                          Precursor_Pcnt3=as.numeric(Precursor_Pcnt3[names(peakcar_ug)]),
                          Precursor_CD8_Pcnt3=as.numeric(Precursor_CD8_Pcnt3[names(peakcar_ug)]),
                          B_Cell_Recovery_Time=as.numeric(b_cell_recovery_time[names(peakcar_ug)]),
                          PeakCAR_ug=as.numeric(peakcar_ug),
                          Wk1CAR_ug=as.numeric(wk1car_ug[names(peakcar_ug)]),
                          Wk2CAR_ug=as.numeric(wk2car_ug[names(peakcar_ug)]),
                          Wk3CAR_ug=as.numeric(wk3car_ug[names(peakcar_ug)]),
                          PeakCAR_ml=as.numeric(peakcar_ml),
                          Wk1CAR_ml=as.numeric(wk1car_ml[names(peakcar_ug)]),
                          Wk2CAR_ml=as.numeric(wk2car_ml[names(peakcar_ug)]),
                          Wk3CAR_ml=as.numeric(wk3car_ml[names(peakcar_ug)]),
                          Tumor_Burden=as.numeric(Tumor_Burden[names(peakcar_ug)]),
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    for(a in factor1_list) {
      for(b in factor2_list) {
        x <- as.numeric(plot_df[,a])
        y <- as.numeric(plot_df[,b])
        Cor <- round(cor(x, y, method = "spearman", use = "complete.obs"), 2)
        Cor_PV <- round(cor.test(x, y, method = "spearman", use = "complete.obs")$p.value, 2)
        
        test_df$Variable1[cnt] <- a
        test_df$Variable2[cnt] <- b
        test_df$Threshold[cnt] <- thresh
        test_df$Cor[cnt] <- Cor
        test_df$PVal[cnt] <- Cor_PV
        
        cnt <- cnt + 1
      }
    }
    
    test_df$Adj.Pval[which(test_df$Threshold == thresh)] <- p.adjust(test_df$PVal[which(test_df$Threshold == thresh)], method = "BH")
    
    ### Multiple correlation plots
    pdf(file = paste0("./results/New3/Manuscript/Module_Score_Correlations/Original_Cors_pairwise_plot_", de_gene_num, "_", thresh, ".pdf"), width = 20, height = 15)
    pairs(data=plot_df,
          ~TIGIT_Cell_Num + TIGIT_CD8_Cell_Num + Precursor_Pcnt + Precursor_CD8_Pcnt + Precursor_Pcnt2 +
            Precursor_CD8_Pcnt2 + Precursor_Pcnt3 + Precursor_CD8_Pcnt3 + B_Cell_Recovery_Time + PeakCAR_ug + PeakCAR_ml +
            Wk1CAR_ug + Wk2CAR_ug + Wk3CAR_ug + Wk1CAR_ml + Wk2CAR_ml + Wk3CAR_ml + Tumor_Burden,
          upper.panel = panel.lm, lower.panel = panel.cor, diag.panel = panel.hist,
          font.labels = 2, pch = 19)
    dev.off()
  }
  
  test_df <- test_df[order(test_df$Adj.Pval),]
  test_df$DE_Gene_Num <- de_gene_num
  
  test_set[[test_cnt]] <- test_df
  test_cnt <- test_cnt + 1
}
names(test_set) <- de_gene_num

### unlist test_df and order it
test_df2 <- Reduce(function(d1, d2) rbind(d1, d2), test_set)
test_df2 <- test_df2[order(test_df2$Adj.Pval),]


### just look at the 3 genes - TIGIT, SELL, CD27

### add module scores
JCC_Seurat_Obj <- AddModuleScore(JCC_Seurat_Obj,
                                 features = list(c("TIGIT")),
                                 nbin = nbins,
                                 ctrl = nctrl_genes,
                                 seed = seed,
                                 name="TIGIT_Module")
JCC_Seurat_Obj <- AddModuleScore(JCC_Seurat_Obj,
                                 features = list(c("SELL", "CD27")),
                                 nbin = nbins,
                                 ctrl = nctrl_genes,
                                 seed = seed,
                                 name="SELL_CD27_Module")

### pairwise factor set
factor1_list <- c("TIGIT_Pcnt", "TIGIT_CD8_Pcnt",
                  "SELL_CD27_Pcnt", "SELL_CD27_CD8_Pcnt", "Precursor_Pcnt3", "Precursor_CD8_Pcnt3")
factor2_list <- c("B_Cell_Recovery_Time", "PeakCAR_ug", "Wk1CAR_ug", "Wk2CAR_ug", "Wk3CAR_ug",
                  "PeakCAR_ml", "Wk1CAR_ml", "Wk2CAR_ml", "Wk3CAR_ml", "Tumor_Burden")
test_df3 <- data.frame(Variable1=rep("", length(factor1_list)*length(factor2_list)*length(module_score_cut_off)),
                       Variable2=rep("", length(factor1_list)*length(factor2_list)*length(module_score_cut_off)),
                       Threshold=rep("", length(factor1_list)*length(factor2_list)*length(module_score_cut_off)),
                       Cor=NA,
                       PVal=NA,
                       Adj.Pval=NA,
                       stringsAsFactors = FALSE, check.names = FALSE)

cnt <- 1
for(thresh in module_score_cut_off) {
  
  ### get the percentage
  TIGIT_Pcnt <- rep(0, length(unique(JCC_Seurat_Obj$px)))
  names(TIGIT_Pcnt) <- unique(JCC_Seurat_Obj$px)
  TIGIT_CD8_Pcnt <- rep(0, length(unique(JCC_Seurat_Obj$px)))
  names(TIGIT_CD8_Pcnt) <- unique(JCC_Seurat_Obj$px)
  SELL_CD27_Pcnt <- rep(0, length(unique(JCC_Seurat_Obj$px)))
  names(SELL_CD27_Pcnt) <- unique(JCC_Seurat_Obj$px)
  SELL_CD27_CD8_Pcnt <- rep(0, length(unique(JCC_Seurat_Obj$px)))
  names(SELL_CD27_CD8_Pcnt) <- unique(JCC_Seurat_Obj$px)
  Precursor_Pcnt3 <- rep(0, length(unique(JCC_Seurat_Obj$px)))
  names(Precursor_Pcnt3) <- unique(JCC_Seurat_Obj$px)
  Precursor_CD8_Pcnt3 <- rep(0, length(unique(JCC_Seurat_Obj$px)))
  names(Precursor_CD8_Pcnt3) <- unique(JCC_Seurat_Obj$px)
  
  for(px in unique(JCC_Seurat_Obj$px)) {
    TIGIT_Pcnt[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                 which(JCC_Seurat_Obj$time2 == "GMP")),
                                       which(JCC_Seurat_Obj$TIGIT_Module1 > thresh))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                               which(JCC_Seurat_Obj$time2 == "GMP")))
    
    TIGIT_CD8_Pcnt[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                     which(JCC_Seurat_Obj$time2 == "GMP")),
                                           intersect(which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"),
                                                     which(JCC_Seurat_Obj$TIGIT_Module1 > thresh)))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                              intersect(which(JCC_Seurat_Obj$time2 == "GMP"),
                                                                                                                                        which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"))))
    
    SELL_CD27_Pcnt[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                     which(JCC_Seurat_Obj$time2 == "GMP")),
                                           which(JCC_Seurat_Obj$SELL_CD27_Module1 < thresh))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                       which(JCC_Seurat_Obj$time2 == "GMP")))
    
    SELL_CD27_CD8_Pcnt[px] <- length(intersect(intersect(which(JCC_Seurat_Obj$px == px),
                                                         which(JCC_Seurat_Obj$time2 == "GMP")),
                                               intersect(which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"),
                                                         which(JCC_Seurat_Obj$SELL_CD27_Module1 < thresh)))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == px),
                                                                                                                                      intersect(which(JCC_Seurat_Obj$time2 == "GMP"),
                                                                                                                                                which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8"))))
    
    Precursor_Pcnt3[px] <- (TIGIT_Pcnt[px] + SELL_CD27_Pcnt[px]) / 2
    Precursor_CD8_Pcnt3[px] <- (TIGIT_CD8_Pcnt[px] + SELL_CD27_CD8_Pcnt[px]) / 2
    
  }
  
  plot_df <- data.frame(Patient=names(peakcar_ug),
                        TIGIT_Pcnt=as.numeric(TIGIT_Pcnt[names(peakcar_ug)]),
                        TIGIT_CD8_Pcnt=as.numeric(TIGIT_CD8_Pcnt[names(peakcar_ug)]),
                        SELL_CD27_Pcnt=as.numeric(SELL_CD27_Pcnt[names(peakcar_ug)]),
                        SELL_CD27_CD8_Pcnt=as.numeric(SELL_CD27_CD8_Pcnt[names(peakcar_ug)]),
                        Precursor_Pcnt3=as.numeric(Precursor_Pcnt3[names(peakcar_ug)]),
                        Precursor_CD8_Pcnt3=as.numeric(Precursor_CD8_Pcnt3[names(peakcar_ug)]),
                        B_Cell_Recovery_Time=as.numeric(b_cell_recovery_time[names(peakcar_ug)]),
                        PeakCAR_ug=as.numeric(peakcar_ug),
                        Wk1CAR_ug=as.numeric(wk1car_ug[names(peakcar_ug)]),
                        Wk2CAR_ug=as.numeric(wk2car_ug[names(peakcar_ug)]),
                        Wk3CAR_ug=as.numeric(wk3car_ug[names(peakcar_ug)]),
                        PeakCAR_ml=as.numeric(peakcar_ml),
                        Wk1CAR_ml=as.numeric(wk1car_ml[names(peakcar_ug)]),
                        Wk2CAR_ml=as.numeric(wk2car_ml[names(peakcar_ug)]),
                        Wk3CAR_ml=as.numeric(wk3car_ml[names(peakcar_ug)]),
                        Tumor_Burden=as.numeric(Tumor_Burden[names(peakcar_ug)]),
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  for(a in factor1_list) {
    for(b in factor2_list) {
      x <- as.numeric(plot_df[,a])
      y <- as.numeric(plot_df[,b])
      Cor <- round(cor(x, y, method = "spearman", use = "complete.obs"), 2)
      Cor_PV <- round(cor.test(x, y, method = "spearman", use = "complete.obs")$p.value, 2)
      
      test_df3$Variable1[cnt] <- a
      test_df3$Variable2[cnt] <- b
      test_df3$Threshold[cnt] <- thresh
      test_df3$Cor[cnt] <- Cor
      test_df3$PVal[cnt] <- Cor_PV
      
      cnt <- cnt + 1
    }
  }
  
  test_df3$Adj.Pval[which(test_df3$Threshold == thresh)] <- p.adjust(test_df3$PVal[which(test_df3$Threshold == thresh)], method = "BH")
}
test_df3 <- test_df3[order(test_df3$Adj.Pval),]

### write out
saveRDS(test_df2, file = "./Original_Cors.rds")
saveRDS(test_df3, file = "./Only_3_Cors.rds")
write.xlsx2(test_df2, file = "./Original_Cors.xlsx", row.names = FALSE)
write.xlsx2(test_df3, file = "./Only_3_Cors.xlsx", row.names = FALSE)

### distribution
pdf(file = paste0("./Original_Cors_positive_distribution.pdf"), width = 10, height = 10)
plot(density(JCC_Seurat_Obj$PM_positive1[which(JCC_Seurat_Obj$time2 == "GMP")]),
     main = "Precursor_Module_Score_Positive")
dev.off()
pdf(file = paste0("./Original_Cors_negative_distribution.pdf"), width = 10, height = 10)
plot(density(JCC_Seurat_Obj$PM_negative1[which(JCC_Seurat_Obj$time2 == "GMP")]),
     main = "Precursor_Module_Score_Negative")
dev.off()

pdf(file = paste0("./Only_3_Cors_positive_distribution.pdf"), width = 10, height = 10)
plot(density(JCC_Seurat_Obj$TIGIT_Module1[which(JCC_Seurat_Obj$time2 == "GMP")]),
     main = "Precursor_Module_Score_Positive")
dev.off()
pdf(file = paste0("./Only_3_Cors_negative_distribution.pdf"), width = 10, height = 10)
plot(density(JCC_Seurat_Obj$SELL_CD27_Module1[which(JCC_Seurat_Obj$time2 == "GMP")]),
     main = "Precursor_Module_Score_Negative")
dev.off()


### subset GMP cells only and calculate the module scores
### is there a difference with calculating with all the cells?
### compare the distribution to the result from PI only?

### subset the object into GMP only and PI only
GMP_Seurat_Obj <- subset(JCC_Seurat_Obj,
                         cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$GMP == "GMP")])
PI_Seurat_Obj <- subset(JCC_Seurat_Obj,
                        cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$GMP == "PI")])

### calculate module scores
JCC_Seurat_Obj <- AddModuleScore(JCC_Seurat_Obj,
                                 features = list(rownames(de_result)[which(de_result$avg_log2FC > 0)][1:20]),
                                 nbin = nbins,
                                 ctrl = nctrl_genes,
                                 seed = seed,
                                 name="Precursor_Pcnt")

GMP_Seurat_Obj <- AddModuleScore(GMP_Seurat_Obj,
                                 features = list(rownames(de_result)[which(de_result$avg_log2FC > 0)][1:20]),
                                 nbin = nbins,
                                 ctrl = nctrl_genes,
                                 seed = seed,
                                 name="Precursor_Pcnt")

PI_Seurat_Obj <- AddModuleScore(PI_Seurat_Obj,
                                features = list(rownames(de_result)[which(de_result$avg_log2FC > 0)][1:20]),
                                nbin = nbins,
                                ctrl = nctrl_genes,
                                seed = seed,
                                name="Precursor_Pcnt")

pdf(file = paste0("./results/New3/Manuscript/Module_Score_Tests/GMP_From_All_Module_Score_Positive_Distribution.pdf"), width = 10, height = 10)
plot(density(JCC_Seurat_Obj$Precursor_Pcnt1[which(JCC_Seurat_Obj$GMP == "GMP")]),
     main = "GMP_From_All_Precursor_Module_Score_Positive")
dev.off()

pdf(file = paste0("./results/New3/Manuscript/Module_Score_Tests/GMP_Only_Module_Score_Positive_Distribution.pdf"), width = 10, height = 10)
plot(density(GMP_Seurat_Obj$Precursor_Pcnt1[which(GMP_Seurat_Obj$GMP == "GMP")]),
     main = "GMP_ONLY_Precursor_Module_Score_Positive")
dev.off()

pdf(file = paste0("./results/New3/Manuscript/Module_Score_Tests/PI_Only_Module_Score_Positive_Distribution.pdf"), width = 10, height = 10)
plot(density(PI_Seurat_Obj$Precursor_Pcnt1[which(PI_Seurat_Obj$GMP == "PI")]),
     main = "PI_ONLY_Precursor_Module_Score_Positive")
dev.off()

a <- JCC_Seurat_Obj$Precursor_Pcnt1[which(JCC_Seurat_Obj$GMP == "GMP")]
b <- GMP_Seurat_Obj$Precursor_Pcnt1[which(GMP_Seurat_Obj$GMP == "GMP")]
identical(names(a), names(b))
cor(a, b)
identical(a, b)
mean(a)
mean(b)

### why length(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
###                      which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% cluster38_pi_subsister_clones)))
### and length(which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "PI_Subsisters_In_Cluster_3_And_8"))
### are different? 321 vs 291
### ANSWER: In the first one, those clones can be also in other PI clusters but the second one is only in PI cluster3&8


### GMP effector precursor proportion correlates with the proportion of effectors post infusion
PI_Effector_Pcnt <- sapply(unique(JCC_Seurat_Obj$px), function(x) {
  return(length(intersect(which(JCC_Seurat_Obj$px == x),
                          which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "PI_Subsisters_In_Cluster_3_And_8"))) * 100 / length(intersect(which(JCC_Seurat_Obj$px == x),
                                                                                                                                                    which(JCC_Seurat_Obj$GMP == "PI"))))
})

### calculate module score
JCC_Seurat_Obj <- AddModuleScore(JCC_Seurat_Obj,
                                 features = list(rownames(de_result)[which(de_result$avg_log2FC > 0)][1:20]),
                                 nbin = nbins,
                                 ctrl = nctrl_genes,
                                 seed = seed,
                                 name="Precursor_Pcnt")
JCC_Seurat_Obj <- AddModuleScore(JCC_Seurat_Obj,
                                 features = list(rownames(de_result)[which(de_result$avg_log2FC < 0)][1:20]),
                                 nbin = nbins,
                                 ctrl = nctrl_genes,
                                 seed = seed,
                                 name="Precursor_Pcnt2")

### positive
Avg_Pos_Module_Score <- sapply(unique(JCC_Seurat_Obj$px), function(x) {
  return(mean(JCC_Seurat_Obj$Precursor_Pcnt1[intersect(which(JCC_Seurat_Obj$px == x),
                                                       which(JCC_Seurat_Obj$GMP == "GMP"))]))
})
Avg_Neg_Module_Score <- sapply(unique(JCC_Seurat_Obj$px), function(x) {
  return(mean(JCC_Seurat_Obj$Precursor_Pcnt21[intersect(which(JCC_Seurat_Obj$px == x),
                                                        which(JCC_Seurat_Obj$GMP == "GMP"))]))
})

### draw the correlation plot
plot_df <- data.frame(Patient=names(PI_Effector_Pcnt),
                      PI_Effector_Pcnt=PI_Effector_Pcnt,
                      Avg_Pos_Module_Score=Avg_Pos_Module_Score,
                      Avg_Neg_Module_Score=Avg_Neg_Module_Score,
                      stringsAsFactors = FALSE, check.names = FALSE)
### positive
p_cor <- round(cor(plot_df$PI_Effector_Pcnt,
                   plot_df$Avg_Pos_Module_Score, method = "spearman", use = "complete.obs"), 2)
pv <- round(cor.test(plot_df$PI_Effector_Pcnt,
                     plot_df$Avg_Pos_Module_Score, method = "spearman", use = "complete.obs")$p.value, 2)
p <- ggplot(data = plot_df, aes(x=PI_Effector_Pcnt, y=Avg_Pos_Module_Score)) +
  geom_point(col = "#487A8F", size = 8) +
  labs(title = paste0("Spearman Correlation:", p_cor),
       subtitle = paste0("P.Val:", pv)) +
  xlab("Average_Positive_Module_Score") +
  ylab("PI_Effector_Precursor_Pcnt") +
  geom_label_repel(aes(label = Patient),
                   size = 5,
                   col = "#3B3B53",
                   segment.color = "#3B3B53") +
  # geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
  theme_classic(base_size = 40) +
  theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 30, color = "black", face = "bold"),
        plot.subtitle = element_text(hjust = 0, vjust = 0.5, size = 25, color = "black", face = "bold"),
        axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
        axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
ggsave(file = paste0("./results/New3/Manuscript/Cor_Avg_Pos_MS_vs_PI_Effect_Prcs_Pcnt.pdf"), plot = p, width = 12, height = 10, dpi = 400)

### negative
p_cor <- round(cor(plot_df$PI_Effector_Pcnt,
                   plot_df$Avg_Neg_Module_Score, method = "spearman", use = "complete.obs"), 2)
pv <- round(cor.test(plot_df$PI_Effector_Pcnt,
                     plot_df$Avg_Neg_Module_Score, method = "spearman", use = "complete.obs")$p.value, 2)
p <- ggplot(data = plot_df, aes(x=PI_Effector_Pcnt, y=Avg_Neg_Module_Score)) +
  geom_point(col = "#487A8F", size = 8) +
  labs(title = paste0("Spearman Correlation:", p_cor),
       subtitle = paste0("P.Val:", pv)) +
  xlab("Average_Negative_Module_Score") +
  ylab("PI_Effector_Precursor_Pcnt") +
  geom_label_repel(aes(label = Patient),
                   size = 5,
                   col = "#3B3B53",
                   segment.color = "#3B3B53") +
  # geom_smooth(method = lm, color="#AA4C26", se=TRUE) +
  theme_classic(base_size = 40) +
  theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 30, color = "black", face = "bold"),
        plot.subtitle = element_text(hjust = 0, vjust = 0.5, size = 25, color = "black", face = "bold"),
        axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
        axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
ggsave(file = paste0("./results/New3/Manuscript/Cor_Avg_Neg_MS_vs_PI_Effect_Prcs_Pcnt.pdf"), plot = p, width = 12, height = 10, dpi = 400)

### look at the correlations of average values
plot_df <- data.frame(Patient=names(peakcar_ug),
                      Avg_Pos_Module_Score=Avg_Pos_Module_Score[names(peakcar_ug)],
                      Avg_Neg_Module_Score=Avg_Neg_Module_Score[names(peakcar_ug)],
                      B_Cell_Recovery_Time=as.numeric(b_cell_recovery_time[names(peakcar_ug)]),
                      PeakCAR_ug=as.numeric(peakcar_ug),
                      Wk1CAR_ug=as.numeric(wk1car_ug[names(peakcar_ug)]),
                      Wk2CAR_ug=as.numeric(wk2car_ug[names(peakcar_ug)]),
                      Wk3CAR_ug=as.numeric(wk3car_ug[names(peakcar_ug)]),
                      PeakCAR_ml=as.numeric(peakcar_ml),
                      Wk1CAR_ml=as.numeric(wk1car_ml[names(peakcar_ug)]),
                      Wk2CAR_ml=as.numeric(wk2car_ml[names(peakcar_ug)]),
                      Wk3CAR_ml=as.numeric(wk3car_ml[names(peakcar_ug)]),
                      Tumor_Burden=as.numeric(Tumor_Burden[names(peakcar_ug)]),
                      stringsAsFactors = FALSE, check.names = FALSE)

pdf(file = paste0("./results/New3/Manuscript/Module_Score_Correlations/Original_Cors_pairwise_plot_20_Avg.pdf"), width = 20, height = 15)
pairs(data=plot_df,
      ~Avg_Pos_Module_Score + Avg_Neg_Module_Score + B_Cell_Recovery_Time + PeakCAR_ug + PeakCAR_ml +
        Wk1CAR_ug + Wk2CAR_ug + Wk3CAR_ug + Wk1CAR_ml + Wk2CAR_ml + Wk3CAR_ml + Tumor_Burden,
      upper.panel = panel.lm, lower.panel = panel.cor, diag.panel = panel.hist,
      font.labels = 2, pch = 19)
dev.off()

### correlation with PI (Wk2/Peak weak) and average module scores
peak_time <- rep(NA, length(unique(JCC_Seurat_Obj$px)))
names(peak_time) <- unique(JCC_Seurat_Obj$px)
peak_time["SJCAR19-01"] <- "Wk2"
peak_time["SJCAR19-02"] <- "Wk2"
peak_time["SJCAR19-03"] <- "Wk1"
peak_time["SJCAR19-04"] <- "Wk2"
peak_time["SJCAR19-05"] <- "Wk1"
peak_time["SJCAR19-06"] <- "Wk2"
peak_time["SJCAR19-07"] <- "Wk3"
peak_time["SJCAR19-08"] <- "Wk2"
peak_time["SJCAR19-09"] <- "Wk2"
peak_time["SJCAR19-10"] <- "Wk1"
peak_time["SJCAR19-11"] <- "Wk2"
peak_time["SJCAR19-13"] <- "Wk2"
peak_time["SJCAR19-14"] <- "Wk1"
peak_time["SJCAR19-15"] <- "Wk1"

wk2_precursor_pcnt <- sapply(unique(JCC_Seurat_Obj$px), function(x) {
  return(length(intersect(intersect(which(JCC_Seurat_Obj$px == x),
                                    which(JCC_Seurat_Obj$time2 == "Wk2")),
                          union(which(JCC_Seurat_Obj$AllSeuratClusters == "3"),
                                which(JCC_Seurat_Obj$AllSeuratClusters == "8")))) / length(intersect(which(JCC_Seurat_Obj$px == x),
                                                                                                     which(JCC_Seurat_Obj$time2 == "Wk2"))))
})

peaktime_precursor_pcnt <- sapply(unique(JCC_Seurat_Obj$px), function(x) {
  if(is.na(peak_time[x])) {
    return(NA)
  } else {
    return(length(intersect(intersect(which(JCC_Seurat_Obj$px == x),
                                      which(JCC_Seurat_Obj$time2 == peak_time[x])),
                            union(which(JCC_Seurat_Obj$AllSeuratClusters == "3"),
                                  which(JCC_Seurat_Obj$AllSeuratClusters == "8")))) / length(intersect(which(JCC_Seurat_Obj$px == x),
                                                                                                       which(JCC_Seurat_Obj$time2 == peak_time[x]))))
  }
})

plot_df <- data.frame(Patient=names(peakcar_ug),
                      Avg_Pos_Module_Score=Avg_Pos_Module_Score[names(peakcar_ug)],
                      Avg_Neg_Module_Score=Avg_Neg_Module_Score[names(peakcar_ug)],
                      wk2_precursor_pcnt=wk2_precursor_pcnt[names(peakcar_ug)],
                      peaktime_precursor_pcnt=peaktime_precursor_pcnt[names(peakcar_ug)],
                      B_Cell_Recovery_Time=as.numeric(b_cell_recovery_time[names(peakcar_ug)]),
                      PeakCAR_ug=as.numeric(peakcar_ug),
                      Wk1CAR_ug=as.numeric(wk1car_ug[names(peakcar_ug)]),
                      Wk2CAR_ug=as.numeric(wk2car_ug[names(peakcar_ug)]),
                      Wk3CAR_ug=as.numeric(wk3car_ug[names(peakcar_ug)]),
                      PeakCAR_ml=as.numeric(peakcar_ml),
                      Wk1CAR_ml=as.numeric(wk1car_ml[names(peakcar_ug)]),
                      Wk2CAR_ml=as.numeric(wk2car_ml[names(peakcar_ug)]),
                      Wk3CAR_ml=as.numeric(wk3car_ml[names(peakcar_ug)]),
                      Tumor_Burden=as.numeric(Tumor_Burden[names(peakcar_ug)]),
                      stringsAsFactors = FALSE, check.names = FALSE)

pdf(file = paste0("./results/New3/Manuscript/Module_Score_Correlations/Wk2_Peaktime_effector_correlations.pdf"), width = 20, height = 15)
pairs(data=plot_df,
      ~Avg_Pos_Module_Score + Avg_Neg_Module_Score + wk2_precursor_pcnt + peaktime_precursor_pcnt +
        B_Cell_Recovery_Time + PeakCAR_ug + PeakCAR_ml +
        Wk1CAR_ug + Wk2CAR_ug + Wk3CAR_ug + Wk1CAR_ml + Wk2CAR_ml + Wk3CAR_ml + Tumor_Burden,
      upper.panel = panel.lm, lower.panel = panel.cor, diag.panel = panel.hist,
      font.labels = 2, pch = 19)
dev.off()



### If we split our pre-cursor effectors into those that first appear in cluster 3 versus those that first appear in cluster 8, what do the numbers look like?
### if we split the effector precursors into those that first end up in cluster 3 vs those that first end up in cluster 8. if we have sufficient lineages, we might want to look at DEGs between those two subsets of precursors

### set clone types
only_cluster3_pi_subsister_clones <- setdiff(cluster3_pi_subsister_clones, cluster8_pi_subsister_clones)
only_cluster8_pi_subsister_clones <- setdiff(cluster8_pi_subsister_clones, cluster3_pi_subsister_clones)
shared_cluster38_subsister_clones <- intersect(cluster3_pi_subsister_clones, cluster8_pi_subsister_clones)
names(shared_cluster38_subsister_clones) <- shared_cluster38_subsister_clones

for(clone in names(shared_cluster38_subsister_clones)) {
  shared_cluster38_subsister_clones[clone] <- NA
  
  cluster_list <- JCC_Seurat_Obj$AllSeuratClusters[which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta == clone)]
  time_list <- JCC_Seurat_Obj$time2[which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta == clone)]
  time_list <- factor(time_list, levels = unique(JCC_Seurat_Obj$time2))
  
  cluster_time_df <- data.frame(Cluster=cluster_list,
                                Time=time_list,
                                stringsAsFactors = FALSE, check.names = FALSE)
  cluster_time_df <- cluster_time_df[order(cluster_time_df$Time),]
  
  
  after_gmp_earliest_time <- NULL
  for(i in 1:nrow(cluster_time_df)) {
    if((cluster_time_df$Time[i] != "GMP") && ((as.character(cluster_time_df$Cluster[i]) == "3") || (as.character(cluster_time_df$Cluster[i]) == "8"))) {
      after_gmp_earliest_time <- as.character(time_list[i])
      break
    }
  }
  
  after_gmp_earliest_time_idx <- which(as.character(time_list) == after_gmp_earliest_time)
  
  if(length(after_gmp_earliest_time_idx) == 1) {
    shared_cluster38_subsister_clones[clone] <- as.character(cluster_list[names(time_list)[after_gmp_earliest_time_idx]])
  }
}

JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First <- NA
JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% only_cluster3_pi_subsister_clones))] <- "GMP_Precursor_Effector_Cluster3_First"
JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% only_cluster8_pi_subsister_clones))] <- "GMP_Precursor_Effector_Cluster8_First"
JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% names(shared_cluster38_subsister_clones)[which(shared_cluster38_subsister_clones == "3")]))] <- "GMP_Precursor_Effector_Cluster3_First"
JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% names(shared_cluster38_subsister_clones)[which(shared_cluster38_subsister_clones == "8")]))] <- "GMP_Precursor_Effector_Cluster8_First"
JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% names(shared_cluster38_subsister_clones)[which(is.na(shared_cluster38_subsister_clones))]))] <- "GMP_Precursor_Effector_Both_First"

### print out some numbers
print(length(which(JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First == "GMP_Precursor_Effector_Cluster3_First")))
print(length(which(JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First == "GMP_Precursor_Effector_Cluster8_First")))
print(length(which(JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First == "GMP_Precursor_Effector_Both_First")))

print(length(only_cluster3_pi_subsister_clones) + length(which(shared_cluster38_subsister_clones == "3")))
print(length(only_cluster8_pi_subsister_clones) + length(which(shared_cluster38_subsister_clones == "8")))
print(length(which(is.na(shared_cluster38_subsister_clones))))

sapply(unique(JCC_Seurat_Obj$px), function(x) {length(intersect(which(JCC_Seurat_Obj$px == x),
                                                                which(JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First == "GMP_Precursor_Effector_Cluster3_First")))})
sapply(unique(JCC_Seurat_Obj$px), function(x) {length(intersect(which(JCC_Seurat_Obj$px == x),
                                                                which(JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First == "GMP_Precursor_Effector_Cluster8_First")))})
sapply(unique(JCC_Seurat_Obj$px), function(x) {length(intersect(which(JCC_Seurat_Obj$px == x),
                                                                which(JCC_Seurat_Obj$GMP_Precursor_Effector_Cluster38_First == "GMP_Precursor_Effector_Both_First")))})


### DE genes between GMP precursor cluster3 first vs cluster8 first
JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                           cells = rownames(JCC_Seurat_Obj@meta.data),
                           value = JCC_Seurat_Obj@meta.data$GMP_Precursor_Effector_Cluster38_First)
de_result <- FindMarkers(JCC_Seurat_Obj,
                         ident.1 = "GMP_Precursor_Effector_Cluster3_First",
                         ident.2 = "GMP_Precursor_Effector_Cluster8_First",
                         min.pct = 0,
                         logfc.threshold = 0,
                         test.use = "wilcox")

### write out the DE result
write.xlsx2(data.frame(Gene=rownames(de_result),
                       de_result,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = paste0("./results/New3/Manuscript/GMP_Precursor_Effector_Cluster3_First_vs_Cluster8_First.xlsx"),
            sheetName = "GMP_Precursor_Effector_Cluster3_First_vs_Cluster8_First", row.names = FALSE)










