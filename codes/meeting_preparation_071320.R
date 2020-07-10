### This is just to create some useful plots for the immuno-core meeting presentation

Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj"

outputDir="./etc/071320/"
dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)

### load libraries
if(!require(Seurat, quietly = TRUE)) {
  install.packages("Seurat")
  require(Seurat, quietly = TRUE)
}
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}
if(!require(ggplot2, quietly = TRUE)) {
  install.packages("ggplot2")
  require(ggplot2, quietly = TRUE)
}

### load the Seurat object and save the object name
tmp_env <- new.env()
load(Seurat_RObj_path, tmp_env)
obj_name <- ls(tmp_env)
assign("Seurat_Obj", get(obj_name, envir = tmp_env))
rm(tmp_env)
gc()

### rownames in the meta.data should be in the same order as colnames in the counts
Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]

### set unique patients and time points
px <- unique(Seurat_Obj@meta.data$Px)
tp <- levels(Seurat_Obj@meta.data$TimeF)

### make tables for all the cells (GEX & TCR) - the number of cells
all_gex_mat <- matrix(0, length(px), length(tp))
colnames(all_gex_mat) <- tp
rownames(all_gex_mat) <- px
all_tcr_mat <- matrix(0, length(px), length(tp))
colnames(all_tcr_mat) <- tp
rownames(all_tcr_mat) <- px

### fill out the tables
for(p in px) {
  for(t in tp) {
    all_gex_mat[p, t] <- length(intersect(which(Seurat_Obj@meta.data$Px == p),
                                          which(Seurat_Obj@meta.data$Time == t)))
  }
}
for(p in px) {
  for(t in tp) {
    all_tcr_mat[p, t] <- length(intersect(intersect(which(Seurat_Obj@meta.data$Px == p),
                                                    which(Seurat_Obj@meta.data$Time == t)),
                                          which(!is.na(Seurat_Obj@meta.data$cdr3_aa))))
  }
}

### make tables for CAR+ cells (GEX & TCR) - the number of cells
car_gex_mat <- matrix(0, length(px), length(tp))
colnames(car_gex_mat) <- tp
rownames(car_gex_mat) <- px
car_tcr_mat <- matrix(0, length(px), length(tp))
colnames(car_tcr_mat) <- tp
rownames(car_tcr_mat) <- px

### fill out the tables
for(p in px) {
  for(t in tp) {
    car_gex_mat[p, t] <- length(intersect(intersect(which(Seurat_Obj@meta.data$Px == p),
                                                    which(Seurat_Obj@meta.data$Time == t)),
                                          which(Seurat_Obj@meta.data$CAR == "CARpos")))
  }
}
for(p in px) {
  for(t in tp) {
    car_tcr_mat[p, t] <- length(intersect(intersect(intersect(which(Seurat_Obj@meta.data$Px == p),
                                                              which(Seurat_Obj@meta.data$Time == t)),
                                                    which(!is.na(Seurat_Obj@meta.data$cdr3_aa))),
                                          which(Seurat_Obj@meta.data$CAR == "CARpos")))
  }
}

### write out the tables
write.xlsx2(data.frame(all_gex_mat, check.names = FALSE), file = paste0(outputDir, "all_gex_table.xlsx"))
write.xlsx2(data.frame(all_tcr_mat, check.names = FALSE), file = paste0(outputDir, "all_tcr_table.xlsx"))
write.xlsx2(data.frame(car_gex_mat, check.names = FALSE), file = paste0(outputDir, "car_gex_table.xlsx"))
write.xlsx2(data.frame(car_tcr_mat, check.names = FALSE), file = paste0(outputDir, "car_tcr_table.xlsx"))

### make a heatmap for CAR+ GEX & TCR






### how many cells are in the TCR but not in the GEX?






