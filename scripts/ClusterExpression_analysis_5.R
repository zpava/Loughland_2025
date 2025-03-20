
##Read in data
immune.combined <- readRDS( "data/IBSM_monoSCT_integbyDayDonor_220823.rds")

###FIGURE 1a.
##Clustering

#Checking new clustering

##Changing clusters names
immune.combined@meta.data$celltype2 <- NA
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-COQ7_11213"))] <- "cMonocyte_COQ7"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-IL_02614"))] <- "cMonocyte_IL1B"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-ISG15_10"))] <- "cMonocyte_ISG15"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-S100A8_79"))] <- "cMonocyte_S100A8"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "DCs1_15"))] <- "DCs1"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "DCs2_511"))] <- "DCs2"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "intMonocyte_3"))] <- "int_Monocyte"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "ncMonocyte-CD16_48"))] <- "nc_Monocyte"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "pDCs_16"))] <- "pDCs"

##Fixing levels of some variables
celltype_levels <- c("cMonocyte_IL1B","cMonocyte_S100A8",
                     "cMonocyte_COQ7","cMonocyte_ISG15",
                     "int_Monocyte", "nc_Monocyte", 
                     "DCs1", "DCs2", "pDCs")

immune.combined@meta.data$celltype2 <- factor(x = immune.combined@meta.data$celltype2, levels = celltype_levels)

Idents(immune.combined) <- "celltype2"
Idents(immune.combined) <- "integrated_snn_res.0.3"

Tfh_colors_list <- c("#9e0142", "#f4a582", "#5e4fa2", "#66c2a5","#74add1",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")



DefaultAssay(immune.combined) <- "SCT"

immune.combined <- SCTransform(immune.combined, vars.to.regress = "percent.mt", verbose = FALSE)
immune.combined <- RunPCA(immune.combined, features = VariableFeatures(object = immune.combined), npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
ElbowPlot(immune.combined, ndims = 30)## we'll go with 30.


##How classical markers are driven by donors RETN is one of them!
Idents(immune.combined) <- "day_donor"
pdf("output_for_paper/clusters_expression/Dotplot_monomarkers_per_daydonor_Res0.3_Feb25.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", 
        features=c("CD14", "FCGR3A", "CD3E","CD4", "CCL3", "RETN", "CD8A", 
                   "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB","CCL5", 
                   "FOXP3", "CTLA4", "HLA-DRA", "PPBP", "CLEC4C", "CLIC3","IL1B", "CD1C"),
        cols=c("white", "#9e0142")) +
        theme_bw(base_size=14) + 
        theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1)) + 
  theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
dev.off()


##Creating cluster day variable
immune.combined@meta.data$cluster_day <- paste0(immune.combined@meta.data$celltype2, "_", immune.combined@meta.data$day)

Idents(object = immune.combined) <- immune.combined@meta.data$'cluster_day'
pdf("output_for_paper/clusters_expression/Dotplot_monomarkers_per_cluster_day_Feb25.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", 
        features=c("CD14", "FCGR3A", "CD3E","CD4", "CCL3", "RETN", "CD8A", 
                   "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB","CCL5", 
                   "FOXP3", "CTLA4", "HLA-DRA", "PPBP", "CLEC4C", "CLIC3","IL1B", "CD1C"),
        cols=c("#74add1","#9e0142",  "#b2abd2"), split.by = "day") +
  theme_bw(base_size=14) + 
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1)) + 
  theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())+ coord_flip()
dev.off()

Idents(object = immune.combined) <- immune.combined@meta.data$'celltype2'
pdf("output_for_paper/clusters_expression/Dotplot_monomarkers_per_cluster_day_seq_Feb25.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", 
        features=c("CD14","IL1B","LRRK2","COQ7", "ISG15", "FCGR3A", "CDKN1C", "CX3CR1",
                   "CLEC9A","CADM1", "XCR1",
                   "CD1C","FCER1A","CLEC10A",
                   "CLIC3","NRP1","GZMB"),
        cols=c("#74add1","#9e0142",  "#b2abd2"), split.by = "day") +
  theme_bw(base_size=14) + 
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1)) + 
  theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())+ coord_flip()
dev.off()
###Differential expression analysis Clusters were identified via cannonical expression of the previous markers
##Diff expression among classical monocytes
Idents(object = immune.combined) <- immune.combined@meta.data$'celltype2'
cMonocyte_IL1B <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_IL1B", ident.2 = c("cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_IL1B$cluster <- "cMonocyte_IL1B"
cMonocyte_IL1B$Gene <- row.names(cMonocyte_IL1B)
head(cMonocyte_IL1B)

cMonocyte_S100A8 <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_S100A8", ident.2 = c("cMonocyte_IL1B", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_S100A8$cluster <- "cMonocyte_S100A8"
cMonocyte_S100A8$Gene <- row.names(cMonocyte_S100A8)
head(cMonocyte_S100A8)

cMonocyte_COQ7 <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_COQ7", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_COQ7$cluster <- "cMonocyte_COQ7"
cMonocyte_COQ7$Gene <- row.names(cMonocyte_COQ7)
head(cMonocyte_COQ7)

cMonocyte_ISG15 <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_ISG15", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_ISG15$cluster <- "cMonocyte_ISG15"
cMonocyte_ISG15$Gene <- row.names(cMonocyte_ISG15)
head(cMonocyte_ISG15)

int_Monocyte <-  FindMarkers(immune.combined, ident.1 = "int_Monocyte", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
int_Monocyte$cluster <- "int_Monocyte"
int_Monocyte$Gene <- row.names(int_Monocyte)
head(int_Monocyte)

nc_Monocyte <-  FindMarkers(immune.combined, ident.1 = "nc_Monocyte", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
nc_Monocyte$cluster <- "nc_Monocyte"
nc_Monocyte$Gene <- row.names(nc_Monocyte)
head(nc_Monocyte)

combined
##Already known identity by classical markers. 
cluster2Con.markers <- FindConservedMarkers(immune.combined, ident.1 = 2, ident.2 = c(0,1,3,4,5,6,7,8,9), min.pct = 0.25, grouping.var = "day")
cluster4Con.markers <- FindConservedMarkers(immune.combined, ident.1 = 4, ident.2 = c(0,1,2,3,5,6,7,8,9), min.pct = 0.25, grouping.var = "day")
cluster7Con.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, ident.2 = c(0,1,2,3,4,5,6,8,9), min.pct = 0.25, grouping.var = "day")
cluster9Con.markers <- FindConservedMarkers(immune.combined, ident.1 = 9, ident.2 = c(0,1,2,3,4,5,6,7,8), min.pct = 0.25, grouping.var = "day")


cluster0Con.markers$cluster <- "CD14-IL1B-Mono"
cluster1Con.markers$cluster <- "CD14-HNRNPU-Mono"
cluster2Con.markers$cluster <-  "CD16 Mono"
cluster3Con.markers$cluster <- "CD14-RETN-Mono"
cluster4Con.markers$cluster <- "DCs"
cluster5Con.markers$cluster <- "CD14-RPS26 Mono"
cluster6Con.markers$cluster <- "CD14-FOLR3 Mono"
cluster7Con.markers$cluster <- "DCs1"
cluster9Con.markers$cluster <- "pDC"


C0_10.markers <- rbind(cluster0Con.markers,cluster1Con.markers,cluster2Con.markers,cluster3Con.markers,cluster4Con.markers,cluster5Con.markers,cluster6Con.markers,cluster7Con.markers,cluster9Con.markers)
rownames(C0_10.markers)
C0_10.markers <-setDT(C0_10.markers, keep.rownames = "genes")
head(C0_10.markers)

C0_10.markers_p <- C0_10.markers %>% 
  group_by(cluster) %>%
  arrange(day0_p_val_adj) %>%
  slice_head(n=10)

C0_10.markers_p <- C0_10.markers_p [ , c(1,19, 8, 11, 9, 10, 7, 13, 16, 14, 15, 12, 3, 6, 4, 5, 2, 17, 18)]

write.csv(C0_10.markers_p, "/Users/zuleip/Documents/Seurat_analysis_test/IBSM_40/data/DF_conMarkers_int_0.3_040822.csv")

##There is a test for discriminatory markers AUC, which quantifies how accurate a marker is.
cluster0.markers <- FindMarkers(immune.combined, ident.1 = 0, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)

cluster1.markers <- FindMarkers(immune.combined, ident.1 = 1, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers, n = 5)

cluster2.markers <- FindMarkers(immune.combined, ident.1 = 2, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster2.markers, n = 5)

cluster3.markers <- FindMarkers(immune.combined, ident.1 = 3, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster3.markers, n = 5)

cluster4.markers <- FindMarkers(immune.combined, ident.1 = 4, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster4.markers, n = 5)

cluster5.markers <- FindMarkers(immune.combined, ident.1 = 5, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster5.markers, n = 5)

cluster6.markers <- FindMarkers(immune.combined, ident.1 = 6, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster6.markers, n = 5)

cluster7.markers <- FindMarkers(immune.combined, ident.1 = 7, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster7.markers, n = 5)

cluster9.markers <- FindMarkers(immune.combined, ident.1 = 9, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster9.markers, n = 5)


### Adding cluster name for merging

cluster0.markers$cluster <- "CD14-IL1B-Mono"
cluster1.markers$cluster <- "CD14-HNRNPU-Mono"
cluster2.markers$cluster <-  "CD16 Mono"
cluster3.markers$cluster <- "CD14-RETN-Mono"
cluster4.markers$cluster <- "DCs"
cluster5.markers$cluster <- "CD14-RPS26 Mono"
cluster6.markers$cluster <- "CD14-FOLR3 Mono"
cluster7.markers$cluster <- "DCs1"
cluster9.markers$cluster <- "pDC"


AUC_10.markers <- rbind(cluster0.markers,cluster1.markers,cluster2.markers,cluster3.markers,cluster4.markers,cluster5.markers,cluster6.markers,cluster7.markers,cluster9.markers)
rownames(AUC_10.markers)
AUC_10.markers <-setDT(AUC_10.markers, keep.rownames = "genes")
head(AUC_10.markers)

AUC_10.markers_p <- AUC_10.markers %>% 
  group_by(cluster) %>%
  arrange(desc(power)) %>%
  slice_head(n=10)

AUC_10.markers_p <- AUC_10.markers_p [ , c(1,8,4, 2:3, 5:7)]

write.csv(AUC_10.markers_p, "/Users/zuleip/Documents/Seurat_analysis_test/IBSM_40/data/DF_AUCMarkers_int_0.3_040822.csv")
