##Plotting with seurat object
##Author Zuly P
##Date 14 Nov 2023

library(glmGamPoi)
library(tidyverse)
library(Seurat) 
library(patchwork)
library(sctransform)
library(limma)
library(data.table)
library(ape)
library(SC3)
library(clustree)
library(MAST)
library(knitr)

##Read in data
immune.combined <- readRDS( "data/IBSM_monoSCT_integbyDayDonor_220823.rds")

###FIGURE 1a.
##Clustering

#Checking new clustering
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
#immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)

##Changing clusters names
immune.combined@meta.data$celltype2 <- NA
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-COQ7_11213"))] <- "cMonocyte-COQ7"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-IL_02614"))] <- "cMonocyte-IL1B"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-ISG15_10"))] <- "cMonocyte-ISG15"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-S100A8_79"))] <- "cMonocyte-S100A8"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "DCs1_15"))] <- "DCs1"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "DCs2_511"))] <- "DCs2"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "intMonocyte_3"))] <- "int-Monocyte"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "ncMonocyte-CD16_48"))] <- "nc-Monocyte"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "pDCs_16"))] <- "pDCs"

##Fixing levels of some variables
celltype_levels <- c("cMonocyte-IL1B","cMonocyte-S100A8",
                     "cMonocyte-COQ7","cMonocyte-ISG15",
                     "int-Monocyte", "nc-Monocyte", 
                     "DCs1", "DCs2", "pDCs")

immune.combined@meta.data$celltype2 <- factor(x = immune.combined@meta.data$celltype2, levels = celltype_levels)

Idents(immune.combined) <- "celltype2"
Idents(immune.combined) <- "integrated_snn_res.0.3"

Tfh_colors_list <- c("#9e0142", "#f4a582", "#5e4fa2", "#66c2a5","#74add1",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")

###FIGURE 1a Dimplot with final clusters####.

pdf("output_for_paper/Fig1a_Dimplot_celltype_Res0.3_July24.pdf", h=15, w=10)
DimPlot(immune.combined, label = FALSE, pt.size = 2, group.by = "celltype2", cols = Tfh_colors_list) + NoLegend() 
dev.off()

pdf("output_for_paper/Dimplot_celltype_Res0.3_July24withlabels.pdf", h=15, w=10)
DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.3",  cols = Tfh_colors_list) + NoLegend() 
dev.off()

pdf("output_for_paper/Dimplot_celltype_Res0.3_July24withlabels.pdf", h=15, w=10)
DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "celltype2",cols = Tfh_colors_list) + NoLegend() 
dev.off()

pdf("output_for_paper/PCA_celltype_Res0.3_July24withlabels.pdf", h=15, w=10)
DimPlot(immune.combined, reduction= 'pca',label = TRUE, pt.size = 2, group.by = "celltype2",cols = Tfh_colors_list) + NoLegend() 
dev.off()

####Just classical monocytes

# Subset based on PCA_2 using the subset function
#immune.combined_subset  <- subset(immune.combined, cells = rownames(immune.combined@reductions$pca@cell.embeddings[immune.combined@reductions$pca@cell.embeddings[, 2] < 8, ]))

# Subset based on celltypes using the subset function
cMonocytes <- c("int-Monocyte", "nc-Monocyte", "cMonocyte-IL1B", "cMonocyte-S100A8", "cMonocyte-COQ7")
immune.combined_subset2  <- subset(immune.combined, idents = cMonocytes)

# Dimplot with only classical monocytes
pdf("output_for_paper/PCA_cMonocytes_Res0.3_July24withlabels.pdf", h=15, w=10)
DimPlot(immune.combined_subset2, reduction= 'pca', label = TRUE, pt.size = 2, group.by = "celltype2",cols = Tfh_colors_list) + NoLegend() 
dev.off()

# Dimplot with only classical monocytes
pdf("output_for_paper/umap_cMonocytes_Res0.3_July24withlabels.pdf", h=15, w=10)
DimPlot(immune.combined_subset2, reduction= 'umap', label = TRUE, pt.size = 2, group.by = "celltype2",cols = Tfh_colors_list) + NoLegend() 
dev.off()


###FIGURE 1b Dotplot_monomarkers_Res0.8_March24.####

##Cluster id
DefaultAssay(immune.combined) <- "SCT"

#immune.combined <- SCTransform(immune.combined, vars.to.regress = "percent.mt", verbose = FALSE)
immune.combined <- RunPCA(immune.combined, features = VariableFeatures(object = immune.combined), npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
Idents(immune.combined) <- "celltype2"
table(immune.combined$celltype2)

pdf("output_for_paper/Fig1b_Dotplot_monomarkers_Res0.8_July24.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", features=c("CD14","IL1B","S100A8","COQ7", "ISG15", "FCGR3A", "CDKN1C", "CX3CR1",
                                                   "CLEC9A","CADM1", "XCR1",
                                                   "CD1C","FCER1A","CD33", 
                                                   "CLEC4C","NRP1","IL3RA"),
        cols=c("#fee090","#9e0142"))+theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
dev.off()

pdf("output_for_paper/Fig1b_Dotplot_monomarkers_Res0.8_July24B.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", features=c("CD14","IL1B","LRRK2","COQ7", "ISG15", "FCGR3A", "CDKN1C", "CX3CR1",
                                                   "CLEC9A","CADM1", "XCR1",
                                                   "CD1C","FCER1A","CLEC10A",
                                                   "CLIC3","NRP1","GZMB"),
        cols=c("#fee090","#9e0142"))+theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
dev.off()

####Figure 1e Barplot cellnumbers per donor per different days####

Fig1EMeta <- immune.combined@meta.data
colnames(Fig1EMeta)
Fig1EMeta_Freq <- as.data.frame(table(Fig1EMeta$donor_id, Fig1EMeta$celltype2, Fig1EMeta$day))
#Renaming columns
colnames(Fig1EMeta_Freq)[1] ="Donor"
colnames(Fig1EMeta_Freq)[2] ="Celltype"
colnames(Fig1EMeta_Freq)[3] ="Day"
colnames(Fig1EMeta_Freq)[4] ="Count"
head(Fig1EMeta_Freq)

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "int-Monocyte", "nc-Monocyte", "DCs1", "DCs2", "pDCs")
##Creating day factors
Fig1EMeta_Freq$Celltype <- factor(Fig1EMeta_Freq$Celltype, levels = celltype_fact)

##Creating day factors
day_fact<- c("day0","day8", "day36")
##Creating day factors
Fig1EMeta_Freq$Day <- factor(Fig1EMeta_Freq$Day, levels = day_fact)

##Creating donors factors
donor_fact<- c("donor0","donor1", "donor2", "donor3")
##Creating day factors
Fig1EMeta_Freq$Donor <- factor(Fig1EMeta_Freq$Donor, levels = donor_fact)

# Calculate the proportions
Fig1EMeta_Prop <- Fig1EMeta_Freq %>%
  group_by(Donor, Day) %>%
  mutate(Proportion = Count / sum(Count)*100) %>%
  ungroup()

#Figure1E-Barplot with cell numbers per donor
pop_clus <- ggplot(data = Fig1EMeta_Prop, aes(x = Day, y = Proportion, fill = Celltype)) +
  geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = Tfh_colors_list) +
  scale_y_continuous(limits = c(0, 100)) +
  facet_grid(~ Donor, scales = "free") +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(size = 12, color = "black")
  ) 

pop_clus

ggsave("output_for_paper/Fig1E_Barplot_NoCells_per_DonorDay_16724.pdf", plot = pop_clus,
      height = 8, width = 16, dpi = 300)

#SupFigure Barplot with cell numbers per day
SupFig1EMeta_Freq <- as.data.frame(table(Fig1EMeta$celltype2, Fig1EMeta$day))
#Renaming columns
colnames(SupFig1EMeta_Freq)[1] ="Celltype"
colnames(SupFig1EMeta_Freq)[2] ="Day"
colnames(SupFig1EMeta_Freq)[3] ="Count"
head(SupFig1EMeta_Freq)

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "int-Monocyte", "nc-Monocyte", "DCs1", "DCs2", "pDCs")
##Creating day factors
SupFig1EMeta_Freq$Celltype <- factor(SupFig1EMeta_Freq$Celltype, levels = celltype_fact)

##Creating day factors
day_fact<- c("day0","day8", "day36")
##Creating day factors
SupFig1EMeta_Freq$Day <- factor(SupFig1EMeta_Freq$Day, levels = day_fact)
# Calculate the proportions
SupFig1EMeta_Freq <- SupFig1EMeta_Freq %>%
  group_by(Day) %>%
  mutate(Proportion = Count / sum(Count)*100) %>%
  ungroup()

pop_clus2 <- ggplot(data = SupFig1EMeta_Freq, aes(x = Day, y = Proportion, fill = Celltype)) +
  geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = Tfh_colors_list) +
  scale_y_continuous(limits = c(0, 100)) +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(size = 12, color = "black")
  ) 

pop_clus2

ggsave("output_for_paper/SupFigure_Barplot_ClusterProp_per_Day_290125.pdf", plot = pop_clus,
       height = 8, width = 16, dpi = 300)


####NOT USED###
# Idents(immune.combined) <- "subtype"
# 
# pdf("graphs/paper/Dimplot_subtype_Res0.8_Nov23.pdf", h=15, w=10)
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "celltype") + NoLegend() 
# dev.off()

# ##clustree
# 
# # ClusTree of reference editors. Each level of the tree (from top to bottom) corresponds to the k used, and each node on that level corresponds to a cluster of that k. The edges (arrows) of the tree represent the editors that "move" from one cluster in level k to another cluster in level k+1. The legend of the graph displays 4 scales, from top to bottom: (1) the transparency level of arrows (in_prop) shows the proportion editors from one group that end up in another group. (2) The arrow colour (count) shows the number of editors that "move" from one cluster to another, (3) the node size is proportional to the number
# pdf("graphs/paper/clustree_SCTintDayDonor_Nov23.pdf", h=20, w=10)
# clustree(immune.combined, prefix = "integrated_snn_res.", layout = "sugiyama")
# dev.off()
# 
# pdf("graphs/paper/clustreeStability_SCTintDayDonor_Nov23.pdf", h=20, w=10)
# clustree(immune.combined, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
# dev.off()
# ##Based on the tree we will explore resolution 0.8 (17 clusters) 
# 
# 
# ##Cluster id
# DefaultAssay(immune.combined) <- "SCT"
# 
# #immune.combined <- SCTransform(immune.combined, vars.to.regress = "percent.mt", verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, features = VariableFeatures(object = immune.combined), npcs = 30, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
# 
# table(immune.combined$subtype)

# ##Fixing levels of some variables
# subtype_levels <- c("cMonocyte-IL1B_0","cMonocyte-COQ7_1",
#                     "cMonocyte-CDKN1A_2","intMonocyte_3",
#                     "ncMonocyte-CD16_4", "DCs2_5", 
#                     "cMonocyte-DUSP1_6", "cMonocyte-MALAT_7",
#                     "ncMonocyte-CD16_8", "cMonocyte-S100A8_9",
#                     "cMonocyte-ISG15_10","DCs2-RPS18_11",
#                     "cMonocyte-_12","cMonocyte-HLA_13",
#                     "cMonocyte-IL1RN_14", "DCs1_15", "pDCs_16")
# immune.combined@meta.data$subtype <- factor(x = immune.combined@meta.data$subtype, levels = subtype_levels)
# 
# 
# 
# pdf("graphs/paper/Dotplot_monomarkers_Res0.8_Nov23.pdf", h=15, w=10)
# DotPlot(immune.combined, assay = "SCT", features=c("CD14", "FCGR3A", "CD3E","CD4", "CCL3", "RETN", "CD8A", "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB","CCL5", 
#                                                    "FOXP3", "CTLA4", "HLA-DRA", "PPBP", "CLEC4C", "CLIC3","IL1B", "CD1C"),
#         cols=c("white", "red3"))+theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ()) + scale_fill_manual(values=Tfh_colors_list)
# dev.off()
# 
# pdf("graphs/paper/Dotplot_monomarkers2_Res0.8_Nov23.pdf", h=15, w=10)
# DotPlot(immune.combined, assay = "SCT", features=c("CD14", "FCGR3A", "CX3CR1", "CDKN1C", "CD1C", "FCER1A", "CLEC4C", "NRP1","XCR1", "CLEC9A","CADM1","IL3RA"),
#         cols=c("white", "red3"))+theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
# dev.off()

# ##subtype plots
# subtypes_freq <- as.data.frame(table(immune.combined@meta.data$subtype, immune.combined@meta.data$day))
# colnames(subtypes_freq)[1] ="Subtype"
# colnames(subtypes_freq)[2] ="Day"
# colnames(subtypes_freq)[3] ="Frequency"
# kable(subtypes_freq)
# table(immune.combined@meta.data$day)
# # day0  day8 day36 
# # 6640  7032  1963 
# subtype_levels2 <- c("cMonocyte-IL1B_0","cMonocyte-CDKN1A_2",
#                      "cMonocyte-DUSP1_6","ncMonocyte-CD16_8",
#                      "cMonocyte-COQ7_1","intMonocyte_3",
#                     "ncMonocyte-CD16_4","DCs2-RPS18_11", 
#                     "cMonocyte-HLA_13","cMonocyte-IL1RN_14",
#                     "DCs2_5", "cMonocyte-MALAT_7",
#                      "cMonocyte-S100A8_9", "cMonocyte-ISG15_10",
#                     "cMonocyte-_12", "pDCs_16", "DCs1_15")
# 
# subtypes_freq$Subtype <- factor(subtypes_freq$Subtype, levels = subtype_levels2)
# 
# subtypes_freq = subtypes_freq %>% 
#   group_by(Day) %>%
#   mutate(total= sum(Frequency)) %>%
#   mutate(`Proportion of clusters` = round(Frequency/total*100, 1)) %>%
#   mutate(high_day=case_when(Subtype=="cMonocyte-IL1B_0"~"Day0",
#                             Subtype=="cMonocyte-CDKN1A_2"~"Day0",
#                             Subtype=="cMonocyte-DUSP1_6"~"Day0",
#                             Subtype=="ncMonocyte-CD16_8"~"Day0",
#                             Subtype=="cMonocyte-COQ7_1"~"Day8",
#                             Subtype=="intMonocyte_3"~"Day8",
#                             Subtype=="ncMonocyte-CD16_4"~"Day8",
#                             Subtype=="DCs2-RPS18_11"~"Day8", 
#                             Subtype=="cMonocyte-HLA_13"~"Day8",
#                             Subtype=="cMonocyte-IL1RN_14"~"Day8",
#                             Subtype=="DCs2_5"~"Day36", 
#                             Subtype=="cMonocyte-MALAT_7"~"Day36",
#                             Subtype=="cMonocyte-S100A8_9"~"Day36", 
#                             Subtype=="cMonocyte-ISG15_10"~"Day36",
#                             Subtype=="cMonocyte-_12"~"Day36", 
#                             Subtype=="pDCs_16"~"Day36", 
#                             Subtype=="DCs1_15"~"Day36"))
# 
# subtypes_freq$high_day <- factor(subtypes_freq$high_day, levels = c("Day0", "Day8", "Day36"))
# 
# head(subtypes_freq)
# 
# 
# #using adj_frequency
# pop_clus <- ggplot(data=subtypes_freq, aes(x=Subtype, y=`Proportion of clusters`, fill=Day)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   scale_y_continuous(limits = c(0,20))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#                                    panel.grid.major.y = element_blank(),
#                                    panel.grid.minor.y = element_blank()) + 
#   facet_grid(~high_day, scales = "free")
# 
# pop_clus
# 
# ggsave("graphs/paper/Cluster_proportion_per_Day.pdf", plot = pop_clus,
#        height = 8, width = 16, dpi = 300)


####SUPPLEMENTARY PLOTS####

#### Sup_Fig_1A_vlnplot_prepostQC####
renv::init()

Tfh_colors_list <- c("#9e0142", "#f4a582", "#5e4fa2", "#66c2a5","#74add1",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")

ibsm_mono <- readRDS( "data/IBSM_mono_prefiltered_210823.rds")

##Changing variable names
colnames(metadata)
colnames(metadata)[2]<- "nUMI"
colnames(metadata)[3]<- "nGene"
colnames(metadata)[8]<- "mitoRatio"

##A. pre and post filter quality variables

colnames(ibsm_mono@meta.data)
colnames(ibsm_mono@meta.data)[2]<- "nUMI"
colnames(ibsm_mono@meta.data)[3]<- "nGene"
colnames(ibsm_mono@meta.data)[8]<- "mitoRatio"

ibsm_mono_f <- subset(ibsm_mono, subset = nGene > 200 & nUMI < 7000 & mitoRatio < 20 & log10GenesPerUMI > 0.8) 

a <- VlnPlot(ibsm_mono, features = c("nUMI", "nGene", "mitoRatio"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id", pt.size = -1, ncol = 4, flip = TRUE)

b <- VlnPlot(ibsm_mono_f, features = c("nUMI", "nGene", "mitoRatio"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id",  pt.size = -1, ncol = 4, flip = TRUE) 

(a/b)
ggsave("output_for_paper/Sup_Fig_1A_vlnplot_QCprepostQC_v2a.pdf", plot =(a/b),
       height = 16, width = 12, dpi = 300)



##vertical option f<-(a-b)  + plot_layout(guides='collect')

a <- VlnPlot(ibsm_mono, features = c("nUMI"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id", pt.size = -1, flip = TRUE)

b <- VlnPlot(ibsm_mono_f, features = c("nUMI"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id",  pt.size = -1, flip = TRUE) 

c <- VlnPlot(ibsm_mono, features = c("nGene"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id", pt.size = -1,  flip = TRUE)

d <- VlnPlot(ibsm_mono_f, features = c("nGene"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id",  pt.size = -1,  flip = TRUE) 

e <- VlnPlot(ibsm_mono, features = c("mitoRatio"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id", pt.size = -1,  flip = TRUE)

f <- VlnPlot(ibsm_mono_f, features = c("mitoRatio"), cols = Tfh_colors_list, group.by = "day", split.by = "donor_id",  pt.size = -1,flip = TRUE) 

final <- (a+b+c+d+e+f) + plot_layout(guides='collect', ncol = 2, widths = c(2,2))

ggsave("output_for_paper/Sup_Fig_1A_vlnplot_QCprepostQC_v2b.pdf", plot = final,
       height = 20, width = 16, dpi = 300)

####Supp_Fig_1B, barplot with cell numbers/donor and mito Ratio ####
head(ibsm_mono@meta.data)
metadata <- ibsm_mono@meta.data
donor_day_lev <- c("day0_donor0", "day8_donor0", "day36_donor0", "day0_donor1", "day8_donor1", "day36_donor1",
                   "day0_donor2",  "day8_donor2","day36_donor2", "day0_donor3", "day8_donor3", "day36_donor3")
metadata$day_donor <- factor(metadata$day_donor, levels = donor_day_lev)
metadata %>%
  ggplot(aes(x=day_donor, fill=day_donor)) +
  geom_bar() +
  scale_fill_manual(values=Tfh_colors_list)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none")+
  ggtitle("Number of cells per donor and day") -> aa

ibsm_mono_f@meta.data$day_donor <- factor(ibsm_mono_f@meta.data$day_donor, levels = donor_day_lev)

Idents(object = ibsm_mono_f) <- "donor_id"
bb <- VlnPlot(ibsm_mono_f, features = c("mitoRatio"), fill.by = "ID", group.by = "day_donor", pt.size = -1, ncol = ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + scale_fill_manual(values=Tfh_colors_list)+
  theme(legend.position = "none")

ff <- bb/aa
 
ggsave("output_for_paper/Sup_Fig_1B_barplot_ncell_mito.pdf", plot = (aa/bb),
       height = 10, width = 6, dpi = 300)


####Figure 1C, dimplot before and after integration ####

##Read in data
ibsm.mono.F <- readRDS(file = "data/IBSM_mono_filtered_210823.rds")
##Before integration
#### NORMALISATION #########
#adjust the limit for allowable object sizes within R 
options(future.globals.maxSize = 4000 * 1024^2)

ibsm.mono.F <- SCTransform(ibsm.mono.F, vars.to.regress = "percent.mt", verbose = TRUE)


#### DIMENSIONALITY REDUCTION
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
ibsm.mono.F   <- RunPCA(ibsm.mono.F  , npcs = 30, verbose = FALSE)
ibsm.mono.F   <- RunUMAP(ibsm.mono.F  , reduction = "pca", dims = 1:30)
ibsm.mono.F   <- FindNeighbors(ibsm.mono.F  , reduction = "pca", dims = 1:30)

##Plotting
DefaultAssay(ibsm.mono.F) <- "SCT"
Idents(ibsm.mono.F) <- "day_donor"
ddd <- DimPlot(ibsm.mono.F, label = FALSE, cols= Tfh_colors_list, pt.size = 2, group.by = "day_donor") +  scale_fill_manual(values=Tfh_colors_list)+
  theme(legend.position = "none")
ddd
##After integration

immune.combined <- readRDS( "data/IBSM_monoSCT_integbyDayDonor_220823.rds")

DefaultAssay(immune.combined) <- "integrated"
Idents(immune.combined) <- "day_donor"

eee <- DimPlot(immune.combined, label = FALSE, cols= Tfh_colors_list, pt.size = 2, group.by = "day_donor") + scale_fill_manual(values=Tfh_colors_list)
ggsave("output_for_paper/Sup_Fig_1C_dimplot_prepostIntegration_V2.pdf", plot = (ddd+eee),
       height = 6, width = 10, dpi = 300)
eee

##Violin plots with differential markers



