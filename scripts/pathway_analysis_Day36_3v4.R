##This script performs pathway analysis for 
##each cluster found during single cell data analysis.
##Date: 25-Jan 2024

library(pathfindR)
library(pathfindR.data)
library(knitr)
library(tidyverse)
library(purrr)
library(clusterProfiler)
library(here)

##workflow from https://bioc.ism.ac.jp/packages/3.11/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html#supported-organisms
### Day 36 pathway analysis
#Excludes D8VSd36
set.seed(1234)

psbDE <- read.csv(file=here("output_for_paper/edger/psbulk_sigGens__all_DayClusters.csv"))

colnames(psbDE)
#Renaming columns
colnames(psbDE)[2] ="celltype"
colnames(psbDE)[4] ="genes"
colnames(psbDE)[5] ="Log2FC"
colnames(psbDE)[9] ="adj_Pval"

colnames(psbDE)
#Create a new variable adding cluster + day if not present already.
psbDE$cluster_day <- NA
psbDE$cluster_day <- paste0(psbDE$celltype,"_",psbDE$day)

##Selecting classical monocytes day 36.
psbDE <- psbDE %>%
  select(-X) %>%
  filter(cluster_day %in% c("cMonocyte-COQ7_day36", "cMonocyte-IL_day36", "cMonocyte-ISG15_day36", "cMonocyte-S100A8_day36"))

unique(psbDE$cluster_day)
colnames(psbDE)
##This forloop is:
##selecting the columns required for input ("genes", "Log2FC""Pval", "cluster_day")
##creating lists of the sets for splitting,
##removing sets (cluster_day) that don't contain enriched pathways

psbDE %>% select("genes", "Log2FC", "adj_Pval", "celltype", "day")%>% 
  mutate(cluster_day = paste0(celltype,"_",day)) %>%
  select("genes", "Log2FC", "adj_Pval", "cluster_day")-> psbDEup

head(psbDEup)

#Creating a list per set or cluster_day that only has gene, log2FC and adj_pval
psbDESplit = split(psbDEup[,c(1:3)], f = psbDEup$cluster_day)
#Getting the clusters_day names
clust_names <- names(psbDESplit)
clust_names
#Creating the output for the for loop
pathway_outcome = data.frame(matrix(NA, nrow= 1, ncol = 10))
#Adding column names
colnames(pathway_outcome) <- c("ID", "Term_Description", "Fold_Enrichment", "occurrence",
                               "support", "lowest_p", "highest_p", "Up_regulated", 
                               "Down_regulated", "cluster_names")
head(pathway_outcome)

##KEGG is default for other geneset databases “KEGG”, “Reactome”, “BioCarta”,
##“GO-All”, “GO-BP”, “GO-CC”, “GO-MF”, “cell_markers”, “mmu_KEGG”, “Custom”

##running pathfindR for each cluster

for (i in seq_along(psbDESplit)) {
  print(clust_names[i])
  tmp = try(run_pathfindR(psbDESplit[[i]], gene_sets = 'Reactome', adj_method = "fdr",
                          enrichment_threshold = 0.01))
  tmp$cluster_names = try(clust_names[i])
  pathway_outcome = rbind(pathway_outcome, tmp)
}

#Formatting output
colnames(pathway_outcome)
head(pathway_outcome)
pathway_outcome = pathway_outcome[-1,]
head(pathway_outcome)
##Renaming columns
colnames(pathway_outcome)[10] = "Cluster_name"
##checking out
clustEnriTerms <- unique(pathway_outcome$Cluster)
clustEnriTerms

##write results for reactome
#write.csv(pathway_outcome, file=here("data/pathway_analysis/reactome/pathways_D8D36_cMonos_180924.csv"))
write.csv(pathway_outcome, file=here("data/pathway_analysis/reactome/pathways_D36_cMonos_021024.csv"))

#enrichement table
#Creating a list per set or cluster_day that only has gene, log2FC and adj_pval
pathway_outcomeSplit = split(pathway_outcome[,c(1,8,9)], f = pathway_outcome$Cluster)
head(pathway_outcomeSplit)

#experimental matrix
# Creating a list per cluster and donor
counts0 <- read.csv(file=here("data/pathway_analysis/monoIBSM_Averagemat_LOGNORMALISED_donordaycluster_cMonsDay36.csv"), sep = "," )
#counts0 <- read.csv(file=here("output_for_paper/scpa/monoIBSM_Averagemat_LOGNORMALISED_donordaycluster_150424.csv"), sep = "," )

colnames(counts0)
counts0 <- counts0%>%select(-X)
rownames(counts0) <- counts0$donor_day_cluster

counts1 <- separate(counts0, col = donor_day_cluster, into = c("donor", "day","cluster", "subtype"), sep="-", remove=FALSE)

counts1$celltype <- paste0(counts1$cluster,"_",counts1$subtype)

dim(counts1)
counts2 <- counts1 %>%
  filter(day== "day0"|day=="day36") %>%
  filter(cluster == "cMonocyte") %>%
  select(-donor_day_cluster, -donor, -day, -cluster, -subtype)


##split the exp matrix by cluster

countsT<-t(counts2)
countsT2 <- as.data.frame(countsT)%>%
  dplyr::select("donor0-day0-cMonocyte-COQ7", "donor0-day0-cMonocyte-IL1B", "donor0-day0-cMonocyte-ISG15", "donor0-day0-cMonocyte-S100A8", 
         "donor1-day0-cMonocyte-COQ7", "donor1-day0-cMonocyte-IL1B", "donor1-day0-cMonocyte-ISG15", "donor1-day0-cMonocyte-S100A8",
         "donor2-day0-cMonocyte-COQ7", "donor2-day0-cMonocyte-IL1B", "donor2-day0-cMonocyte-ISG15", "donor2-day0-cMonocyte-S100A8",
         "donor3-day0-cMonocyte-COQ7", "donor3-day0-cMonocyte-IL1B", "donor3-day0-cMonocyte-ISG15", "donor3-day0-cMonocyte-S100A8",
         "donor0-day36-cMonocyte-COQ7", "donor0-day36-cMonocyte-IL1B", "donor0-day36-cMonocyte-ISG15", "donor0-day36-cMonocyte-S100A8",
         "donor1-day36-cMonocyte-COQ7", "donor1-day36-cMonocyte-IL1B", "donor1-day36-cMonocyte-ISG15", "donor1-day36-cMonocyte-S100A8", 
         "donor2-day36-cMonocyte-COQ7", "donor2-day36-cMonocyte-IL1B", "donor2-day36-cMonocyte-ISG15", "donor2-day36-cMonocyte-S100A8",
         "donor3-day36-cMonocyte-COQ7", "donor3-day36-cMonocyte-IL1B", "donor3-day36-cMonocyte-ISG15", "donor3-day36-cMonocyte-S100A8")

head(countsT2)

countsT2_1 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("COQ7")))
countsT2_2 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("IL1B")))
countsT2_3 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("ISG15")))
countsT2_4 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("S100A8")))

countsTSplit <- list(countsT2_1,countsT2_2,countsT2_3, countsT2_4)

countsT3 <-as.matrix(countsT2)

Cluster_name = c("cMonocyte-COQ7_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-S100A8_day36")
##Change names to get unique clusters
pathway_outcome <- read.csv(here("data/pathway_analysis/reactome/pathways_D36_cMonos_021024.csv"))
pathway_1 <- pathway_outcome%>%
  filter(str_detect(Cluster_name, "ISG15"))%>%
  select(ID, Term_Description, Up_regulated, Down_regulated, Cluster_name)
head(pathway_1)  

colnames(pathway_1)[5]<- "Cluster"
unique(pathway_1$Cluster)
head(pathway_1)

pathways_Z <- score_terms(pathway_1, countsT2_4, use_description = TRUE)

a <- plot_scores(pathways_Z, low ="#5e4fa2", mid = "white", high = "#9e0142" )
a

#ggsave(here("data/pathway_analysis/reactome/plots/zcosre_persample_COQ7.pdf"), plot = a, device = "pdf", height = 16, width = 20)
#ggsave(here("data/pathway_analysis/reactome/plots/zcosre_persample_IL1B.pdf"), plot = a, device = "pdf", height = 16, width = 20)
#ggsave(here("data/pathway_analysis/reactome/plots/zcosre_persample_ISG15.pdf"), plot = a, device = "pdf", height = 16, width = 20)
#ggsave(here("data/pathway_analysis/reactome/plots/zcosre_persample_S100A8.pdf"), plot = a, device = "pdf", height = 16, width = 20)

##write results for reactome
write.csv(pathways_Z, file=here("data/pathway_analysis/reactome/pathways_D36_cMonos_S100A8_Zscore_180924.csv"))
write.csv(pathways_Z, file=here("data/pathway_analysis/reactome/pathways_D36_cMonos_COQ7_Zscore_180924.csv"))
write.csv(pathways_Z, file=here("data/pathway_analysis/reactome/pathways_D36_cMonos_IL_Zscore_180924.csv"))
write.csv(pathways_Z, file=here("data/pathway_analysis/reactome/pathways_D36_cMonos_ISG15_Zscore_180924.csv"))


####Calculating overall z-score from the adj p values edgeR workflow####

source(here("scripts/helper_functions.R"))

##Filtering more relevant enriched terms
summary(pathway_outcome$occurrence)
summary(pathway_outcome$Fold_Enrichment)
summary(pathway_outcome$support)

pathway_outcome0 <- pathway_outcome %>% 
  filter (occurrence >5 & Fold_Enrichment>7)
dim(pathway_outcome0)
table(pathway_outcome0$Cluster)
##input dataset !!Do not change colnames for the love of GoD

psbDE <- read.csv(file=here("output_for_paper/edger/psbulk_sigGens__all_DayClusters.csv"))

colnames(psbDE)
#Renaming columns
colnames(psbDE)[2] ="celltype"
colnames(psbDE)[4] ="genes"
colnames(psbDE)[5] ="Log2FC"
colnames(psbDE)[9] ="adj_Pval"

colnames(psbDE)
#Create a new variable adding cluster + day if not present already.
psbDE$cluster_day <- NA
psbDE$cluster_day <- paste0(psbDE$celltype,"_",psbDE$day)

##Selecting classical monocytes day 36.
psbDEup <- psbDE %>%
  select(-X) %>%
  filter(cluster_day %in% c("cMonocyte-COQ7_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-S100A8_day36"))

unique(psbDE$cluster_day)
colnames(psbDE)

psbDEup %>% select("genes", "Log2FC", "adj_Pval", "celltype", "day")%>% 
  mutate(cluster_day = paste0(celltype,"_",day)) %>%
  select("genes", "Log2FC", "adj_Pval", "cluster_day")-> psbDEup


##Calculating z score based on adj p value with helper function, the direction is based on logfc
head(psbDEup)
psbDEup$z_score <- calc_zscore_from_adjp_logfc(psbDEup$adj_Pval, psbDEup$Log2FC)
##Filtering z-scores less than -2 and more than 2
psbDEup_f <- psbDEup %>%
  filter(z_score >= 2 | z_score <=-2)
summary(psbDEup$z_score)

##Calculating z score on the adj p value of the enrichment analysis (NOT the DEA
head(pathway_outcome)
##Filtering z-scores less than -2 and more than 2
psbDEup_f <- psbDEup %>%
  filter(z_score >= 2 | z_score <=-2)
summary(psbDEup$z_score)

##Using heatmap function
# Define your color scheme
low_color <- "#5e4fa2"
mid_color <- "#74add1"
high_color <- "#9e0142"

#Define sample list
sample_groups = c("cMonocyte-S100A8_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-COQ7_day36" )
names(sample_groups) <- c("cMonocyte-S100A8_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-COQ7_day36" )

# Create the heatmap
heatmap <- create_pathfindr_heatmap(pathway_outcome0, psbDEup_f, sample_groups, low_color, mid_color, high_color)
#head(heatmap@matrix)

# Save the heatmap as a PDF
pdf(here("data/pathway_analysis/reactome/plots/pathfindr_z_score_Filt_Day36_heatmap.pdf"), width = 12, height = 20)
draw(heatmap, heatmap_legend_side = "left")
dev.off()

#Extract information
Path_Av_Zcore <- as.data.frame(heatmap@matrix)
write.csv(Path_Av_Zcore, file=here("data/pathway_analysis/reactome/heatmap_pathways_z_score_Filt_Day36_060225.csv"))

##Picking the top 25 for each group
# Split by cluster_day within the dataframe
psbDEup_f_outcome = data.frame(matrix(NA, nrow= 1, ncol = 5))
#Adding column names
colnames(psbDEup_f_outcome) <- c("genes", " Log2FC", "adj_Pval", "cluster_day",
                                 "z_score")

psbDEup_f_Split <- split(psbDEup_f, f = psbDEup_f$cluster_day)

# Iterate over each dataframe in psbDEup_f_Split
psbDEup_f_top25 <- lapply(names(psbDEup_f_Split), function(cluster_name) {
  # Get the dataframe for the current cluster
  df <- psbDEup_f_Split[[cluster_name]]
  
  # Arrange by adj_Pval and select the top 25 rows
  tmp <- df %>%
    dplyr::arrange(adj_Pval) %>%
    head(25)
  
  # Add the cluster_name to the dataframe
  tmp$cluster_names <- cluster_name
  
  # Store the result in the list
  psbDEup_f_outcome <<- append(psbDEup_f_outcome, list(tmp))
})

# Combine the list of dataframes into a single dataframe
psbDEup_f_outcome <- do.call(rbind, psbDEup_f_outcome)
head(psbDEup_f_outcome)

psbDEup_f_top25_df <- psbDEup_f_outcome[-c(1:5),]
head(psbDEup_f_top25_df)
# Create the heatmap
heatmap0 <- create_pathfindr_heatmap(pathway_outcome0, psbDEup_f_top25_df, sample_groups, low_color, mid_color, high_color)
# Save the heatmap as a PDF
pdf(here("data/pathway_analysis/reactome/plots/pathfindr_z_scoreFDR_Filt_Day36_heatmap.pdf"), width = 12, height = 20)
draw(heatmap0, heatmap_legend_side = "left")
dev.off()

#Extract information
Path_Av_Zcore <- as.data.frame(heatmap@matrix)
write.csv(Path_Av_Zcore, file=here("data/pathway_analysis/reactome/heatmap_pathways_z_score_Filt_Day36_060225.csv"))


##Modyfing heatmap
heatmap2 <- heatmap
#Extract information
Path_Av_Zcore <- as.data.frame(heatmap2@matrix)
#Change order of values
Path_Av_Zcore <- Path_Av_Zcore[c(149,1:148),]
#convert to matrix
reorder <- as.matrix(Path_Av_Zcore)
#Replace in object
heatmap2@matrix<- reorder

##Changing rownames in heatmap
path_row_names <- as.data.frame(heatmap2@row_names_param$labels)
#Change order of value
path_row_names2 <- path_row_names[c(149,1:148),]                               
#convert to character
reorder2 <- as.character(path_row_names2)
#Replace in object
heatmap2@row_names_param$labels<- reorder2

heatmap2@row_names_param$side <- "right"

heatmap2@column_names_param$rot <- as.numeric(89)

draw(heatmap2)

write.csv(Path_Av_Zcore, file=here("data/pathway_analysis/reactome/heatmap_pathways_D36_cMonos_AvZscore_180924.csv"))

####Plotting from Day0 vs Day8####

counts2 <- counts1 %>%
  filter(day== "day0"|day=="day8") %>%
  filter(cluster == "cMonocyte") %>%
  select(-donor_day_cluster, -donor, -day, -cluster, -subtype)

counts3 <- counts1 %>%
  select(-donor_day_cluster, -donor, -day, -cluster, -subtype, -celltype)


##split the exp matrix by cluster

countsT<-t(counts3)
countsT2 <- as.data.frame(countsT)%>%
  dplyr::select("donor0-day0-cMonocyte-COQ7", "donor0-day0-cMonocyte-IL1B", "donor0-day0-cMonocyte-ISG15", "donor0-day0-cMonocyte-S100A8", 
                "donor1-day0-cMonocyte-COQ7", "donor1-day0-cMonocyte-IL1B", "donor1-day0-cMonocyte-ISG15", "donor1-day0-cMonocyte-S100A8",
                "donor2-day0-cMonocyte-COQ7", "donor2-day0-cMonocyte-IL1B", "donor2-day0-cMonocyte-ISG15", "donor2-day0-cMonocyte-S100A8",
                "donor3-day0-cMonocyte-COQ7", "donor3-day0-cMonocyte-IL1B", "donor3-day0-cMonocyte-ISG15", "donor3-day0-cMonocyte-S100A8",
                "donor0-day36-cMonocyte-COQ7", "donor0-day36-cMonocyte-IL1B", "donor0-day36-cMonocyte-ISG15", "donor0-day36-cMonocyte-S100A8",
                "donor1-day36-cMonocyte-COQ7", "donor1-day36-cMonocyte-IL1B", "donor1-day36-cMonocyte-ISG15", "donor1-day36-cMonocyte-S100A8", 
                "donor2-day36-cMonocyte-COQ7", "donor2-day36-cMonocyte-IL1B", "donor2-day36-cMonocyte-ISG15", "donor2-day36-cMonocyte-S100A8",
                "donor3-day36-cMonocyte-COQ7", "donor3-day36-cMonocyte-IL1B", "donor3-day36-cMonocyte-ISG15", "donor3-day8-cMonocyte-S100A8")

head(countsT2)

countsT2_1 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("COQ7")))
countsT2_2 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("IL1B")))
countsT2_3 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("ISG15")))
countsT2_4 <- as.matrix(countsT2%>%
                          dplyr::select(ends_with("S100A8")))

countsTSplit <- list(countsT2_1,countsT2_2,countsT2_3, countsT2_4)

##test
pathway_1 <- pathway_outcome%>%
  filter(str_detect(Cluster_name, "COQ7"))%>%
  select(ID, Term_Description, Up_regulated, Down_regulated)


pathways_Z <- score_terms(pathway_1, countsT2_1, use_description = TRUE)

a <- plot_scores(pathways_Z)

##Change the filtering to get each cluster.
ggsave(here("data/pathway_analysis/reactome/plots/zcosre_persample.pdf"), plot = a, device = "pdf", height = 16, width = 20)
##write results for reactome
write.csv(pathways_Z, file=here("data/pathway_analysis/reactome/zcosre_persample_COQ7.csv"))


##Dotplot of enriched terms
colnames(pathway_outcome)[10]<- "Cluster"
colnames(pathway_outcome)
#Renaming column

##High clusters day 0
pathway_outcome0 <- pathway_outcome %>% 
  #filter(Cluster %in% day0clusters) %>%
  filter (occurrence >8 & Fold_Enrichment>8)
pathway_outcome0 <- factor(pathway_outcome0$Cluster, levels = )
p0 <- enrichment_chart(pathway_outcome0, top_terms = 5, plot_by_cluster = TRUE)

p0
ggsave(plot = p0, filename = "data/pathway_analysis/reactome/plots/dotplot_Patwhay_Day836_190924.pdf", 
       width = 300, height= 350, units = "mm", device = "pdf") 

##Clustering the enrichded terms

# change agglomeration method (default="average") for hierarchical clustering
clustered_df <- cluster_enriched_terms(pathway_outcome, clu_method="average", plot_hmap=TRUE, plot_dend=TRUE,
                                       use_description=TRUE)

max(clustered_df$Cluster)
cluster_Sum <- as.data.frame(table(clustered_df$Cluster))
colnames(cluster_Sum)[1] <- "Cluster"

max(cluster_Sum$Freq)
cluster_Sum2 <- cluster_Sum%>%
  filter(Freq>9)
keep <- as.vector(cluster_Sum2$Cluster)
keep
##subsetting clusters with 10 or more terms per cluster
clustered_df_sub <- clustered_df[clustered_df$Cluster%in%keep,]
clustered_df_sub <- clustered_df_sub %>%
  filter(Status == "Representative")
##write results for reactome
write.csv(clustered_df_sub, file=here("data/pathway_analysis/reactome/Clustered_Pathways_D36_cMonos_180924.csv"))

unique(clustered_df_sub$Term_Description)

##plotting clustered terms with more than 10 pathways

kappa_mx <- create_kappa_matrix(pathway_outcome)

clus_obj <- hierarchical_term_clustering(
  kappa_mat=kappa_mx,
  enrichment_res=pathway_outcome,
  num_clusters = NULL,
  use_description = FALSE,
  clu_method = "average",
  plot_hmap = FALSE,
  plot_dend = TRUE
)

##Trying to select clusters from the clus object
clus_obj_sub <- clus_obj[clus_obj%in%keep]
ID_sub <- names(clus_obj_sub)
pathway_outcome_sub <- pathway_outcome[pathway_outcome$ID%in%ID_sub,]
##selecting columns
kappa_mx_sub0 <- kappa_mx[,ID_sub]
kappa_mx_sub <- kappa_mx_sub0[ID_sub,]

#Subset clusters
cluster_graph_vis(clus_obj_sub, 
                  kappa_mat = kappa_mx_sub,
                  enrichment_res = pathway_outcome_sub, kappa_threshold = 0.4,vertex.label.cex = 0.9, vertex.size.scaling=10)

ggsave(plot = , filename = "data/pathway_analysis/reactome/plots/clusters_subEnrichPath_190924.pdf", 
       width = 600, height= 650, units = "mm", device = "pdf") 


##all 250 clusters
cluster_graph_vis(clus_obj, 
                  kappa_mat = kappa_mx,
                  enrichment_res = pathway_outcome, kappa_threshold = 0.4)

ggsave(plot = all, filename = "data/pathway_analysis/reactome/plots/clusters250_allEnrichPath_190924.pdf", 
       width = 600, height= 650, units = "mm", device = "pdf") 



term_names <- c("Chk1/Chk2(Cds1) mediated inactivation of Cyclin B:Cdk1 complex",
                "SARS-CoV-2 targets host intracellular signalling and regulatory pathways",
                "Activation of BAD and translocation to mitochondria",
                "Mitochondrial tRNA aminoacylation",
                "Metalloprotease DUBs",
                "DNA Damage Recognition in GG-NER","Alpha-protein kinase 1 signaling pathway",
                "ERKs are inactivated",
                "ERK/MAPK targets",
                "Regulation of RUNX1 Expression and Activity",
                "MAPK3 (ERK1) activation",
                "VLDLR internalisation and degradation",
                "RAF-independent MAPK1/3 activation",
                "Modulation by Mtb of host immune system",
                "Prevention of phagosomal-lysosomal fusion",
                "TICAM1, RIP1-mediated IKK complex recruitment",
                "PERK regulates gene expression",
                "ATF4 activates genes in response to endoplasmic reticulum  stress",
                "Senescence-Associated Secretory Phenotype (SASP)",
                "CaMK IV-mediated phosphorylation of CREB",
                "CREB1 phosphorylation through the activation of Adenylate Cyclase",
                "PKA-mediated phosphorylation of CREB",
                "MyD88 cascade initiated on plasma membrane",
                "Toll Like Receptor 10 (TLR10) Cascade",
                "Toll Like Receptor 5 (TLR5) Cascade",
                "TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation",
                "MyD88 dependent cascade initiated on endosome",
                "Toll Like Receptor 7/8 (TLR7/8) Cascade",
                "Toll Like Receptor 3 (TLR3) Cascade",
                "Toll Like Receptor 9 (TLR9) Cascade",
                "MyD88-independent TLR4 cascade",
                "TRIF(TICAM1)-mediated TLR4 signaling",
                "MyD88:MAL(TIRAP) cascade initiated on plasma membrane",
                "Toll Like Receptor 2 (TLR2) Cascade",
                "Toll Like Receptor TLR1:TLR2 Cascade",
                "Toll Like Receptor TLR6:TLR2 Cascade",
                "MAP kinase activation",
                "Interleukin-17 signaling",
                "Transcriptional Regulation by MECP2",
                "AKT phosphorylates targets in the nucleus",
                "Estrogen-dependent nuclear events downstream of ESR-membrane signaling",
                "DARPP-32 events",
                "TAK1-dependent IKK and NF-kappa-B activation",
                "Nucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways",
                "NOD1/2 Signaling Pathway",
                "FLT3 Signaling",
                "Regulation of FOXO transcriptional activity by acetylation",
                "Regulation of localization of FOXO transcription factors",
                "FOXO-mediated transcription of cell death genes",
                "RIP-mediated NFkB activation via ZBP1",
                "ZBP1(DAI) mediated induction of type I IFNs",
                "TRAF6 mediated NF-kB activation",
                "Diseases associated with the TLR signaling cascade",
                "Diseases of Immune System",
                "Regulated proteolysis of p75NTR",
                "The NLRP3 inflammasome",
                "Negative regulation of MAPK pathway",
                "Aberrant regulation of mitotic G1/S transition in cancer due to RB1 defects",
                "Defective binding of RB1 mutants to E2F1,(E2F2, E2F3)",
                "Interleukin-10 signaling",
                "GAB1 signalosome",
                "ERBB2 Activates PTK6 Signaling",
                "GRB2 events in EGFR signaling",
                "SHC1 events in EGFR signaling",
                "ERBB2 Regulates Cell Motility",
                "GRB2 events in ERBB2 signaling",
                "PI3K events in ERBB2 signaling",
                "Signaling by ERBB2 TMD/JMD mutants",
                "SHC1 events in ERBB2 signaling",
                "Signaling by ERBB2 KD Mutants",
                "Signaling by EGFR in Cancer",
                "Signaling by ERBB2 in Cancer",
                "Downregulation of ERBB2 signaling",
                "EGFR downregulation",
                "Signaling by EGFR",
                "FOXO-mediated transcription",
                "Constitutive Signaling by NOTCH1 HD+PEST Domain Mutants",
                "Constitutive Signaling by NOTCH1 PEST Domain Mutants",
                "Signaling by NOTCH1 HD+PEST Domain Mutants in Cancer",
                "Signaling by NOTCH1 PEST Domain Mutants in Cancer",
                "Signaling by NOTCH1 in Cancer",
                "Signaling by NOTCH1",
                "Interleukin-1 signaling",
                "NF-kB is activated and signals survival",
                "p75NTR signals via NF-kB",
                "SUMOylation of immune response proteins",
                "Tristetraprolin (TTP, ZFP36) binds and destabilizes mRNA",
                "Negative regulation of FLT3",
                "Regulation of KIT signaling",
                "Amino acid transport across the plasma membrane",
                "Mitophagy",
                "Transcriptional activity of SMAD2/SMAD3:SMAD4 heterotrimer",
                "IRAK2 mediated activation of TAK1 complex",
                "IRAK2 mediated activation of TAK1 complex upon TLR7/8 or 9 stimulation",
                "TRAF6-mediated induction of TAK1 complex within TLR4 complex",
                "Constitutive Signaling by NOTCH1 HD Domain Mutants",
                "Signaling by NOTCH1 HD Domain Mutants in Cancer",
                "PECAM1 interactions",
                "Signal regulatory protein family interactions",
                "p38MAPK events",
                "InlA-mediated entry of Listeria monocytogenes into host cells",
                "DCC mediated attractive signaling",
                "VEGFR2 mediated cell proliferation",
                "GP1b-IX-V activation signalling",
                "Spry regulation of FGF signaling",
                "MET activates PTK2 signaling",
                "Signaling by BMP",
                "SARS-CoV-1 activates/modulates innate immune responses",
                "Downregulation of SMAD2/3:SMAD4 transcriptional activity",
                "TNF receptor superfamily (TNFSF) members mediating non-canonical NF-kB pathway",
                "Regulation of gene expression by Hypoxia-inducible Factor",
                "Signaling by ALK",
                "STAT3 nuclear events downstream of ALK signaling",
                "TP53 regulates transcription of additional cell cycle genes whose exact role in the p53 pathway remain uncertain",
                "Gastrin-CREB signalling pathway via PKC and MAPK",
                "CASP8 activity is inhibited",
                "Dimerization of procaspase-8",
                "Regulation by c-FLIP",
                "Caspase activation via Death Receptors in the presence of ligand",
                "SARS-CoV-1 targets host intracellular signalling and regulatory pathways",
                "SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription",
                "Downregulation of TGF-beta receptor signaling",
                "Toll-like Receptor Cascades",
                "MECP2 regulates neuronal receptors and channels",
                "Regulation of MECP2 expression and activity",
                "Signaling by NOTCH2",
                "Transcriptional regulation of granulopoiesis",
                "Constitutive Signaling by AKT1 E17K in Cancer",
                "Endosomal/Vacuolar pathway",
                "NOTCH2 intracellular domain regulates transcription",
                "NOTCH1 Intracellular Domain Regulates Transcription",
                "NOTCH4 Intracellular Domain Regulates Transcription",
                "RUNX3 regulates NOTCH signaling",
                "Regulation of gene expression in late stage (branching morphogenesis) pancreatic bud precursor cells",
                "NOTCH3 Intracellular Domain Regulates Transcription",
                "RUNX2 regulates osteoblast differentiation",
                "FOXO-mediated transcription of cell cycle genes",
                "Formation of Senescence-Associated Heterochromatin Foci (SAHF)",
                "AKT phosphorylates targets in the cytosol",
                "Aberrant regulation of mitotic cell cycle due to RB1 defects",
                "Diseases of mitotic cell cycle",
                "Cyclin D associated events in G1",
                "G1 Phase",
                "TP53 Regulates Transcription of Cell Cycle Genes",
                "MAP3K8 (TPL2)-dependent MAPK1/3 activation",
                "Response of EIF2AK1 (HRI) to heme deficiency",
                "Nuclear events stimulated by ALK signaling in cancer",
                "Transcriptional Regulation by VENTX",
                "Reduction of cytosolic Ca++ levels")


#### Calculating z score from enrichement adjusted p value


pathway_outcome0 = pathway_outcome %>%
  filter (occurrence >8 & Fold_Enrichment>8) %>%
  select(-final_dir, -z.score, -dir, -Down_Reg_No, -Up_Reg_No)

write.csv(pat)
##Using heatmap function
# Define your color scheme
low_color <- "#5e4fa2"
mid_color <- "#74add1"
high_color <- "#9e0142"

#Define sample list
sample_groups = c("cMonocyte-S100A8_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-COQ7_day36" )
names(sample_groups) <- c("cMonocyte-S100A8_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-COQ7_day36" )

# Create the heatmap
heatmap <- create_pathfindr_heatmap2(pathway_outcome0, sample_groups, low_color, mid_color, high_color)
heatmap
head(heatmap@matrix)


# Save the heatmap as a PDF
pdf(here("data/pathway_analysis/reactome/plots/pathfindr_z_score_Filt_Day36_heatmap_v2fix.pdf"), width = 12, height = 20)
draw(heatmap, heatmap_legend_side = "left")
dev.off()



###Selecting only the 25 more significant terms

pathway_outcome0 = pathway_outcome %>%
  filter (occurrence >8 & Fold_Enrichment>8) %>%
  select(-final_dir, -z.score, -dir, -Down_Reg_No, -Up_Reg_No)

pathway_f_Split <- split(pathway_outcome0, f = pathway_outcome0$Cluster_name)
# Split by cluster_day within the dataframe
pathway_f_outcome = data.frame(matrix(NA, nrow= 1, ncol = 10))
#Adding column names
colnames(pathway_f_outcome) <- c("ID", "Term_Description", "Fold_Enrichment", "occurrence",
                                 "support", "lowest_p", "highest_p", "Up_regulated",
                                 "Down_regulated", "Cluster_name")


pathway_f_top25 <- lapply(names(pathway_f_Split), function(Cluster_name) {
  # Get the dataframe for the current cluster
  df <- pathway_f_Split[[Cluster_name]]
  
  # Arrange by adj_Pval (lowest_p) and select the top 25 rows
  tmp <- df %>%
    dplyr::arrange(lowest_p) %>%
    head(25)
  
  # Add the cluster_name to the dataframe
  tmp$Cluster_name <- Cluster_name
  
  # Store the result in the list
  pathway_f_outcome <<- append(pathway_f_outcome, list(tmp))
})

# Combine the list of dataframes into a single dataframe
pathway_f_outcomes <- do.call(rbind, pathway_f_outcome)
head(pathway_f_outcomes)


pathway_f_outcome <- pathway_f_outcomes[-c(1:10),]
head(pathway_f_outcome)

##Using heatmap function
# Define your color scheme
low_color <- "#5e4fa2"
mid_color <- "#74add1"
high_color <- "#9e0142"

#Define sample list
sample_groups = c("cMonocyte-S100A8_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-COQ7_day36" )
names(sample_groups) <- c("cMonocyte-S100A8_day36","cMonocyte-IL_day36","cMonocyte-ISG15_day36","cMonocyte-COQ7_day36" )

# Create the heatmap
heatmap2 <- create_pathfindr_heatmap2(pathway_f_outcome, sample_groups, low_color, mid_color, high_color)
heatmap2
#head(heatmap@matrix)

# Save the heatmap as a PDF
pdf(here("data/pathway_analysis/reactome/plots/pathfindr_z_score_Filt_Day36_heatmap_25f_v2fix.pdf"), width = 12, height = 20)
draw(heatmap, heatmap_legend_side = "left")
dev.off()

##writing down results with z scores
# Create direction
pathway_outcome <- pathway_outcome %>%
  mutate(
    dir = case_when(
      Down_regulated != "" & Up_regulated == "" ~ -1,
      Down_regulated == "" & Up_regulated != "" ~ 1,
      TRUE ~ 0
    ),
    Down_Reg_No = ifelse(Down_regulated == "", 0, str_count(Down_regulated, ",") + 1),
    Up_Reg_No = ifelse(Up_regulated == "", 0, str_count(Up_regulated, ",") + 1),
    final_dir = case_when(
      dir == 0 & Down_Reg_No > Up_Reg_No ~ -1,
      dir == 0 & Down_Reg_No < Up_Reg_No ~ 1,
      TRUE ~ dir
    )
  )

# Calculate Z-scores
pathway_outcome$z.score <- calc_zscore_from_adjp_logfc(pathway_outcome$lowest_p, pathway_outcome$final_dir)

write.csv(pathway_outcome, here("data/pathway_analysis/reactome/all_pathways_D36_cMonos_Zscore_lowp.csv"))
