##

library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(UpSetR)
library(here)
#Plots from pseudobulk DE analysis. 


####Figure1c_Barplot showing Prop of up and downregulated genes per celltype per day __Updated####
psbDE <- read.csv(here(file="output_for_paper/edger/psbulk_sigGens__all_DayClusters.csv"))

colnames(psbDE)
#Renaming columns
colnames(psbDE)[2] ="celltype"
colnames(psbDE)[4] ="genes"
colnames(psbDE)[5] ="Log2FC"
colnames(psbDE)[9] ="Pval"

head(psbDE)

##removing first column
psbDE <- psbDE[,-1]
## checking clusters name
unique(psbDE$celltype)
##fixing names
psbDE$celltype <- ifelse(psbDE$celltype == "cMonocyte-IL", "cMonocyte-IL1B", psbDE$celltype)

psbDE$day <- ifelse(psbDE$day == "day8", "Day0 vs Day8", psbDE$day)
psbDE$day <- ifelse(psbDE$day == "day36", "Day0 vs Day36", psbDE$day)
psbDE$day <- ifelse(psbDE$day == "day836", "Day8 vs Day36", psbDE$day)
##creating day factors
psbDE$day <- factor(psbDE$day, levels = c("Day0 vs Day8", "Day0 vs Day36", "Day8 vs Day36" ))

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16", "DCs1", "DCs2", "pDCs" )
##Creating day factors
psbDE$celltype <- factor(psbDE$celltype, levels = celltype_fact)

#Creating dataset with direction 
psbDE$Dir <- NA
psbDE <- psbDE%>%
  mutate(
    Dir =case_when(
      Log2FC>0 ~ "Up",
      Log2FC<0 ~ "Down",
      TRUE ~ NA_character_ 
    )
  )
table(psbDE$Dir)
##Frequency table
celltypes_freq <- as.data.frame(table(psbDE$celltype, psbDE$day, psbDE$Dir))
colnames(celltypes_freq)[1] ="Celltype"
colnames(celltypes_freq)[2] ="Day"
colnames(celltypes_freq)[3] ="Direction"
colnames(celltypes_freq)[4] ="No_of_DEGs"

##Calculating the proportion of Up and Down regulated genes per day per celltype

# Calculate the total number of DEGs for each cell type and day
celltypes_freq <- celltypes_freq %>%
  group_by(Celltype, Day) %>%
  mutate(Total_DEGs = sum(No_of_DEGs)) %>%
  ungroup()

# Calculate the proportion of up and down regulated genes
celltypes_freq <- celltypes_freq %>%
  mutate(Proportion = No_of_DEGs / Total_DEGs * 100)

#Removing one comparison and changing day input
celltypes_freq<- celltypes_freq %>%
  filter(Day != "Day8 vs Day36")

head(celltypes_freq)
celltypes_freq$Day <- factor(celltypes_freq$Day, levels = c("Day0 vs Day8", "Day0 vs Day36"))
celltypes_freq$Direction <- factor(celltypes_freq$Direction, levels = c("Up", "Down"))
celltypes_freq$Celltype <- factor(celltypes_freq$Celltype, levels = celltype_fact)
celltypes_freq 

Tfh_colors_list <- c("#9e0142", "#74add1","#f4a582", "#5e4fa2", "#66c2a5",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")


pop_clus <- ggplot(data=celltypes_freq, aes(x=Celltype, y=No_of_DEGs, fill=Celltype)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values=Tfh_colors_list)+
  scale_y_continuous(limits = c(0,400))+
  facet_grid(Day ~Direction, scales = "free") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(
          size = 12, color = "black"))
  

pop_clus

ggsave("output_for_paper/Fig1c_Barplot_DEGs_per_Daycomp_July2024.pdf", plot = pop_clus,
       height = 8, width = 16, dpi = 300)

####Fig1d1_upsetplots_perday_Which genes are shared among clusters. (3plots)####

####day0 vs day8####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$day=="Day0 vs Day8",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1 & upset1$day== "Day0 vs Day8"] <- "up"
upset1$dir[upset1$Log2FC< -1 & upset1$day== "Day0 vs Day8"] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$celltype,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16", "DCs1", "DCs2", "pDCs" )
##Creating day factors
upset1$celltype <- factor(upset1$celltype, levels = celltype_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_D0vsD8_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 40,
      nsets = 9,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()

###Day 0 vs Day 36
upset1 <- psbDE1[psbDE1$day=="Day0 vs Day36",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1 & upset1$day== "Day0 vs Day36"] <- "up"
upset1$dir[upset1$Log2FC< -1 & upset1$day== "Day0 vs Day36"] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$celltype,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16", "DCs1", "DCs2", "pDCs" )
##Creating day factors
upset1$celltype <- factor(upset1$celltype, levels = celltype_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_D0vsD36_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 40,
      nsets = 9,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()


###Day 8 vs Day 36
upset1 <- psbDE1[psbDE1$day=="Day8 vs Day36",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1 & upset1$day== "Day8 vs Day36"] <- "up"
upset1$dir[upset1$Log2FC< -1 & upset1$day== "Day8 vs Day36"] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$celltype,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16", "DCs1", "DCs2", "pDCs" )
##Creating day factors
upset1$celltype <- factor(upset1$celltype, levels = celltype_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_D8vsD36_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 40,
      nsets = 9,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()


####Fig1d2_upsetplots_percluster_Which genes are shared among days (9plots)####

####cMonocyte-IL1B8####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="cMonocyte-IL1B",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_cMonocyte-IL1B_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()


####cMonocyte-COQ7####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="cMonocyte-COQ7",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_cMonocyte-COQ7_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()

####cMonocyte-ISG15####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="cMonocyte-ISG15",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_cMonocyte-ISG15_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()

####cMonocyte-S100A8####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="cMonocyte-S100A8",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_cMonocyte-S100A8_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()

####intMonocyte####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="intMonocyte",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_intMonocyte_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()

####ncMonocyte-CD16####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="ncMonocyte-CD16",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_ncMonocyte-CD16_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()
####DCs1####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="DCs1",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_DCs1_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()

####DCs2####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="DCs2",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating celltype factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_DCs2_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 60,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()

####pDCs####
#split by celltype
psbDE1 <- psbDE %>% select(celltype, genes, day, Log2FC, Pval)
upset1 <- psbDE1[psbDE1$celltype=="pDCs",]

head(upset1)

##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$Log2FC>1] <- "up"
upset1$dir[upset1$Log2FC< -1] <- "dw"
tapply(upset1$Log2FC, upset1$dir, summary)
upset1 <- upset1[upset1$dir=="up",]
#upset1$cluster_day_dir <- NA
upset1$cluster_dir <- paste0(upset1$day,"_",upset1$dir)
colnames(upset1)

##Creating cell type factors
day_fact<- c("Day0 vs Day8","Day0 vs Day36", "Day8 vs Day36")
##Creating day factors
upset1$day <- factor(upset1$day, levels = day_fact)


##creating lists of the sets
#upset2 <- upset1[,c(1,7)]
upset2Split = split(upset1[,2], f = upset1$cluster_dir)
head(upset2Split)

pdf("output_for_paper/UpsetR_pDCs_Upgenes_orderbyFreq.pdf", h=16, w =28)
upset(fromList(upset2Split),
      nintersects = 10,
      nsets = 6,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = c(3,3,3,3,3,3),
      point.size = 6,
      line.size = 1)
dev.off()
#c("3","3","3","3","3","3")
#  intersection size tick labels, set size title, set size tick labels, 
#  set names, numbers above bars)
#text.scale=c(intersection size title,
#  intersection size tick labels, set size title, set size tick labels, 
#  set names, numbers above bars)


####Figure1e1. Dotplot with DEGs for classical monocytes ####

##Calling only the clustered terms
pathsD8 <- read.csv("data/pathway_analysis/reactome/Clustered_Pathways_D8_cMonos_031024.csv")
pathsD36 <- read.csv("data/pathway_analysis/reactome/Clustered_Pathways_D36_cMonos_180924.csv")

pathsD8a <- pathsD8%>%
  select(-X)
pathsD36a <- pathsD36%>%
  select(-X)

pathsD8_list <-unique(pathsD8$Term_Description)
pathsD36_list <-unique(pathsD36$Term_Description)  

all_genes_D8 <- paste(pathsD8a$Up_regulated, pathsD8a$Down_regulated, sep = ",")
all_genes_D36 <- paste(pathsD36a$Up_regulated, pathsD36a$Down_regulated, sep = ",")

allgenesD8D36 <-list(all_genes_D8, all_genes_D36 ) 

##Selecting the genes based on clustered pathway analysis 

# AP <- c("CREB1", "CTSL", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
#         "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HSPA5", "HSPA8", "TAP1", "HSPA1A", "HSPA1B", "PSME2")
# 
# Apopt<- c("APAF1", "BIRC2", "GADD45B", "AIFM1", "BCL2A1", 
#          "CASP6", "DDIT3", "EIF2AK3", "NFKBIA", "NFKB1", "PARP1", "TNF", "TNFSF10","TUBA4A")
# 
# Chemo <- c("CCL2","CCR2", "CXCL2", "CXCL8","CXCL9","CXCL10", "CXCL16","CXCR4","CX3CR1",
#            "FOXO3","NCF1","NFKB1","NFKBIA","PIK3CG","SHC1","SRC","STAT1", "STAT2")
# 
# VEGF <- c("PLCG2", "PTGS2", "SRC", "VEGFA")
# NOD <- c("GABARAPL1","GBP1","GBP2","GBP4","GBP5","IL6","IL1B","NLRP3","TRIP6","CARD6","NLRP12")
# ##Important pathways in D8 vs D36
# ##Same across clusters for Day8
# MHClassI_Presentation <- c("FCGR1A", "PSMB10", "PSMB9", "PSME2", "PSMB8", "PSME2")
# Interferon_Sign<- c("FCGR1A", "GBP1", "GBP2", "GBP4", "GBP5", "IRF1", 
#                           "STAT1","IFI35", "PSMB8",  "UBE2L6",
#                           "FCGR1B", "HLA-DQA2", "IRF2", "ISG20",
#                           "IFIT2", "IFIT3", "IFITM3", "ISG15", "MT2A", "PML","XAF1","STAT2")
# Allgenes <- c("CREB1", "CTSL", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
# "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HSPA5", "HSPA8", "TAP1", "HSPA1A", "HSPA1B", "PSME2",
# "APAF1", "BIRC2", "GADD45B", "AIFM1", "BCL2A1",
# "CASP6", "DDIT3", "EIF2AK3", "NFKB1", "PARP1", "TNF", "TNFSF10","TUBA4A",
# "CCL2","CCR2", "CXCL2", "CXCL8","CXCL9","CXCL10", "CXCL16","CXCR4","CX3CR1",
# "FOXO3","NCF1","NFKBIA","PIK3CG","SHC1","STAT1", "STAT2",
# "PTGS2", "SRC","PLCG2", "VEGFA",
# "GABARAPL1","GBP1","GBP2","GBP4","GBP5","IL6","IL1B","NLRP3","TRIP6","CARD6","NLRP12")

##Unique per cluster on Day8
Interleukin_1_signaling <- c("PSMB10", "PSMB8", "PSMB9"," PSME2", "IL1R1", "IL1R2")
CrossPresentation_of_exogenous_antigens <- c("FCGR1A", "PSMB10", "PSMB8", "PSMB9", "PSME2")
Interleukin_10 <- c("CCL2", "CCL4", "CXCL10")
FCGR3Amediated_IL10_synthesis <- c("FCGR1A", "FCGR3A","PLCG2")
Platelet_activation <- c("PLCG2", "THBS1", "VEGFA", "SERPING1")

##Unique per cluster on Day36
TGFR_SMADs <- c("SMAD3", "SMAD7", "UBC", "PPP1R15A", "SKI", "SKIL")
Activated_NOTCH1 <- c("JAG1", "UBC") ##only ILB1
RIPK1_mediated_regulated_necrosis <- c("TNFSF10",	"BIRC2", "UBC")
Toll_Like_Receptor3 <- c("BIRC2", "DUSP4", "DUSP6", "IRAK2", "MAP3K8", "MAPK7", "NFKBIA", "UBC", "TIFA")

Allgenes <- c("PSMB10", "PSMB8", "PSMB9"," PSME2", "IL1R1", "IL1R2",
              "FCGR1A","CCL2", "CCL4", "CXCL10", "FCGR3A","PLCG2","THBS1", "VEGFA", "SERPING1",
              "SMAD3", "SMAD7", "UBC", "PPP1R15A", "SKI", "SKIL",
              "JAG1", "UBC",
              "TNFSF10",	"BIRC2", "UBC",
              "DUSP4", "DUSP6", "IRAK2", "MAP3K8", "MAPK7", "NFKBIA", "UBC", "TIFA")
# Define Paramaters for dotplot
psbDE <- read.csv(file="D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/output_for_paper/edger/psbulk_sigGens__all_DayClusters.csv")
colnames(psbDE)
head(psbDE)
#Renaming columns
colnames(psbDE)[2] ="celltype"
colnames(psbDE)[4] ="genes"
colnames(psbDE)[5] ="Log2FC"
colnames(psbDE)[9] ="Pval"
psbDE <- psbDE[,-1]
##fixing names
psbDE$celltype <- ifelse(psbDE$celltype == "cMonocyte-IL", "cMonocyte-IL1B", psbDE$celltype)
##REmoving anything that is not monocytes
psbDE <- psbDE%>%filter(celltype=="cMonocyte-IL1B"|
                        celltype=="cMonocyte-S100A8"|
                        celltype=="cMonocyte-COQ7"|
                        celltype=="cMonocyte-ISG15"|
                        celltype=="intMonocyte"|
                        celltype=="ncMonocyte-CD16")
##Removing D8Day36 comparison
psbDE<- psbDE%>%filter(day!="day836")
##creating day factors
psbDE$day <- factor(psbDE$day, levels = c("day8","day36"))

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16")
##Creating day factors
psbDE$celltype <- factor(psbDE$celltype, levels = celltype_fact)

##Nicks plot
psbDE$log_p <- -log10(psbDE$Pval)
psbDE$log_p_cat[psbDE$log_p <2] <- "<2"
psbDE$log_p_cat[psbDE$log_p >2 & psbDE$log_p<=10] <- "2-10"
psbDE$log_p_cat[psbDE$log_p >10 ] <- ">10"
psbDE$log_p_cat <- factor(psbDE$log_p_cat, levels = c("<2","2-10", ">10"))

##Selecting all genes from the 5pathwasy of interest in only monocytes
DEGdfAllgenes <- psbDE[psbDE$genes %in% Allgenes,]
unique(DEGdfAllgenes$celltype)
DEGdf60 <- DEGdfAllgenes
#DEGdf60 <- DEGdfAllgenes[DEGdfAllgenes$celltype!="DCs1",]
#DEGdf60 <- DEGdf60[DEGdf60$celltype!="DCs2",]
#DEGdf60 <- DEGdf60[DEGdf60$celltype!="pDCs",]

unique(DEGdf60$celltype)
Tfh_colors_list <- c("#9e0142", "#f4a582", "#5e4fa2", "#66c2a5","#74add1",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")


##Creating celltype factors
celltype_DEGdf<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16" )
##Creating day factors
DEGdf60$celltype <- factor(DEGdf60$celltype, levels = celltype_DEGdf)
DEGdf60$celltype_day <- paste0(DEGdf60$celltype, "_", DEGdf60$day)
cell_day_lev <- c("cMonocyte-COQ7_day8", "cMonocyte-COQ7_day36",
                  "cMonocyte-IL1B_day8", "cMonocyte-IL1B_day36",
                  "cMonocyte-ISG15_day8", "cMonocyte-ISG15_day36",
                  "cMonocyte-S100A8_day8", "cMonocyte-S100A8_day36",
                  "intMonocyte_day8","intMonocyte_day36", "ncMonocyte-CD16_day8","ncMonocyte-CD16_day36")
##Michelle's suggestion
cell_day_lev2 <- c("cMonocyte-COQ7_day8","cMonocyte-IL1B_day8",
                   "cMonocyte-ISG15_day8","cMonocyte-S100A8_day8",
                   "intMonocyte_day8","ncMonocyte-CD16_day8",
                   "cMonocyte-COQ7_day36", "cMonocyte-IL1B_day36",
                   "cMonocyte-ISG15_day36", "cMonocyte-S100A8_day36",
                  "intMonocyte_day36", "ncMonocyte-CD16_day36")
DEGdf60$celltype_day <- factor(DEGdf60$celltype_day, levels = cell_day_lev2)
DEGdf60$genes<- factor(DEGdf60$genes, levels = Allgenes)
# Create dotplot from selected genes
DEGdf_dotplot <- ggplot(DEGdf60, aes(reorder(genes, celltype_day), celltype_day))+
  geom_point(aes(colour=Log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
DEGdf_dotplot
ggsave(DEGdf_dotplot, filename = "output_for_paper/Fig1e_pbDEG_Mono_Allgens_Sept24_NewPAth.pdf", h = 10,w =10, dpi = 400,limitsize = FALSE)





##Selecting AP genes and only monocytes
DEGdfAP <- psbDE[psbDE$genes %in% AP,]
unique(DEGdfAP$celltype)
DEGdf60 <- DEGdfAP[DEGdfAP$celltype!="DCs1",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="DCs2",]

unique(DEGdf$celltype)
##Creating celltype factors
celltype_DEGdf<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16" )
##Creating day factors
DEGdf60$celltype <- factor(DEGdf60$celltype, levels = celltype_DEGdf)
DEGdf60$celltype_day <- paste0(DEGdf60$celltype, "_", DEGdf60$day)
cell_day_lev <- c("cMonocyte-COQ7_day8", "cMonocyte-COQ7_day836", "cMonocyte-COQ7_day36",
                  "cMonocyte-IL1B_day8", "cMonocyte-IL1B_day836", "cMonocyte-IL1B_day36",
                  "cMonocyte-ISG15_day8", "cMonocyte-ISG15_day836", "cMonocyte-ISG15_day36",
                  "cMonocyte-S100A8_day8", "cMonocyte-S100A8_day836", "cMonocyte-S100A8_day36",
                  "intMonocyte_day8","intMonocyte_day36", "ncMonocyte-CD16_day8","ncMonocyte-CD16_day836")
DEGdf60$celltype_day <- factor(DEGdf60$celltype_day, levels = cell_day_lev)
  
# Create dotplot from selected genes
DEGdf_dotplot <- ggplot(DEGdf60, aes(reorder(genes, celltype_day), celltype_day))+
  geom_point(aes(colour=Log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="blue", mid = "white", high = "green", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
DEGdf_dotplot
ggsave(filename = "output_for_paper/pbDEG_Mono_APGen.pdf", h = 5,w =10, limitsize = FALSE)

##Selecting Apopt genes and only monocytes
DEGdfApopt <- psbDE[psbDE$genes %in% Apopt,]
unique(DEGdfApopt$celltype)
DEGdf60 <- DEGdfApopt[DEGdfApopt$celltype!="DCs1",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="DCs2",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="pDCs",]

unique(DEGdf60$celltype)
##Creating celltype factors
celltype_DEGdf<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16" )
##Creating day factors
DEGdf60$celltype <- factor(DEGdf60$celltype, levels = celltype_DEGdf)
DEGdf60$celltype_day <- paste0(DEGdf60$celltype, "_", DEGdf60$day)
cell_day_lev <- c("cMonocyte-COQ7_day8", "cMonocyte-COQ7_day836", "cMonocyte-COQ7_day36",
                  "cMonocyte-IL1B_day8", "cMonocyte-IL1B_day836", "cMonocyte-IL1B_day36",
                  "cMonocyte-ISG15_day8", "cMonocyte-ISG15_day836", "cMonocyte-ISG15_day36",
                  "cMonocyte-S100A8_day8", "cMonocyte-S100A8_day836", "cMonocyte-S100A8_day36",
                  "intMonocyte_day8","intMonocyte_day36", "ncMonocyte-CD16_day8","ncMonocyte-CD16_day836")
DEGdf60$celltype_day <- factor(DEGdf60$celltype_day, levels = cell_day_lev)

# Create dotplot from selected genes
DEGdf_dotplot <- ggplot(DEGdf60, aes(reorder(genes, celltype_day), celltype_day))+
  geom_point(aes(colour=Log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="blue", mid = "white", high = "green", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
DEGdf_dotplot
ggsave(filename = "output_for_paper/pbDEG_Mono_Apopt.pdf", h = 5,w =10, limitsize = FALSE)

##Selecting Chemo genes and only monocytes
DEGdfChemo <- psbDE[psbDE$genes %in% Chemo,]
unique(DEGdfChemo$celltype)
DEGdf60 <- DEGdfChemo[DEGdfChemo$celltype!="DCs1",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="DCs2",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="pDCs",]

unique(DEGdf60$celltype)
##Creating celltype factors
celltype_DEGdf<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16" )
##Creating day factors
DEGdf60$celltype <- factor(DEGdf60$celltype, levels = celltype_DEGdf)
DEGdf60$celltype_day <- paste0(DEGdf60$celltype, "_", DEGdf60$day)
cell_day_lev <- c("cMonocyte-COQ7_day8", "cMonocyte-COQ7_day836", "cMonocyte-COQ7_day36",
                  "cMonocyte-IL1B_day8", "cMonocyte-IL1B_day836", "cMonocyte-IL1B_day36",
                  "cMonocyte-ISG15_day8", "cMonocyte-ISG15_day836", "cMonocyte-ISG15_day36",
                  "cMonocyte-S100A8_day8", "cMonocyte-S100A8_day836", "cMonocyte-S100A8_day36",
                  "intMonocyte_day8","intMonocyte_day36", "ncMonocyte-CD16_day8","ncMonocyte-CD16_day836")
DEGdf60$celltype_day <- factor(DEGdf60$celltype_day, levels = cell_day_lev)

# Create dotplot from selected genes
DEGdf_dotplot <- ggplot(DEGdf60, aes(reorder(genes, celltype_day), celltype_day))+
  geom_point(aes(colour=Log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="blue", mid = "white", high = "green", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
DEGdf_dotplot
ggsave(filename = "output_for_paper/pbDEG_Mono_Chemo.pdf", h = 5,w =10, limitsize = FALSE)

##Selecting VEGF genes and only monocytes
DEGdfVEGF <- psbDE[psbDE$genes %in% VEGF,]
unique(DEGdfVEGF$celltype)
DEGdf60 <- DEGdfVEGF[DEGdfVEGF$celltype!="DCs1",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="DCs2",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="pDCs",]

unique(DEGdf60$celltype)
##Creating celltype factors
celltype_DEGdf<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16" )
##Creating day factors
DEGdf60$celltype <- factor(DEGdf60$celltype, levels = celltype_DEGdf)
DEGdf60$celltype_day <- paste0(DEGdf60$celltype, "_", DEGdf60$day)
cell_day_lev <- c("cMonocyte-COQ7_day8", "cMonocyte-COQ7_day836", "cMonocyte-COQ7_day36",
                  "cMonocyte-IL1B_day8", "cMonocyte-IL1B_day836", "cMonocyte-IL1B_day36",
                  "cMonocyte-ISG15_day8", "cMonocyte-ISG15_day836", "cMonocyte-ISG15_day36",
                  "cMonocyte-S100A8_day8", "cMonocyte-S100A8_day836", "cMonocyte-S100A8_day36",
                  "intMonocyte_day8","intMonocyte_day36", "ncMonocyte-CD16_day8","ncMonocyte-CD16_day836")
DEGdf60$celltype_day <- factor(DEGdf60$celltype_day, levels = cell_day_lev)

# Create dotplot from selected genes
DEGdf_dotplot <- ggplot(DEGdf60, aes(reorder(genes, celltype_day), celltype_day))+
  geom_point(aes(colour=Log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="blue", mid = "white", high = "green", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
DEGdf_dotplot
ggsave(filename = "output_for_paper/pbDEG_Mono_VEGF.pdf", h = 5,w =10, limitsize = FALSE)

##Selecting NOD genes and only monocytes
DEGdfNOD <- psbDE[psbDE$genes %in% NOD,]
unique(DEGdfNOD$celltype)
DEGdf60 <- DEGdfNOD[DEGdfNOD$celltype!="DCs1",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="DCs2",]
DEGdf60 <- DEGdf60[DEGdf60$celltype!="pDCs",]

unique(DEGdf60$celltype)
##Creating celltype factors
celltype_DEGdf<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16" )
##Creating day factors
DEGdf60$celltype <- factor(DEGdf60$celltype, levels = celltype_DEGdf)
DEGdf60$celltype_day <- paste0(DEGdf60$celltype, "_", DEGdf60$day)
cell_day_lev <- c("cMonocyte-COQ7_day8", "cMonocyte-COQ7_day836", "cMonocyte-COQ7_day36",
                  "cMonocyte-IL1B_day8", "cMonocyte-IL1B_day836", "cMonocyte-IL1B_day36",
                  "cMonocyte-ISG15_day8", "cMonocyte-ISG15_day836", "cMonocyte-ISG15_day36",
                  "cMonocyte-S100A8_day8", "cMonocyte-S100A8_day836", "cMonocyte-S100A8_day36",
                  "intMonocyte_day8","intMonocyte_day36", "ncMonocyte-CD16_day8","ncMonocyte-CD16_day836")
DEGdf60$celltype_day <- factor(DEGdf60$celltype_day, levels = cell_day_lev)

# Create dotplot from selected genes
DEGdf_dotplot <- ggplot(DEGdf60, aes(reorder(genes, celltype_day), celltype_day))+
  geom_point(aes(colour=Log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="blue", mid = "white", high = "green", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
DEGdf_dotplot
ggsave(filename = "output_for_paper/pbDEG_Mono_NOD.pdf", h = 5,w =10, limitsize = FALSE)
























