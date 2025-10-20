library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(readr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot)
library(celldex)
library(SingleR)
library(data.table)
library(scGate)
library(ggpubr)
library(scCustomize)
library(viridis)
library(ggalluvial)
library(SeuratExtend)
library(gridExtra)
library(ggrepel)

#uploading ovarian cancer samples from the GSE180661 dataset ( OV009-CD45P-LEFT-OVARY-r0; OV026-CD45P-RIGHT-OVARY-r0;OV067-CD45P-RIGHT-OVARY-AND-TUBE-r0;OV082-CD45P-LEFT-r0;OV112-CD45P-LEFT-r0)
#uploading lung cancer samples from the GSE131907 dataset (BRONCHO_58; EBUS_06; EBUS_49;LUNG_T28; LUNG_T31)
#Low-quality cells were excluded from each sample individually based on the following thresholds: number of detected features (nFeature) ranging between 200 and 700, and mitochondrial gene content exceeding 15% (see method in manuscript) 
#genes that were not common to both datasets were removed from the datasets
#The ovarian and lung datasets were merged and further analyzed
load("OV.LUNG_data.rda",verbose=T)


OV.LUNG<- NormalizeData(OV.LUNG)
OV.LUNG<- FindVariableFeatures(OV.LUNG)
OV.LUNG<- ScaleData(OV.LUNG)
OV.LUNG<- RunPCA(OV.LUNG)

# visualize the results of a standard analysis without integration

OV.LUNG<- FindNeighbors(OV.LUNG, dims = 1:30, reduction = "pca")
OV.LUNG<- FindClusters(OV.LUNG, resolution = 2, cluster.name = "unintegrated_clusters")
OV.LUNG<- RunUMAP(OV.LUNG, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
OV.LUNG<- IntegrateLayers(
  object = OV.LUNG, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
OV.LUNG[["RNA"]] <- JoinLayers(OV.LUNG[["RNA"]])

# RPCA
OV.LUNG<- FindNeighbors(OV.LUNG, reduction = "integrated.rpca", dims = 1:30)
OV.LUNG<- FindClusters(OV.LUNG, resolution = 2, cluster.name = "rpca_clusters")

OV.LUNG<- RunUMAP(OV.LUNG, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(
  OV.LUNG,
  reduction = "umap.rpca",
  group.by = c("rpca_clusters"),
  combine = FALSE, label.size = 2
)


#Cluster annotation using the CellDex and SingleR packages
sample<- as.SingleCellExperiment(OV.LUNG)
IGD<-BlueprintEncodeData()
pred.IGD<- SingleR(test = sample, ref =IGD, assay.type.test=1,labels = IGD$label.main)
OV.LUNG <- AddMetaData( object=OV.LUNG,metadata=make.names(sub( " [(].*$","",pred.IGD$pruned.labels )),col.name="BlueprintEncodeData")
IGD<-HumanPrimaryCellAtlasData()
pred.IGD<- SingleR(test = sample, ref =IGD, assay.type.test=1,labels = IGD$label.main)
OV.LUNG <- AddMetaData( object=OV.LUNG,metadata=make.names(sub( " [(].*$","",pred.IGD$pruned.labels )),col.name="HumanPrimaryCellAtlasData")
scGate_models_DB <- get_scGateDB()
scGate_models_DB
OV.LUNG<- scGate(OV.LUNG,reduction="umap.rpca", model = scGate_models_DB$human$TME_HiRes$NK)


Idents(OV.LUNG)="BlueprintEncodeData"

OV.LUNG@meta.data=OV.LUNG@meta.data%>%mutate(Population=case_when(BlueprintEncodeData=="CD8..T.cells"~"CD8T",
BlueprintEncodeData=="CD4..T.cells"~"CD4T",
BlueprintEncodeData=="B.cells"~"B.cells",
BlueprintEncodeData=="NK.cells"~"NK.cells",
BlueprintEncodeData=="DC"~"DC",
BlueprintEncodeData=="Macrophages"~"Macrophages",
BlueprintEncodeData=="Monocytes"~"Monocytes",
BlueprintEncodeData=="Neutrophils"~"Neutrophils",
BlueprintEncodeData=="Endothelial.cells"~"Endothelial.cells",
BlueprintEncodeData=="Epithelial.cells"~"Epithelial.cells",
BlueprintEncodeData=="Fibroblasts"~"Fibroblasts",
BlueprintEncodeData=="Adipocytes"~"Adipocytes",
BlueprintEncodeData=="Astrocytes"~"Others",
BlueprintEncodeData=="Eosinophils"~"Eosinophils",
BlueprintEncodeData=="Erythrocytes"~"Others",
BlueprintEncodeData=="HSC"~"Others",
BlueprintEncodeData=="Chondrocytes"~"Others",
BlueprintEncodeData=="Keratinocytes"~"Others",
BlueprintEncodeData=="Melanocytes"~"Others",
BlueprintEncodeData=="Mesangial.cells"~"Others",
BlueprintEncodeData=="Myocytes"~"Others",
BlueprintEncodeData=="Pericytes"~"Others",
BlueprintEncodeData=="Skeletal.muscle"~"Others",
BlueprintEncodeData=="Smooth.muscle"~"Others",
BlueprintEncodeData=="NA."~"Others"))

Idents(OV.LUNG)<-"Population"
OV.LUNG@meta.data$Population=factor(OV.LUNG@meta.data$Population,levels=c(
"CD8T",
"CD4T",
"B.cells",
"NK.cells",
"DC",
"Macrophages",
"Monocytes",
"Neutrophils",
"Endothelial.cells",
"Epithelial.cells",
"Fibroblasts",
"Adipocytes",
"Eosinophils",
"Others"))
OV.LUNG@meta.data$Group1=factor(OV.LUNG@meta.data$Group1,levels=c(
"OV","LUAD"))

ClusterDistrBar(OV.LUNG$Group1, OV.LUNG$Population,flip = F)+theme(axis.text.x = element_text(angle = 90))

Idents(OV.LUNG)<-"Population"
OV.LUNG@meta.data$Population=factor(OV.LUNG@meta.data$Population,levels=c(
"CD8T",
"CD4T",
"B.cells",
"NK.cells",
"DC",
"Macrophages",
"Monocytes",
"Neutrophils",
"Endothelial.cells",
"Epithelial.cells",
"Fibroblasts",
"Adipocytes",
"Eosinophils",
"Others"))

DimPlot2(OV.LUNG,reduction="umap.rpca",split.by="Group1",theme = NoLegend())

#----------------selection of NK cluster-------------
Idents(OV.LUNG)="rpca_clusters"
nk.subset=subset(OV.LUNG,idents=c("8","12"))
Idents(nk.subset)="BlueprintEncodeData"
nk.subset=subset(nk.subset,ident="NK.cells")
Idents(nk.subset)="is.pure"
nk.subset=subset(nk.subset,ident="Pure")
Idents(nk.subset)="HumanPrimaryCellAtlasData"
nk.subset=subset(nk.subset,ident="NK_cell")

DimPlot(nk.subset,reduction="umap.rpca")

nk.subset[["RNA"]] <- split(nk.subset[["RNA"]], f = nk.subset$Sample)

nk.subset<- NormalizeData(nk.subset)
nk.subset<- FindVariableFeatures(nk.subset)
nk.subset<- ScaleData(nk.subset)
nk.subset<- RunPCA(nk.subset)

ElbowPlot(nk.subset)

nk.subset<- FindNeighbors(nk.subset, dims = 1:15, reduction = "pca")
nk.subset<- FindClusters(nk.subset,resolution=0.8, cluster.name = "unintegrated_clusters")
nk.subset<- RunUMAP(nk.subset, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")


nk.subset<- IntegrateLayers(
  object = nk.subset, method = RPCAIntegration,dims = 1:15,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  k.weight=20,
  verbose = FALSE
)

nk.subset[["RNA"]] <- JoinLayers(nk.subset[["RNA"]])

# RPCA
nk.subset<- FindNeighbors(nk.subset, reduction = "integrated.rpca", dims = 1:15)
nk.subset<- FindClusters(nk.subset,resolution=0.8, cluster.name = "rpca_clusters")
nk.subset<- RunUMAP(nk.subset, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")

nk.subset@meta.data[["Cohort"]]=nk.subset@meta.data$Group1

#laod Meatastatic data

load("merged.seurat.rda", verbose=T)

#Primary tumor data

Idents(OV.LUNG)="Group1"
OV=subset(OV.LUNG,idents="OV")
counts <- GetAssayData(OV, slot = "counts")
OV.prim <- CreateSeuratObject(counts = counts)
OV.prim[["percent.mt"]] <- PercentageFeatureSet(OV.prim, pattern = "^MT-")

# Extract the last character of orig.ident
last_char <- substr(rownames(OV.prim@meta.data), nchar(rownames(OV.prim@meta.data)), nchar(rownames(OV.prim@meta.data)))

# Assign based on condition
OV.prim@meta.data$orig.ident <- ifelse(last_char == "1", "OV1",
                               ifelse(last_char == "2", "OV2",
                               ifelse(last_char == "3", "OV3",
                               ifelse(last_char == "4", "OV4",
                               ifelse(last_char == "5", "OV5", NA)))))

#compare gene symbols between OV.prim and Metastatic cohort

DF=all(rownames(merged_seurat[["RNA"]]$counts)%in% rownames(OV.prim[["RNA"]]$counts))
# Find common features between the two objects
intersect=intersect(rownames(merged_seurat[["RNA"]]$counts),rownames(OV.prim[["RNA"]]$counts))
PRIM.MET=merge(merged_seurat[intersect, ], OV.prim[intersect, ])


PRIM.MET[["RNA"]] <- split(PRIM.MET[["RNA"]], f = PRIM.MET$orig.ident)

PRIM.MET<- NormalizeData(PRIM.MET)
PRIM.MET<- FindVariableFeatures(PRIM.MET)
PRIM.MET<- ScaleData(PRIM.MET)
PRIM.MET<- RunPCA(PRIM.MET)

ElbowPlot(PRIM.MET)

PRIM.MET<- FindNeighbors(PRIM.MET, dims = 1:20, reduction = "pca")
PRIM.MET<- FindClusters(PRIM.MET, cluster.name = "unintegrated_clusters")
PRIM.MET<- RunUMAP(PRIM.MET, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")


PRIM.MET<- IntegrateLayers(
  object = PRIM.MET, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

PRIM.MET[["RNA"]] <- JoinLayers(PRIM.MET[["RNA"]])

# Harmony
PRIM.MET<- FindNeighbors(PRIM.MET, reduction = "harmony", dims = 1:20)
PRIM.MET<- FindClusters(PRIM.MET, cluster.name = "harmony_clusters")
PRIM.MET<- RunUMAP(PRIM.MET, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

Idents(PRIM.MET)<-"orig.ident"
PRIM.MET@meta.data=PRIM.MET@meta.data%>%mutate(Cohort=case_when(orig.ident=="S1"~"omentum_met_distant",
orig.ident=="S2"~"omentum_met_peritumoral",
orig.ident=="S3"~"omentum_met_tumor",
orig.ident=="S4"~"omentum_met_peritumoral",
orig.ident=="S5"~"omentum_met_distant",
orig.ident=="S6"~"omentum_met_peritumoral",
orig.ident=="S7"~"omentum_met_tumor",
orig.ident=="S8"~"omentum_met_distant",
orig.ident=="S9"~"omentum_met_peritumoral",
orig.ident=="S10"~"omentum_met_tumor",
orig.ident=="S11"~"omentum_met_distant",
orig.ident=="S12"~"omentum_met_peritumoral",
orig.ident=="S13"~"omentum_met_tumor",
orig.ident=="S14"~"omentum_met_distant",
orig.ident=="S15"~"omentum_met_peritumoral",
orig.ident=="S16"~"omentum_met_tumor",
orig.ident=="S17"~"omentum_met_distant",
orig.ident=="S18"~"omentum_met_peritumoral",
orig.ident=="S19"~"omentum_met_tumor",
orig.ident=="S20"~"benign_omentum",
orig.ident=="S21"~"benign_omentum",
orig.ident=="S22"~"benign_omentum",
orig.ident=="S23"~"benign_omentum",
orig.ident=="S24"~"benign_omentum",
orig.ident=="S25"~"benign_omentum",
orig.ident=="S26"~"benign_omentum",
orig.ident=="S27"~"benign_omentum",
orig.ident=="S28"~"benign_omentum",
orig.ident=="S29"~"benign_omentum",
orig.ident=="S30"~"benign_omentum",
orig.ident=="S31"~"benign_omentum",
orig.ident=="S32"~"benign_omentum",
orig.ident=="S33"~"benign_omentum",
orig.ident=="S34"~"benign_omentum",
orig.ident=="S35"~"benign_omentum",
orig.ident=="S36"~"benign_omentum",
orig.ident=="OV1"~"Primary_tumor",
orig.ident=="OV2"~"Primary_tumor",
orig.ident=="OV3"~"Primary_tumor",
orig.ident=="OV4"~"Primary_tumor",
orig.ident=="OV5"~"Primary_tumor"))


# CellDex
sample<- as.SingleCellExperiment(PRIM.MET)
IGD<-celldex::BlueprintEncodeData()
pred.IGD<- SingleR(test = sample, ref =IGD, assay.type.test=1,labels = IGD$label.main)
PRIM.MET<- AddMetaData( object=PRIM.MET,metadata=make.names(sub( " [(].*$","",pred.IGD$pruned.labels )),col.name="BlueprintEncodeData" )

Idents(PRIM.MET)="BlueprintEncodeData"
DimPlot(PRIM.MET,reduction="umap.harmony",label=T)


Idents(PRIM.MET)<-"BlueprintEncodeData"
PRIM.MET@meta.data=PRIM.MET@meta.data%>%mutate(Population=case_when(BlueprintEncodeData=="CD8..T.cells"~"CD8T",
BlueprintEncodeData=="CD4..T.cells"~"CD4T",
BlueprintEncodeData=="B.cells"~"B.cells",
BlueprintEncodeData=="NK.cells"~"NK.cells",
BlueprintEncodeData=="DC"~"DC",
BlueprintEncodeData=="Macrophages"~"Macrophages",
BlueprintEncodeData=="Monocytes"~"Monocytes",
BlueprintEncodeData=="Neutrophils"~"Neutrophils",
BlueprintEncodeData=="Endothelial.cells"~"Endothelial.cells",
BlueprintEncodeData=="Epithelial.cells"~"Epithelial.cells",
BlueprintEncodeData=="Fibroblasts"~"Fibroblasts",
BlueprintEncodeData=="Adipocytes"~"Adipocytes",
BlueprintEncodeData=="Astrocytes"~"Others",
BlueprintEncodeData=="Eosinophils"~"Eosinophils",
BlueprintEncodeData=="Erythrocytes"~"Others",
BlueprintEncodeData=="HSC"~"Others",
BlueprintEncodeData=="Chondrocytes"~"Others",
BlueprintEncodeData=="Keratinocytes"~"Others",
BlueprintEncodeData=="Melanocytes"~"Others",
BlueprintEncodeData=="Mesangial.cells"~"Others",
BlueprintEncodeData=="Myocytes"~"Others",
BlueprintEncodeData=="Pericytes"~"Others",
BlueprintEncodeData=="Skeletal.muscle"~"Others",
BlueprintEncodeData=="Smooth.muscle"~"Others",
BlueprintEncodeData=="NA."~"Others"))

Idents(PRIM.MET)<-"harmony_clusters"
NK=subset(PRIM.MET,idents="6")
Idents(NK)="BlueprintEncodeData"
NK=subset(NK,idents="NK.cells")

#Re-clustering of NK subset

NK[["RNA"]] <- split(NK[["RNA"]], f = NK$orig.ident)

NK<- NormalizeData(NK)
NK<- FindVariableFeatures(NK)
NK<- ScaleData(NK)
NK<- RunPCA(NK)

ElbowPlot(NK)

NK<- FindNeighbors(NK, dims = 1:20, reduction = "pca")
NK<- FindClusters(NK,resolution=0.3, cluster.name = "unintegrated_clusters")
NK<- RunUMAP(NK, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")


NK<- IntegrateLayers(
  object = NK, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

NK[["RNA"]] <- JoinLayers(NK[["RNA"]])

# Harmony
NK<- FindNeighbors(NK, reduction = "harmony", dims = 1:20)
NK<- FindClusters(NK,resolution=0.3, cluster.name = "harmony_clusters")
NK<- RunUMAP(NK, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")


#compare gene symbols between OV_LUNG and Metastatic cohort

DF=all(rownames(nk.subset[["RNA"]]$counts)%in% rownames(NK[["RNA"]]$counts))
# Find common features between the two objects
intersect=intersect(rownames(nk.subset[["RNA"]]$counts),rownames(NK[["RNA"]]$counts))
ALL=merge(nk.subset[intersect, ], NK[intersect, ])

ALL@meta.data <- ALL@meta.data[, !(colnames(ALL@meta.data) %in% c("unintegrated_clusters", "rpca_clusters","harmony_clusters","NK_UCell",
"NK.group","Group","seurat_clusters","cca_clusters","Group1","HumanPrimaryCellAtlasData","is.pure",
"NK.cells","Group2"))]

ALL[["RNA"]] <- split(ALL[["RNA"]], f = ALL$Sample)

ALL<- NormalizeData(ALL)
ALL<- FindVariableFeatures(ALL)
ALL<- ScaleData(ALL)
ALL<- RunPCA(ALL)

ElbowPlot(ALL)

ALL<- FindNeighbors(ALL, dims = 1:10, reduction = "pca")
ALL<- FindClusters(ALL,resolution = 0.3, cluster.name = "unintegrated_clusters")
ALL<- RunUMAP(ALL, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")


ALL<- IntegrateLayers(
  object = ALL, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE,dims = 1:10
)

ALL[["RNA"]] <- JoinLayers(ALL[["RNA"]])

# Harmony
ALL<- FindNeighbors(ALL, reduction = "harmony", dims = 1:15)
ALL<- FindClusters(ALL, cluster.name = "harmony_clusters")
ALL<- RunUMAP(ALL, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")

signatures <- list(
NK1=c("SPON2","NKG7","FCER1G","PRF1","GZMB","CLIC3","MYOM2","FGFBP2","AKR1C3","IGFBP7","LAIR2","KLRB1","CX3CR1","CHST2","CD247","CD38","ADGRG1",
"CCL4","CTSD","CD160"),
NK2=c("GZMK","XCL1","XCL2","TPT1","TCF7","KLRC1","CMC1","GPR183","SELL","CD44","DUSP2","ZFP36L2","AREG","PIK3R1","IL7R","COTL1","IFITM3","IL2RB","FOS",
"LTB"),
NK3=c("HBA1",	"CD3G",	"CD3D",	"DSTN",	"ZBTB38",
"CD2",	"PPDPF",	"PRDM1",	"LGALS1",	"S100A4",	"ITGB1",
"IL32",	"VIM",	"LINC01871",	"S100A6",	"PTMS",	"GZMH",
"CD3E",	"CCL5",	"KLRC2"))

library(UCell)

ALL<- AddModuleScore_UCell(ALL,features=signatures, name=NULL)

Idents(NK)<-"harmony_clusters"
ALL@meta.data=ALL@meta.data%>%mutate(NK.pop=case_when(harmony_clusters=="1"~"NK1",
harmony_clusters=="10"~"NK1",
harmony_clusters=="6"~"NK2",
harmony_clusters=="2"~"NK2",
harmony_clusters=="8"~"NK2",
harmony_clusters=="5"~"NK2",
harmony_clusters=="4"~"NK2",
harmony_clusters=="0"~"NK3",
harmony_clusters=="11"~"NK3",
harmony_clusters=="3"~"NK3",
harmony_clusters=="9"~"NK3",
harmony_clusters=="14"~"Unknown",
harmony_clusters=="15"~"Unknown",
harmony_clusters=="13"~"Unknown",
harmony_clusters=="7"~"Unknown",
harmony_clusters=="12"~"Unknown",
harmony_clusters=="16"~"Unknown"))

#remove "unknown"
Idents(ALL)<-"NK.pop"
ALL.clean=subset(ALL,idents="Unknown",invert=T)

DimPlot(ALL.clean,reduction="umap.harmony")

counts=ALL.clean@meta.data%>%group_by(Cohort,NK.pop)%>%summarize(no=n())%>% mutate(percent = 100 * no / sum(no))

counts$Cohort <- factor(counts$Cohort, levels = c("LUAD","OV","benign_omentum","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))


p1=ggplot(counts, aes(x = Cohort, y = percent, fill = NK.pop)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(percent, 1)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  labs(y = "Percentage", x = "Group", fill = "NK Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))


p2=ggplot(counts, aes(x = Cohort, y = no, fill = NK.pop)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(no, 1)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  labs(y = "No", x = "Group", fill = "NK Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))

p1+p2

#------------------DotPlot Fig.1J----------------------
genes <- c("SPON2","NKG7","FCER1G","PRF1","GZMB","CLIC3","MYOM2","FGFBP2","AKR1C3","IGFBP7","LAIR2","KLRB1","CX3CR1","CHST2","CD247","CD38","ADGRG1",
"CCL4","CTSD","CD160","GZMK","XCL1","XCL2","TPT1","TCF7","KLRC1","CMC1","GPR183","SELL","CD44","DUSP2","ZFP36L2","AREG","PIK3R1","IL7R","COTL1","IFITM3","IL2RB","FOS",
"LTB","KLRC2","CCL5","CD3E","GZMH","PTMS","S100A6","LINC01871","VIM","IL32","ITGB1","S100A4","LGALS1","PRDM1","PPDPF","CD2","ZBTB38","DSTN","CD3D","CD3G","HBA1")

Idents(nk.subset)="NK.pop"
levels(nk.subset) <- c("NK1", "NK2","NK3")

pdf("DotPlot_NK.pop.pdf",width=5,height=12)
DotPlot(nk.subset, features = genes)+scale_colour_gradientn(
    colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
    breaks = -2:2
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


dev.off()

Idents(nk.subset)="Group1"
levels(nk.subset) <- c("OV", "LUAD")

pdf("DotPlot_Cohort.pdf",width=5,height=12)
DotPlot(nk.subset, features = genes)+scale_colour_gradientn(
    colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
    breaks = -2:2
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


dev.off()



#-----------------------Heatmap Fig.1K--------------------------------------

d= c("IFNG", "CCL4","CCL4L2","CCL5","IL32","IL16","CCL3","XCL2","XCL1","FLT3LG",
"KLRK1", "KLRC2","CD160","NCR3","SLAMF6","SLAMF7",
"KLRC1","CD300A","TIGIT","KIR3DL2","KIR3DL1","KIR2DL3","KIR2DL1","KLRB1","SIGLEC7","HAVCR2",
"TGFBR1","IL12RB1","IL2RG","TGFBR3","TGFBR2","IL10RA","IL18R1","IL18RAP","IL2RB","IL10RB",
"FASLG","GSDMD","GZMB","NKG7","PRF1","GZMA","GZMH","GZMK","TNFSF10")

nk.subset@meta.data[["Group.f"]]=paste0(nk.subset@meta.data$Group1,"_",nk.subset@meta.data$NK.pop)

Idents(nk.subset)="Group.f"
levels(nk.subset) <- c("OV_NK1", "LUAD_NK1","OV_NK2","LUAD_NK2","OV_NK3","LUAD_NK3")

cal <- CalcStats(nk.subset, features = d, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)


library(circlize)
library(ComplexHeatmap)
library(pheatmap)

row_annot=read.table(file="cal_row_split.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

# Clean group names in case of spaces
row_annot$Groups <- trimws(row_annot$Groups)


# Desired block order
block_levels <- c(
  "Cytotoxicity",
  "Cytokines",
  "Cytokines_receptors",
  "Activatory_receptors",
  "Inhibitory_receptors"
)

Rowsplit <- factor(row_annot$Groups, levels = block_levels)

# Cluster within each block but keep block order
order_within_blocks <- unlist(lapply(block_levels, function(block) {
  rows_in_block <- which(Rowsplit == block)
  if (length(rows_in_block) > 1) {
    # distance + hierarchical clustering inside this block
    dist_mat <- dist(mat[rows_in_block, ])
    hc <- hclust(dist_mat, method = "ward.D2")
    rows_in_block[hc$order]
  } else {
    rows_in_block
  }
}))

mat_clustered <- mat[order_within_blocks, ]
Rowsplit_clustered <- Rowsplit[order_within_blocks]

# Plot without re-clustering rows

pdf("Fig1K.pdf",width=6,height=12)

ComplexHeatmap::pheatmap(
  mat_clustered,
  row_split = Rowsplit_clustered,
  cluster_rows = FALSE, # we already clustered manually
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)
dev.off()


write.table(cal,file="cal_Fig1K.txt",sep="\t",col.names=NA)

#selection of normalize gene expression

gene_list <- c("IFNG", "CCL4","CCL4L2","CCL5","IL32","IL16","CCL3","XCL2","XCL1","FLT3LG",
"KLRK1", "KLRC2","CD160","NCR3","SLAMF6","SLAMF7",
"KLRC1","CD300A","TIGIT","KIR3DL2","KIR3DL1","KIR2DL3","KIR2DL1","KLRB1","SIGLEC7","HAVCR2",
"TGFBR1","IL12RB1","IL2RG","TGFBR3","TGFBR2","IL10RA","IL18R1","IL18RAP","IL2RB","IL10RB",
"FASLG","GSDMD","GZMB","NKG7","PRF1","GZMA","GZMH","GZMK","TNFSF10")

data_subset <- GetAssayData(nk.subset, slot = "data")[gene_list, , drop = FALSE]
data=data.frame(t(data_subset))
data$Group=paste0(nk.subset@meta.data$Group.f)

write.table(data,file="normilze_data_Fig1K.txt",sep="\t",col.names=NA)

#--------------------------------Fig1L--Shell---------------
library(clusterProfiler)
library(org.Hs.eg.db)

Idents(nk.subset)="NK.pop"
markers <- FindAllMarkers(nk.subset, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25,
                          min.pct = 0.25)

markers_entrez <- bitr(markers$gene,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
markers <- merge(markers, markers_entrez, by.x = "gene", by.y = "SYMBOL")

go_results <- lapply(unique(markers$cluster), function(clust) {
  genes <- markers %>% filter(cluster == clust) %>% pull(ENTREZID)
  
  ego <- enrichGO(gene = genes,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",        # Can be "MF", "CC" or "BP"
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  ego@result$cluster <- clust
  return(ego@result)
})

go_df <- do.call(rbind, go_results)

#upload data FoldEnrichment

FE=fread("FoldEnrichment.txt")

ggplot(FE, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = as.factor(cluster))) +
  geom_bar(stat = "identity",position="dodge") +
  coord_polar(start = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),   # remove x-axis text
    axis.text.y = element_blank(),   # remove y-axis text
    axis.title.x = element_blank(),  # remove x-axis title
    axis.title.y = element_blank(),  # remove y-axis title
    legend.title = element_blank(),  # remove legend title
    legend.text = element_blank()    # remove legend text
  )  +scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
))+
  labs(
    fill = "Cluster"
  )

#upload data p.adjust

PA=fread("p.adjust.txt")

ggplot(PA, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = as.factor(cluster))) +
  geom_bar(stat = "identity",position="dodge") +
  coord_polar(start = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text( vjust = 0.5, hjust=1, size = 8),
    
    legend.position = "right"
  ) +scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
))+
  labs(
    fill = "Cluster"
  )

data=read.table(file="top10_go_df_update.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

ggplot(data, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = as.factor(cluster))) +
  geom_bar(stat = "identity",position="dodge") +
  coord_polar(start = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text( vjust = 0.5, hjust=1, size = 8),
    
    legend.position = "right"
  ) +scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))+
  labs(
    fill = "Cluster",
    title = "Top 10 GO Terms per Cluster (Circular Barplot)",
    subtitle = "-log10 adjusted p-values"
  )

#without text

ggplot(data, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = as.factor(cluster))) +
  geom_bar(stat = "identity",position="dodge") +
  coord_polar(start = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # hides x-axis text
    axis.text.y = element_blank(),
    
    legend.position = "right"
  ) +
scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))+
  labs(
    fill = "Cluster",
    title = NULL,                   # removes title
    subtitle = NULL
  )


#----------------CellChat NK+CD8-----------------
#Figure 1L 
library(CellChat)
library(ProjecTILs)

load("CD8.subset.rda")
CD8.subset[["RNA3"]] <- as(object = CD8.subset[["RNA"]], Class = "Assay")
ref <- load.reference.map(ref = "CD8T_human_ref_v1.rds")
palette <- ref@misc$atlas.palette
DimPlot(ref, label = T)
DefaultAssay(CD8.subset) <- "RNA3"

sample.list <- SplitObject(CD8.subset, split.by = "Group1")
names(sample.list)

sample.class <- list()

sample.class$OV <- ProjecTILs.classifier(query = sample.list$OV, ref = ref,
    filter.cells = F)

sample.class$LUAD <- ProjecTILs.classifier(query = sample.list$LUAD, ref = ref,
    filter.cells = F)

a <- DimPlot(sample.class$OV,reduction = "umap.rpca", group.by = "functional.cluster", cols = palette) +
    theme(aspect.ratio = 1) + ggtitle("OV")
b <- DimPlot(sample.class$LUAD,reduction = "umap.rpca", group.by = "functional.cluster", cols = palette) +
    theme(aspect.ratio = 1) + ggtitle("LUAD")
a | b

c <- plot.statepred.composition(ref = ref, query = sample.class$OV, metric = "Percent") +
    ggtitle("OV") + ylim(0, 60)
d <- plot.statepred.composition(ref = ref, query = sample.class$LUAD, metric = "Percent") +
    ggtitle("LUAD") + ylim(0, 60)

obj.projected <- Run.ProjecTILs(query = sample.list, ref = ref, filter.cells = FALSE)

e <- plot.projection(ref = ref, query = obj.projected$OV, linesize = 0.5) +
    theme(aspect.ratio = 1) + ggtitle("OV")
f <- plot.projection(ref = ref, query = obj.projected$LUAD, linesize = 0.5) +
    theme(aspect.ratio = 1) + ggtitle("LUAD")
library(gridExtra)
grid.arrange(a,b,c,d,e,f, nrow = 3)

# Alluvial plot of ProjecTILs.classifier data
library(scales)

OV=sample.class$OV@meta.data%>%filter(functional.cluster !="CD8.MAIT")%>%group_by(functional.cluster)%>%summarize(no=n())%>%mutate(Percentage= 100*no/sum(no))%>%mutate(Cohort= paste("OV"))
LUAD=sample.class$LUAD@meta.data%>%filter(functional.cluster !="CD8.MAIT")%>%group_by(functional.cluster)%>%summarize(no=n())%>%mutate(Percentage= 100*no/sum(no))%>%mutate(Cohort= paste("LUAD"))

data=data.frame(rbind(OV,LUAD))

data$Cohort <- factor(data$Cohort, levels=c("OV", "LUAD"))
data$functional.cluster <- factor(data$functional.cluster, levels=c("CD8.CM", "CD8.EM","CD8.TEMRA","CD8.TPEX","CD8.TEX","CD8.NaiveLike"))

show_col(hue_pal()(6))

ggplot(data,
       aes(x = Cohort, stratum = functional.cluster, alluvium = functional.cluster,y=Percentage,
           fill = functional.cluster, label = functional.cluster)) +
  scale_fill_manual(values = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  theme(legend.position = "bottom")+theme_classic()


CD8=merge(sample.class$OV,y=sample.class$LUAD,add.cell.ids=c("OV","LUAD"))
load("NK.subset_new1.rda",verbose=T)
CD8.NK=merge(CD8,y=nk.subset)

CD8.NK$Group_final <- case_when(
  grepl("CD8.TEX",CD8.NK@meta.data$functional.cluster) ~ "CD8.TEX",
  grepl("CD8.CM", CD8.NK@meta.data$functional.cluster) ~ "CD8.CM",
  grepl("CD8.TPEX", CD8.NK@meta.data$functional.cluster) ~ "CD8.TPEX",
  grepl("CD8.TEMRA", CD8.NK@meta.data$functional.cluster) ~ "CD8.TEMRA",
  grepl("CD8.EM", CD8.NK@meta.data$functional.cluster) ~ "CD8.EM",
  grepl("CD8.NaiveLike", CD8.NK@meta.data$functional.cluster) ~ "CD8.NaiveLike",
  grepl("CD8.MAIT", CD8.NK@meta.data$functional.cluster) ~ "CD8.MAIT",
  grepl("NK1", CD8.NK@meta.data$NK.pop) ~ "NK1",
  grepl("NK2", CD8.NK@meta.data$NK.pop) ~ "NK2",
  grepl("NK3", CD8.NK@meta.data$NK.pop) ~ "NK3",
    TRUE ~ NA_character_
)

#removing unused columns
CD8.NK@meta.data <- CD8.NK@meta.data[, !(colnames(CD8.NK@meta.data) %in% c("STAGE1", "STAGE", "my_selection","Group","NK_UCell","is.pure","Sample","NK.cells","NK.group","NK1","NK2","NK3","combined_clusters","NK.pop"))]

Idents(CD8.NK)="Group_final"
CD8.NK=subset(CD8.NK,subset =Group_final != "CD8.MAIT")
save(CD8.NK,file="CD8.NK.rda")
data.list <- SplitObject(CD8.NK, split.by = "Group1")

OV=data.list$OV

LUAD=data.list$LUAD
#----------------------OV-------------
Idents(OV)="Group_final"
OV[["RNA3"]] <- as(object = OV[["RNA"]], Class = "Assay")

data.input = OV[["RNA3"]]$data
meta = data.frame( labels=OV@meta.data$Group_final,row.names=colnames(data.input) )

cellchat.OV<- createCellChat( object=data.input,meta=meta,group.by="labels" )
cellchat.OV<- addMeta(cellchat.OV, meta = meta)
cellchat.OV<- setIdent(cellchat.OV, ident.use = "labels") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat.OV@idents))
CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)
#  use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). CellChatDB.use <- subsetDB(CellChatDB)
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat.OV@DB <- CellChatDB.use
cellchat.OV<- subsetData(cellchat.OV)
options(future.globals.maxSize = 16000 * 1024^2)
cellchat.OV<- identifyOverExpressedGenes(cellchat.OV)
cellchat.OV<- identifyOverExpressedInteractions(cellchat.OV)
cellchat.OV<- projectData(cellchat.OV, PPI.human)
cellchat.OV<- computeCommunProb(cellchat.OV, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.OV<- filterCommunication(cellchat.OV, min.cells = 50)
cellchat.OV<- computeCommunProbPathway(cellchat.OV)
cellchat.OV<- aggregateNet(cellchat.OV)

mat <- cellchat.OV@net$weight
pdf( "OV--circlePerCellPopulation.pdf",width=20,height=25 )
 par(mfrow = c(4,5), xpd=TRUE)
 for ( i in 1:nrow(mat) ) {
  print(i)
  mat2      <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
 }
dev.off()

#----------------------------------CellChat_LUNG cancer --------------------
Idents(LUAD)="Group_final"
LUAD[["RNA3"]] <- as(object = LUAD[["RNA"]], Class = "Assay")

data.input = LUAD[["RNA3"]]$data
meta = data.frame( labels=LUAD@meta.data$Group_final,row.names=colnames(data.input) )

cellchat.LUAD<- createCellChat( object=data.input,meta=meta,group.by="labels" )
cellchat.LUAD<- addMeta(cellchat.LUAD, meta = meta)
cellchat.LUAD<- setIdent(cellchat.LUAD, ident.use = "labels") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat.LUAD@idents))
CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)
#  use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). CellChatDB.use <- subsetDB(CellChatDB)
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat.LUAD@DB <- CellChatDB.use
cellchat.LUAD<- subsetData(cellchat.LUAD)
options(future.globals.maxSize = 16000 * 1024^2)
cellchat.LUAD<- identifyOverExpressedGenes(cellchat.LUAD)
cellchat.LUAD<- identifyOverExpressedInteractions(cellchat.LUAD)
cellchat.LUAD<- projectData(cellchat.LUAD, PPI.human)
cellchat.LUAD<- computeCommunProb(cellchat.LUAD, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.LUAD<- filterCommunication(cellchat.LUAD, min.cells = 50)
cellchat.LUAD<- computeCommunProbPathway(cellchat.LUAD)
cellchat.LUAD<- aggregateNet(cellchat.LUAD)

mat <- cellchat.LUAD@net$weight
pdf( "LUAD--circlePerCellPopulation.pdf",width=20,height=25 )
 par(mfrow = c(4,5), xpd=TRUE)
 for ( i in 1:nrow(mat) ) {
  print(i)
  mat2      <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
 }
dev.off()

#counts Dotplot


OV=data.frame(cellchat.OV@net$count)
Numbers.OV=OV[c(7:9),-c(7,8,9)]

LUNG=data.frame(cellchat.LUAD@net$count)
Numbers.LUNG=LUNG[c(7:9),-c(7,8,9)]

combined <- bind_rows(Numbers.OV, Numbers.LUNG)
rownames(combined)=c("NK1_OV","NK2_OV","NK3_OV","NK1_LUNG","NK2_LUNG","NK3_LUNG")
combined$Cohort <- c("OV", "OV","OV", "LUNG","LUNG","LUNG")
combined$cells <- c("NK1", "NK2","NK3", "NK1","NK2","NK3")
combined=as.data.table(combined)
inp.melt.dt <-
 melt.data.table(
  data          = combined,
  id.vars       = c("Cohort","cells"),
  measure.vars  = c("CD8.CM","CD8.EM","CD8.NaiveLike","CD8.TEMRA","CD8.TEX","CD8.TPEX"),
  variable.name = "Population"
 )

write.table(inp.melt.dt,file="inp.melt.dt.txt",sep="\t",col.names=NA)

inp.melt.dt$Cohort_Cell <- paste(inp.melt.dt$Cohort, inp.melt.dt$cells, sep = "_")

inp.melt.dt$Cohort_Cell=factor(inp.melt.dt$Cohort_Cell,levels=c("OV_NK1","LUNG_NK1","OV_NK2","LUNG_NK2","OV_NK3","LUNG_NK3"))
inp.melt.dt$Population=factor(inp.melt.dt$Population,levels=rev(c("CD8.NaiveLike","CD8.TEX","CD8.TPEX","CD8.TEMRA","CD8.EM","CD8.CM")))


pdf("counts.pdf",width=6,height=4)
ggplot(inp.melt.dt, aes(x= Cohort_Cell, y=Population, size=value, group=Cohort)) + geom_point(aes(colour=value)) + 
scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu")))+
theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+coord_flip()
dev.off()




load("ALL.clean.rda")
Idents(ALL.clean)="Cohort"
nk.subset=subset(ALL.clean,ident=c("OV","LUAD"))

nk.subset[["RNA"]] <- split(nk.subset[["RNA"]], f = nk.subset$Sample)

nk.subset<- NormalizeData(nk.subset)
nk.subset<- FindVariableFeatures(nk.subset)
nk.subset<- ScaleData(nk.subset)
nk.subset<- RunPCA(nk.subset)

ElbowPlot(nk.subset)

nk.subset<- FindNeighbors(nk.subset, dims = 1:15, reduction = "pca")
nk.subset<- FindClusters(nk.subset,resolution=0.3, cluster.name = "unintegrated_clusters")
nk.subset<- RunUMAP(nk.subset, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

options(future.globals.maxSize = 20000 * 1024^4)

nk.subset<- IntegrateLayers(
  object = nk.subset, method = RPCAIntegration,dims = 1:15,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  k.weight=15,k.anchor = 15,k.filter = 15,k.score = 15,
  verbose = FALSE
)

nk.subset[["RNA"]] <- JoinLayers(nk.subset[["RNA"]])

# RPCA
nk.subset<- FindNeighbors(nk.subset, reduction = "integrated.rpca", dims = 1:15)
nk.subset<- FindClusters(nk.subset,resolution=0.3, cluster.name = "rpca_clusters")
nk.subset<- RunUMAP(nk.subset, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")

Idents(nk.subset)="rpca_clusters"
DimPlot(nk.subset,reduction="umap.rpca")


counts=nk.subset@meta.data%>%group_by(Cohort,NK.pop)%>%summarize(no=n())%>% mutate(percent = 100 * no / sum(no))

counts$Cohort <- factor(counts$Cohort, levels = c("OV","LUAD"))


p1=ggplot(counts, aes(x = Cohort, y = percent, fill = NK.pop)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(percent, 1)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  labs(y = "Percentage", x = "Group", fill = "NK Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))

p1

#----------------------Supplementary files Fig1-------------------
library(UCell)
load("NK.subset_new1.rda")

signatures <- list(
NK1=c("SPON2","NKG7","FCER1G","PRF1","GZMB","CLIC3","MYOM2","FGFBP2","AKR1C3","IGFBP7","LAIR2","KLRB1","CX3CR1","CHST2","CD247","CD38","ADGRG1",
"CCL4","CTSD","CD160"),
NK2=c("GZMK","XCL1","XCL2","TPT1","TCF7","KLRC1","CMC1","GPR183","SELL","CD44","DUSP2","ZFP36L2","AREG","PIK3R1","IL7R","COTL1","IFITM3","IL2RB","FOS",
"LTB"),
NK3=c("HBA1",	"CD3G",	"CD3D",	"DSTN",	"ZBTB38",
"CD2",	"PPDPF",	"PRDM1",	"LGALS1",	"S100A4",	"ITGB1",
"IL32",	"VIM",	"LINC01871",	"S100A6",	"PTMS",	"GZMH",
"CD3E",	"CCL5",	"KLRC2"))


nk.subset<- AddModuleScore_UCell(nk.subset,features=signatures, name=NULL)

Idents(nk.subset)="NK3"

FeaturePlot(nk.subset,features="NK3",reduction="umap.rpca")& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


library(SCpubr)
library(Nebulosa)
Nebulosa::plot_density(nk.subset,reduction="umap.rpca",
                            features = "NK3")
#--------------------BarPlot-------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(ggpubr)
library(purrr)

load("NK.subset_new1.rda")
# 1. Calculate counts and percentages per sample
replicate_counts <- nk.subset@meta.data %>%
  group_by(Sample, Group1, NK.pop) %>%
  summarise(no = n(), .groups = "drop") %>%
  group_by(Sample, Group1) %>%
  mutate(percent = 100 * no / sum(no)) %>%
  ungroup()

# 2. Run t-tests per NK.pop comparing Group1
ttest_results <- replicate_counts %>%
  group_by(NK.pop) %>%
  summarise(
    wilcox_test = list(wilcox.test(percent ~ Group1, data = cur_data())),
    .groups = "drop"
  ) %>%
  mutate(tidy = map(wilcox_test, broom::tidy)) %>%
  unnest(tidy)

print(ttest_results)


replicate_counts$Group1 <- factor(replicate_counts$Group1, levels = c("LUAD", "OV"))

# 4. Plot boxplots with significance
p <- ggplot(replicate_counts, aes(x = Group1, y = percent, fill = Group1)) +
  geom_boxplot() +
  facet_wrap(~ NK.pop, scales = "free_y") +
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. Add significance labels
# For each facet, add label above max boxplot value
library(ggplot2)
library(dplyr)

sig_labels <- ttest_results %>%
  select(NK.pop, p.value) %>%
  mutate(label = paste0("p = ", signif(p.value, digits = 3)))  # round p-value to 3 significant digits

sig_labels <- ttest_results %>%
  select(NK.pop, p.value) %>%
  mutate(label = paste0("p = ", signif(p.value, digits = 3))) %>%
  left_join(
    replicate_counts %>%
      group_by(NK.pop) %>%
      summarise(y_pos = max(percent) * 1.05),
    by = "NK.pop"
  )

# Add p-value labels to plot
p + geom_text(data = sig_labels, aes(x = 1.5, y = y_pos, label = label),
              inherit.aes = FALSE, size = 5)


#-----------------------BarPlot of gene expression-------------------------

load("NK.subset_new1.rda",verbose=T)
nk.subset@meta.data[["GROUP"]]=paste0(nk.subset@meta.data$Group1,"_",nk.subset@meta.data$NK.pop)

#Cytokines

gene_list <- c("IFNG",	"CCL4",	"CCL4L2",	"CCL5",	"IL32",	"IL16",	"CCL3",	"XCL2",	"XCL1",	"FLT3LG")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "OV_NK1"   # First group
group2 <- "LUAD_NK1" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "GROUP"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK1_Cytokines.txt",sep="\t",col.names=NA)


markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

write.table(markers ,file="NK1_Cytokines.txt",sep="\t",col.names=NA)

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK1_cytokines.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()
write.table(inp.melt.dt,file="inp.melt.dt.txt",sep="\t",col.names=NA)


# Activatory receptors

gene_list <- c("KLRK1",	"KLRC2",	"CD160",	"NCR3",	"SLAMF6",	"SLAMF7")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "OV_NK1"   # First group
group2 <- "LUAD_NK1" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "GROUP"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK1_Activatory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK1_Activatory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()


# Inhibitory receptors

gene_list <- c("KLRC1",	"CD300A",	"TIGIT",	"KIR3DL2",	"KIR3DL1",	"KIR2DL3",	"KIR2DL1",	"KLRB1",	"SIGLEC7",	"HAVCR2")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "OV_NK1"   # First group
group2 <- "LUAD_NK1" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "GROUP"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK1_Inhibitory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK1_Inhibitory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

# Cytotoxicity

gene_list <- c("NKG7","GZMH","FASLG","PRF1","GZMA","GZMK","GSDMD","GZMB","TNFSF10")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "OV_NK1"   # First group
group2 <- "LUAD_NK1" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "GROUP"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",
  logfc.threshold = 0,p.adjust.method = "BH",
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK1_Cytotoxicity.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK1_Cytotoxicity.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()


# Cytokine receptors

gene_list <- c("IL10RA","IL10RB","IL12RB1","IL18R1","IL18RAP","IL2RB","TGFBR2","TGFBR1","IL2RG","TGFBR3")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "OV_NK1"   # First group
group2 <- "LUAD_NK1" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "GROUP"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",
  logfc.threshold = 0,p.adjust.method = "BH",
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK1_Cytokine_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK1_Cytokine_receptors.pdf",width=5,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

#-------------------- Cohort-----------------------------------------

#Cytotoxicity

gene_list <- c("NKG7","GZMH","FASLG","PRF1","GZMA","GZMK","GSDMD","GZMB","TNFSF10")  # Replace with your gene list
group_col <- "Group1"  # Name of metadata column with group info
group1 <- "OV"   # First group
group2 <- "LUAD" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "Group1"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Cytotoxicity.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Cytotoxicity.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

#Cytokines recepotors

gene_list <- c("IL10RA","IL10RB","IL12RB1","IL18R1","IL18RAP","IL2RB","TGFBR2","TGFBR1","IL2RG","TGFBR3")  # Replace with your gene list
group_col <- "Group1"  # Name of metadata column with group info
group1 <- "OV"   # First group
group2 <- "LUAD" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "Group1"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Cytokines_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Cytokines_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()




#Inhibitory recepotors

gene_list <- c("KLRC1",	"CD300A",	"TIGIT",	"KIR3DL2",	"KIR3DL1",	"KIR2DL3",	"KIR2DL1",	"KLRB1",	"SIGLEC7",	"HAVCR2")  # Replace with your gene list
group_col <- "Group1"  # Name of metadata column with group info
group1 <- "OV"   # First group
group2 <- "LUAD" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "Group1"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Inhibitory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Inhibitory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()


#Cytokines

gene_list <- c("IFNG",	"CCL4",	"CCL4L2",	"CCL5",	"IL32",	"IL16",	"CCL3",	"XCL2",	"XCL1",	"FLT3LG")  # Replace with your gene list
group_col <- "Group1"  # Name of metadata column with group info
group1 <- "OV"   # First group
group2 <- "LUAD" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "Group1"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Cytokines.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Cytokines.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

#Activatory receptors

gene_list <- c("KLRK1","SLAMF6","CD160","SLAMF7","KLRC2","NCR3")  # Replace with your gene list
group_col <- "Group1"  # Name of metadata column with group info
group1 <- "OV"   # First group
group2 <- "LUAD" # Second group

genes_present <- gene_list[gene_list %in% rownames(nk.subset)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(nk.subset) <- "Group1"
markers <- FindMarkers(
  object =nk.subset,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Activatory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Activatory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()


#-----------------------------Fig2---------------------------------------
#------------------------------------------------------------------------

load("ALL.clean.rda")
Idents(ALL.clean)="Cohort"
PRIM.MET=subset(ALL.clean,ident=c("OV","omentum_met_distant","omentum_met_tumor","omentum_met_peritumoral","benign_omentum"))
PRIM.MET@meta.data$first20 <- substr(rownames(PRIM.MET@meta.data), 1, 20)
#laod NK without CD8T
#C:\Users\hensler\OneDrive - PPF\BackUp_disk_D\Projects\NK_paper\Revize\NatComm\VER2\Fig2\OV_Francis_cohort\test
load("NK.ver1.rda",verbose=T)
NK@meta.data$first20 <- substr(rownames(NK@meta.data), 1, 20)
common_cells <- intersect(NK@meta.data$first20,PRIM.MET@meta.data$first20)
prim_labels <- setNames(PRIM.MET@meta.data$NK.pop, PRIM.MET@meta.data$first20)

save(NK,file="NK.subset_new.rda")

Idents(NK)="NK.pop"
NK=subset(NK,ident=c("NK1","NK2","NK3"))
NK@meta.data$NK.pop <- factor(
  NK@meta.data$NK.pop,
  levels = c("NK1","NK2","NK3")
)
DimPlot(NK,reduction="umap.harmony")+
   scale_color_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
))

save(NK,file="NK.subset_new1.rda")

counts=NK@meta.data%>%group_by(Cohort,NK.pop)%>%summarize(no=n())%>% mutate(percent = 100 * no / sum(no))

counts$Cohort <- factor(counts$Cohort, levels = c("benign_omentum","Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))


p1=ggplot(counts, aes(x = Cohort, y = percent, fill = NK.pop)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(percent, 1)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  labs(y = "Percentage", x = "Group", fill = "NK Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))


p2=ggplot(counts, aes(x = Cohort, y = no, fill = NK.pop)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(no, 1)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  labs(y = "No", x = "Group", fill = "NK Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))

p1+p2

#Fig2D

d= c("IFNG", "CCL4","CCL4L2","CCL5","IL32","IL16","CCL3","XCL2","XCL1","FLT3LG",
"KLRK1", "KLRC2","CD160","NCR3","SLAMF6","SLAMF7",
"KLRC1","CD300A","TIGIT","KIR3DL2","KIR3DL1","KIR2DL3","KIR2DL1","KLRB1","SIGLEC7","HAVCR2",
"TGFBR1","IL12RB1","IL2RG","TGFBR3","TGFBR2","IL10RA","IL18R1","IL18RAP","IL2RB","IL10RB",
"FASLG","GSDMD","GZMB","NKG7","PRF1","GZMA","GZMH","GZMK","TNFSF10")

#NK
Idents(NK)="NK.pop"
Idents(NK)="Cohort"
levels(NK) <- c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK, features = d, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

write.table(cal,file="cal_Fig2D.txt",sep="\t",col.names=NA)

#selection of normalize gene expression

gene_list <- c("IFNG", "CCL4","CCL4L2","CCL5","IL32","IL16","CCL3","XCL2","XCL1","FLT3LG",
"KLRK1", "KLRC2","CD160","NCR3","SLAMF6","SLAMF7",
"KLRC1","CD300A","TIGIT","KIR3DL2","KIR3DL1","KIR2DL3","KIR2DL1","KLRB1","SIGLEC7","HAVCR2",
"TGFBR1","IL12RB1","IL2RG","TGFBR3","TGFBR2","IL10RA","IL18R1","IL18RAP","IL2RB","IL10RB",
"FASLG","GSDMD","GZMB","NKG7","PRF1","GZMA","GZMH","GZMK","TNFSF10")

data_subset <- GetAssayData(NK, slot = "data")[gene_list, , drop = FALSE]
data=data.frame(t(data_subset))
data$Cohort=paste0(NK@meta.data$Cohort)

write.table(data,file="normilze_data_Fig2D.txt",sep="\t",col.names=NA)

row_annot=read.table(file="cal_row_split.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

# Clean group names in case of spaces
row_annot$Groups <- trimws(row_annot$Groups)

# Desired block order
block_levels <- c(
  "Cytotoxicity",
  "Cytokines",
  "Cytokines_receptors",
  "Activatory_receptors",
  "Inhibitory_receptors"
)

Rowsplit <- factor(row_annot$Groups, levels = block_levels)

# Cluster within each block but keep block order
order_within_blocks <- unlist(lapply(block_levels, function(block) {
  rows_in_block <- which(Rowsplit == block)
  if (length(rows_in_block) > 1) {
    # distance + hierarchical clustering inside this block
    dist_mat <- dist(mat[rows_in_block, ])
    hc <- hclust(dist_mat, method = "ward.D2")
    rows_in_block[hc$order]
  } else {
    rows_in_block
  }
}))

mat_clustered <- mat[order_within_blocks, ]
Rowsplit_clustered <- Rowsplit[order_within_blocks]

# Plot without re-clustering rows
pdf("Heatmap_all.pdf",width=6,height=12)
ComplexHeatmap::pheatmap(
  mat_clustered,
  row_split = Rowsplit_clustered,
  cluster_rows = FALSE, # we already clustered manually
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

#NK1 subset
Idents(NK)="NK.pop"
NK1=subset(NK,idents="NK1")

Idents(NK1)="Cohort"
levels(NK1) <- c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK1, features = d, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

# Clean group names in case of spaces
row_annot$Groups <- trimws(row_annot$Groups)

# Desired block order
block_levels <- c(
  "Cytotoxicity",
  "Cytokines",
  "Cytokines_receptors",
  "Activatory_receptors",
  "Inhibitory_receptors"
)

Rowsplit <- factor(row_annot$Groups, levels = block_levels)

# Cluster within each block but keep block order
order_within_blocks <- unlist(lapply(block_levels, function(block) {
  rows_in_block <- which(Rowsplit == block)
  if (length(rows_in_block) > 1) {
    # distance + hierarchical clustering inside this block
    dist_mat <- dist(mat[rows_in_block, ])
    hc <- hclust(dist_mat, method = "ward.D2")
    rows_in_block[hc$order]
  } else {
    rows_in_block
  }
}))

mat_clustered <- mat[order_within_blocks, ]
Rowsplit_clustered <- Rowsplit[order_within_blocks]

# Plot without re-clustering rows

pdf("Heatmap_NK1.pdf",width=6,height=12)

ComplexHeatmap::pheatmap(
  mat_clustered,
  row_split = Rowsplit_clustered,
  cluster_rows = FALSE, # we already clustered manually
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

#NK2 subset

Idents(NK)="NK.pop"
NK2=subset(NK,idents="NK2")

Idents(NK2)="Cohort"
levels(NK2) <- c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK2, features = d, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

row_annot=read.table(file="cal_row_split.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
# Clean group names in case of spaces
row_annot$Groups <- trimws(row_annot$Groups)

# Desired block order
block_levels <- c(
  "Cytotoxicity",
  "Cytokines",
  "Cytokines_receptors",
  "Activatory_receptors",
  "Inhibitory_receptors"
)

Rowsplit <- factor(row_annot$Groups, levels = block_levels)

# Cluster within each block but keep block order
order_within_blocks <- unlist(lapply(block_levels, function(block) {
  rows_in_block <- which(Rowsplit == block)
  if (length(rows_in_block) > 1) {
    # distance + hierarchical clustering inside this block
    dist_mat <- dist(mat[rows_in_block, ])
    hc <- hclust(dist_mat, method = "ward.D2")
    rows_in_block[hc$order]
  } else {
    rows_in_block
  }
}))

mat_clustered <- mat[order_within_blocks, ]
Rowsplit_clustered <- Rowsplit[order_within_blocks]

# Plot without re-clustering rows
pdf("Heatmap_NK2.pdf",width=6,height=12)
ComplexHeatmap::pheatmap(
  mat_clustered,
  row_split = Rowsplit_clustered,
  cluster_rows = FALSE, # we already clustered manually
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

#NK3 subset

Idents(NK)="NK.pop"
NK3=subset(NK,idents="NK3")

Idents(NK3)="Cohort"
levels(NK3) <- c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK3, features = d, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

row_annot=read.table(file="cal_row_split.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
# Clean group names in case of spaces
row_annot$Groups <- trimws(row_annot$Groups)

# Desired block order
block_levels <- c(
  "Cytotoxicity",
  "Cytokines",
  "Cytokines_receptors",
  "Activatory_receptors",
  "Inhibitory_receptors"
)

Rowsplit <- factor(row_annot$Groups, levels = block_levels)

# Cluster within each block but keep block order
order_within_blocks <- unlist(lapply(block_levels, function(block) {
  rows_in_block <- which(Rowsplit == block)
  if (length(rows_in_block) > 1) {
    # distance + hierarchical clustering inside this block
    dist_mat <- dist(mat[rows_in_block, ])
    hc <- hclust(dist_mat, method = "ward.D2")
    rows_in_block[hc$order]
  } else {
    rows_in_block
  }
}))

mat_clustered <- mat[order_within_blocks, ]
Rowsplit_clustered <- Rowsplit[order_within_blocks]

# Plot without re-clustering rows
pdf("Heatmap_NK3.pdf",width=6,height=12)
ComplexHeatmap::pheatmap(
  mat_clustered,
  row_split = Rowsplit_clustered,
  cluster_rows = FALSE, # we already clustered manually
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

Idents(NK)="harmony_clusters"
DimPlot2(NK,reduction="umap.harmony",theme = theme_umap_arrows())

#------------------DotPlot Fig.2E----------------------
genes <- c("SPON2","NKG7","FCER1G","PRF1","GZMB","CLIC3","MYOM2","FGFBP2","AKR1C3","IGFBP7","LAIR2","KLRB1","CX3CR1","CHST2","CD247","CD38","ADGRG1",
"CCL4","CTSD","CD160","GZMK","XCL1","XCL2","TPT1","TCF7","KLRC1","CMC1","GPR183","SELL","CD44","DUSP2","ZFP36L2","AREG","PIK3R1","IL7R","COTL1","IFITM3","IL2RB","FOS",
"LTB","KLRC2","CCL5","CD3E","GZMH","PTMS","S100A6","LINC01871","VIM","IL32","ITGB1","S100A4","LGALS1","PRDM1","PPDPF","CD2","ZBTB38","DSTN","CD3D","CD3G","HBA1")

Idents(NK)="NK.pop"
levels(NK) <- c("NK1", "NK2","NK3")

pdf("DotPlot_NK.pop.pdf",width=5,height=12)
DotPlot(NK, features = genes)+scale_colour_gradientn(
    colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
    breaks = -2:2
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


dev.off()

Idents(NK)="Cohort"
levels(NK) <- c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

pdf("DotPlot_Cohort.pdf",width=5,height=12)
DotPlot(NK, features = genes)+scale_colour_gradientn(
    colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
    breaks = -2:2
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


dev.off()

#heatmap
Genes <- rev(c("SPON2","NKG7","FCER1G","PRF1","GZMB","CLIC3","MYOM2","FGFBP2","AKR1C3","IGFBP7","LAIR2","KLRB1","CX3CR1","CHST2","CD247","CD38","ADGRG1",
"CCL4","CTSD","CD160","GZMK","XCL1","XCL2","TPT1","TCF7","KLRC1","CMC1","GPR183","SELL","CD44","DUSP2","ZFP36L2","AREG","PIK3R1","IL7R","COTL1","IFITM3","IL2RB","FOS",
"LTB","KLRC2","CCL5","CD3E","GZMH","PTMS","S100A6","LINC01871","VIM","IL32","ITGB1","S100A4","LGALS1","PRDM1","PPDPF","CD2","ZBTB38","DSTN","CD3D","CD3G","HBA1"))


cal <- CalcStats(NK, features = Genes, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)
pdf("Heatmap_NK.pop.pdf",width=6,height=12)
ComplexHeatmap::pheatmap(
  mat,
  cluster_rows = FALSE, # we already clustered manually
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()


#-------------------LinePlot Fig2D-----------------

signatures <- list(
  Inhibitory_receptors = c("SIGLEC7","CD300A","KIR3DL2","TIGIT","KLRC1","HAVCR2","KIR3DL1","KLRB1","KIR2DL3","KIR2DL1"),
  Cytokines_receptors  = c("IL10RA","IL10RB","IL12RB1","IL18R1","IL18RAP","IL2RB","TGFBR2","TGFBR1","IL2RG","TGFBR3"),
  Cytotoxicity         = c("NKG7","GZMH","FASLG","PRF1","GZMA","GZMK","GSDMD","GZMB","TNFSF10"),
  Activatory_receptors = c("KLRK1","SLAMF6","CD160","SLAMF7","KLRC2","NCR3"),
  Cytokines            = c("FLT3LG","CCL4","CCL3","XCL2","XCL1","IFNG","CCL4L2","IL16","CCL5","IL32")
)

# Calculate mean expression per cell for each signature
for (sig_name in names(signatures)) {
  genes <- signatures[[sig_name]]
  genes_present <- genes[genes %in% rownames(NK)]
  
  expr_mat <- GetAssayData(NK, assay = "RNA", slot = "data")[genes_present, , drop = FALSE]
  NK[[sig_name]] <- Matrix::colMeans(expr_mat)
}

# Gather into long format
df <- FetchData(NK, vars = c("Cohort", "NK.pop", names(signatures))) %>%
  pivot_longer(
    cols = -c(Cohort, NK.pop),
    names_to = "Signature",
    values_to = "MeanExpression"
  )

df$Cohort <- factor(df$Cohort, levels = c("benign_omentum", "Primary_tumor", "omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))
df$Signature <- factor(df$Signature, levels =c(
  "Cytotoxicity",
  "Cytokines",
  "Cytokines_receptors",
  "Activatory_receptors",
  "Inhibitory_receptors"))

# Line plot with points
ggplot(df, aes(x = Cohort, y = MeanExpression, color = NK.pop, group = NK.pop)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2, shape = 16) +  # solid circle
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~Signature, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mean Expression", x = "Cohort", color = "NK Population") +
  scale_color_manual(values = c(
    "NK1" = "#99002F",
    "NK2" = "#006498",
    "NK3" = "#006100"
  ))

g <- c("SIGLEC7","CD300A","KIR3DL2","TIGIT","KLRC1","HAVCR2","KIR3DL1","KLRB1","KIR2DL3","KIR2DL1",
  "IL10RA","IL10RB","IL12RB1","IL18R1","IL18RAP","IL2RB","TGFBR2","TGFBR1","IL2RG","TGFBR3",
  "NKG7","GZMH","FASLG","PRF1","GZMA","GZMK","GSDMD","GZMB","TNFSF10",
  "KLRK1","SLAMF6","CD160","SLAMF7","KLRC2","NCR3",
  "FLT3LG","CCL4","CCL3","XCL2","XCL1","IFNG","CCL4L2","IL16","CCL5","IL32")

data_subset <- GetAssayData(NK, slot = "data")[g, , drop = FALSE]
data=data.frame(t(data_subset))
data$Cohort=NK@meta.data$Cohort
data$NK.pop=NK@meta.data$NK.pop
write.table(data,file="normilze_data_Fig1K.txt",sep="\t",col.names=NA)

#-------------------LinePlot Fig2D all -----------------

signatures <- list(
  Inhibitory_receptors = c("SIGLEC7","CD300A","KIR3DL2","TIGIT","KLRC1","HAVCR2","KIR3DL1","KLRB1","KIR2DL3","KIR2DL1"),
  Cytotoxicity         = c("NKG7","GZMH","FASLG","PRF1","GZMA","GZMK","GSDMD","GZMB","TNFSF10"),
  Cytokines            = c("FLT3LG","CCL4","CCL3","XCL2","XCL1","IFNG","CCL4L2","IL16","CCL5","IL32")
)

# Calculate mean expression per cell for each signature
for (sig_name in names(signatures)) {
  genes <- signatures[[sig_name]]
  genes_present <- genes[genes %in% rownames(NK)]
  
  expr_mat <- GetAssayData(NK, assay = "RNA", slot = "data")[genes_present, , drop = FALSE]
  NK[[sig_name]] <- Matrix::colMeans(expr_mat)
}

# Gather into long format
df <- FetchData(NK, vars = c("Cohort", names(signatures))) %>%
  pivot_longer(
    cols = -c(Cohort),
    names_to = "Signature",
    values_to = "MeanExpression"
  )

df$Cohort <- factor(df$Cohort, levels = c("benign_omentum", "Primary_tumor", "omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))
df$Signature <- factor(df$Signature, levels =c(
  "Cytotoxicity",
  "Cytokines",
  "Inhibitory_receptors"))

# Line plot with points
ggplot(df, aes(x = Cohort, y = MeanExpression, color = Signature, group = Signature)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2, shape = 16) +  # solid circle
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~Signature, scales = "free_y") +
  scale_color_manual(values = c("#d62f27", "#fdad60", "#4574b3")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mean Expression", x = "Cohort", color = "Signature")


#-----------------------BarPlot of gene expression-------------------------

load("NK.subset_new1.rda",verbose=T)
NK@meta.data[["GROUP"]]=paste0(NK@meta.data$Cohort,"_",NK@meta.data$NK.pop)

#Cytokines

gene_list <- c("IFNG",	"CCL4",	"CCL4L2",	"CCL5",	"IL32",	"IL16",	"CCL3",	"XCL2",	"XCL1",	"FLT3LG")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral_NK3"   # First group
group2 <- "omentum_met_tumor_NK3" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "GROUP"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK3_cytokines.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK3_cytokines.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()



# Activatory receptors

gene_list <- c("KLRK1",	"KLRC2",	"CD160",	"NCR3",	"SLAMF6",	"SLAMF7")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral_NK3"   # First group
group2 <- "omentum_met_tumor_NK3" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "GROUP"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK3_Activatory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK3_Activatory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()


# Inhibitory receptors

gene_list <- c("KLRC1",	"CD300A",	"TIGIT",	"KIR3DL2",	"KIR3DL1",	"KIR2DL3",	"KIR2DL1",	"KLRB1",	"SIGLEC7",	"HAVCR2")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral_NK3"   # First group
group2 <- "omentum_met_tumor_NK3" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "GROUP"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK3_Inhibitory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK3_Inhibitory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

# Cytotoxicity

gene_list <- c("NKG7","GZMH","FASLG","PRF1","GZMA","GZMK","GSDMD","GZMB","TNFSF10")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral_NK3"   # First group
group2 <- "omentum_met_tumor_NK3" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "GROUP"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK3_Cytotoxicity.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK3_Cytotoxicity.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()


# Cytokine receptors

gene_list <- c("IL10RA","IL10RB","IL12RB1","IL18R1","IL18RAP","IL2RB","TGFBR2","TGFBR1","IL2RG","TGFBR3")  # Replace with your gene list
group_col <- "GROUP"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral_NK3"   # First group
group2 <- "omentum_met_tumor_NK3" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "GROUP"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="NK3_Cytokine_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("NK3_Cytokine_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

#-------------------- Cohort-----------------------------------------

#Cytotoxicity

gene_list <- c("NKG7","GZMH","FASLG","PRF1","GZMA","GZMK","GSDMD","GZMB","TNFSF10")  # Replace with your gene list
group_col <- "Cohort"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral"   # First group
group2 <- "omentum_met_tumor" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "Cohort"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Cytotoxicity.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Cytotoxicity.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

#Cytokines recepotors

gene_list <- c("IL10RA","IL10RB","IL12RB1","IL18R1","IL18RAP","IL2RB","TGFBR2","TGFBR1","IL2RG","TGFBR3")  # Replace with your gene list
group_col <- "Cohort"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral"   # First group
group2 <- "omentum_met_tumor" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "Cohort"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Cytokines_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Cytokines_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()




#Inhibitory recepotors

gene_list <- c("KLRC1",	"CD300A",	"TIGIT",	"KIR3DL2",	"KIR3DL1",	"KIR2DL3",	"KIR2DL1",	"KLRB1",	"SIGLEC7",	"HAVCR2")  # Replace with your gene list
group_col <- "Cohort"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral"   # First group
group2 <- "omentum_met_tumor" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "Cohort"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Inhibitory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Inhibitory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()


#Cytokines

gene_list <- c("IFNG",	"CCL4",	"CCL4L2",	"CCL5",	"IL32",	"IL16",	"CCL3",	"XCL2",	"XCL1",	"FLT3LG")  # Replace with your gene list
group_col <- "Cohort"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral"   # First group
group2 <- "omentum_met_tumor" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "Cohort"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Cytokines.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Cytokines.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

#Activatory receptors

gene_list <- c("KLRK1","SLAMF6","CD160","SLAMF7","KLRC2","NCR3")  # Replace with your gene list
group_col <- "Cohort"  # Name of metadata column with group info
group1 <- "omentum_met_peritumoral"   # First group
group2 <- "omentum_met_tumor" # Second group

genes_present <- gene_list[gene_list %in% rownames(NK)]
if (length(genes_present) == 0) stop("None of the genes in gene_list are in the Seurat object.")

# --- STEP 2: Extract expression ---
# Seurat's FindMarkers can run Wilcoxon test directly
Idents(NK) <- "Cohort"
markers <- FindMarkers(
  object =NK,
  ident.1 = group1,
  ident.2 = group2,
  features = genes_present,
  test.use = "wilcox",p.adjust.method = "BH",
  logfc.threshold = 0,
  min.pct = 0
)

markers$gene <- rownames(markers)

write.table(markers ,file="Activatory_receptors.txt",sep="\t",col.names=NA)

markers <- markers %>%
  mutate(
    neg_log10_p = -log10(p_val_adj),  # adjusted p-value
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0  ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0 ~ "Down",
      TRUE ~ "NS"
    )
  )

markers <- markers %>%
  arrange(p_val_adj) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

pdf("Activatory_receptors.pdf",width=7,height=3)
ggplot(markers, aes(x = gene, y = avg_log2FC, fill = significance)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  coord_flip() +  # rotate for readability
  labs(
    title = paste("Wilcoxon test:", group1, "vs", group2),
    x = NULL,
    y = "Log2 Fold Change")+
theme(
    axis.text.y = element_text(color = "black")
  )
dev.off()

#---------------------View NK subset-----------------------
Idents(NK)="NK1"

FeaturePlot(NK,features="NK1",reduction="umap.harmony")& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


library(SCpubr)
library(Nebulosa)
Nebulosa::plot_density(NK,reduction="umap.harmony",
                            features = "NK1")



#---------NK ligands-------------------------

load("NK.subset_new1.rda",verbose=T)
NK.ligands=read.table(file="NK_ligands.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

d= NK.ligands$Gene

Idents(nk.subset)="Group1"
levels(nk.subset) <- c("OV", "LUAD")

cal <- CalcStats(nk.subset, features = d, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK_Ligands.pdf",width=4,height=12)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

#---------------Tumor.cells----------------------------------

load("Tumor.cells.rda",verbose=T)
Epi=subset(Tumor.cells,idents=c("Epithelial.cells"))
Idents(Epi)="Cohort"
Epi$Group <- case_when(
  grepl("OVT1",Epi@meta.data$Cohort) ~ "OV",
  grepl("OVT2",Epi@meta.data$Cohort) ~ "OV",
  grepl("OVT3",Epi@meta.data$Cohort) ~ "OV",
  grepl("OVT4",Epi@meta.data$Cohort) ~ "OV",
  grepl("LUNGT1",Epi@meta.data$Cohort) ~ "LUAD",
  grepl("LUNGT2",Epi@meta.data$Cohort) ~ "LUAD",
  grepl("LUNGT3",Epi@meta.data$Cohort) ~ "LUAD",
  grepl("LUNGT4",Epi@meta.data$Cohort) ~ "LUAD",
  grepl("LUNGT5",Epi@meta.data$Cohort) ~ "LUAD",
      TRUE ~ NA_character_
)

Tumor.ligands=read.table(file="Tumor.cells_ligands.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=Tumor.ligands$Gene

Idents(Epi)="Group"
levels(Epi) <- c("OV", "LUAD")

cal <- CalcStats(Epi, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("Tumor.cells_ligands.pdf",width=4,height=6)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()


NK_counts=GetAssayData(object = nk.subset, assay = "RNA", slot = "counts")
T_counts=GetAssayData(object = Epi, assay = "RNA", slot = "counts")
NK <- CreateSeuratObject(counts = NK_counts)
NK[["Group"]]=paste0(nk.subset@meta.data$Group1)
NK[["Sample"]]=paste0(nk.subset@meta.data$Sample)
NK[["NK.pop"]]=paste0(nk.subset@meta.data$NK.pop)

T <- CreateSeuratObject(counts = T_counts)

T[["Group"]]=paste0(Epi@meta.data$Group)
T[["Sample"]]=paste0(Epi@meta.data$Cohort)

#-----------combine NK+Tumor cells (Epithelial.cells)---------------

NKT=merge(NK,T)
save(NKT,file="NKT.rda")

NKT[["RNA"]] <- split(NKT[["RNA"]], f = NKT$Sample)

NKT<- NormalizeData(NKT)
NKT<- FindVariableFeatures(NKT)
NKT<- ScaleData(NKT)
NKT<- RunPCA(NKT)

ElbowPlot(NKT)

NKT<- FindNeighbors(NKT, dims = 1:15, reduction = "pca")
NKT<- FindClusters(NKT, cluster.name = "unintegrated_clusters")
NKT<- RunUMAP(NKT, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

options(future.globals.maxSize = 20000 * 1024^4)

NKT<- IntegrateLayers(
  object = NKT, method = RPCAIntegration,dims = 1:15,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  k.weight=15,k.anchor = 15,k.filter = 15,k.score = 15,
  verbose = FALSE
)

NKT[["RNA"]] <- JoinLayers(NKT[["RNA"]])

# RPCA
NKT<- FindNeighbors(NKT, reduction = "integrated.rpca", dims = 1:15)
NKT<- FindClusters(NKT, cluster.name = "rpca_clusters")
NKT<- RunUMAP(NKT, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")

#-------------------------------Ligand/Receptor file------------

LigRec=read.table(file="Lig.Rec.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=LigRec$Gene

NKT[["G"]]=paste0(NKT@meta.data$Group,"_",NKT@meta.data$NK.pop)

NKT$Group_F <- case_when(
  grepl("OVT1",NKT@meta.data$Sample) ~ "OV_T",
  grepl("OVT2",NKT@meta.data$Sample) ~ "OV_T",
  grepl("OVT3",NKT@meta.data$Sample) ~ "OV_T",
  grepl("OVT4",NKT@meta.data$Sample) ~ "OV_T",
  grepl("LUNGT1",NKT@meta.data$Sample) ~ "LUAD_T",
  grepl("LUNGT2",NKT@meta.data$Sample) ~ "LUAD_T",
  grepl("LUNGT3",NKT@meta.data$Sample) ~ "LUAD_T",
  grepl("LUNGT4",NKT@meta.data$Sample) ~ "LUAD_T",
  grepl("LUNGT5",NKT@meta.data$Sample) ~ "LUAD_T",
  grepl("OV_NK1",NKT@meta.data$G) ~ "OV_NK1",
  grepl("OV_NK2",NKT@meta.data$G) ~ "OV_NK2",
  grepl("OV_NK3",NKT@meta.data$G) ~ "OV_NK3",
  grepl("LUAD_NK1",NKT@meta.data$G) ~ "LUAD_NK1",
  grepl("LUAD_NK2",NKT@meta.data$G) ~ "LUAD_NK2",
  grepl("LUAD_NK3",NKT@meta.data$G) ~ "LUAD_NK3",
      TRUE ~ NA_character_
)

save(NKT,file="NKT_ver1.rda")

Idents(NKT)="Group_F"
levels(NKT) <- c("OV_T", "OV_NK1","OV_NK2","OV_NK3","LUAD_T","LUAD_NK1","LUAD_NK2","LUAD_NK3")

cal <- CalcStats(NKT, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("LigRec.pdf",width=4,height=6)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()


#-------------------------Correlation-------------------

LigRec=read.table(file="Lig.Rec.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=LigRec$Gene

expr <- FetchData(NKT, vars = Gene)
expr$group <- NKT@meta.data$Group_F

#---------------Compute correlation matrices per group--------------
gene1="TNFSF10"
gene2="TNFRSF10B"

cor_results <- expr %>%
  group_by(group) %>%
  summarise(
    cor_test = list(cor.test(.data[[gene1]], .data[[gene2]], method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    cor  = sapply(cor_test, function(x) x$estimate),
    pval = sapply(cor_test, function(x) x$p.value)
  ) %>%
  mutate(
    pval_adj = p.adjust(pval, method = "BH")   # FDR adjustment
  ) %>%
  select(group, cor, pval, pval_adj)

print(cor_results)

df <- expr%>% left_join(cor_results, by = "group")


#---------------EPI_OV+Francis---------------------------------------
#Get tumor OV tumor cells

setwd("C:/Users/hensler/OneDrive - PPF/BackUp_disk_D/Projects/NK_paper/Revize/NatComm/VER3/Fig1/NK_Ligands/Tumor_cells")
load("NKT.rda")
Idents(NKT)="Sample"
OV=subset(NKT,idents=c("OVT1","OVT2","OVT3","OVT4"))
OV[["Cohort"]]=paste0(OV@meta.data$Group)

#Get NK cells from OV PRIMARY +MET

setwd("C:/Users/hensler/OneDrive - PPF/BackUp_disk_D/Projects/NK_paper/Revize/NatComm/VER3/Fig2")
load("NK.subset_new1.rda")
Idents(NK)="Cohort"
NK[["Sample"]]=paste0(NK@meta.data$orig.ident)

#Combine OV tumor+NK

NKOV=merge(NK,OV)
NKOV[["RNA"]] <- split(NKOV[["RNA"]], f = NKOV$Sample)

NKOV<- NormalizeData(NKOV)
NKOV<- FindVariableFeatures(NKOV)
NKOV<- ScaleData(NKOV)
NKOV<- RunPCA(NKOV)

ElbowPlot(NKOV)

NKOV<- FindNeighbors(NKOV, dims = 1:15, reduction = "pca")
NKOV<- FindClusters(NKOV, cluster.name = "unintegrated_clusters")
NKOV<- RunUMAP(NKOV, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

options(future.globals.maxSize = 20000 * 1024^4)

NKOV<- IntegrateLayers(
  object = NKOV, method = HarmonyIntegration,dims = 1:15,
  orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE
)

NKOV[["RNA"]] <- JoinLayers(NKOV[["RNA"]])

# Harmony
NKOV<- FindNeighbors(NKOV, reduction = "harmony", dims = 1:15)
NKOV<- FindClusters(NKOV, cluster.name = "harmony_clusters")
NKOV<- RunUMAP(NKOV, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")

DimPlot(NKOV,reduction="umap.harmony")

#-------------------------------Ligand/Receptor file------------
setwd( "C:/Users/hensler/OneDrive - PPF/BackUp_disk_D/Projects/NK_paper/Revize/NatComm/VER3/Fig2/NK_ligands")
LigRec=read.table(file="Lig.Rec.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=LigRec$Gene

Idents(NKOV)="Cohort"
levels(NKOV) <- c("OV","benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NKOV, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("LigRec.pdf",width=4,height=8)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

#NK1
Idents(NKOV)="NK.pop"
NK1=subset(NKOV,idents=c("NK2","NK3"),invert=T)
Idents(NK1)="Cohort"
levels(NK1) <- c("OV","benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK1, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK1_LigRec.pdf",width=4,height=8)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

#NK2
Idents(NKOV)="NK.pop"
NK2=subset(NKOV,idents=c("NK1","NK3"),invert=T)
Idents(NK2)="Cohort"
levels(NK2) <- c("OV","benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK2, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK2_LigRec.pdf",width=4,height=8)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

#NK3
Idents(NKOV)="NK.pop"
NK3=subset(NKOV,idents=c("NK1","NK2"),invert=T)
Idents(NK3)="Cohort"
levels(NK3) <- c("OV","benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK3, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK3_LigRec.pdf",width=4,height=8)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()
#----NK Ligand/Receptor file ---------------------
load("NK.subset_new1.rda")
Idents(NKOV)="Cohort"

NK.ligands=read.table(file="NK_ligands.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=NK.ligands$Gene

Idents(NK)="Cohort"
levels(NK) <- c("benign_omentum","Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK_Ligands.pdf",width=4,height=9)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()


#----NK1sub Ligand/Receptor file ---------------------
Idents(NK)="NK.pop"
NK1=subset(NK,idents="NK1")

NK.ligands=read.table(file="NK_ligands.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=NK.ligands$Gene

Idents(NK1)="Cohort"
levels(NK1) <- c("benign_omentum","Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK1, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK1_Ligands.pdf",width=4,height=9)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

#----NK2sub Ligand/Receptor file ---------------------
Idents(NK)="NK.pop"
NK2=subset(NK,idents="NK2")

NK.ligands=read.table(file="NK_ligands.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=NK.ligands$Gene

Idents(NK2)="Cohort"
levels(NK2) <- c("benign_omentum","Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK2, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK2_Ligands.pdf",width=4,height=9)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()


#----NK3sub Ligand/Receptor file ---------------------
Idents(NK)="NK.pop"
NK3=subset(NK,idents="NK3")

NK.ligands=read.table(file="NK_ligands.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)

Gene=NK.ligands$Gene

Idents(NK3)="Cohort"
levels(NK3) <- c("benign_omentum","Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor")

cal <- CalcStats(NK3, features = Gene, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("NK3_Ligands.pdf",width=4,height=9)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -1.5:1.5,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

#-------------Ligands on Tumor cells, epithelial cells-----------
setwd("C:/Users/hensler/OneDrive - PPF/BackUp_disk_D/Projects/NK_paper/Revize/NatComm")
load("merged.rda",verbose=T)
#harmony_clusters 10 annotated as a Epithelial cells

Idents(merged)="harmony_clusters"
EPI=subset(merged,idents="10")
View(EPI@meta.data)
EPI[["Sample"]]=paste0(EPI@meta.data$orig.ident)
EPI[["Cohort"]]=paste0(EPI@meta.data$Group)
Idents(EPI)="BlueprintEncodeData"
EPI=subset(epi_Francis,idents="Epithelial.cells")


#upload epithelial cells from OV primary tumor
setwd("C:/Users/hensler/OneDrive - PPF/BackUp_disk_D/Projects/NK_paper/Revize/NatComm/VER3/Fig2/NK_ligands")
load("NKOV.rda")
Idents(NKOV)="Cohort"
epi=subset(NKOV,idents="OV")
colnames(epi@meta.data)[colnames(epi@meta.data) == "Group"] <- "Cohort"
epi[["Sample"]]=paste0(epi@meta.data$Sample)
epi[["Cohort"]]=paste0(epi@meta.data$Cohort)

intersect=intersect(rownames(EPI[["RNA"]]$counts),rownames(epi[["RNA"]]$counts))
OVT.MET=merge(epi[intersect, ], EPI[intersect, ])
OVT.MET[["orig.ident"]]=paste0(OVT.MET@meta.data$Sample)

OVT.MET<- NormalizeData(OVT.MET)
OVT.MET<- FindVariableFeatures(OVT.MET)
OVT.MET<- ScaleData(OVT.MET)
OVT.MET<- RunPCA(OVT.MET)

ElbowPlot(OVT.MET)

OVT.MET<- FindNeighbors(OVT.MET, dims = 1:15, reduction = "pca")
OVT.MET<- FindClusters(OVT.MET, cluster.name = "unintegrated_clusters")
OVT.MET<- RunUMAP(OVT.MET, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

options(future.globals.maxSize = 20000 * 1024^4)

OVT.MET<- IntegrateLayers(
  object = OVT.MET, method = HarmonyIntegration,dims = 1:15,
  orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE,group.by = "orig.ident"
)

OVT.MET[["RNA"]] <- JoinLayers(OVT.MET[["RNA"]])

# Harmony
OVT.MET<- FindNeighbors(OVT.MET, reduction = "harmony", dims = 1:15)
OVT.MET<- FindClusters(OVT.MET, cluster.name = "harmony_clusters")
OVT.MET<- RunUMAP(OVT.MET, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")

save(OVT.MET,file="OVT.MET.rda")

DimPlot(OVT.MET,reduction="umap.harmony",split.by="Cohort")

#------------------SLAMF Receptors----------------

setwd("C:/Users/hensler/OneDrive - PPF/BackUp_disk_D/Projects/NK_paper/Revize/NatComm/VER3/Fig1")
load("NK.subset_new1.rda",verbose=T)

SLAMF <- c("SLAMF1","CD48","LY9","CD244","CD84","SLAMF6","SLAMF7","SLAMF8","SLAMF9")
nk.subset@meta.data[["Group.f"]]=paste0(nk.subset@meta.data$Group1,"_",nk.subset@meta.data$NK.pop)
nk.subset@meta.data$Group.f=factor(nk.subset@meta.data$Group.f,levels=c("OV_NK1","LUAD_NK1","OV_NK2","LUAD_NK2","OV_NK3","LUAD_NK3"))

Idents(nk.subset)="Group.f"

cal <- CalcStats(nk.subset, features = SLAMF, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)

pdf("SLAMF_genes.pdf",width=4,height=3)
ComplexHeatmap::pheatmap(
  mat,
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

write.table(cal,file="SLAMPF_OV_LUNG_cal.data.txt",sep="\t",col.names=NA)

DATA <- subset(nk.subset, features = SLAMF)
norm.DATA <- SeuratObject::GetAssayData(DATA, layer = "data")
meta_cols <- c("Sample", "NK.pop") 
meta_data <- DATA@meta.data[, meta_cols, drop = FALSE]
gene_matrix_t <- t(as.matrix(norm.DATA))
final_matrix <- cbind(meta_data, gene_matrix_t)

write.table(final_matrix,file="SLAMF_OV_LUNG_norm.data.txt",sep="\t",col.names=NA)


#Francis cohort
setwd("C:/Users/hensler/OneDrive - PPF/BackUp_disk_D/Projects/NK_paper/Revize/NatComm/VER3/Fig2")

load("NK.subset_new1.rda",verbose=T)
#
Idents(NK)="Cohort"
NK@meta.data$Cohort <-factor(NK@meta.data$Cohort,level=c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))

cal <- CalcStats(NK, features = SLAMF, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

write.table(cal,file="PRIM.Francis_SLAMF_cal.data.txt",sep="\t",col.names=NA)

library(circlize)
library(ComplexHeatmap)
library(pheatmap)

set.seed(123456)
pdf("PRIM.Francis_SLAMF.pdf",width=4,height=4)
ComplexHeatmap::pheatmap(
  mat,
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

DATA <- subset(NK, features = SLAMF)
norm.DATA <- SeuratObject::GetAssayData(DATA, layer = "data")
meta_cols <- c("orig.ident", "Cohort", "NK.pop") 
meta_data <- DATA@meta.data[, meta_cols, drop = FALSE]
gene_matrix_t <- t(as.matrix(norm.DATA))
final_matrix <- cbind(meta_data, gene_matrix_t)

write.table(final_matrix,file="PRIM.Francis_SLAMF_norm.data.txt",sep="\t",col.names=NA)

#------------NK1.sub------------------------
Idents(NK)="NK.pop"
NK1=subset(NK,idents="NK1")

Idents(NK1)="Cohort"
NK1@meta.data$Cohort <-factor(NK1@meta.data$Cohort,level=c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))


cal <- CalcStats(NK1, features = SLAMF, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

write.table(cal,file="NK1_PRIM.Francis_SLAMF_cal.data.txt",sep="\t",col.names=NA)

library(circlize)
library(ComplexHeatmap)
library(pheatmap)

set.seed(123456)
pdf("NK1_PRIM.Francis_SLAMF.pdf",width=4,height=4)
ComplexHeatmap::pheatmap(
  mat,
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

#------------NK2.sub------------------------
Idents(NK)="NK.pop"
NK2=subset(NK,idents="NK2")

Idents(NK2)="Cohort"
NK2@meta.data$Cohort <-factor(NK2@meta.data$Cohort,level=c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))


cal <- CalcStats(NK2, features = SLAMF, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

write.table(cal,file="NK2_PRIM.Francis_SLAMF_cal.data.txt",sep="\t",col.names=NA)

library(circlize)
library(ComplexHeatmap)
library(pheatmap)

set.seed(123456)
pdf("NK2_PRIM.Francis_SLAMF.pdf",width=4,height=4)
ComplexHeatmap::pheatmap(
  mat,
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()

#------------NK3.sub------------------------
Idents(NK)="NK.pop"
NK3=subset(NK,idents="NK3")

Idents(NK3)="Cohort"
NK3@meta.data$Cohort <-factor(NK3@meta.data$Cohort,level=c("benign_omentum", "Primary_tumor","omentum_met_distant","omentum_met_peritumoral","omentum_met_tumor"))


cal <- CalcStats(NK3, features = SLAMF, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

write.table(cal,file="NK3_PRIM.Francis_SLAMF_cal.data.txt",sep="\t",col.names=NA)

library(circlize)
library(ComplexHeatmap)
library(pheatmap)

set.seed(123456)
pdf("NK3_PRIM.Francis_SLAMF.pdf",width=4,height=4)
ComplexHeatmap::pheatmap(
  mat,
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  legend_breaks = -1.5:1.5,
  show_row_dend = FALSE,
  show_parent_dend_line = FALSE,
  row_title = NULL,
  color = rev(brewer.pal(n = 11, name = "RdYlBu"))
)

dev.off()


