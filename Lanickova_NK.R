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
library(ProjecTILs)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(wordcloud2)
library(wordcloud)
library(SeuratExtend)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)

#uploading ovarian cancer samples from the GSE180661 dataset ( OV009-CD45P-LEFT-OVARY-r0; OV026-CD45P-RIGHT-OVARY-r0;OV067-CD45P-RIGHT-OVARY-AND-TUBE-r0;OV082-CD45P-LEFT-r0;OV112-CD45P-LEFT-r0)
#uploading lung cancer samples from the GSE131907 dataset (BRONCHO_58; EBUS_06; EBUS_49;LUNG_T28; LUNG_T31)
#Low-quality cells were excluded from each sample individually based on the following thresholds: number of detected features (nFeature) ranging between 200 and 700, and mitochondrial gene content exceeding 15% (see method in manuscript) 
#genes that were not common to both datasets were removed from the datasets
#The ovarian and lung datasets were merged and further analyzed
OV.LUNG=merge(OV, LUNG)


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

#Selection of Natural killer subset
Idents(OV.LUNG)="BlueprintEncodeData"
nk.subset=subset(OV.LUNG,subset = NK.cells == "NK.cells")

# pseudobulk DE analysis
#NK1 and NK2 labels were added to the metadata, column was named as "Samples" 

pseudo <- AggregateExpression(nk.subset, assays = "RNA", return.seurat = T, group.by = c("Samples"))
bulk <- FindMarkers(object = pseudo, 
                         ident.1 = "NK1", 
                         ident.2 = "NK2",
                         test.use = "DESeq2")

#DotPlot Fig.1J; The most significant genes between NK1, NK2

Genes=rev(c("SRGAP3",	"CD44",	"CASK",	"ITGA1",	"CLNK",
"LDLRAD4",	"PDE4B",	"COTL1",	"XCL1",	"CD96",	"ATP8B4",
"LDB2",	"SAMSN1",	"KRT81",	"ZBTB20",	"SFMBT2",	"PDE7A",
"RBPJ",	"ZNF331",	"GZMK","NCAM1","NCR1",	"FGFBP2",	"SPON2",
"FCGR3A",	"NKG7",	"PRF1",	"GZMH",	"CST7",
"PLAC8",	"CYBA",	"GZMB",	"CX3CR1",	"S100A4",	"KLRF1",
"ITGB2",	"PFN1",	"S1PR5",	"PRSS23",	"IGFBP7",	"ARL4C",
"AKR1C3"))

DotPlot2(nk.subset, features = Genes, cols = c("blue", "red"))+coord_flip()

#Heatmap Fig.1K

d= c("IFNG", "CCL4","CCL4L2","CCL5","IL32","IL16","CCL3","XCL2","XCL1","FLT3LG",
"KLRK1", "KLRC2","CD160","NCR3","SLAMF6","SLAMF7",
"KLRC1","CD300A","TIGIT","KIR3DL2","KIR3DL1","KIR2DL3","KIR2DL1","KLRB1","SIGLEC7","HAVCR2",
"TGFBR1","IL12RB1","IL2RG","TGFBR3","TGFBR2","IL10RA","IL18R1","IL18RAP","IL2RB","IL10RB",
"FASLG","GSDMD","GZMB","NKG7","PRF1","GZMA","GZMH","GZMK","TNFSF10")

cal <- CalcStats(nk.subset, features = d, method = "zscore")
mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("Heatmap.pdf",width=3.6,height=10)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,legend_breaks = -2:2,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

# A column indicating OV or LUNG cohort was added to the metadata

CD8.subset=subset(OV.LUNG,idents=c("CD8..T.cells"))

CD8.subset[["RNA3"]] <- as(object = CD8.subset[["RNA"]], Class = "Assay")
ref <- load.reference.map(ref = "CD8T_human_ref_v1.rds")
palette <- ref@misc$atlas.palette
DimPlot(ref, label = T)
DefaultAssay(CD8.subset) <- "RNA3"

sample.list <- SplitObject(CD8.subset, split.by = "cohort")
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

OV=sample.class$OV@meta.data%>%filter(functional.cluster !="CD8.MAIT")%>%group_by(functional.cluster)%>%summarize(no=n())%>%mutate(Percentage= 100*no/sum(no))%>%mutate(Cohort= paste("OV"))
LUAD=sample.class$LUAD@meta.data%>%filter(functional.cluster !="CD8.MAIT")%>%group_by(functional.cluster)%>%summarize(no=n())%>%mutate(Percentage= 100*no/sum(no))%>%mutate(Cohort= paste("LUAD"))

data=data.frame(rbind(OV,LUAD))

data$Cohort <- factor(data$Cohort, levels=c("OV", "LUAD"))
data$functional.cluster <- factor(data$functional.cluster, levels=c("CD8.CM", "CD8.EM","CD8.TEMRA","CD8.TPEX","CD8.TEX","CD8.NaiveLike"))

ggplot(data,
       aes(x = Cohort, stratum = functional.cluster, alluvium = functional.cluster,y=Percentage,
           fill = functional.cluster, label = functional.cluster)) +
  scale_fill_manual(values = c("#99FF00","#FF3300","#FFFF00","#66CCCC","#663333","#FF9933","#00FF99","pink","#006600","#0099CC")) +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  theme(legend.position = "bottom")+theme_classic()

# Fig.1M; CellChat of CD8, NK cells

CD8=merge(sample.class$OV,y=sample.class$LUAD,add.cell.ids=c("OV","LUAD"))
CD8.NK=merge(CD8,y=nk.subset)

#The "CellChat.group" column was added to the metadata to reflect information regarding the corresponding cells from TillPred analyses and NK cells

data.list <- SplitObject(CD8.NK, split.by = "cohort")
OV=data.list$OV
LUAD=data.list$LUAD


#----------------------------------CellChat_Ovarina cancer --------------------
Idents(OV)="CellChat.group"
OV[["RNA3"]] <- as(object = OV[["RNA"]], Class = "Assay")

data.input = OV[["RNA3"]]$data
meta = data.frame( labels=OV@meta.data$CellChat.group,row.names=colnames(data.input) )

cellchat.OV<- createCellChat( object=data.input,meta=meta,group.by="labels" )
cellchat.OV<- addMeta(cellchat.OV, meta = meta)
cellchat.OV<- setIdent(cellchat.OV, ident.use = "labels") # set "labels" as default cell identity
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


ptm = Sys.time()
groupSize <- as.numeric(table(cellchat.OV@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.OV@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.OV@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat.OV@net$weight
pdf( "OV--circlePerCellPopulation.pdf",width=20,height=25 )
 par(mfrow = c(4,5), xpd=TRUE)
# pepap # tmp.name <- rownames(mat)[1]
 for ( i in 1:nrow(mat) ) {
  print(i)
  mat2      <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
# pepap #  mat2[1, ] <- mat[i, ]
# pepap #  mat2[ ,1] <- mat[ ,i]
# pepap #  rownames(mat2)[1] <- rownames(mat)[i]
# pepap #  colnames(mat2)[1] <- colnames(mat)[i]
# pepap #  rownames(mat2)[i] <- tmp.name
# pepap #  colnames(mat2)[i] <- tmp.name
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
 }
dev.off()

#----------------------------------CellChat_LUNG cancer --------------------
Idents(LUAD)="CellChat.group"
LUAD[["RNA3"]] <- as(object = LUAD[["RNA"]], Class = "Assay")

data.input = LUAD[["RNA3"]]$data
meta = data.frame( labels=LUAD@meta.data$CellChat.group,row.names=colnames(data.input) )

cellchat.LUAD<- createCellChat( object=data.input,meta=meta,group.by="labels" )
cellchat.LUAD<- addMeta(cellchat.LUAD, meta = meta)
cellchat.LUAD<- setIdent(cellchat.LUAD, ident.use = "labels") # set "labels" as default cell identity
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



ptm = Sys.time()
groupSize <- as.numeric(table(cellchat.LUAD@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.LUAD@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.LUAD@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat.LUAD@net$weight
pdf( "LUAD--circlePerCellPopulation.pdf",width=20,height=25 )
 par(mfrow = c(4,5), xpd=TRUE)
# pepap # tmp.name <- rownames(mat)[1]
 for ( i in 1:nrow(mat) ) {
  print(i)
  mat2      <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
# pepap #  mat2[1, ] <- mat[i, ]
# pepap #  mat2[ ,1] <- mat[ ,i]
# pepap #  rownames(mat2)[1] <- rownames(mat)[i]
# pepap #  colnames(mat2)[1] <- colnames(mat)[i]
# pepap #  rownames(mat2)[i] <- tmp.name
# pepap #  colnames(mat2)[i] <- tmp.name
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
 }
dev.off()



#--------------------------------------------------------------------------------------
#						Figure 2
#--------------------------------------------------------------------------------------

#combining all processed samples (n=36)

merged<- NormalizeData(merged_seurat)
merged<- FindVariableFeatures(merged)
merged<- ScaleData(merged)
merged<- RunPCA(merged)

# visualize the results of a standard analysis without integration

merged<- FindNeighbors(merged, dims = 1:30, reduction = "pca")
merged<- FindClusters(merged, cluster.name = "unintegrated_clusters")
merged<- RunUMAP(merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

options(future.globals.maxSize = 20000 * 1024^4)

merged<- IntegrateLayers(
  object = merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

# Harmony
merged<- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged<- FindClusters(merged, cluster.name = "harmony_clusters")
merged<- RunUMAP(merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

Idents(merged)<-"orig.ident"
merged@meta.data=merged@meta.data%>%mutate(Group=case_when(orig.ident=="S1"~"omentum_met_distant",
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
orig.ident=="S31"~"benign_omentum",
orig.ident=="S32"~"benign_omentum",
orig.ident=="S33"~"benign_omentum",
orig.ident=="S34"~"benign_omentum",
orig.ident=="S35"~"benign_omentum",
orig.ident=="S36"~"benign_omentum"))

# CellDex
sample<- as.SingleCellExperiment(merged)
IGD<-celldex::BlueprintEncodeData()
pred.IGD<- SingleR(test = sample, ref =IGD, assay.type.test=1,labels = IGD$label.main)
merged <- AddMetaData( object=merged,metadata=make.names(sub( " [(].*$","",pred.IGD$pruned.labels )),col.name="BlueprintEncodeData" )

ClusterDistrBar(merged$Group, merged$BlueprintEncodeData)


#--------------------------------------NK.cells-----------------------------------------------------
Idents(merged)="BlueprintEncodeData"
NK=subset(merged,idents=c("NK.cells"))

NK[["RNA"]] <- split(NK[["RNA"]], f = NK$orig.ident)

NK<- NormalizeData(NK)
NK<- FindVariableFeatures(NK)
NK<- ScaleData(NK)
NK<- RunPCA(NK)

NK<- FindNeighbors(NK, dims = 1:30, reduction = "pca")
NK<- FindClusters(NK, cluster.name = "unintegrated_clusters")
NK<- RunUMAP(NK, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

options(future.globals.maxSize = 20000 * 1024^4)

NK<- IntegrateLayers(
  object = NK, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

NK[["RNA"]] <- JoinLayers(NK[["RNA"]])

# Harmony
NK<- FindNeighbors(NK, reduction = "harmony", dims = 1:30)
NK<- FindClusters(NK, cluster.name = "harmony_clusters")
NK<- RunUMAP(NK, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

#Fig.2C

Genes=rev(c("SRGAP3",	"CD44",	"CASK",	"ITGA1",	"CLNK",
"LDLRAD4",	"PDE4B",	"COTL1",	"XCL1",	"CD96",	"ATP8B4",
"LDB2",	"SAMSN1",	"KRT81",	"ZBTB20",	"SFMBT2",	"PDE7A",
"RBPJ",	"ZNF331",	"GZMK","NCAM1","NCR1",	"FGFBP2",	"SPON2",
"FCGR3A",	"NKG7",	"PRF1",	"GZMH",	"CST7",
"PLAC8",	"CYBA",	"GZMB",	"CX3CR1",	"S100A4",	"KLRF1",
"ITGB2",	"PFN1",	"S1PR5",	"PRSS23",	"IGFBP7",	"ARL4C",
"AKR1C3"))

genes.zscore <- CalcStats(
  NK,
  features = Genes,
  group.by = "Group")

mat = as.matrix(genes.zscore)
colnames(mat)=colnames(mat)
library(circlize)
library(ComplexHeatmap)
library(pheatmap)

set.seed(123456)
pdf("Heatmap.pdf",width=3,height=10)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,cluster_rows=F,legend_breaks = -2:2,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()

#Fig.2D

d= c("IFNG", "CCL4","CCL4L2","CCL5","IL32","IL16","CCL3","XCL2","XCL1","FLT3LG",
"KLRK1", "KLRC2","CD160","NCR3","SLAMF6","SLAMF7",
"KLRC1","CD300A","TIGIT","KIR3DL2","KIR3DL1","KIR2DL3","KIR2DL1","KLRB1","SIGLEC7","HAVCR2",
"TGFBR1","IL12RB1","IL2RG","TGFBR3","TGFBR2","IL10RA","IL18R1","IL18RAP","IL2RB","IL10RB",
"FASLG","GSDMD","GZMB","NKG7","PRF1","GZMA","GZMH","GZMK","TNFSF10")

Idents(NK)="Group"
levels(NK) <- c("benign_omentum", "omentum_met_distant", "omentum_met_peritumoral","omentum_met_tumor")
cal <- CalcStats(NK, features = d, method = "zscore")

mat = as.matrix(cal)
colnames(mat)=colnames(mat)

set.seed(123456)
pdf("Heatmap_of_selected_genes.pdf",width=3.6,height=10)
ComplexHeatmap::pheatmap(mat,cluster_cols=F,legend_breaks = -2:2,
clustering_method="ward.D2",show_parent_dend_line=F,row_title = NULL,
    color = rev(brewer.pal(n = 11, name = "RdYlBu")))

dev.off()




