#load libraries
library(harmony)
library(Seurat)
library(Matrix)
library(SeuratDisk)
library(tidyverse)
library(reshape2)

#load Helper functions 
#source("helper_functions.R")

#Create Importer function
Importer <- function(pathway,id, TenX=TRUE, performNormalisation=TRUE, performScaling = FALSE,performVariableGeneDetection=TRUE) {
  if (TenX) {
    Matrix <- Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  seuratObject =CreateSeuratObject(counts = Matrix, project = id, min.cells = 5)
  seuratObject$sample <- id
  tmp<-unlist(strsplit(id,split = "-"))
  seuratObject$condition <- paste0(tmp[1:length(tmp)-1],collapse = "-")
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > 200)
  if (performNormalisation==TRUE) {
    seuratObject<-NormalizeData(object = seuratObject,verbose = FALSE)
  }
  if(performVariableGeneDetection){
    seuratObject<-FindVariableFeatures(object = seuratObject, do.plot = FALSE, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  if (performScaling==TRUE) {
    seuratObject<-ScaleData(object = seuratObject)
  }
  cat("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n")
  return(seuratObject)
}

#Load Object
#Download the file for the Heart Cell ATLAS here: https://cellgeni.cog.sanger.ac.uk/heartcellatlas/data/global_raw.h5ad

#convert and import the data
Convert("global.h5ad","h5Seurat")
data_hca <- LoadH5Seurat("/global.h5seurat")

#Filter nuclei and septum to match the diseased patients by region
names(data_hca@meta.data)[6] <- c("orig.ident")
data_hca <- data_hca[,data_hca$cell_source %in% c("Harvard-Nuclei","Sanger-Nuclei")]
data_hca <- data_hca[,data_hca$region=="SP"]
data_hca$dataset <- factor(rep("hca",ncol(data_hca)))


#import data for AS patients (plase download here: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11268/samples/)

#!!! Adjust the Paths accordingly
Sample.Paths <- c("/media/Helios_scStorage/Luka/103837/cellranger/103837-001-001/outs/filtered_feature_bc_matrix/",
                  "/media/Helios_scStorage/Luka/103837/cellranger/103837-001-003/outs/filtered_feature_bc_matrix/",
                  "/media/Helios_scStorage/Luka/103837/cellranger/103837-001-004/outs/filtered_feature_bc_matrix/",
                  "/media/Helios_scStorage/Luka/Frankfurt_Human_Biopsies/2nd_Try_2019_08_20/103799-001-001/outs/filtered_feature_bc_matrix/",
                  "/media/Helios_scStorage/Luka/103837/cellranger/103837-001-002/outs/filtered_feature_bc_matrix/")
Samplenames <- c("3011","3015","3018","3005","3013")

SeuratObjectList <- list()
for (i in 1:length(Sample.Paths)) {
  SeuratObjectList[[i]]<-Importer(pathway = Sample.Paths[i],id = Samplenames[i])
}

#integrate Data
combined.anchors <- FindIntegrationAnchors(object.list = SeuratObjectList, dims = 1:10, k.filter = 150)
combined <- IntegrateData(anchorset = combined.anchors, dims = 1:10)
combined <- ScaleData(object = combined, verbose = FALSE)

#Filter data
combined_filtered <- subset(combined, subset = percent.mito < 10 & nFeature_RNA < 2500 & nFeature_RNA > 200)
combined_filtered$dataset <- factor(rep("aort_sten",ncol(combined_filtered)))
combined_filtered$region <- factor(rep("SP",ncol(combined_filtered)))

#combine datasets
Object <- merge(data_hca,combined_filtered,c("hca","aort_sten"))

Object <- NormalizeData(Object) 
Object <- FindVariableFeatures(Object)
Object <- ScaleData(Object)
Object <- RunPCA(Object)

#Run Harmony integration
Object <- RunHarmony(Object, "dataset", plot_convergence = T)
Object <- RunUMAP(Object, reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.35)

saveRDS(object = Object, file = "~/seurat_object_hca_as_harmonized_AS_SP_nuc_refined_cells.rds")
Object <- readRDS("~/seurat_object_hca_as_harmonized_AS_SP_nuc_refined_cells.rds")
Object@meta.data$dataset <- factor(Object@meta.data$dataset, levels = c("hca", "aort_sten"))




#Figure 1

#Figure 1a
Idents(Object) <- "seurat_clusters"
DimPlot(Object, label = T)

#Figure 1b
CM_Object <- subset(Object, subset = seurat_clusters %in% c(0,1,3,6))
VlnPlot(CM_Object, group.by = "dataset", features = c("MYH6","MYH7"))
Idents(CM_Object) <- "dataset"
sig <- FindAllMarkers(CM_Object, test.use = "bimod", only.pos = TRUE, features = c("MYH6","MYH7"))

#Figure 1c
hyp_genes <- c("NPPA","NPPB","MYH7","MYH7B","XIRP2","PFKP","CMYA5","ACTA1","TNNI3","ANKRD1")
Idents(CM_Object) <- "dataset"
hyp_df <- data.frame(Barcode = row.names(CM_Object@meta.data),
                     Condition = CM_Object@meta.data$dataset,
                     NPPA = CM_Object@assays$RNA@data["NPPA",],
                     NPPB = CM_Object@assays$RNA@data["NPPB",],
                     MYH7 = CM_Object@assays$RNA@data["MYH7",],
                     MYH7B = CM_Object@assays$RNA@data["MYH7B",],
                     XIRP2 = CM_Object@assays$RNA@data["XIRP2",],
                     PFKP = CM_Object@assays$RNA@data["PFKP",],
                     CMYA5 = CM_Object@assays$RNA@data["CMYA5",],
                     ACTA1 = CM_Object@assays$RNA@data["ACTA1",],
                     TNNI3 = CM_Object@assays$RNA@data["TNNI3",],
                     ANKRD1 = CM_Object@assays$RNA@data["ANKRD1",])
hyp_df_score <- data.frame(hyp_df,
                           Disease_Score = rowSums(hyp_df[3:12]))

Disease_Score <- hyp_df_score$Disease_Score
names(Disease_Score) <- hyp_df_score$Barcode
CM_Object <- AddMetaData(object = CM_Object, metadata = Disease_Score, col.name = "Disease_Score")
VlnPlot(CM_Object, group.by = "dataset", features = "Disease_Score")
sig <- t.test(subset(hyp_df_score, subset = Condition == "hca")$Disease_Score,
              subset(hyp_df_score, subset = Condition == "aort_sten")$Disease_Score)

#Figure 1d
library(Seurat)
CM_Object <- FindVariableFeatures(CM_Object)
CM_Object <- ScaleData(object = CM_Object)
CM_Object <- RunPCA(object = CM_Object)
CM_Object <- RunUMAP(object = CM_Object, reduction = "harmony", dims = 1:12)
CM_Object <- FindNeighbors(object = CM_Object, reduction = "harmony", dims = 1:12)
CM_Object <- FindClusters(CM_Object, resolution = 0.3)
Idents(CM_Object) <- "RNA_snn_res.0.3"
library(viridis)
DimPlot(CM_Object, label = T, cols = c(viridis_pal()(6)))

#Figure 1e
Idents(CM_Object) <- "RNA_snn_res.0.3"
cm_clusters_markers <- FindAllMarkers(CM_Object)

#Figure 1f
perc <- NULL
cl <- NULL
cond <- NULL
for (i in levels(CM_Object@meta.data$RNA_snn_res.0.3)){
  for (ii in levels(CM_Object@meta.data$dataset)){
    perc <- c(perc, length(row.names(subset(CM_Object, subset = dataset == ii & RNA_snn_res.0.3 == i)@meta.data))/
                length(row.names(subset(CM_Object, subset = dataset == ii)@meta.data)))
    cl <- c(cl, i)
    cond <- c(cond, ii)
  }}
cl_repr_df <- data.frame(Cluster = cl,
                         Condition = cond,
                         Percentage = perc)
cl_repr_df$Condition <- factor(cl_repr_df$Condition, levels = c("hca", "aort_sten"))
ggplot(cl_repr_df, aes(x=Condition, y=Percentage, fill=Cluster))+
  geom_bar(stat="identity")


#Figure 1g
VlnPlot(CM_Object, group.by = "dataset", features = c("ERBB4","TBX20","GATA4","ADRA1A","NLGN1","FGF12"))
Idents(CM_Object) <- "dataset"
sig <- FindAllMarkers(CM_Object, test.use = "bimod", only.pos = TRUE, features = c("ERBB4","TBX20","GATA4","ADRA1A","NLGN1","FGF12"))

#Figure 1h
VlnPlot(CM_Object, group.by = "dataset", features = c("VEGFA","VEGFB"))
Idents(CM_Object) <- "dataset"
sig <- FindAllMarkers(CM_Object, test.use = "bimod", only.pos = TRUE, features = c("VEGFA","VEGFB"))

#Figure 1i
sample_ids <- levels(as.factor(CM_Object@meta.data$sample))
ind_means <- NULL
sem <- NULL
cond <- NULL
for (i in sample_ids) {
  ind_means <- c(ind_means, mean(expm1(subset(CM_Object, subset = sample == i)@assays$RNA@data["VEGFA",])))
  sem <- c(sem, sd(expm1(subset(CM_Object, subset = sample == i)@assays$RNA@data["VEGFA",]))/
             sqrt(length(row.names(subset(CM_Object, subset = sample == i)@meta.data))))
  cond <- c(cond, as.character(subset(CM_Object, subset = sample == i)@meta.data$dataset[1]))
}
data_ind <- data.frame(ID = sample_ids, VEGFA_mean = ind_means, SEM = sem, Condition = cond)
data_ind$Condition <- factor(data_ind$Condition, levels = c("hca", "aort_sten"))
ggplot(data_ind, aes(x=Condition, y=VEGFA_mean, fill = Condition))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values = c("lightblue","coral2"))+
  theme_classic()+
  xlab("")+
  ylab("Mean Expression per Patient")+
  theme(axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 18))
t.test(data_ind$VEGFA_mean ~ data_ind$Condition)
sample_ids <- levels(as.factor(CM_Object@meta.data$sample))
ind_means <- NULL
sem <- NULL
cond <- NULL
for (i in sample_ids) {
  ind_means <- c(ind_means, mean(expm1(subset(CM_Object, subset = sample == i)@assays$RNA@data["VEGFB",])))
  sem <- c(sem, sd(expm1(subset(CM_Object, subset = sample == i)@assays$RNA@data["VEGFB",]))/
             sqrt(length(row.names(subset(CM_Object, subset = sample == i)@meta.data))))
  cond <- c(cond, as.character(subset(CM_Object, subset = sample == i)@meta.data$dataset[1]))
}
data_ind <- data.frame(ID = sample_ids, VEGFB_mean = ind_means, SEM = sem, Condition = cond)
data_ind$Condition <- factor(data_ind$Condition, levels = c("hca", "aort_sten"))
ggplot(data_ind, aes(x=Condition, y=VEGFB_mean, fill = Condition))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values = c("lightblue","coral2"))+
  theme_classic()+
  xlab("")+
  ylab("Mean Expression per Patient")+
  theme(axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 18))
t.test(data_ind$VEGFB_mean ~ data_ind$Condition)




#Figure 2

#Figure 2a
Idents(Object)<-"cell_type"
seurat_obj.noLow <- subset(Object, idents = c("CM","PC","FB","EC","Myeloid","SMC","Lympoid","Neuro"))

Idents(seurat_obj.noLow)<-"dataset"
seurat_obj.noLow.AS<-subset(seurat_obj.noLow, idents = "aort_sten")
seurat_obj.noLow.HCA<-subset(seurat_obj.noLow, idents = "hca")

#output matrix and Cellname and Celltype
write.table(x = data.frame(seurat_obj.noLow.AS@assays$RNA@data), file = "noLowQCClusters/AS/AS.counts.noLowQC.txt", col.names = NA, sep = "\t", dec = ".")
write.table(x = data.frame(seurat_obj.noLow.HCA@assays$RNA@data), file = "noLowQCClusters/HCA/HCA.counts.noLow.txt", col.names = NA, sep = "\t", dec = ".")

meta.AS<-data.frame(Cell=names(seurat_obj.noLow.AS$cell_type), celltype=seurat_obj.noLow.AS$cell_type)
write.table(meta.AS, "noLowQCClusters/AS/AS.meta.noLow.txt", quote = FALSE, row.names = FALSE, sep = "\t")

meta.HCA<-data.frame(Cell=colnames(tmp.noLow), celltype=seurat_obj.noLow.HCA$cell_type)
write.table(meta.HCA, "noLowQCClusters/HCA/HCA.meta.noLow.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#!!! To run cellphone DB, the following commands were run in bash
  #AS
    #cellphonedb method statistical_analysis AS.meta.txt AS.counts.txt --counts-data=gene_name --iterations=100 --threads=20
    #cd AS
    #cellphonedb plot dot_plot
    #cellphonedb plot heatmap_plot AS.meta.txt
  
  #HCA
    #cd ../HCA
    #cellphonedb method statistical_analysis HCA.meta.txt HCA.counts.txt --counts-data=gene_name --iterations=100 --threads=20
    #cd HCA
    #cellphonedb plot dot_plot
    #cellphonedb plot heatmap_plot HCA.meta.txt

#Figure 2b
sig_healthy_inter <- read.table(file = "~/interaction_count.txt")
sig_as_inter <- read.table(file = "~/interaction_count.txt")
total_int_count <- data.frame(Condition = c("Healthy", "AS"),
                              Interactions = c(sum(sig_healthy_inter$all_sum),
                                               sum(sig_as_inter$all_sum)))
total_int_count$Condition <- factor(total_int_count$Condition, levels = c("Healthy", "AS"))
ggplot(total_int_count, aes(x = Condition, y = Interactions))+
  geom_bar(stat = "identity")

#Figure 2c
cm_healthy_inter <- as.numeric(sig_healthy_inter["CM",])
cm_as_inter <- as.numeric(sig_as_inter["CM",])
cm_int_count <- data.frame(Condition = c("Healthy", "AS"),
                           Interactions = c(cm_healthy_inter, 
                                            cm_as_inter))
cm_int_count$Condition <- factor(cm_int_count$Condition, levels = c("Healthy", "AS"))
ggplot(cm_int_count, aes(x = Condition, y = Interactions))+
  geom_bar(stat = "identity")

#Figure 2d
total_int <- data.frame(Celltype = c(rep(c("EC", "PC", "FB", "SMC", "NLC", "MC", "CM", "LC"), times = 2)),
                        Interactions = c(as.numeric(sig_healthy_inter["EC",]),
                                    as.numeric(sig_healthy_inter["PC",]),
                                    as.numeric(sig_healthy_inter["FB",]),
                                    as.numeric(sig_healthy_inter["SMC",]),
                                    as.numeric(sig_healthy_inter["Neuro",]),
                                    as.numeric(sig_healthy_inter["Myeloid",]),
                                    as.numeric(sig_healthy_inter["CM",]),
                                    as.numeric(sig_healthy_inter["Lympoid",]),
                                    as.numeric(sig_as_inter["EC",]),
                                    as.numeric(sig_as_inter["PC",]),
                                    as.numeric(sig_as_inter["FB",]),
                                    as.numeric(sig_as_inter["SMC",]),
                                    as.numeric(sig_as_inter["Neuro",]),
                                    as.numeric(sig_as_inter["Myeloid",]),
                                    as.numeric(sig_as_inter["CM",]),
                                    as.numeric(sig_as_inter["Lympoid",])), 
                        Condition = c(rep("Healthy", times = 8),
                                      rep("AS", times = 8)))
total_int$Condition <- factor(total_int$Condition, levels = c("Healthy", "AS"))
ggplot(total_int, aes(x = Celltype, y = Interactions, fill = Condition))+
  geom_bar(stat = "identity", position = position_dodge2())

#Figure 2e
healthy <-  read.table(file = "~/For R_Healthy.csv", header = T, dec = ",", sep = ";")
as <-  read.table(file = "~/For R_AS.csv", header = T, dec = ",", sep = ";")
ggplot(healthy, aes(x = reorder(interacting_pair, -CM.EC), y = CM.EC))+
  geom_bar(position = "dodge", stat = "identity", color = "black")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(colour = "black", size = 18))+
  xlab("")+
  ylab("Interaction Score")
ggplot(as, aes(x = reorder(interacting_pair, -CM.EC), y = CM.EC))+
  geom_bar(position = "dodge", stat = "identity", color = "black")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(colour = "black", size = 18))+
  xlab("")+
  ylab("Interaction Score")

#Figure 2f
healthy <-  read.table(file = "~/Documents/PhD/Cell Atlas - Healthy vs AS/Plots for Interactions CMs/FB For R_Healthy.csv", header = T, dec = ",", sep = ";")
as <-  read.table(file = "~/Documents/PhD/Cell Atlas - Healthy vs AS/Plots for Interactions CMs/FB For R_AS.csv", header = T, dec = ",", sep = ";")
ggplot(healthy, aes(x = reorder(interacting_pair, -CM.FB), y = CM.FB))+
  geom_bar(position = "dodge", stat = "identity", color = "black")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(colour = "black", size = 18))+
  xlab("")+
  ylab("Interaction Score")
ggplot(as, aes(x = reorder(interacting_pair, -CM_FB), y = CM_FB))+
  geom_bar(position = "dodge", stat = "identity", color = "black")+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(colour = "black", size = 18))+
  xlab("")+
  ylab("Interaction Score")




#Figure 3

#Figure 3a
Objects <- c(CM_Object)
EPHs <- c("EPHA1","EPHA2","EPHA3","EPHA4","EPHA5","EPHA6","EPHA7","EPHA8","EPHA10","EPHB1","EPHB2","EPHB3","EPHB4","EPHB6")
means <- NULL
a <- NULL
g <- NULL
for (j in EPHs) {
  for (i in Objects) {
    a <- mean(expm1(i@assays$RNA@data[j,]))
    means <- c(means, a)
    g <- c(g, j)
  }
}
df_new_barplots <- data.frame(Gene = g, Mean = means)
df_new_barplots <- df_new_barplots[order(-df_new_barplots$Mean),]
ggplot(df_new_barplots, aes(y = Mean, x = reorder(Gene, Mean)))+
  geom_bar(stat = "identity", color = "black", fill = viridis_pal()(14))+
  theme_classic()+
  coord_flip()+
  theme(axis.text = element_text(colour = "black", size = 16))

#Figure 3b
DoHeatmap(CM_Object, features = EPHs, slot = "data", group.by = "dataset", group.colors = c("blue","red"))+
  scale_fill_gradientn(colors = c("beige", "red", "darkred"))

#Figure 3c
CM <- CM_Object
means_Healthy <- c(mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA1",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA2",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA3",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA4",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA5",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA6",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA7",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA8",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHA10",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHB1",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHB2",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHB3",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHB4",])),
                   mean(expm1(subset(CM, subset = dataset == "hca")@assays$RNA@data["EPHB6",])))
means_Hypertr <- c(mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA1",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA2",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA3",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA4",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA5",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA6",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA7",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA8",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHA10",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHB1",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHB2",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHB3",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHB4",])),
                   mean(expm1(subset(CM, subset = dataset == "aort_sten")@assays$RNA@data["EPHB6",])))
df_EPHs <- data.frame(Celltype = rep(c("CM"), times = 14),
                      Gene = g,
                      Expression = means,
                      Healthy = means_Healthy,
                      Hypertr = means_Hypertr)
diff <- NULL
for (i in EPHs){
  diff <- c(diff, subset(df_EPHs, subset = Gene == i)$Hypertr - subset(df_EPHs, subset = Gene == i)$Healthy)
}
df_EPHs <- data.frame(df_EPHs, 
                      Diff = diff)
ggplot(df_EPHs, aes(x = reorder(Gene, -Diff), y = Diff))+
  geom_bar(stat = "identity", color = "black", fill = "darkred")+
  theme_classic()+
  coord_flip()+
  theme(axis.text = element_text(colour = "black", size = 16))

#Figure 3d
FeaturePlot(Object, features = "EPHB1", order = TRUE)

#Figure 3e
CM <- subset(Object, subset = seurat_clusters %in% c(0,1,3,6))
FB <- subset(Object, subset = seurat_clusters %in% c(4,11,12,16))
EC <- subset(Object, subset = seurat_clusters %in% c(5))
MC <- subset(Object, subset = seurat_clusters %in% c(9))
LC <- subset(Object, subset = seurat_clusters %in% c(13,18))
SMC <- subset(Object, subset = seurat_clusters %in% c(10))
PC <- subset(Object, subset = seurat_clusters %in% c(2))
NLC <- subset(Object, subset = seurat_clusters %in% c(15))
Objects <- c(CM,EC,FB,SMC,PC,MC,LC,NLC)
means <- NULL
a <- NULL
for (i in Objects) {
  a <- mean(expm1(i@assays$RNA@data["EPHB1",]))
  means <- c(means, a)
}

df_EPHB1 <- data.frame(Gene = rep("EPHB1", times = 8),
                       Celltype = c("CM","EC","FB","SMC","PC","MC","LC","NLC"),
                       Mean_EPHB1 = means)
df_EPHB1 <- df_EPHB1[order(-df_EPHB1$Mean_EPHB1),]
ggplot(df_EPHB1, aes(y = Mean_EPHB1, x = reorder(Celltype, Mean_EPHB1)))+
  geom_bar(stat = "identity", color = "black", fill = viridis_pal()(8))+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 16))

#Figure 3f
VlnPlot(CM_Object, group.by = "dataset", features = "EPHB1")
Idents(CM_Object) <- "dataset"
sig <- FindAllMarkers(CM_Object, test.use = "bimod", only.pos = TRUE, features = c("EPHB1"))

#Figure 3g
VlnPlot(CM_Object, group.by = "sample", features = "EPHB1")

#Figure 3h
sample_ids <- levels(as.factor(CM_Object@meta.data$sample))
ind_means <- NULL
sem <- NULL
cond <- NULL
for (i in sample_ids) {
  ind_means <- c(ind_means, mean(expm1(subset(CM_Object, subset = sample == i)@assays$RNA@data["EPHB1",])))
  sem <- c(sem, sd(expm1(subset(CM_Object, subset = sample == i)@assays$RNA@data["EPHB1",]))/
             sqrt(length(row.names(subset(CM_Object, subset = sample == i)@meta.data))))
  cond <- c(cond, as.character(subset(CM_Object, subset = sample == i)@meta.data$dataset[1]))
}
data_ind <- data.frame(ID = sample_ids, EPHB1_mean = ind_means, SEM = sem, Condition = cond)
data_ind$Condition <- factor(data_ind$Condition, levels = c("hca", "aort_sten"))
ggplot(data_ind, aes(x=Condition, y=EPHB1_mean, fill = Condition))+
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values = c("lightblue","coral2"))+
  theme_classic()+
  xlab("")+
  ylab("Mean Expression per Patient")+
  theme(axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 18))
t.test(data_ind$EPHB1_mean ~ data_ind$Condition)




#Figure 4

#Figure 4a
FeaturePlot(Object, features = "EFNB2")

#Figure 4b
CM <- subset(Object, subset = seurat_clusters %in% c(0,1,3,6))
FB <- subset(Object, subset = seurat_clusters %in% c(4,11,12,16))
EC <- subset(Object, subset = seurat_clusters %in% c(5))
MC <- subset(Object, subset = seurat_clusters %in% c(9))
LC <- subset(Object, subset = seurat_clusters %in% c(13,18))
SMC <- subset(Object, subset = seurat_clusters %in% c(10))
PC <- subset(Object, subset = seurat_clusters %in% c(2))
NLC <- subset(Object, subset = seurat_clusters %in% c(15))
Objects <- c(CM,EC,FB,SMC,PC,MC,LC,NLC)
means <- NULL
a <- NULL
for (i in Objects) {
  a <- mean(expm1(i@assays$RNA@data["EFNB2",]))
  means <- c(means, a)
}

df_EFNB2 <- data.frame(Gene = rep("EFNB2", times = 8),
                       Celltype = c("CM","EC","FB","SMC","PC","MC","LC","NLC"),
                       Mean_EFNB2 = means)
df_EFNB2 <- df_EFNB2[order(-df_EFNB2$Mean_EFNB2),]
ggplot(df_EFNB2, aes(y = Mean_EFNB2, x = reorder(Celltype, Mean_EFNB2)))+
  geom_bar(stat = "identity", color = "black", fill = viridis_pal()(8))+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 16))







