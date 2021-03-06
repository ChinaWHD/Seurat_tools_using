�ҵ���ʶ������о�AD���ں�����С����ϸ�����׺ʹ��׷���Ĺ���̽Ѱ���������߷�Ӧ������
����gse103334

#install.packages("Seurat")
#install.packages("outliers")
#install.packages("pbmcapply")
#install.packages("doFuture")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("singscore")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")

#install.packages("devtools")
#library(devtools)
#devtools::install_github('dviraran/SingleR')
########
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("SingleR")



###################################04.����ǰ�ڴ����ͽ���###################################
#��ȡ����
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)

setwd("C:\\Users\\lexb4\\Desktop\\scRNA\\04-07.Seurat")             #���ù���Ŀ¼

#��ȡ�ļ��������ظ�����ȡ��ֵ
rt=read.table("geneMatrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#������ת��ΪSeurat���󣬲������ݽ��й���
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
#ʹ��PercentageFeatureSet�����������������İٷֱ�
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6)           #�����������С����ͼ
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #�����ݽ��й���

#������ȵ�����Ի�ͼ
pdf(file="04.featureCor.pdf",width=10,height=6)              #����������������ͼ
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#�����ݽ��б�׼��
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#��ȡ��Щ��ϸ�������ϵ���ϴ�Ļ���
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#�����������ͼ
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)              #���������������ͼ
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()





###################################05.PCA���ɷַ���###################################
##PCA����
pbmc=ScaleData(pbmc)                     #PCA��ά֮ǰ�ı�׼Ԥ��������
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA����

#����ÿ��PCA�ɷֵ���ػ���
pdf(file="05.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#���ɷַ���ͼ��
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#���ɷַ�����ͼ
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#ÿ��PC��pֵ�ֲ��;��ȷֲ�
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()





###################################06.TSNE���������marker����###################################
##TSNE�������
pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #�����ڽӾ���
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #��ϸ������,�Ż���׼ģ�黯
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)                      #TSNE����
pdf(file="06.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 2, label = TRUE)    #TSNE���ӻ�
dev.off()
write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)

##Ѱ�Ҳ�����������
logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="06.markers.xls",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#����marker�ڸ���cluster����ͼ
pdf(file="06.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

#����marker��С����ͼ
pdf(file="06.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = c("IGLL5", "MBOAT1"))
dev.off()

#����marker�ڸ���cluster��ɢ��ͼ
pdf(file="06.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = c("IGLL5", "MBOAT1"),cols = c("green", "red"))
dev.off()

#����marker�ڸ���cluster������ͼ
pdf(file="06.markerBubble.pdf",width=12,height=6)
cluster10Marker=c("MBOAT1", "NFIB", "TRPS1", "SOX4", "CNN3", "PIM2", "MZB1", "MS4A1", "ELK2AP", "IGLL5")
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()





###################################07.ע��ϸ������###################################
library(SingleR)
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
singler = CreateSinglerObject(counts, annot = ann, "pbmc", min.genes = 0,
  species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
  reduce.file.size = T, numCores = 1)
singler$seurat = pbmc
singler$meta.data$xy = pbmc@reductions$tsne@cell.embeddings
clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
write.table(clusterAnn,file="07.clusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(singler$other,file="07.cellAnn.txt",quote=F,sep="\t",col.names=F)

#׼��monocle������Ҫ���ļ�
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="07.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="07.monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="07.monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(sig.markers,file="07.monocleMarkers.txt",sep="\t",row.names=F,quote=F)

