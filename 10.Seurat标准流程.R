library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- NormalizeData(pbmc)#（同上一步骤相同，用于标准化数据）

#识别高度可变的特征（特征选择）
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#缩放数据 
#移动每个基因的表达，以使整个细胞的平均表达为0
#缩放每个基因的表达，从而使细胞之间的方差为1
#此步骤在下游分析中具有相等的权重，因此高表达的基因不会占主导地位
#结果存储在 pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#耗时太长，可以使用ScaleData()是对先前标识的可变特征执行默认缩放（默认为2,000）
pbmc <- ScaleData(pbmc)

#与Seurat v2中一样，可以“消退”与（例如）细胞周期阶段或线粒体污染相关的异质性。
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

#执行线性尺寸缩减
#默认情况下，仅将先前确定的变量特征用作输入，但是features如果您希望选择其他子集，则可以使用arguments进行定义。
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#瑟拉提供可视化细胞和定义PCA，包括功能的几种有用的方法VizDimReduction()，DimPlot()和DimHeatmap()
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#可视化前两个PC
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

#DimHeatmap()可以轻松探索数据集中异质性的主要来源，并且在尝试确定要包括哪些PC以便进行下游分析。
#细胞和基因均根据其PCA分数排序。设置cells会在频谱的两端绘制“极端”细胞，从而极大地加快了大型数据集的绘制速度。
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)#绘制前15个PC热图

#确定数据集的“维数”
#为了克服scRNA-seq数据的任何单个基因中的广泛技术噪声，Seurat根据其PCA分数对细胞进行聚类，
#每个PC本质上代表一种“元特征”，该特征将跨相关特征集的信息进行组合。
#在Macosko等人中，我们实施了受JackStraw程序启发的重采样测试。我们随机置换数据的一部分（默认为1％），
#然后重新运行PCA，以构建特征分数的“零分布”，然后重复此过程。我们将“重要”的PC识别为具有丰富的低p值功能的PC。
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time（花费时间较长）
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#另一种启发式方法会生成“肘形图”：基于每个分量（ElbowPlot()函数）所解释的方差百分比对主成分进行排序。
ElbowPlot(pbmc)#显示重要的PC

#聚集细胞 
#首先基于PCA空间中的欧式距离构造一个KNN图，
#并根据两个细胞在当地社区的共享重叠(Jaccard相似性)来优化它们之间的边权值。
#使用该FindNeighbors()功能执行此步骤
#为了对细胞进行聚类，我们接下来将应用模块化优化技术（例如Louvain算法（默认）和SLM）
#将细胞迭代地分组在一起，以优化标准模块化功能。该FindClusters()函数实现此过程，
#并包含一个分辨率参数，该参数设置下游集群的“粒度”，其值越大，导致集群的数量越多。
#可以使用Idents()函数找到簇。
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#运行非线性降维（UMAP / tSNE），以可视化和探索这些数据集。
#作为UMAP和tSNE的输入，我们建议使用相同的PC作为聚类分析的输入。
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

#可以在此时保存该对象，以便可以轻松地将其重新装入，而不必重新运行上面执行的计算密集型步骤
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

#查找差异表达的特征（簇生物标志物）
#Seurat可以帮助找到通过差异表达定义聚类的标记。
#默认情况下，ident.1与所有其他细胞相比，它识别单个群集的正向和负向标记。
#FindAllMarkers()自动执行所有集群的此过程，但是您也可以测试集群组之间的相互关系，或针对所有细胞进行测试
#min.pct参数要求在两组细胞的任何一组中以最小百分比检测基因
#thresh.test参数要求在两组细胞之间平均表达一定程度的基因。
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#Seurat可以使用test.use参数设置一些针对差异表达的测试，
#例如，ROC测试返回任何单个标记的“分类能力”（范围从0-随机，到1-完美）。
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#提供了几种可视化标记表达的工具。
#VlnPlot()（显示整个集群中的表达概率分布），以及FeaturePlot()（在tSNE或PCA图上可视化特征表达）是我们最常用的可视化。
#我们还建议探索RidgePlot()，CellScatter()和DotPlot()作为查看数据集的其他方法。
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

#DoHeatmap()生成给定单元和特征的表达式热图。在这种情况下，我们将为每个聚类绘制前20个标记（如果小于20，则绘制所有标记）。
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#将细胞类型标识分配给集群
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()














































































