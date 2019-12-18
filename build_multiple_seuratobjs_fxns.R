# May we never copy and paste again.

packageVersion('Seurat')
#[1] "3.0.2"
packageVersion('SingleCellExperiment')
#[1] "1.6.0"

# read in Seurat Object, calc meta features of blood, mito, ribo ####

read_create_meta = function(dirvar,projvar){
  mat = Read10X(dirvar)
  seuratobj = CreateSeuratObject(counts = mat, project = projvar)
  seuratobj[["percent.blood"]] = PercentageFeatureSet(object = seuratobj, pattern = "^Hb")
  seuratobj[["percent.mito"]] = PercentageFeatureSet(object = seuratobj, pattern = "^mt-")
  seuratobj[["percent.ribo"]] = PercentageFeatureSet(object = seuratobj, pattern = c("^Rps"))
  seuratobj
}

# QC Pre and Post plotting functions ####
preplots = function(dir,seuratobj){
  print(VlnPlot(seuratobj, features = c(colnames(seuratobj@meta.data)[2:length(colnames(seuratobj@meta.data))]), ncol = 5)) 
  png(paste0(as.character(dir),Sys.Date(),"-prefilter-violin-",deparse(substitute(seuratobj)),".png"), width = 1500, height = 600)
  print(VlnPlot(seuratobj, features = c(colnames(seuratobj@meta.data)[2:length(colnames(seuratobj@meta.data))]), ncol = 5))
  dev.off()
  plot1a = FeatureScatter(object = seuratobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  
  plot1b = FeatureScatter(object = seuratobj, feature1 = "nCount_RNA", feature2 = "percent.blood")
  plot1c = FeatureScatter(object = seuratobj, feature1 = "nCount_RNA", feature2 = "percent.mito")
  print(CombinePlots(list(plot1a, plot1b, plot1c),ncol = 3))
  png(paste0(as.character(dir),Sys.Date(),"-prefilter-scatterplots-",deparse(substitute(seuratobj)),".png"), width = 1300, height = 700)
  print(CombinePlots(list(plot1a, plot1b, plot1c),ncol = 3))
  dev.off()
}

postplots = function(dir,seuratobj){
  print(VlnPlot(seuratobj, features = c(colnames(seuratobj@meta.data)[2:length(colnames(seuratobj@meta.data))]), ncol = 5)) 
  png(paste0(as.character(dir),Sys.Date(),"-postfilter-violin-",deparse(substitute(seuratobj)),".png"), width = 1500, height = 600)
  print(VlnPlot(seuratobj, features = c(colnames(seuratobj@meta.data)[2:length(colnames(seuratobj@meta.data))]), ncol = 5))
  dev.off()
  plot1a = FeatureScatter(object = seuratobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  
  plot1b = FeatureScatter(object = seuratobj, feature1 = "nCount_RNA", feature2 = "percent.blood")
  plot1c = FeatureScatter(object = seuratobj, feature1 = "nCount_RNA", feature2 = "percent.mito")
  print(CombinePlots(list(plot1a, plot1b, plot1c),ncol = 3))
  png(paste0(as.character(dir),Sys.Date(),"-postfilter-scatterplots-",deparse(substitute(seuratobj)),".png"), width = 1300, height = 700)
  print(CombinePlots(list(plot1a, plot1b, plot1c),ncol = 3))
  dev.off()
}

# Normalize, scale (scales all genes as written), find variable genes in seurat object ####

norm_scale_var = function(seuratobj,nfeat){
  seuratobj = NormalizeData(seuratobj)
  seuratobj = ScaleData(seuratobj, features = rownames(seuratobj))
  seuratobj = FindVariableFeatures(seuratobj, nfeatures = nfeat)
  seuratobj
}

# Convert multiple seurat objects to SCE objects ####
seurat_to_scefxn = function(seurobj,genelist) {
  sce = SingleCellExperiment(
    assays = list(counts = seurobj@assays$RNA@counts[genelist,],
                  logcounts  = seurobj@assays$RNA@data[genelist,],
                  scale.data = seurobj@assays$RNA@scale.data[genelist,]),
    colData = seurobj@meta.data[,1:3]
  )
  sce
}

# Dimplot wrapper which creates unique files names and also names the plot based on input parameters
dpwrapper = function(fileloc = fileloc,
                     seuratobj = seuratobj,
                     reduct = reduct,
                     groupmeta = groupmeta,
                     title = paste0(deparse(substitute(seuratobj))," ",reduct," ",groupmeta),
                     w = 1000,
                     h = 1000){
  plota = DimPlot(object = seuratobj,
                  reduction = reduct,
                  group.by = groupmeta,
                  pt.size = 2) + ggtitle(title)
  print(plota)
  png(paste0(fileloc,Sys.Date(),"-",deparse(substitute(seuratobj)),"-",reduct,"-",groupmeta,".png"), width = w, height = h)
  print(plota)
  dev.off()
}

# dpwrapper(fileloc = './directory/here/', #specified directory to place the file in
#           seuratobj = seuratobj, #specifies which seurat object to use
#           reduct = "umap", #could be umap, pca, tsne, etc
#           groupmeta = "orig.ident") #could be anything stored in metadata