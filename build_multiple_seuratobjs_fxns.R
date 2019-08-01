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
preplots = function(VAR){
  print(VlnPlot(VAR, features = c(colnames(VAR@meta.data)[2:length(colnames(VAR@meta.data))]), ncol = 5)) 
  png(paste0(Sys.Date(),"-prefilter-violin-",deparse(substitute(VAR)),".png"), width = 1500, height = 600)
  print(VlnPlot(VAR, features = c(colnames(VAR@meta.data)[2:length(colnames(VAR@meta.data))]), ncol = 5))
  dev.off()
  plot1a = FeatureScatter(object = VAR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  
  plot1b = FeatureScatter(object = VAR, feature1 = "nCount_RNA", feature2 = "percent.blood")
  plot1c = FeatureScatter(object = VAR, feature1 = "nCount_RNA", feature2 = "percent.mito")
  print(CombinePlots(list(plot1a, plot1b, plot1c),ncol = 3))
  png(paste0(Sys.Date(),"-prefilter-scatterplots-",deparse(substitute(VAR)),".png"), width = 1300, height = 700)
  print(CombinePlots(list(plot1a, plot1b, plot1c),ncol = 3))
  dev.off()
}

postplots = function(VAR){
  print(VlnPlot(VAR, features = c(colnames(VAR@meta.data)[2:length(colnames(VAR@meta.data))]), ncol = 5))
  png(paste0(Sys.Date(),"-postfilter-violin-",deparse(substitute(VAR)),".png"), width = 1500, height = 600)
  print(VlnPlot(VAR, features = c(colnames(VAR@meta.data)[2:length(colnames(VAR@meta.data))]), ncol = 5))
  dev.off()
  plot1a = FeatureScatter(object = VAR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  plot1b = FeatureScatter(object = VAR, feature1 = "nCount_RNA", feature2 = "percent.blood")
  plot1c = FeatureScatter(object = VAR, feature1 = "nCount_RNA", feature2 = "percent.mito")
  print(CombinePlots(list(plot1a, plot1b, plot1c),ncol = 3))
  png(paste0(Sys.Date(),"-postfilter-scatterplots-",deparse(substitute(VAR)),".png"), width = 1300, height = 700)
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