#assumes normal pre-processing of Seurat object up "Cluster the cells" step of https://satijalab.org/seurat/pbmc3k_tutorial.html

seurat_object_var <- FindClusters(object = seurat_object_var, 
                             reduction.type = "pca", 
                             dims.use = 1:20, #must be the same number specified in RunTSNE()
                             resolution = seq(0.1,1,0.1), 
                             print.output = 0, 
                             save.SNN = TRUE,
                             force.recalc = TRUE) #calculate clusters for varying resolutions; seq() can be varied depending on the dataset

#create a dataframe of the resolution / cluster ID information from the meta.data slot of your Seurat object
res_dataframe = sapply(grep("^res",colnames(seurat_object_var@meta.data), value=T), function(x) length(unique(seurat_object_var@meta.data[,x])))
res_dataframe = as.data.frame(res_dataframe)
res_dataframe = tibble::rownames_to_column(res_dataframe, "res")
length(rownames(res_dataframe)) #add this value to the interior for loop of the nested pair below

#before running the below, amend x in 1:__ value
setwd("~/") #dir to output all these pngs
for (perp in seq(30,110,20)){ #sets "perp" to each of the values specified in seq()
  seurat_object_var <- RunTSNE(object = seurat_object_var, dims.use = 1:20, do.fast = TRUE, perplexity = perp) #alters perp for each of the different previously computed res's; ensure dims.use = the variable used in FindClusters()
  for (x in 1:10) {
   png(paste0("tSNE-perp_",perp,"-",res_dataframe$res[x],".png"), width = 700, height = 600) #creates a blank file titled with the perp and resolution number
   TSNEPlot(object = D1_G1_e11_5, group.by = res_dataframe$res[x], do.label = T, pt.size=1, plot.title = paste0("tSNE: perp ", perp, ", ", res_dataframe$res[x])) #plots the specific tsne with specific res/perp
   dev.off() #saves the photo
 }
}

seurat_object_var = SetAllIdent(seurat_object_var, id = "res.0.3") #use to switch between the res options previously computed
