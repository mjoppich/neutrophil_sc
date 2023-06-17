

source("/mnt/input/own/SeqMappings/Pekayvaz/21068_M12_macs/functions.R")
source("/mnt/extproj/projekte/pekayvaz/scripts/functions.R")




files <- Sys.glob("../calls/*/outs/raw*/*.mtx.gz")
files = grep(x=files, pattern="old", value=T, invert=T)
files = grep(x=files, pattern="sepsis", value=T)

allfiles.raw = list()
allABs.raw = list()

for (file in files)
{
  samplename = str_split(dirname(file), "/")[[1]][3]
  foldername = dirname(file)
  
  print(paste(samplename, foldername))
  
  h5file = Read10X(foldername,unique.features = TRUE)

  if (is.null(names(h5file)))
  {
      print(paste("WITHOUT AB", samplename))
    allfiles.raw[[samplename]] = h5file
  } else {
      print(paste("WITH AB", samplename))
    allfiles.raw[[samplename]] = h5file$`Gene Expression`
    allABs.raw[[samplename]] = h5file$`Antibody Capture`
  }

  print(paste(samplename, nrow(allfiles.raw[[samplename]]), "x", ncol(allfiles.raw[[samplename]]), "genes x cells"))
}



statDF = data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("RawAll", "RawFiltered"))))
objlist = list()

for (x in names(allfiles.raw))
{

    matrix = allfiles.raw[[x]]
    origSize = ncol(matrix)
    
    filteredAll = dim(matrix)[2]
    filteredObj = makeSeuratObj(matrix, x, patternList.human)
    rawFilteredSize = ncol(filteredObj)
    
    
    filteredObj <- NormalizeData(filteredObj, verbose = FALSE)
    filteredObj <- FindVariableFeatures(filteredObj, verbose = FALSE)
    
    objlist[[x]] = filteredObj

    print(x)
    print(filteredObj)
    
    statDF <- rbind(statDF, data.frame("RawAll"=origSize, "RawFiltered"=rawFilteredSize))
    
}


for (name in names(objlist))
{
    print(name)
  #p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  #save_plot(p, paste(name, "violins_qc", sep="_"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste(name, "scatter_ncount_mt", sep="_"), fig.width=10, fig.height=6)

  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste(name, "scatter_ncount_rp", sep="_"), fig.width=10, fig.height=6)
}



objlist.raw = objlist

objlist <- lapply(X = objlist.raw, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  print(paste("Seurat obj project", obj@project.name))
  print(obj)
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 100)
  obj <- subset(obj, subset = percent.mt < 15)
  print(obj)
  
  return(obj)
})


for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste(name, "prefiltered_violins_qc", sep="_"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = cowplot::plot_grid(plotlist=p, ncol=3)
  save_plot(p, paste(name, "prefiltered_violins_detail_qc", sep="_"), fig.width=18, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste(name, "prefiltered_scatter_ncount_mt", sep="_"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste(name, "prefiltered_scatter_ncount_rp", sep="_"), fig.width=10, fig.height=6)
}


objlist <- lapply(X = objlist.raw, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  print(paste("Seurat obj project", obj@project.name))
  print(obj)
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 100)
  obj <- subset(obj, subset = percent.mt < 15)
  print(obj)
  
  return(obj)
})



for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste(name, "filtered_violins_qc", sep="_"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = cowplot::plot_grid(plotlist=p, ncol=3)
  save_plot(p, paste(name, "filtered_violins_detail_qc", sep="_"), fig.width=18, fig.height=6)
  
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste(name, "filtered_scatter_ncount_mt", sep="_"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste(name, "filtered_scatter_ncount_rp", sep="_"), fig.width=10, fig.height=6)
}


relevantHTOs = list()
relevantHTOs[["samples_CATHR"]] = c("HAB10")


htoObjList = list()

for (name in names(allABs.raw))
{

    relHTOs = relevantHTOs[[name]]
    print("Relevant HTOs:")
    print(relHTOs)

  cellsGEX = length(colnames(objlist[[name]]))
  cellsAB = length(colnames(allABs.raw[[name]]))
  
  gexnames = substring(colnames(objlist[[name]]), str_length(name)+2)
  print(head(gexnames))
  cellsJoint = intersect(gexnames, colnames(allABs.raw[[name]]))
  cellsABsub = allABs.raw[[name]][, cellsJoint]
  cellsABsub = cellsABsub[relHTOs,]
  
  colnames(cellsABsub) = paste(objlist[[name]]@project.name, colnames(cellsABsub), sep="_")
  print(head(colnames(cellsABsub)))
  print(paste(cellsGEX, cellsAB, length(cellsJoint), ncol(cellsABsub)))
  
  # Normalize RNA data with log normalization
  xobj <- NormalizeData(objlist[[name]])
  # Find and scale variable features
  xobj <- FindVariableFeatures(xobj, selection.method = "mean.var.plot")
  xobj <- ScaleData(xobj, features = VariableFeatures(xobj))
    
  
  xobj[["HTO"]] <- CreateAssayObject(counts = cellsABsub)
  xobj <- NormalizeData(xobj, assay = "HTO", normalization.method = "CLR")
  xobj <- HTODemux(xobj, assay = "HTO", positive.quantile = 0.99)
  
  print(table(xobj$HTO_classification.global))
  print(table(xobj$HTO_classification))

  htoObjList[[name]] = xobj
    
}


for (name in names(htoObjList))
{
  p=VlnPlot(htoObjList[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "HTO_classification.global")
  save_plot(p, paste(name, "hto_violins_qc", sep="_"), fig.width=10, fig.height=6)

  allHTOFeatures = rownames(htoObjList[[name]][["HTO"]])
  figHeight = round(length(allHTOFeatures)*0.5 * 6)
  r=RidgePlot(htoObjList[[name]], assay = "HTO", features = allHTOFeatures, ncol = 2)
  save_plot(r, paste(name, "hto_ridge_qc", sep="_"), fig.width=10, fig.height=figHeight)

  allHTOFeatures = rownames(htoObjList[[name]][["HTO"]])
  figHeight = round(length(unique(htoObjList[[name]]$HTO_classification))*0.5 * 4)
  r=RidgePlot(htoObjList[[name]], assay = "HTO", features = allHTOFeatures, group.by="HTO_classification", ncol = 2)
  save_plot(r, paste(name, "hto_ridge_detail_qc", sep="_"), fig.width=10, fig.height=figHeight)

  
  print(name)
  print(table(htoObjList[[name]]$HTO_classification))

}


htoFiltered <- lapply(X = htoObjList, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  #obj <- subset(obj, subset = HTO_classification == "Blood" | HTO_classification == "Thrombus")
    print(paste("Seurat obj project", obj@project.name))
    print(obj)
    print(table(obj$HTO_classification))
    obj = subset(obj, subset=HTO_classification.global == "Singlet")
    print(table(obj$HTO_classification))
  return(obj)
})


#
#
#
# This is meant to add ABA to all samples
#
#
rownametransform = list(
  "CD10"="CD10",
  "CD15"="CD15",
  "CD16"="CD16",
  "CD3"="CD3",
  "CD4.1"="CD4",
  "SIGLEC10.1"="SIGLEC10",
  "CD41"="CD41",
  "CD56"="CD56",
  "CD14.1"="CD14",
  "CD19.1"="CD19"
  )

citeFiltered <- lapply(X = htoFiltered, FUN = function(x) {


    print(paste("Seurat obj project", x@project.name))
    print(x)

    origSample = x@project.name
    abAssayMat <- allABs.raw[[origSample]]
    colnames(abAssayMat) = paste(origSample, colnames(abAssayMat), sep="_")

    gexnames = colnames(x)
    cellsJoint = intersect(gexnames, colnames(abAssayMat))


    abAssayMat = abAssayMat[,cellsJoint]
    abAssayMat = abAssayMat[grep("HAB", rownames(abAssayMat), invert=T, value=T), ]

    abAssayMat = abAssayMat[names(rownametransform), ]
    rownames(abAssayMat) = as.character(rownametransform)

    abAssay=CreateAssayObject(counts = as.matrix(abAssayMat))


    x[["ABA"]] <- abAssay
    x <- NormalizeData(x, normalization.method = "CLR", margin = 2, assay = "ABA")
    return(x)
})




prepareFinalList = function(finalList)
{



print("cells per experiment")
print(mapply(sum, lapply(finalList, function(x) {dim(x)[2]})))
print("total cells")
print(sum(mapply(sum, lapply(finalList, function(x) {dim(x)[2]}))))


objlist = list()
for (objname in names(finalList))
{

    x = finalList[[objname]]
    print(objname)
    print(paste("Seurat obj project", x@project.name))
    
    Project(x) = objname
    print(paste("Seurat obj project", x@project.name))

    DefaultAssay(x) = "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)

    x$library = objname

    objlist[[objname]] = x
}


features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 4000)
objlist <- lapply(X = objlist, FUN = function(x) {

    print(paste("Seurat obj project", x@project.name))
    print(x)

    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes

    x <- CellCycleScoring(
      x,
      g2m.features = g2m.genes,
      s.features = s.genes)


    x <- ScaleData(x, features = features, verbose = FALSE, assay="RNA", vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score"))
    x <- RunPCA(x, features = features, verbose = FALSE, assay="RNA")
    #x <- suppressWarnings(SCTransform(x, verbose = FALSE,vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score")))


    x$project = x@project.name

    return(x)
})

return(objlist)

}



analyseFinalList = function(objlist,intname, do.cite=TRUE)
{



if (do.cite)
{
  #
  # integrate based on ADT/ABA/CITE/HASHTAGS
  #
  objSamples = objlist

  objSamples = lapply(objSamples, function(x) {
    DefaultAssay(x) <- 'ABA'
    VariableFeatures(x) <- rownames(x) # all HTOs
    x = ScaleData(x, assay="ABA")
    x <- RunPCA(x, features = rownames(x), verbose = FALSE, approx=FALSE, npcs=10, reduction.name="pca",  assay="ABA")
    print(x@reductions$pca)

    return(x)
  })

  print(objSamples)

  features_adt <- rownames(objSamples[[1]][['ABA']])
  objlist.anchors.adt <- FindIntegrationAnchors(object.list = objSamples, assay=rep("ABA", length(objSamples)), normalization.method = "LogNormalize",
                                                anchor.features = features_adt, dims = 1:10, reduction = 'cca', k.filter = NA)
  adt.list.integrated <- IntegrateData(new.assay.name = "integrated_adt", anchorset = objlist.anchors.adt, normalization.method = "LogNormalize", dims=1:9)

  print("ADT integration done")
}


#
# integrate based on RNA/GEX assay
#
objSamples = objlist
print(objSamples)

objSamples = lapply(objSamples, function(x) {
  DefaultAssay(x) <- 'RNA'
#  x <- RunPCA(x, features = features, verbose = FALSE, reduction.name="pca",  assay="RNA")
#  DefaultAssay(x) <- 'SCT'
#  print(x@reductions$pca)
  return(x)
})
print("GEX integration features")
print(objSamples)

features_gex <- SelectIntegrationFeatures(object.list = objSamples, nfeatures = 4000)#, assay=rep("RNA", length(objSamples)))
#objSamples <- PrepSCTIntegration(object.list = objSamples, anchor.features = features_gex)
objlist.anchors <- FindIntegrationAnchors(object.list = objSamples,  reduction = "cca", dims = 1:50, anchor.features = features_gex) #normalization.method = "SCT",
obj.list.integrated <- IntegrateData(new.assay.name = "integrated_gex", anchorset = objlist.anchors, dims = 1:50, verbose=T) #normalization.method = "SCT",
print("GEX integration done")

#
# integrated GEX viz
#
obj.list.integrated = ScaleData(obj.list.integrated, assay="integrated_gex")
obj.list.integrated <- RunPCA(obj.list.integrated, npcs = 50, reduction.name="igpca", assay="integrated_gex")
obj.list.integrated <- RunUMAP(obj.list.integrated, reduction = "igpca", dims = 1:50, reduction.key = "UMAPig_",)
p=DimPlot(obj.list.integrated, group.by="orig_project", reduction="umap")
save_plot(p, paste(intname, "wnn_ig_dimplot", sep="/"), 8, 6)

p=DimPlot(obj.list.integrated, group.by="orig_project", reduction="igpca")
save_plot(p, paste(intname, "wnn_pca_ig_dimplot", sep="/"), 8, 6)

obj.list.gex_adt = NULL

if (do.cite)
{
  #
  # integrated ADT viz
  #
  obj.list.integrated[["integrated_adt"]] = adt.list.integrated[["integrated_adt"]]
  obj.list.integrated = ScaleData(obj.list.integrated, assay="integrated_adt")
  obj.list.integrated <- RunPCA(obj.list.integrated, features = rownames(adt.list.integrated[['ABA']]), verbose = FALSE, approx=FALSE, npcs=10, reduction.name="iapca",  assay="integrated_adt")
  obj.list.integrated <- RunUMAP(obj.list.integrated, reduction = "iapca", dims = 1:10, reduction.key = "UMAPia_",)

  p=DimPlot(obj.list.integrated, group.by="orig_project", reduction="umap")
  save_plot(p, paste(intname, "wnn_ia_dimplot", sep="/"), 8, 6)


  #
  # multi modal neighbors
  #

  obj.list.gex_adt <- FindMultiModalNeighbors(
    obj.list.integrated, reduction.list = list("igpca", "iapca"), 
    dims.list = list(1:50, 1:10), modality.weight.name = "RNA.weight", prune.SNN=1/20
  )
  #
  # multi modal viz
  #
  obj.list.gex_adt <- RunUMAP(obj.list.gex_adt, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  obj.list.gex_adt <- FindClusters(obj.list.gex_adt, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

  p <- DimPlot(obj.list.gex_adt, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5)
  save_plot(p, paste(intname, "wnn_ig_ia_cluster_dimplot", sep="/"), 8, 6)

  p <- DimPlot(obj.list.gex_adt, group.by="orig_project", reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5)
  save_plot(p, paste(intname, "wnn_ig_ia_project_dimplot", sep="/"), 8, 6)

  p=DimPlot(obj.list.gex_adt, reduction="umap")
  save_plot(p, paste(intname, "wnn_ig_wnncluster_dimplot", sep="/"), 8, 6)

  obj.list.integrated$wnn_clusters = Idents(obj.list.gex_adt)
}




return(list("integrated"=obj.list.integrated, "multimodal"=obj.list.gex_adt))

}


finalList_library = splitObjListByLibrary(citeFiltered)
finalList_library = prepareFinalList(finalList_library)
integratedList_library = analyseFinalList(finalList_library, "libint/wnn")

finalList_patient = splitObjListByHTO(citeFiltered)
finalList_patient = prepareFinalList(finalList_patient)
integratedList_patient = analyseFinalList(finalList_patient, "patint/wnn")


finalList_patient_thrombus = splitObjListByHTO(citeFiltered)

obj.thr3.blood = readRDS("../extdata/sample_pat3.obj.integrated.blood.Rds")
obj.thr5.blood = readRDS("../extdata/sample_pat45.obj.integrated.blood.Rds")
obj.thr6.blood = readRDS("../extdata/sample_pat6.obj.integrated.blood.Rds")
#obj.thr12.blood = readRDS("../extdata/thrombus12.blood.Rds")
#obj.thr12 = readRDS("../extdata/thrombus12.Rds")
#obj.thr12.blood = subset(obj.thr12, HTO_classification =="Blood")
#cellnames.thr1 = grep("^21053_0001_AB", colnames(obj.thr12.blood), value=T)
#cellnames.thr2 = grep("^21053_0003_AB", colnames(obj.thr12.blood), value=T)
#obj.thr1.blood = subset(obj.thr12.blood, cells=cellnames.thr1)
#obj.thr2.blood = subset(obj.thr12.blood, cells=cellnames.thr2)
#saveRDS(obj.thr1.blood, "../extdata/thrombus1.blood.Rds")
#saveRDS(obj.thr2.blood, "../extdata/thrombus2.blood.Rds")

obj.thr1.blood = readRDS("../extdata/thrombus1.blood.Rds")
obj.thr2.blood = readRDS("../extdata/thrombus2.blood.Rds")

print(unique(obj.thr45.blood$patient))

obj.thr1.blood$orig_project = "thrombus1"
obj.thr2.blood$orig_project = "thrombus2"
obj.thr3.blood$orig_project = "thrombus3"
obj.thr5.blood$orig_project = "thrombus5"
obj.thr6.blood$orig_project = "thrombus6"

finalList_patient_thrombus[["thrombus1"]] = obj.thr1.blood
finalList_patient_thrombus[["thrombus2"]] = obj.thr2.blood
finalList_patient_thrombus[["thrombus3"]] = obj.thr3.blood
finalList_patient_thrombus[["thrombus5"]] = obj.thr5.blood
finalList_patient_thrombus[["thrombus6"]] = obj.thr6.blood

finalList_patient_thrombus = prepareFinalList(finalList_patient_thrombus)

integratedList_patient_thrombus = analyseFinalList(finalList_patient_thrombus, "patint_thr/wnn", do.cite=FALSE)


save.image("../all.integrated.cca.RData", compress=FALSE)

source("/mnt/input/own/SeqMappings/Pekayvaz/21068_M12_macs/functions.R")
load("../all.integrated.cca.RData")


preprocessObj = function(obj.in, useAssay, inname, reduction="umap")
{


DefaultAssay(obj.in) <- useAssay

if (reduction == "umap")
{
  # Run the standard workflow for visualization and clustering
  obj.in <- ScaleData(obj.in, verbose = FALSE)
  obj.in <- RunPCA(obj.in, npcs = 50, verbose = FALSE)
  obj.in <- RunUMAP(obj.in, reduction = "pca", dims = 1:50)
  obj.in <- FindNeighbors(obj.in, reduction = "pca", dims = 1:50)
  obj.in <- FindClusters(obj.in, resolution = 0.5)
}

obj.in$idents = Idents(obj.in)

p=DimPlot(obj.in, pt.size = 0.001, label=T, reduction = reduction)
save_plot(p, paste(inname, "dimplot_umap", sep="/"), fig.width=12, fig.height=8)

numProjects = length(unique(obj.in$orig_project))
numRows = ceiling(numProjects/2)


p=DimPlot(obj.in, pt.size = 0.001, label=T, split.by="orig_project", reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_project", sep="/"), fig.width=24, fig.height=8*numRows)

obj.in$libraryHTO = paste(obj.in$library, obj.in$HTO_classification, sep="_")

numLibHTO = length(unique(obj.in$libraryHTO))
numRows = ceiling(numLibHTO/2)

p=DimPlot(obj.in, split.by="libraryHTO", pt.size = 0.001, label=T, reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_libraryHTO", sep="/"), fig.width=12, fig.height=8*numRows)


cellList = colnames(obj.in)
featVec <- vector(mode="character", length=length(cellList))

if ("sepsis1_HAB1" %in% obj.in$library)
{
    featVec = obj.in$library
} else {
    featVec = obj.in$libraryHTO
}


featVec[featVec == "sepsis1_HAB1"] = "control"
featVec[featVec == "sepsis1_HAB2"] = "control"
featVec[featVec == "sepsis1_HAB3"] = "control"
featVec[featVec == "sepsis1_HAB4"] = "control"

featVec[featVec == "sepsis1_HAB5"] = "sepsis"
featVec[featVec == "sepsis1_HAB6"] = "sepsis"
featVec[featVec == "sepsis1_HAB7"] = "sepsis"
featVec[featVec == "sepsis1_HAB8"] = "sepsis"

featVec[featVec == "sepsis2_HAB1"] = "sepsis"
featVec[featVec == "sepsis2_HAB2"] = "sepsis"
featVec[featVec == "sepsis2_HAB3"] = "sepsis"
featVec[featVec == "sepsis2_HAB4"] = "sepsis"

featVec[featVec == "sepsis2_HAB5"] = "control"
featVec[featVec == "sepsis2_HAB6"] = "control"
featVec[featVec == "sepsis2_HAB7"] = "control"
featVec[featVec == "sepsis2_HAB8"] = "control"
featVec[featVec == "sepsis2_HAB9"] = "control"

featVec[featVec == "thrombus1"] = "thrombus"
featVec[featVec == "thrombus2"] = "thrombus"
featVec[featVec == "thrombus3"] = "thrombus"
featVec[featVec == "thrombus5"] = "thrombus"
featVec[featVec == "thrombus6"] = "thrombus"

obj.in$condition=featVec

numConditions = length(unique(obj.in$condition))
numRows = ceiling(numConditions/2)

p=DimPlot(obj.in, pt.size = 0.001, label=T, split.by="condition", reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_condition", sep="/"), fig.width=16, fig.height=8*numRows)


cellList = colnames(obj.in)
featVec <- vector(mode="character", length=length(cellList))

if ("sepsis1_HAB1" %in% obj.in$library)
{
    featVec = obj.in$library
} else {
    featVec = obj.in$libraryHTO
}


featVec[featVec == "sepsis1_HAB1"] = "K1"
featVec[featVec == "sepsis1_HAB2"] = "K2"
featVec[featVec == "sepsis1_HAB3"] = "K3"
featVec[featVec == "sepsis1_HAB4"] = "K4"

featVec[featVec == "sepsis1_HAB5"] = "S4"
featVec[featVec == "sepsis1_HAB6"] = "S5"
featVec[featVec == "sepsis1_HAB7"] = "S6"
featVec[featVec == "sepsis1_HAB8"] = "S7"

featVec[featVec == "sepsis2_HAB1"] = "S8"
featVec[featVec == "sepsis2_HAB2"] = "S9"
featVec[featVec == "sepsis2_HAB3"] = "S10"
featVec[featVec == "sepsis2_HAB4"] = "S11"

featVec[featVec == "sepsis2_HAB5"] = "K5"
featVec[featVec == "sepsis2_HAB6"] = "K6"
featVec[featVec == "sepsis2_HAB7"] = "K7"
featVec[featVec == "sepsis2_HAB8"] = "K8"
featVec[featVec == "sepsis2_HAB9"] = "K9"

featVec[featVec == "thrombus1"] = "thr1"
featVec[featVec == "thrombus2"] = "thr2"
featVec[featVec == "thrombus3"] = "thr3"
featVec[featVec == "thrombus5"] = "thr5"
featVec[featVec == "thrombus6"] = "thr6"

obj.in$patient=featVec

numPatients = length(unique(obj.in$patient))
numRows = ceiling(numPatients/2)


p=DimPlot(obj.in, pt.size = 0.001, label=T, split.by="patient", reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_patient", sep="/"), fig.width=16, fig.height=8*numRows)

return(obj.in)

}


obj.integrated.library = preprocessObj(integratedList_library$integrated, "integrated_gex", "libint")
obj.integrated.librarywnn = preprocessObj(integratedList_library$multimodal, "integrated_gex", "libint_wnn", 'wnn.umap')

obj.integrated.patient = preprocessObj(integratedList_patient$integrated, "integrated_gex", "patint")
obj.integrated.patientwnn = preprocessObj(integratedList_patient$multimodal, "integrated_gex", "patint_wnn", 'wnn.umap')

obj.integrated.patient_thr = preprocessObj(integratedList_patient_thrombus$integrated, "integrated_gex", "patint_thr")



save.image("../all.thr.integrated.cca.RData", compress=FALSE)

############
############
############
############ REANALYSIS of patient_thr!
############
############
############


load("../all.thr.integrated.cca.RData")
source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")

preprocessObj = function(obj.in, useAssay, inname, reduction="umap", cluster.dims=50, cluster.k=20, cluster.resolution=1.0)
{

  if (!dir.exists(inname)){
    dir.create(inname)
    print(inname)
    print("Dir created!")

  } else {
    print(inname)
    print("Dir already exists!")
  }

DefaultAssay(obj.in) <- useAssay

if (reduction == "umap")
{
  # Run the standard workflow for visualization and clustering
  obj.in <- ScaleData(obj.in, verbose = FALSE)
  obj.in <- RunPCA(obj.in, npcs = 50, verbose = FALSE)
  obj.in <- RunUMAP(obj.in, reduction = "pca", dims = 1:50)
  obj.in <- FindNeighbors(obj.in, reduction = "pca", dims = 1:cluster.dims, k.param=cluster.k)
  obj.in <- FindClusters(obj.in, resolution = cluster.resolution)

}

obj.in$idents = Idents(obj.in)

p=DimPlot(obj.in, pt.size = 0.001, label=T, reduction = reduction)
save_plot(p, paste(inname, "dimplot_umap", sep="/"), fig.width=12, fig.height=8)

numProjects = length(unique(obj.in$orig_project))
numRows = ceiling(numProjects/2)


p=DimPlot(obj.in, pt.size = 0.001, label=T, split.by="orig_project", reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_project", sep="/"), fig.width=24, fig.height=8*numRows)

obj.in$libraryHTO = paste(obj.in$library, obj.in$HTO_classification, sep="_")

numLibHTO = length(unique(obj.in$libraryHTO))
numRows = ceiling(numLibHTO/2)

p=DimPlot(obj.in, split.by="libraryHTO", pt.size = 0.001, label=T, reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_libraryHTO", sep="/"), fig.width=12, fig.height=8*numRows)


cellList = colnames(obj.in)
featVec <- vector(mode="character", length=length(cellList))

if ("sepsis1_HAB1" %in% obj.in$library)
{
    featVec = obj.in$library
} else {
    featVec = obj.in$libraryHTO
}


featVec[featVec == "sepsis1_HAB1"] = "control"
featVec[featVec == "sepsis1_HAB2"] = "control"
featVec[featVec == "sepsis1_HAB3"] = "control"
featVec[featVec == "sepsis1_HAB4"] = "control"

featVec[featVec == "sepsis1_HAB5"] = "sepsis"
featVec[featVec == "sepsis1_HAB6"] = "sepsis"
featVec[featVec == "sepsis1_HAB7"] = "sepsis"
featVec[featVec == "sepsis1_HAB8"] = "sepsis"

featVec[featVec == "sepsis2_HAB1"] = "sepsis"
featVec[featVec == "sepsis2_HAB2"] = "sepsis"
featVec[featVec == "sepsis2_HAB3"] = "sepsis"
featVec[featVec == "sepsis2_HAB4"] = "sepsis"

featVec[featVec == "sepsis2_HAB5"] = "control"
featVec[featVec == "sepsis2_HAB6"] = "control"
featVec[featVec == "sepsis2_HAB7"] = "control"
featVec[featVec == "sepsis2_HAB8"] = "control"
featVec[featVec == "sepsis2_HAB9"] = "control"

featVec[featVec == "thrombus1"] = "thrombus"
featVec[featVec == "thrombus2"] = "thrombus"
featVec[featVec == "thrombus3"] = "thrombus"
featVec[featVec == "thrombus5"] = "thrombus"
featVec[featVec == "thrombus6"] = "thrombus"

obj.in$condition=featVec

numConditions = length(unique(obj.in$condition))
numRows = ceiling(numConditions/2)

p=DimPlot(obj.in, pt.size = 0.001, label=T, split.by="condition", reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_condition", sep="/"), fig.width=16, fig.height=8*numRows)


cellList = colnames(obj.in)
featVec <- vector(mode="character", length=length(cellList))

if ("sepsis1_HAB1" %in% obj.in$library)
{
    featVec = obj.in$library
} else {
    featVec = obj.in$libraryHTO
}


featVec[featVec == "sepsis1_HAB1"] = "K1"
featVec[featVec == "sepsis1_HAB2"] = "K2"
featVec[featVec == "sepsis1_HAB3"] = "K3"
featVec[featVec == "sepsis1_HAB4"] = "K4"

featVec[featVec == "sepsis1_HAB5"] = "S4"
featVec[featVec == "sepsis1_HAB6"] = "S5"
featVec[featVec == "sepsis1_HAB7"] = "S6"
featVec[featVec == "sepsis1_HAB8"] = "S7"

featVec[featVec == "sepsis2_HAB1"] = "S8"
featVec[featVec == "sepsis2_HAB2"] = "S9"
featVec[featVec == "sepsis2_HAB3"] = "S10"
featVec[featVec == "sepsis2_HAB4"] = "S11"

featVec[featVec == "sepsis2_HAB5"] = "K5"
featVec[featVec == "sepsis2_HAB6"] = "K6"
featVec[featVec == "sepsis2_HAB7"] = "K7"
featVec[featVec == "sepsis2_HAB8"] = "K8"
featVec[featVec == "sepsis2_HAB9"] = "K9"

featVec[featVec == "thrombus1"] = "thr1"
featVec[featVec == "thrombus2"] = "thr2"
featVec[featVec == "thrombus3"] = "thr3"
featVec[featVec == "thrombus5"] = "thr5"
featVec[featVec == "thrombus6"] = "thr6"

obj.in$patient=featVec

numPatients = length(unique(obj.in$patient))
numRows = ceiling(numPatients/2)


p=DimPlot(obj.in, pt.size = 0.001, label=T, split.by="patient", reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_patient", sep="/"), fig.width=16, fig.height=8*numRows)


cellList = colnames(obj.in)
featVec <- vector(mode="character", length=length(cellList))

if ("sepsis1_HAB1" %in% obj.in$library)
{
    featVec = obj.in$library
} else {
    featVec = obj.in$libraryHTO
}


featVec[featVec == "sepsis1_HAB1"] = -1
featVec[featVec == "sepsis1_HAB2"] = -1
featVec[featVec == "sepsis1_HAB3"] = -1
featVec[featVec == "sepsis1_HAB4"] = -1

featVec[featVec == "sepsis1_HAB5"] = 4
featVec[featVec == "sepsis1_HAB6"] = 0
featVec[featVec == "sepsis1_HAB7"] = 2
featVec[featVec == "sepsis1_HAB8"] = 3

featVec[featVec == "sepsis2_HAB1"] = 1
featVec[featVec == "sepsis2_HAB2"] = 2
featVec[featVec == "sepsis2_HAB3"] = 5
featVec[featVec == "sepsis2_HAB4"] = 5

featVec[featVec == "sepsis2_HAB5"] = -1
featVec[featVec == "sepsis2_HAB6"] = -1
featVec[featVec == "sepsis2_HAB7"] = -1
featVec[featVec == "sepsis2_HAB8"] = -1
featVec[featVec == "sepsis2_HAB9"] = -1

featVec[featVec == "thrombus1"] = -1
featVec[featVec == "thrombus2"] = -1
featVec[featVec == "thrombus3"] = -1
featVec[featVec == "thrombus5"] = -1
featVec[featVec == "thrombus6"] = -1

obj.in$patient=featVec

numPatients = length(unique(obj.in$patient))
numRows = ceiling(numPatients/2)


p=DimPlot(obj.in, pt.size = 0.001, label=T, split.by="patient", reduction = reduction, ncol=2)
save_plot(p, paste(inname, "dimplot_umap_patient", sep="/"), fig.width=16, fig.height=8*numRows)






return(obj.in)

}

obj.test = FindNeighbors(obj.integrated.patient_thr, reduction = "igpca", k.param=50, dims = 1:30)
obj.test = FindClusters(obj.test, resolution = 0.9)

obj.integrated.patient_thr = preprocessObj(integratedList_patient_thrombus$integrated, "integrated_gex", "patint_thr", cluster.dims=30, cluster.k=40, cluster.resolution=1.2)


#obj.integrated.patient_thr = FindNeighbors(obj.integrated.patient_thr, reduction = "pca", k.param=40, dims = 1:30)
#obj.integrated.patient_thr = FindClusters(obj.integrated.patient_thr, resolution = 1.2)
#p=DimPlot(obj.integrated.patient_thr, label=T)
#save_plot(p, "dim_test", 12,12)

#
#
##
# some final QC plots
##
#
#

obj.integrated.patient_thr  = makeQCPlots(obj.integrated.patient_thr, "patint_thr")
saveRDS(obj.integrated.patient_thr, "../obj.integrated.patient_thr.final.Rds")

obj.integrated.patient_thr = readRDS("../obj.integrated.patient_thr.final.Rds")

#
#
##
# Now check CITE+GEX side by side
##
#
#

cdGene2RNA = list()
cdGene2RNA[["CD10"]] = c("MME")
cdGene2RNA[["CD14"]] = c("CD14")
cdGene2RNA[["CD15"]] = c("FUT4", "FUT7")
cdGene2RNA[["CD16"]] = c("FCGR3A", "FCGR3B")
cdGene2RNA[["CD3"]] = c("CD3D")
cdGene2RNA[["CD4"]] = c("CD4")
cdGene2RNA[["SIGLEC10"]] = c("SIGLEC10")
cdGene2RNA[["CD41"]] = c("ITGA2B")
cdGene2RNA[["CD56"]] = c("NCAM1")
cdGene2RNA[["CD19"]] = c("CD19")

makeCITEPlots = function(obj.in, outfolder, reduction="umap")
{
  DefaultAssay(obj.in) = "RNA"
  cdGenes = rownames(obj.in[["ABA"]])

  for (cdGene in cdGenes)
  {
    DefaultAssay(obj.in) = "ABA"
    p1 <- FeaturePlot(obj.in, cdGene, cols = c("lightgrey", "darkgreen"), order=T, reduction=reduction) + ggtitle(paste(cdGene, "protein"))

    p2 = VlnPlot(obj.in, cdGene, group.by="idents", pt.size=0.01)
    save_plot(p2, paste(paste(outfolder,"cite_vlnplot", sep="/"), cdGene, sep="_"), fig.width=12, fig.height=4)


    rnagenes = cdGene2RNA[[cdGene]]

    DefaultAssay(obj.in) = "RNA"

    fplots = list("protein"=p1)
    for (rnagene in rnagenes)
    {
      fplots[[rnagene]] =  FeaturePlot(obj.in, rnagene, order=T, reduction=reduction) + ggtitle(paste(rnagene, "RNA"))

      p2 = VlnPlot(obj.in, rnagene, group.by="idents", pt.size=0.01)
      save_plot(p2, paste(paste(outfolder,"gex_vlnplot", sep="/"), rnagene, sep="_"), fig.width=12, fig.height=4)

    }
    save_plot(combine_plot_grid_list(plotlist=fplots, ncol=length(fplots)), paste(paste(outfolder,"cite_fplot", sep="/"), cdGene, sep="_"), fig.width=length(fplots)*10, fig.height=10)

  }

}

makeCITEPlots(obj.integrated.patient_thr, "patint_thr")

#
##
### same plot, divided for ctrl and sepsis
##
#



makeLargeCITEPlots = function(obj.in, outfolder, reduction="umap")
{
  DefaultAssay(obj.in) = "RNA"
  cdGenes = rownames(obj.in[["ABA"]])

  allConditions = unique(obj.in$condition)
  print(allConditions)

  for (cdGene in cdGenes)
  {

    fplots = list()
    for (conditionState in allConditions)
    {

      #
      # CONTROL
      #
      ctrlSubset = subset(obj.in, condition==conditionState)

      DefaultAssay(ctrlSubset) = "ABA"
      p1 <- FeaturePlot(ctrlSubset, cdGene, cols = c("lightgrey", "darkgreen"), order=T, reduction=reduction) + ggtitle(paste(cdGene, conditionState, "protein"))
      fplots[[paste("ABA", cdGene, conditionState, sep="_")]] = p1


      p2 = VlnPlot(ctrlSubset, cdGene, group.by="idents", pt.size=0.01)
      save_plot(p2, paste(paste(outfolder,paste("cite_vlnplot_", conditionState, sep=""), sep="/"), cdGene, sep="_"), fig.width=12, fig.height=4)


      rnagenes = cdGene2RNA[[cdGene]]
      DefaultAssay(ctrlSubset) = "RNA"

      for (rnagene in rnagenes)
      {
        fplots[[paste("RNA", rnagene, conditionState, sep="_")]] =  FeaturePlot(ctrlSubset, rnagene, order=T, reduction=reduction) + ggtitle(paste(rnagene, conditionState, "RNA"))

        p2 = VlnPlot(ctrlSubset, rnagene, group.by="idents", pt.size=0.01)
        save_plot(p2, paste(paste(outfolder,paste("gex_vlnplot_", conditionState, sep=""), sep="/"), rnagene, sep="_"), fig.width=12, fig.height=4)

      }


    }


    print(length(fplots))
    neededColumns = length(fplots) / length(allConditions)
    #
    # FINAL PLOT
    #
    save_plot(combine_plot_grid_list(plotlist=fplots, ncol=neededColumns), paste(paste(outfolder,paste("cite_fplot_", tolower(paste(allConditions, collapse="_")), sep=""), sep="/"), cdGene, sep="_"), fig.width=neededColumns*10, fig.height=20)

  }
}

makeLargeCITEPlots(obj.integrated.patient_thr, "patint_thr")


#
#
##
# Here we continue with just one seurat object!
##
#
#
#

source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")

exprDFOuts = list()

makeDiffAnalysis = function(obj.in, outfolder, reduction="umap", group.by="idents")
{

  if (!dir.exists(outfolder)){
    dir.create(outfolder)
    print(outfolder)
    print("Dir created!")

  } else {
    print(outfolder)
    print("Dir already exists!")
  }

  countByOIAbs = getCellCountDF(obj.in, prefix="", select_by = group.by, group_by="orig.ident", relative=F, show.percent=F, outname=paste(outfolder, "ia_count_by_origident.tsv", sep="/"))
  countByAFAbs = getCellCountDF(obj.in, prefix="", select_by = group.by, group_by="patient", relative=F, show.percent=F, outname=paste(outfolder, "ia_count_by_patient.tsv", sep="/"))
  countByCFAbs = getCellCountDF(obj.in, prefix="", select_by = group.by, group_by="condition", relative=F, show.percent=F, outname=paste(outfolder, "ia_count_by_condition.tsv", sep="/"))
  countByCFAbs = getCellCountDF(obj.in, prefix="", select_by = group.by, group_by="orig_project", relative=F, show.percent=F, outname=paste(outfolder, "ia_count_by_orig_project.tsv", sep="/"))



  DefaultAssay(obj.in) = "RNA"

  deResTT = makeDEResults(obj.in, group.by=group.by, assay="RNA", test="t")
  exprdfTT = getDEXpressionDF(obj.in, deResTT, assay="RNA", group.by=group.by)
  write.table(exprdfTT, paste(outfolder, "expr_test_t.tsv", sep="/"), sep="\t", row.names=F, quote = F)
  write_xlsx(exprdfTT, paste(outfolder, "expr_test_t.xlsx", sep="/"))

  markers.use.tt=subset(exprdfTT ,avg_log2FC>0&p_val_adj<0.01&!startsWith(gene, "MT-")&!startsWith(gene, "RP"))
  finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:20)
  print(finalMarkers.use.tt)

  p_dp_genes_idents = DotPlot(obj.in, features = unique(finalMarkers.use.tt$gene), assay="RNA", group.by=group.by)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))
  save_plot(p_dp_genes_idents, paste(outfolder, "dotplot_idents_genes", sep="/"), 45, 8)

  cluster_mean_expr = get_population_expression_data(obj.in, group="idents", outname=paste(outfolder, "cluster_mean_expr_sepsis", sep="/"))

  return(list("exprdf" = exprdfTT, "meanexpr"=cluster_mean_expr))
}




exprDFOuts[["patint_thr"]] = makeDiffAnalysis(obj.integrated.patient_thr, "patint_thr")

cells.sepsis = cellIDForClusters(obj.integrated.patient_thr, "condition", c("sepsis"))
cells.thrombus = cellIDForClusters(obj.integrated.patient_thr, "condition", c("thrombus"))
cells.control = cellIDForClusters(obj.integrated.patient_thr, "condition", c("control"))

de.sepsis.control = compareCellsByCluster(obj.integrated.patient_thr, cells.sepsis, cells.control, "sepsis", "control", outfolder=paste("patint_thr", "de", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
de.sepsis.thrombus = compareCellsByCluster(obj.integrated.patient_thr, cells.sepsis, cells.thrombus, "sepsis", "thrombus", outfolder=paste("patint_thr", "de", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
de.thrombus.control = compareCellsByCluster(obj.integrated.patient_thr, cells.thrombus, cells.control, "thrombus", "control", outfolder=paste("patint_thr", "de", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))

makeVolcanos(de.sepsis.control, "DE Sepsis vs Control", paste("patint_thr", "de_volcano", "sepsis_ctrl_volcano", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))
makeVolcanos(de.sepsis.thrombus, "DE Sepsis vs Thrombus", paste("patint_thr", "de_volcano", "sepsis_thrombus_volcano", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))
makeVolcanos(de.thrombus.control, "DE Thrombus vs Control", paste("patint_thr", "de_volcano", "thrombus_ctrl_volcano", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))

#
##
###
#### Neutro analysis (patint)
###
##
#

DefaultAssay(obj.integrated.patient_thr) = "RNA"

cluster.granulocyte = c(0,1,2,3,9)
cells.granulocyte = cellIDForClusters(obj.integrated.patient_thr, "idents", cluster.granulocyte)

cells.sepsis = cellIDForClusters(obj.integrated.patient_thr, "condition", c("sepsis"))
cells.control = cellIDForClusters(obj.integrated.patient_thr, "condition", c("control"))
cells.thrombus = cellIDForClusters(obj.integrated.patient_thr, "condition", c("thrombus"))

dir.create("patint_thr_granulocytes")

p=HighlightedDimPlot(obj.integrated.patient_thr, highlightedClusters=cluster.granulocyte, reduction="umap")
save_plot(p, paste("patint_thr_granulocytes", "dimplot_granulocytes", sep="/"), 8, 6)

granulocyteMarkers.sepsis_ctrl = compareClusters(scdata=obj.integrated.patient_thr,
                                  cellsID1=intersect(cells.granulocyte, cells.sepsis),
                                  cellsID2=intersect(cells.granulocyte, cells.control),
                                  prefix= "granulocytes",
                                  suffix1="sepsis",
                                  suffix2="control",
                                  test="t", fcCutoff="0.25", assay="RNA", outfolder="patint_thr_granulocytes/de",
                                  heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))

granulocyteMarkers.thrombus_ctrl = compareClusters(scdata=obj.integrated.patient_thr,
                                  cellsID1=intersect(cells.granulocyte, cells.thrombus),
                                  cellsID2=intersect(cells.granulocyte, cells.control),
                                  prefix= "granulocytes",
                                  suffix1="thrombus",
                                  suffix2="control",
                                  test="t", fcCutoff="0.25", assay="RNA", outfolder="patint_thr_granulocytes/de",
                                  heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))

granulocyteMarkers.sepsis_thrombus = compareClusters(scdata=obj.integrated.patient_thr,
                                  cellsID1=intersect(cells.granulocyte, cells.sepsis),
                                  cellsID2=intersect(cells.granulocyte, cells.thrombus),
                                  prefix= "granulocytes",
                                  suffix1="sepsis",
                                  suffix2="thrombus",
                                  test="t", fcCutoff="0.25", assay="RNA", outfolder="patint_thr_granulocytes/de",
                                  heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))


exprDFOuts[["patint_thr_granulocytes"]] =makeDiffAnalysis(subset(obj.integrated.patient_thr, idents %in% cluster.granulocyte), "patint_thr_granulocytes", 'umap')


de.granclusters.sepsis.control = compareCellsByCluster(subset(obj.integrated.patient_thr, idents %in% cluster.granulocyte), cells.sepsis, cells.control, "sepsis", "control", outfolder=paste("patint_thr_granulocytes", "de", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
de.granclusters.thrombus.control = compareCellsByCluster(subset(obj.integrated.patient_thr, idents %in% cluster.granulocyte), cells.thrombus, cells.control, "thrombus", "control", outfolder=paste("patint_thr_granulocytes", "de", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
de.granclusters.sepsis.thrombus = compareCellsByCluster(subset(obj.integrated.patient_thr, idents %in% cluster.granulocyte), cells.sepsis, cells.thrombus, "sepsis", "thrombus", outfolder=paste("patint_thr_granulocytes", "de", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))


names(de.granclusters.sepsis.control) = paste("sepsis_control_cluster", names(de.granclusters.sepsis.control), sep="")
names(de.granclusters.thrombus.control) = paste("thrombus_control_cluster", names(de.granclusters.thrombus.control), sep="")
names(de.granclusters.sepsis.thrombus) = paste("sepsis_thrombus_cluster", names(de.granclusters.sepsis.thrombus), sep="")

de.granclusters = c(de.granclusters.sepsis.control, de.granclusters.thrombus.control, de.granclusters.sepsis.thrombus)
names(de.granclusters) = str_replace(names(de.granclusters), " ", "_")

de.granclusters[["sepsis_ctrl"]] = granulocyteMarkers.sepsis_ctrl
de.granclusters[["thrombus_ctrl"]] = granulocyteMarkers.thrombus_ctrl
de.granclusters[["sepsis_thrombus"]] = granulocyteMarkers.sepsis_thrombus

makeVolcanos(de.granclusters, "Differential Expression", paste("patint_thr_granulocytes", "de_volcano", "volcano_granclusters", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))


#
##
###
#### Neutro sub (patint)
###
##
#

dir.create("patint_thr_granulocytes_sub")

obj.granulocyte = subset(obj.integrated.patient_thr, cells=cells.granulocyte)

p=DimPlot(obj.granulocyte, pt.size = 0.001, label=T)
save_plot(p, "patint_thr_granulocytes_sub/dimplot_granulocyte", fig.width=12, fig.height=8)

#obj.granulocyte <- ScaleData(obj.granulocyte, verbose = FALSE)
DefaultAssay(obj.granulocyte) = "integrated_gex"
obj.granulocyte <- RunPCA(obj.granulocyte, npcs = 50, verbose = FALSE)
obj.granulocyte <- RunUMAP(obj.granulocyte, reduction = "pca", dims = 1:50)
obj.granulocyte <- FindNeighbors(obj.granulocyte, reduction = "pca", dims = 1:50)
obj.granulocyte <- FindClusters(obj.granulocyte, resolution = 0.5)

obj.granulocyte$idents_subcluster = Idents(obj.granulocyte)

p1=DimPlot(obj.granulocyte, pt.size = 0.001, label=T, group.by="idents_subcluster")
save_plot(p1, "patint_thr_granulocytes_sub/dimplot_subset_granulocyte", fig.width=12, fig.height=8)

p2=DimPlot(obj.granulocyte, pt.size = 0.001, label=T, group.by="idents_subcluster")
save_plot(p2, "patint_thr_granulocytes_sub/dimplot_subset_subcluster_granulocyte", fig.width=12, fig.height=8)

save_plot(p1|p2, "patint_thr_granulocytes_sub/dimplot_subset_subcluster_comparison_granulocyte", fig.width=24, fig.height=8)


p=DimPlot(obj.granulocyte, split.by="condition", label=T, group.by="idents_subcluster")
save_plot(p, "patint_thr_granulocytes_sub/dimplot_subset_granulocyte_bycondition", fig.width=14, fig.height=8)

p=DimPlot(obj.granulocyte, split.by="patient", label=T, ncol=2, group.by="idents_subcluster")
save_plot(p, "patint_thr_granulocytes_sub/dimplot_subset_granulocyte_bypatient", fig.width=14, fig.height=32)

#
## new clustering
#

exprDFOuts[["patint_thr_granulocytes_sub"]] = makeDiffAnalysis(obj.granulocyte, "patint_thr_granulocytes_sub", 'umap', group.by="idents_subcluster")

de.granulocyte.sepsis.control = compareCellsByCluster(obj.granulocyte, cells.sepsis, cells.control, "sepsis", "control", outfolder=paste("patint_thr_granulocytes_sub", "de", sep="/"), group.by="idents_subcluster", heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
de.granulocyte.sepsis.thrombus = compareCellsByCluster(obj.granulocyte, cells.sepsis, cells.thrombus, "sepsis", "thrombus", outfolder=paste("patint_thr_granulocytes_sub", "de", sep="/"), group.by="idents_subcluster", heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
de.granulocyte.thrombus.control = compareCellsByCluster(obj.granulocyte, cells.thrombus, cells.control, "thrombus", "control", outfolder=paste("patint_thr_granulocytes_sub", "de", sep="/"), group.by="idents_subcluster", heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))

names(de.granulocyte.sepsis.control) = paste("sepsis_control_cluster", names(de.granulocyte.sepsis.control), sep="")
names(de.granulocyte.sepsis.thrombus) = paste("sepsis_thrombus_cluster", names(de.granulocyte.sepsis.thrombus), sep="")
names(de.granulocyte.thrombus.control) = paste("thrombus_control_cluster", names(de.granulocyte.thrombus.control), sep="")

de.granulocyte = c(de.granulocyte.sepsis.control, de.granulocyte.sepsis.thrombus, de.granulocyte.thrombus.control)

makeVolcanos(de.granulocyte, "Differential Expression", paste("patint_thr_granulocytes_sub", "de_volcano", "volcano_subclustering", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))

#
##
### COLOR SCHEME
##
#

ucols = c("#000000", "#FF9200", "#FF2600")
#ucols = c("#F47B78","#fac0bf", "#448CCA")
splitColors = as.list(ucols)
names(splitColors) = c("control", "thrombus", "sepsis")

obj.integrated.patient_thr$condition_sorted = factor(obj.integrated.patient_thr$condition, levels=names(splitColors))










obj.goi = subset(obj.integrated.patient_thr, cells=cells.granulocyte)
obj.goi = ScaleData(obj.goi)

genesOfInterest = unique(c("FUT4", "FUT7","CMTM2", "FPR1", "MNDA", "CXCR1", "CXCR2", "CXCL8", "FCGR3A", "FCGR3B", "FCGR2A", "MME", "SELL", "ITGB2", "PECAM1", "C5AR1", "PTPRC", "CXCR4", "ITGAM", "CD177", "IL1R2", "IL7R", "ARG1", "ANXA1", "TOP2A", "CD24", "TLR2", "TLR4", "TLR8", "ITGAX", "CD101",  "S100A6", "S100A8", "S100A9", "S100A11", "S100A12", "SAT1", "SOD2", "MARCKS", "LITAF", "LUCAT1", "ITM2B", "CSF3R", "CSF2RB", "MNDA", "TREM1", "TNFRSF10C", "CEBPD", "CEBPB", "ALPL", "SELPLG", "FPR1", "MMP9", "MPO", "ELANE", "LYZ", "MXD1", "CEACAM8", "CD33", "SERPINA1", "SERPINB1", "CST7", "TREM1", "GCA", "STXBP2", "SORL1", "HIF1A", "ARRB2", "HMGB1", "CLEC2B", "IFITM1", "IFITM2", "MX2", "GBP1", "IFI6", "IFI16", "OASL", "PADI4", "NCF1", "CD83", "PLAC8", "PTGS2", "STAT5A"))
genesOfInterest.avail = genesOfInterest[genesOfInterest %in% rownames(obj.goi)]
genesOfInterest.unavail = genesOfInterest[!genesOfInterest %in% rownames(obj.goi)]
genesOfInterest.unavail

makeGOIAnalysis = function( obj.in, group.by="idents", outfolder="patint_thr_granulocytes_goi", goiName="goi_all", genes.interest=NULL, genes.fplot=T, assay.name="RNA", split.by="idents", cols=NULL)
{

  if (!dir.exists(outfolder)){
    dir.create(outfolder)
    print(outfolder)
    print("Dir created!")

  } else {
    print(outfolder)
    print("Dir already exists!")
  }

  scaleColors = c("#8b0000", "grey", "#008b2b")
  DefaultAssay(obj.in) = assay.name

  p=DoHeatmap(obj.in, genes.interest)+ scale_fill_gradientn(colors = scaleColors)
  save_plot(p, paste(outfolder, paste("hplot_", goiName, sep=""), sep="/"), fig.width=14, fig.height=0.3*length(genes.interest))

  p=DotPlot(obj.in, features=genes.interest, cols=c(scaleColors[1], scaleColors[3])) + coord_flip()# + scale_fill_gradientn(colors = scaleColors)
  save_plot(p, paste(outfolder, paste("dplot_", goiName, sep=""), sep="/"), fig.width=14, fig.height=0.3*length(genes.interest))


  plotElems = list()
  allConditions = unique(obj.in[[group.by]][[group.by]])

  allIdents = unique(obj.in[[split.by]][[split.by]])

  print("Unique conditions")
  print(allConditions)

  for (cond in allConditions)
  {
    cells.condition = cellIDForClusters(obj.in, group.by, c(cond))
    plotElems[[cond]] = list(cells=cells.condition, label=cond)
    p=DoHeatmap(subset(obj.in, cells=cells.condition), genes.interest)+ scale_fill_gradientn(colors = scaleColors)
    save_plot(p, paste(outfolder, paste("hplot_", goiName, "_", cond, sep=""), sep="/"), fig.width=14, fig.height=0.3*length(genes.interest))

  }
  
  
  p=enhancedDotPlot(obj.in, plotElems, featureGenes = genes.interest, group.by=split.by, col.min = -1, col.max = 1, title="", rotate.x=T, cols=scaleColors, abundance.perelem=T)
  save_plot(p, paste(outfolder, paste("sbs_dplot_", goiName, sep=""), sep="/"), fig.height=3*length(allIdents), fig.width=0.2*length(genes.interest)*length(allConditions))

  if (genes.fplot)
  {
    for (tgene in genes.interest)
    {
      p = FeaturePlot(obj.in, tgene, reduction="umap", order=T)
      save_plot(p, paste(outfolder, paste("fplot_",goiName, "_", tgene, sep=""), sep="/"), fig.width=14, fig.height=8)

      splitFeaturePlotName = paste(outfolder, paste("splitfplot_",goiName, "_", tgene, sep=""), sep="/")
      splitFeaturePlot(obj.in, tgene, group.by, paste("GEX", tgene), splitFeaturePlotName)

      p = VlnPlot(obj.in, tgene, pt.size=0.01, split.by=group.by, cols=cols)
      save_plot(p, paste(outfolder, paste("splitvplot_", goiName, "_", tgene, sep=""), sep="/"), fig.width=14, fig.height=4)

      p = VlnPlot(obj.in, tgene, pt.size=0.01)
      save_plot(p, paste(outfolder, paste("vplot_", goiName, "_", tgene, sep=""), sep="/"), fig.width=14, fig.height=4)
    }
  }

}

source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")
makeGOIAnalysis(obj.goi, outfolder="patint_thr_granulocytes_goi", goiName="goi_neutrosubset", assay.name="RNA", genes.interest=genesOfInterest.avail, genes.fplot=T, group.by="condition_sorted", cols=splitColors)

DefaultAssay(obj.integrated.patient_thr) = "RNA"
obj.integrated.patient_thr = ScaleData(obj.integrated.patient_thr)

makeGOIAnalysis(obj.integrated.patient_thr, outfolder="patint_thr_goi", goiName="goi_patientthr", assay.name="RNA", genes.interest=genesOfInterest.avail, genes.fplot=T, group.by="condition_sorted", cols=splitColors)

DefaultAssay(obj.integrated.patient_thr) = "RNA"
p=VlnPlot(obj.integrated.patient_thr, c("HBB", "CTSS", "EIF1"), pt.size=0.001)
save_plot(p, "highly_expressed_genes", fig.width=12, fig.height=3)

cluster_mean_expr = get_population_expression_data(obj.integrated.patient_thr, group="idents", outname=paste(outfolder, "cluster_mean_expr_sepsis", sep="/"), assay="RNA", slot="data")
cluster_mean_expr_granulocyte = get_population_expression_data(obj.granulocyte, group="idents_subcluster", outname=paste("patint_thr_granulocytes_subclustering", "cluster_mean_expr_sepsis", sep="/"), assay="RNA", slot="data")


save.image("../all.thr.integrated.rev.cca.RData", compress=FALSE)
load("../all.thr.integrated.rev.cca.RData")
source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")
#
##
### Manuscript Plots!
##
#
load("../all.thr.integrated.rev.cca.RData")
source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")

dir.create("manuscript")

obj.integrated.patient_thr$cellnames_quick = as.character(obj.integrated.patient_thr$idents)
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(4)] = "Dendritic cells"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(19)] = "Plasmablasts"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(7)] = "RBCs"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(15)] = "Platelets"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(6,10,11,14,16)] = "T cells"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(12)] = "NK cells"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(13)] = "B cells"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(5,8,17,18)] = "Monocytes"
obj.integrated.patient_thr$cellnames_quick[obj.integrated.patient_thr$cellnames_quick %in% c(0,1,2,3,9)] = "Granulocytes"

allCellTypes = c("Granulocytes","Monocytes","NK cells","Plasmablasts","Platelets","RBCs","T cells", "Dendritic cells", "B cells")
obj.integrated.patient_thr$cellnames_quick = fct_relevel(obj.integrated.patient_thr$cellnames_quick, allCellTypes)


p=DimPlot(obj.integrated.patient_thr, pt.size = 0.001, label=T, group.by="idents")
save_plot(p, "manuscript/dimplot", fig.width=12, fig.height=8)

p=DimPlot(obj.integrated.patient_thr, pt.size = 0.001, label=T, group.by="cellnames_quick")
save_plot(p, "manuscript/dimplot_celltype", fig.width=12, fig.height=8)




#
## Where does hash-based approach work better?
#

#nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 100)

p <- FeatureScatter(obj.integrated.patient_thr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="cellnames_quick")
p = p + geom_hline(yintercept=500, linetype='dashed', color=c('red')) # nCountRNA
p = p + geom_vline(xintercept=200, linetype='dashed', color=c('red')) # nFeature_RNA
p = p + xlim(0, 10000) + ylim(0, 2000)
save_plot(p, paste("manuscript", "scatter_ncount_rna_features", sep="/"), fig.width=10, fig.height=6)


p <- VlnPlot(obj.integrated.patient_thr, feature = "nCount_RNA", group.by="idents")
p = p + geom_hline(yintercept=500, linetype='dashed', color=c('red')) # nCountRNA
p = p + ylim(0, 6000)
save_plot(p, paste("manuscript", "vlnplot_ncount_rna", sep="/"), fig.width=10, fig.height=6)

p <- VlnPlot(obj.integrated.patient_thr, feature = "nCount_RNA", group.by="cellnames_quick")
p = p + geom_hline(yintercept=500, linetype='dashed', color=c('red')) # nCountRNA
p = p + ylim(0, 6000)
save_plot(p, paste("manuscript", "vlnplot_ncount_rna_celltypes", sep="/"), fig.width=10, fig.height=6)

#
## in Dotplot mit Cluster defining genes integrieren: CD177, ITGAX, HBB
#

obj.granulocyte

DefaultAssay(obj.granulocyte) = "RNA"


markers.use.tt.granulocytes_sub=subset( exprDFOuts[["patint_thr_granulocytes_sub"]]$exprdf,avg_log2FC>0&p_val_adj<0.01&!startsWith(gene, "MT-")&!startsWith(gene, "RP"))
finalMarkers.use.tt.granulocytes_sub = markers.use.tt.granulocytes_sub %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:20)
finalMarkers.use.tt.granulocytes_sub

p_dp_genes_idents = DotPlot(obj.granulocyte, features = unique(c(finalMarkers.use.tt.granulocytes_sub$gene, "CD177", "ITGAX", "HBB")), assay="RNA", group.by="idents")+ theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))
save_plot(p_dp_genes_idents, paste("manuscript", "granulocyte_sub_dotplot_idents_genes", sep="/"), 45, 4)



markers.use.tt.granulocyte=subset( exprDFOuts[["patint_thr_granulocytes"]]$exprdf,avg_log2FC>0&p_val_adj<0.01&!startsWith(gene, "MT-")&!startsWith(gene, "RP"))
finalMarkers.use.tt.granulocyte = markers.use.tt.granulocyte %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:20)
finalMarkers.use.tt.granulocyte

p_dp_genes_idents = DotPlot(obj.granulocyte, features = unique(c(finalMarkers.use.tt.granulocyte$gene, "CD177", "ITGAX", "HBB")), assay="RNA", group.by="idents")+ theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))
save_plot(p_dp_genes_idents, paste("manuscript", "granulocyte_dotplot_idents_genes", sep="/"), 45, 4)


#
## Heatmaps und Feature Plots bzw. Violin Plots zu den surface marker genes im FACS ergänzen: CD177, CD11b, MMP9, CD16, MME / CD10, CD184, CD182, CD88, IL1R2, CD114, LSEL, CD45, MX1
#

DefaultAssay(obj.integrated.patient_thr) = "RNA"
obj.goi = subset(obj.integrated.patient_thr, cells=cells.granulocyte)
obj.goi = ScaleData(obj.goi)

genesOfInterest = unique(c("CD177", "ITGAM", "MMP9", "FCGR3A", "FCGR3B", "MME", "CXCR4", "CXCR2", "C5AR1", "IL1R2", "CSF3R", "SELL", "PTPRC", "MX1"))
genesOfInterest.avail = genesOfInterest[genesOfInterest %in% rownames(obj.goi)]
genesOfInterest.unavail = genesOfInterest[!genesOfInterest %in% rownames(obj.goi)]
genesOfInterest.unavail


dir.create("manuscript/patint_thr_granulocytes_goi")
p=DimPlot(obj.goi, pt.size = 0.001, label=T, group.by="idents")
save_plot(p, "manuscript/patint_thr_granulocytes_goi/dimplot", fig.width=12, fig.height=8)




#-  Figure 4 + 5: heatmaps + feature plots ctrl vs. sepsis bzw. control vs. stroke relevanter Gene (z.B. CD177, IL1R2, PLAC8) etc.

additionalGOIs = unique(c("FCER1G", "IL1R2", "S100A12", "CST7", "IFITM3", "S100A12", "GBP5", "FCER1G", "TSPO", "C3AR1", "FCGR1A", "CD63", "FCGR1B", "MNDA", "SLPI", "IRAK3", "CYSTM1", "ANXA3", "FKBP5", "SERPINB1", "CD59", "IFIT1", "MX1", "IFIT2", "FCGR3A", "FCGR3B", "FUT4", "MX2", "PADI4", "MPO", "ELANE","CD177", "ITGAM", "MMP9", "FCGR3A", "FCGR3B", "MME", "CXCR4", "CXCR2", "C5AR1", "IL1R2", "CSF3R", "SELL", "PTPRC", "MX1"))
setdiff(additionalGOIs, rownames(obj.goi))

dir.create("manuscript/patint_thr_granulocytes_goi")
dir.create("manuscript/patint_thr_granulocytes_sub_goi")
dir.create("manuscript/patint_thr_goi")

obj.granulocyte$condition_sorted = factor(obj.granulocyte$condition, levels=names(splitColors))


DefaultAssay(obj.granulocyte) = "RNA"
obj.granulocyte = ScaleData(obj.granulocyte)

makeGOIAnalysis(obj.granulocyte, outfolder="manuscript/patint_thr_granulocytes_goi", goiName="goi_granulocytes_patintthr", assay.name="RNA", genes.interest=additionalGOIs, genes.fplot=T, group.by="condition_sorted", split.by="idents", cols=splitColors)
makeGOIAnalysis(obj.granulocyte, outfolder="manuscript/patint_thr_granulocytes_sub_goi", goiName="goi_granulocytes_patintthr", assay.name="RNA", genes.interest=additionalGOIs, genes.fplot=T, group.by="condition_sorted", split.by="idents_subcluster", cols=splitColors)

DefaultAssay(obj.integrated.patient_thr) = "RNA"
obj.integrated.patient_thr = ScaleData(obj.integrated.patient_thr)

dir.create("manuscript/patint_thr_goi")
makeGOIAnalysis(obj.integrated.patient_thr, outfolder="manuscript/patint_thr_goi", goiName="goi_all_patintthr", assay.name="RNA", genes.interest=additionalGOIs, genes.fplot=T, group.by="condition_sorted", cols=splitColors)


#
##
## -  volcano plots ersetzen mit heatmaps der interessanten DEGs (bei cluster-specific reaction to inflammation)
## -  CD177 feature plot (oder violin plot) der subcluster + conditions (für Figure 3)
## -  CD177 feature plot ctrl vs sepsis
##    was wo? welcher datensatz?
#

splitFeaturePlotName = paste("manuscript", paste("splitfplot_all_", "CD177", sep=""), sep="/")
splitFeaturePlot(obj.integrated.patient_thr, "CD177", "condition", paste("GEX", "CD177"), splitFeaturePlotName, limits=c(0,5), mirrorLimits=F, mid="lightgrey")

p = VlnPlot(obj.integrated.patient_thr, "CD177", group.by="idents", pt.size=0.01, split.by="condition_sorted", cols=splitColors)
save_plot(p, paste("manuscript", paste("splitvplot_all_", "CD177", sep=""), sep="/"), fig.width=14, fig.height=4)


save.image("../final_thr_sepsis_analysis.RData", compress=FALSE)



load("../final_thr_sepsis_analysis.RData")
source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")




exprDF_popdata = get_population_expression_data(obj.integrated.patient_thr, "idents", "mean_expression_idents", slot="data")
exprDF_popdata = get_population_expression_data(subset(obj.integrated.patient_thr, condition=="sepsis"), "idents", "mean_expression_sepsis_idents", slot="data")
exprDF_popdata = get_population_expression_data(subset(obj.integrated.patient_thr, condition=="thrombus"), "idents", "mean_expression_thrombus_idents", slot="data")


relgenes = c("CD177", "PECAM1", "CD38")

nicheResults = list()
nicheResults[["9"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(9), use.senders=c(0,1,2,3,5,8, 17, 18), use.condition="condition", use.condition.coi=c("sepsis"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes)
nicheResults[["3"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(3), use.senders=c(0,1,2,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("sepsis"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes)
nicheResults[["2"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(2), use.senders=c(0,1,3,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("sepsis"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes)
nicheResults[["1"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(1), use.senders=c(0,2,3,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("sepsis"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes)
nicheResults[["0"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(0), use.senders=c(1,2,3,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("sepsis"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes)

nicheResultsThr = list()
nicheResultsThr[["9"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(9), use.senders=c(0,1,2,3,5,8, 17, 18), use.condition="condition", use.condition.coi=c("thrombus"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes, outpath="./nichenetthr")
nicheResultsThr[["3"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(3), use.senders=c(0,1,2,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("thrombus"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes, outpath="./nichenetthr")
nicheResultsThr[["2"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(2), use.senders=c(0,1,3,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("thrombus"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes, outpath="./nichenetthr")
nicheResultsThr[["1"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(1), use.senders=c(0,2,3,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("thrombus"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes, outpath="./nichenetthr")
nicheResultsThr[["0"]] = runNicheNet(obj.integrated.patient_thr, use.receivers=c(0), use.senders=c(1,2,3,9,5,8, 17, 18), use.condition="condition", use.condition.coi=c("thrombus"), nichenet.folder="../../nichenet/",number.ligands=50, add.ligands=relgenes, add.receptors=relgenes, add.target=relgenes, outpath="./nichenetthr")


saveRDS(nicheResults, "../obj_int_patintthr_nichenet.Rds")
saveRDS(nicheResultsThr, "../obj_int_patintthr_nichenetthr.Rds")

source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")

relgenes = c("CD177", "PECAM1", "CD38")

x=runNicheNet(obj.integrated.patient_thr, use.receivers=c(1), use.senders=c(0,2,3,9,5,8,17,18), use.condition="condition", use.condition.coi=c("sepsis"), nichenet.folder="../../nichenet/",number.ligands=30, add.ligands=relgenes, add.receptors=relgenes, add.targets=relgenes)

"CD177" %in% colnames(x$active_ligand_target_links)
"CD177" %in% rownames(x$active_ligand_target_links)

#ligand_target_matrix = readRDS("../../nichenet/ligand_target_matrix.rds")

write_cell_barcodes(obj.integrated.patient_thr, "../cell_barcodes/")





library(ComplexHeatmap)

p=Heatmap(mat, row_order = rownames(mat), column_order = colnames(mat), cluster_rows=FALSE, cluster_columns=FALSE)
cowplot::save_plot("test.png", p)

plt <- ggplot(nicheResults[["2"]]$active_ligand_target_links_df,aes( target, ligand,fill=weight))
plt <- plt + geom_tile(colour = "white",  size=1.5) + scale_fill_gradient2(low = "white", high = "red") + scale_x_discrete(
    expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
    position = "top"
  ) + theme_minimal()
cowplot::save_plot("test.png", plt)


save_plot(p, "test", fig.height=12, fig.width=12)

# nohup ~/Rscript.sh ../../../scripts/velocities_step1.R "patint_thr" "../../obj.integrated.patient_thr.final.Rds" "../../loom_folder/"

#
##
### LR analysis by nichenet
##
#


genesetDF = as.data.frame(read.csv("../sepsis_genes.txt", sep="\t"))
genesets = list()

for (col in colnames(genesetDF))
{

  genes = genesetDF[, col]
  genes = genes[genes != ""]

  genesets[[col]] = unlist(lapply(genes, FUN=toupper))

  print(paste(col, length(genesets[[col]])))
  print(setdiff(genesets[[col]], rownames(obj.integrated.patient_thr)))
}


for (gsname in names(genesets))
{

  if (!dir.exists(dirname("genesets")))
  {
    dir.create(dirname("genesets"), recursive = TRUE)
  }

  print(gsname)
  print(genesets[[gsname]])

  print("missing")
  print(setdiff(genesets[[gsname]], rownames(obj.integrated.patient_thr)))

  DefaultAssay(obj.integrated.patient_thr) = "RNA"

  obj.integrated.patient_thr <- AddModuleScore(obj.integrated.patient_thr,
                    features = genesets[[gsname]],
                    name=gsname)

  baseplotname = paste("genesets/", gsname, sep="")

  ffname = paste(gsname, 1, sep="")
  split.by="condition"
  group.by="idents"

  splitFeaturePlot(obj.integrated.patient_thr, feature = ffname, split.by=split.by, title=paste("Module Score", gsname), filename=paste(baseplotname, "split", sep="_"), limits=c(-2, 2), low="#713fdf", high="#df713f",  mid="grey")

  p=plotPValueViolentBoxPlot(obj.integrated.patient_thr, ffname, group.by, split.by, NULL, splitColors)+ylab(paste("Module Score", gsname)) 
  save_plot(p, paste(baseplotname, "violent", sep="_"), fig.height=6, fig.width=20)

  p=plotPValueViolentBoxPlot(subset(obj.integrated.patient_thr, idents %in% c(0,1,2,3,9)), ffname, group.by, split.by, NULL, splitColors)+ylab(paste("Module Score", gsname)) 
  save_plot(p, paste(baseplotname, "neutro_violent", sep="_"), fig.height=6, fig.width=20)

  p=VlnPlot(subset(obj.integrated.patient_thr, idents %in% c(0,1,2,3,9)), ffname, group.by=group.by)+ylab(paste("Module Score", gsname)) 
  save_plot(p, paste(baseplotname, "neutro_violin", sep="_"), fig.height=4, fig.width=10)

  p=plotPValueViolentBoxPlot(subset(obj.integrated.patient_thr, idents %in% c(0,1,2,3,9)), ffname, group.by, split.by=NULL, dsrCols=NULL, onelineLabel=TRUE, yStepIncrease=0.75)+
  ylab(paste("Module Score", gsname))+ylim(c(0,10))
  save_plot(p, paste(baseplotname, "neutro_violent_pooled", sep="_"), fig.height=4, fig.width=10)

  p = FeaturePlot(subset(obj.integrated.patient_thr, idents %in% c(0,1,2,3,9)), features=ffname )+ggtitle(paste("Module Score", gsname))
  save_plot(p, paste(baseplotname, "neutro_featureplot", sep="_"), fig.height=10, fig.width=12)


}

source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")






#
##
### ENRICHMENT!
##
#


saveRDS(obj.granulocyte, "../obj.granulocyte.sub.Rds")
saveRDS(obj.integrated.patient_thr, "../obj.patient_thr.integrated.Rds")

saveRDS(rownames(obj.integrated.patient_thr), "../obj_integrated_patient_thr_universe.Rds")
saveRDS(rownames(obj.granulocyte), "../obj_granulocyte_universe.Rds")

saveRDS(de.granulocyte, "../de.granulocyte_results.Rds")


saveRDS(de.sepsis.control, "../de.patint_thr.sepsis.control.Rds")
saveRDS(de.sepsis.thrombus, "../de.patint_thr.sepsis.thrombus.Rds")
saveRDS(de.thrombus.control, "../de.patint_thr.thrombus.control.Rds")
saveRDS(rownames(obj.integrated.patient_thr), "../de.patint_thr.universe.Rds")

#<infile> <universeFile> <organismName> <outfile> <enrichmentresult>

nohup ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../de.patint_thr.sepsis.control.Rds ../de.patint_thr.universe.Rds human ../de.patint_thr.sepsis.control.gse.Rds ./patint_thr/gse/sepsis_control/ > enrich_sepsis_control.out &
nohup ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../de.patint_thr.sepsis.thrombus.Rds ../de.patint_thr.universe.Rds human ../de.patint_thr.sepsis.thrombus.gse.Rds ./patint_thr/gse/sepsis_thrombus/ > enrich_sepsis_thrombus.out &
nohup ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../de.patint_thr.thrombus.control.Rds ../de.patint_thr.universe.Rds human ../de.patint_thr.thrombus.control.gse.Rds ./patint_thr/gse/thrombus_control/ > enrich_thrombus_control.out &
nohup ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../de.granulocyte_results.Rds ../obj_granulocyte_universe.Rds human ../obj_granulocyte_universe.gse.Rds ./patint_thr_granulocytes/gse/ > enrich_granulocyte_results.out &

enrichment.granulocyte = readRDS("../enrichment_de.granulocyte_results.Rds")
enrichment.granulocyte_subclustering = readRDS("../enrichment_de.granulocyte_subclustering_results.Rds")
enrichment.control.sepsis = readRDS("../enrichment_de.patint_thr.control.sepsis.Rds")
enrichment.thrombus.sepsis = readRDS("../enrichment_de.patint_thr.thrombus.sepsis.Rds")
enrichment.control.thrombus = readRDS("../enrichment_de.patint_thr.control.thrombus.Rds")

#saveRDS(de.granulocyte_subclustering, "../de.granulocyte_subclustering_results.Rds")


cells.sepsis = cellIDForClusters(obj.integrated.patient_thr, "condition", c("sepsis"))
cells.thrombus = cellIDForClusters(obj.integrated.patient_thr, "condition", c("thrombus"))
cells.control = cellIDForClusters(obj.integrated.patient_thr, "condition", c("control"))



dect.sepsis.control = compareCellsByCluster(obj.integrated.patient_thr, cells.sepsis, cells.control, group.by="cellnames_quick", "sepsis", "control", outfolder=paste("patint_thr", "dect", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
dect.sepsis.thrombus = compareCellsByCluster(obj.integrated.patient_thr, cells.sepsis, cells.thrombus, group.by="cellnames_quick", "sepsis", "thrombus", outfolder=paste("patint_thr", "dect", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))
dect.thrombus.control = compareCellsByCluster(obj.integrated.patient_thr, cells.thrombus, cells.control, group.by="cellnames_quick", "thrombus", "control", outfolder=paste("patint_thr", "dect", sep="/"), heatmap.plot=TRUE, heatmap.addgenes=c("CD177"))

makeVolcanos(dect.sepsis.control, "DE Sepsis vs Control", paste("patint_thr", "dect_volcano", "sepsis_ctrl_volcano", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))
makeVolcanos(dect.sepsis.thrombus, "DE Sepsis vs Thrombus", paste("patint_thr", "dect_volcano", "sepsis_thrombus_volcano", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))
makeVolcanos(dect.thrombus.control, "DE Thrombus vs Control", paste("patint_thr", "dect_volcano", "thrombus_ctrl_volcano", sep="/"), turnExpression=F, FCcutoff=0.25, pCutoff = 0.05, highlightGene=c("CD177"))

saveRDS(dect.sepsis.control, "../dect.patint_thr.sepsis.control.Rds")
saveRDS(dect.sepsis.thrombus, "../dect.patint_thr.sepsis.thrombus.Rds")
saveRDS(dect.thrombus.control, "../dect.patint_thr.thrombus.control.Rds")

nohup ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../dect.patint_thr.sepsis.control.Rds ../de.patint_thr.universe.Rds human ../dect.patint_thr.sepsis.control.gse.Rds ./patint_thr/gse_ct/sepsis_control/ > enrich_dect_sepsis_control.out &
nohup ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../dect.patint_thr.sepsis.thrombus.Rds ../de.patint_thr.universe.Rds human ../dect.patint_thr.sepsis.thrombus.gse.Rds ./patint_thr/gse_ct/sepsis_thrombus/ > enrich_dect_sepsis_thrombus.out &
nohup ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../dect.patint_thr.thrombus.control.Rds ../de.patint_thr.universe.Rds human ../dect.patint_thr.thrombus.control.gse.Rds ./patint_thr/gse_ct/thrombus_control/ > enrich_dect_thrombus_control.out &


cp -r ./patint_thr/gse/ /mnt/raidtmp/joppich/lrzsyncshare/pekayvaz_sepsis/seurat_final/patint_thr/gse/
cp -r ./patint_thr/gse_ct/ /mnt/raidtmp/joppich/lrzsyncshare/pekayvaz_sepsis/seurat_final/patint_thr/gse_ct/

cp -r ./patint_thr_granulocytes/gse/ /mnt/raidtmp/joppich/lrzsyncshare/pekayvaz_sepsis/seurat_final/patint_thr_granulocytes/gse/








load("../final_thr_sepsis_analysis.RData")
source("/mnt/extproj/projekte/pekayvaz/scripts/functions.R")


p=DimPlot(subset(obj.integrated.patient_thr, idents %in% cluster.granulocyte), split.by="condition")
save_plot(p, "add_plots/split_neutros_by_condition", 24, 8)


#
##
###
#### SOFA scores
###
##
#

p=DimPlot(obj.integrated.patient_thr, pt.size = 0.001, label=T, group.by="idents")
save_plot(p, "add_plots/dimplot", fig.width=12, fig.height=8)


annotateList.sofa = list(
  list(name=-1, selector="^sepsis1_HAB1"),
  list(name=-1, selector="^sepsis1_HAB2"),
  list(name=-1, selector="^sepsis1_HAB3"),
  list(name=-1, selector="^sepsis1_HAB4"),
  
  list(name=4, selector="^sepsis1_HAB5"),
  list(name=0, selector="^sepsis1_HAB6"),
  list(name=2, selector="^sepsis1_HAB7"),
  list(name=3, selector="^sepsis1_HAB8"),

  list(name=1, selector="^sepsis2_HAB1"),
  list(name=2, selector="^sepsis2_HAB2"),
  list(name=5, selector="^sepsis2_HAB3"),
  list(name=5, selector="^sepsis2_HAB4"),
  
  list(name=-1, selector="^sepsis2_HAB5"),
  list(name=-1, selector="^sepsis2_HAB6"),
  list(name=-1, selector="^sepsis2_HAB7"),
  list(name=-1, selector="^sepsis2_HAB8"),
  list(name=-1, selector="^sepsis2_HAB9"),

  list(name=-1, selector="^thrombus1"),
  list(name=-1, selector="^thrombus2"),
  list(name=-1, selector="^thrombus3"),
  list(name=-1, selector="^thrombus5"),
  list(name=-1, selector="^thrombus6")

)

obj.integrated.patient_thr = annotateByCellnamePattern( obj.integrated.patient_thr, "sofa_score", annotateList.sofa, use.base = "library")



obj.sepsis = subset(obj.integrated.patient_thr, condition %in% c("sepsis", "control")& idents %in% c(0,1,2,3,9))
obj.sepsis@meta.data[ obj.sepsis$sofa_score == -1, "sofa_score"] = 0

exprdf_sepsis_c9 = get_population_expression_data(subset(obj.sepsis, idents %in% c(0,1,2,3)), "patient", NULL, addScores=c("sofa_score"))

for (gene in c("CD177", "CST7", "CXCR1", "CXCR2", "CSF3R", "S100A8", "S100A9", "S100A12", "ANXA3", "CXCR4", "SELL"))
{

  p = scatterAnalysisPlot(exprdf_sepsis_c9, "sofa_score", gene, logg1=FALSE)
  save_plot(p, paste("add_plots/sofa_sepsis_corr_c9", "sofa_score", gene, sep="_"), 8, 8)


}


exprdf_all = get_population_expression_data(obj.integrated.patient_thr, "idents", NULL)
isgGenes = c("ISG15", "RSAD2", "IFIT3", "IFIT1", "IFIT2", "MX1", "IFITM3", "IFI6", "XAF1", "IFI44L", "IFITM1", "MX2", "IFI16", "MT2A", "LY6E", "IFI44", "ISG20", "IFIT5", "IRF1", "IFIH1", "IFITM2", "IFI35", "IRF9", "IRF7", "IRF2", "IFI27L2", "IRF3")
exprdf_isg = exprdf_all[exprdf_all$gene %in% isgGenes,c("mean.0", "mean.1", "mean.2", "mean.3", "mean.9")]

indf = reshape2::melt(log1p(exprdf_isg))
indf$variable = as.character(indf$variable)
colnames(indf) = c("cluster", "expression")

comp_out = pairwise_t_test(as.tibble(indf), expression ~ cluster, paired = FALSE, p.adjust.method = "BH")
comp_out = comp_out[comp_out$p.adj.signif != "ns", ]
comp_out = rstatix::add_xy_position(comp_out, x = "cluster", step.increase=0.1)


p = ggplot(indf, aes(x=cluster, y=expression))+
  geom_violin(aes(fill=cluster)) + geom_boxplot(aes(fill=cluster)) +geom_jitter(aes(fill=cluster))+theme_minimal()
p = p + stat_pvalue_manual(comp_out, hide.ns=TRUE, tip.length = 0)
save_plot(p, "add_plots/isg_expression", fig.width=8, fig.height=6)




flowsetGenes = c('HBB', 'S100A12', 'S100A8', 'S100A9', 'S100A11', 'HBA2', 'FCER1G', 'SAT1', 'FTH1', 'CSF3R', 'PTPRC', 'NAMPT', 'TPT1', 'H3F3B', 'B2M', 'IFITM2', 'SRGN', 'SLC25A37', 'H3F3A', 'S100A6', 'FTL', 'HLA-A', 'HLA-C', 'GCA', 'TMSB10', 'SOD2', 'TXNIP', 'CYBA', 'ITM2B', 'HBA1', 'TYROBP', 'BASP1', 'ATP5F1E', 'GLUL', 'TMSB4X', 'HLA-B', 'EEF1A1', 'SERPINA1', 'NCF1', 'RNF149', 'VIM', 'LCP1', 'MYL6', 'BCL2A1', 'IFITM3', 'G0S2', 'FAU', 'EIF1', 'HLA-E', 'S100A4')
ageingGenes = c("CD177", "CXCR4", "SELL", "MMP9", "MME", "ITGAX", "CD24", "ITGA2", "ICAM1")


obj.integrated.patient_thr <- AddModuleScore(obj.integrated.patient_thr,
                  features = flowsetGenes,
                  name="flowsetGenes")
obj.integrated.patient_thr <- AddModuleScore(obj.integrated.patient_thr,
                  features = ageingGenes,
                  name="ageingGenes")

obj.sepsis_ctrl = subset(obj.integrated.patient_thr, condition %in% c("sepsis", "control"))

exprdf_sepsis_neutros = get_population_expression_data(subset(obj.sepsis_ctrl, idents %in% c(0,1,2,3,9)), "patient", NULL, addScores=c("flowsetGenes", "ageingGenes", "sofa_score"), slot="data")
exprdf_sepsis_neutros9 = get_population_expression_data(subset(obj.sepsis_ctrl, idents %in% c(9)), "patient", NULL, addScores=c("flowsetGenes", "ageingGenes", "sofa_score"), slot="data")

source("/mnt/extproj/projekte/pekayvaz/scripts/functions.R")


p = scatterAnalysisPlot(exprdf_sepsis_neutros, "flowsetGenes", "ageingGenes", logg1=FALSE, logg2=FALSE)
save_plot(p, paste("add_plots/corr_flowset_ageing", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "flowsetGenes", "ageingGenes", logg1=FALSE, logg2=FALSE)
save_plot(p, paste("add_plots/corr_flowset_ageing", "neutros9", sep="_"), 8, 8)

p = scatterAnalysisPlot(exprdf_sepsis_neutros, "sofa_score", "ageingGenes", logg1=FALSE, logg2=FALSE)
save_plot(p, paste("add_plots/corr_sofa_score_ageing", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "sofa_score", "ageingGenes", logg1=FALSE, logg2=FALSE)
save_plot(p, paste("add_plots/corr_sofa_score_ageing", "neutros9", sep="_"), 8, 8)

p = scatterAnalysisPlot(exprdf_sepsis_neutros, "sofa_score", "flowsetGenes", logg1=FALSE, logg2=FALSE)
save_plot(p, paste("add_plots/corr_sofa_score_flowsetGenes", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "sofa_score", "flowsetGenes", logg1=FALSE, logg2=FALSE)
save_plot(p, paste("add_plots/corr_sofa_score_flowsetGenes", "neutros9", sep="_"), 8, 8)


p = scatterAnalysisPlot(exprdf_sepsis_neutros, "CXCR4", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(0.75, 1.5))
save_plot(p, paste("add_plots/corr_CXCR4_flowsetGenes", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "CXCR4", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(0.75, 1.5))
save_plot(p, paste("add_plots/corr_CXCR4_flowsetGenes", "neutros9", sep="_"), 8, 8)


p = scatterAnalysisPlot(exprdf_sepsis_neutros, "MME", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(1.0, 1.5))
save_plot(p, paste("add_plots/corr_MME_flowsetGenes", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "MME", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(1.0, 1.5))
save_plot(p, paste("add_plots/corr_MME_flowsetGenes", "neutros9", sep="_"), 8, 8)

p = scatterAnalysisPlot(exprdf_sepsis_neutros, "SELL", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(1.2, 1.5))
save_plot(p, paste("add_plots/corr_SELL_flowsetGenes", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "SELL", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(1.2, 1.5))
save_plot(p, paste("add_plots/corr_SELL_flowsetGenes", "neutros9", sep="_"), 8, 8)


p = scatterAnalysisPlot(exprdf_sepsis_neutros, "CD177", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(.0, 1.5))
save_plot(p, paste("add_plots/corr_CD177_flowsetGenes", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "CD177", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(.0, 1.5))
save_plot(p, paste("add_plots/corr_CD177_flowsetGenes", "neutros9", sep="_"), 8, 8)

p = scatterAnalysisPlot(exprdf_sepsis_neutros, "CXCR2", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(0.75, 1.5))
save_plot(p, paste("add_plots/corr_CXCR2_flowsetGenes", "neutros", sep="_"), 8, 8)
p = scatterAnalysisPlot(exprdf_sepsis_neutros9, "CXCR2", "flowsetGenes", logg1=TRUE, logg2=FALSE, xlim.plot=c(0.75, 1.5))
save_plot(p, paste("add_plots/corr_CXCR2_flowsetGenes", "neutros9", sep="_"), 8, 8)

p=VlnPlot(subset(obj.sepsis_ctrl, idents %in% c(0,1,2,3,9)), c("CXCR4", "MME", "SELL", "CD177", "CXCR2"))
save_plot(p, paste("add_plots/corr_vplot", "neutros", sep="_"), 8, 4)


avge = AverageExpression(object = subset(obj.sepsis_ctrl, idents %in% c(0,1,2,3,9)), group.by="patient")
avge$RNA[c("CXCR4", "MME", "SELL", "CD177", "CXCR2"),]

exprdf_sepsis_neutros[exprdf_sepsis_neutros$gene=="flowsetGenes",]
exprdf_sepsis_neutros[exprdf_sepsis_neutros$gene=="ageingGenes",]

exprdf_sepsis_neutros[exprdf_sepsis_neutros$gene=="CXCR6",]
exprdf_sepsis_neutros[exprdf_sepsis_neutros$gene=="TMEM240",]



#
##
###
#### Protein Correlation
###
##
#

source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")


library(openxlsx)
protDF = read.xlsx("../proteome_files/neutrophils_proteome.xlsx")


protDF = protDF[protDF$Genes %in% rownames(obj.integrated.patient_thr),]
print(dim(protDF))

deFile1 = data.frame(gene=protDF[, c("Gene.names")], avg_log2FC=protDF[, c("logFC_Sepsis.Block.1.over.ctr")], p_val_adj=protDF[,c("adj.P.Val_Sepsis.Block.1.over.ctr")])
deFile2 = data.frame(gene=protDF[, c("Gene.names")], avg_log2FC=protDF[, c("logFC_Sepsis.Block.2.over.ctr")], p_val_adj=protDF[,c("adj.P.Val_Sepsis.Block.2.over.ctr")])


deFile1 = deFile1 %>%
       mutate(gene = strsplit(as.character(gene), ";")) %>%
       unnest(gene) %>%
       filter(gene != "") %>%
       select(gene, avg_log2FC:p_val_adj)
deFile2 = deFile2 %>%
       mutate(gene = strsplit(as.character(gene), ";")) %>%
       unnest(gene) %>%
       filter(gene != "") %>%
       select(gene, avg_log2FC:p_val_adj)

de.prot = list(block1_ctrl=as.data.frame(deFile1), block2_ctrl=as.data.frame(deFile2))


makeVolcanos(de.prot, "DE Sepsis vs Control", paste("prot_enrichment", "volcano", sep="/"), turnExpression=F, FCcutoff=1, pCutoff = 0.05, highlightGene=c("CD177"), restrict_labels=T)


de.prot.filtered = list()
de.prot.filtered$block1_ctrl = de.prot$block1_ctrl[abs(de.prot$block1_ctrl$avg_log2FC) >= 1,]
de.prot.filtered$block2_ctrl = de.prot$block2_ctrl[abs(de.prot$block2_ctrl$avg_log2FC) >= 1,]

saveRDS(de.prot.filtered, "../de.prot.Rds")


# ~/Rscript422.sh /mnt/input/own/projekte/pekayvaz/scripts/enrichmentAnalysis.R ../de.prot.Rds ../obj_granulocyte_universe.Rds human ../enrichment_de.prot.Rds prot_enrichment/enrichment/

#enrichment.prot = readRDS("../enrichment_de.prot.Rds")
#makeEnrichmentPlots(enrichment.prot, "prot_enrichment/")







cells.sepsis = cellIDForClusters(obj.integrated.patient_thr, "condition", c("sepsis"))
cells.thrombus = cellIDForClusters(obj.integrated.patient_thr, "condition", c("thrombus"))
cells.control = cellIDForClusters(obj.integrated.patient_thr, "condition", c("control"))

cells.sepsis = cellIDForClusters(obj.integrated.patient_thr, "condition", c("sepsis"))
cells.thrombus = cellIDForClusters(obj.integrated.patient_thr, "condition", c("thrombus"))
cells.control = cellIDForClusters(obj.integrated.patient_thr, "condition", c("control"))


cells.neutrophils = cellIDForClusters(obj.integrated.patient_thr, "idents", c(0,1,2,3,9))

prot_sepsis_ctrl = protDF[protDF[["adj.P.Val_Sepsis.Block.1.over.ctr"]] < 0.05,]
prot_logfc_sepsis1_ctrl = prot_sepsis_ctrl$logFC_Sepsis.Block.1.over.ctr
names(prot_logfc_sepsis1_ctrl) = prot_sepsis_ctrl$Genes

prot_sepsis_ctrl = protDF[protDF[["adj.P.Val_Sepsis.Block.2.over.ctr"]] < 0.05,]
prot_logfc_sepsis2_ctrl = prot_sepsis_ctrl$logFC_Sepsis.Block.2.over.ctr
names(prot_logfc_sepsis2_ctrl) = prot_sepsis_ctrl$Genes

block2unique = setdiff(names(prot_logfc_sepsis2_ctrl), names(prot_logfc_sepsis1_ctrl))
prot_sepsis_combined_ctrl = c(prot_logfc_sepsis1_ctrl, prot_logfc_sepsis2_ctrl[block2unique])



rna_sepsis_ctrl = compareClusters(scdata=obj.integrated.patient_thr,
                            cellsID1=intersect(cells.neutrophils, cells.sepsis),
                            cellsID2=intersect(cells.neutrophils, cells.control),
                            prefix= "neutrophils",
                            suffix1="sepsis",
                            suffix2="control",heatmap.plot=TRUE,
                            test="wilcox", fcCutoff=0.1, assay="RNA", outfolder="prot_comparison/rna_de")

makeVolcanos(list(all=rna_sepsis_ctrl), "DE Control vs Sepsis", paste("prot_comparison", "rna_de", "ctrl_sepsis_volcano", sep="/"), turnExpression=F, FCcutoff=0.2, pCutoff = 0.05)


rna_sepsis_ctrl = rna_sepsis_ctrl[rna_sepsis_ctrl$p_val_adj < 0.05,]

rna_logfc_sepsis_ctrl = rna_sepsis_ctrl$avg_log2FC
names(rna_logfc_sepsis_ctrl) = rna_sepsis_ctrl$gene





scatterNamedVectors = function(exprVec1RNA, exprVec2Prot, formula=y ~ x, name1=NULL,name2=NULL)
{

  if (is.null(name1))
  {
    name1 = "RNA"
  }

  if (is.null(name2))
  {
    name2 = "Protein"
  }

  common.names = intersect(names(exprVec1RNA), names(exprVec2Prot))
  exprVec1 = exprVec1RNA[common.names]
  exprVec2 = exprVec2Prot[common.names]

  lowerlimX = plyr::round_any(min(exprVec1), 0.1, f = floor)-1
  upperlimX = plyr::round_any(max(exprVec1), 0.1, f = ceiling)+1

  lowerlimY = plyr::round_any(min(exprVec2), 0.1, f = floor)-1
  upperlimY = plyr::round_any(max(exprVec2), 0.1, f = ceiling)+1


  dflist = list()
  dflist[[name1]] = exprVec1
  dflist[[name2]] = exprVec2

  df=do.call(cbind, dflist)
  df = as.data.frame(df)

  print(paste("Common Genes", length(common.names)))
  print(dim(df))
  print(head(df))

  #geom_text(label=rownames(df))+

  p = ggplot(df, aes(.data[[name1]], .data[[name2]], label=rownames(df))) +
  geom_point() + geom_text_repel(max.overlaps=20, direction="both", force=2) + ggtitle(paste(length(common.names), "measurements")) +
  stat_smooth(formula=formula, method = "lm")+
  stat_cor(r.accuracy = 0.01, label.y = 0.9*upperlimY, size = 4, label.sep='\n') +
  stat_regline_equation(aes(label = after_stat(eq.label)), size = 4, label.y = 0.8*upperlimY)+ theme_classic()+
  geom_vline(xintercept = 0,linetype = "dashed", colour = "red")+geom_hline(yintercept = 0,linetype = "dashed", colour = "red")

  p = p + xlim(lowerlimX, upperlimX)
  p = p + ylim(lowerlimY, upperlimY)

  return(p)
}






p = scatterNamedVectors(rna_logfc_sepsis_ctrl, prot_logfc_sepsis1_ctrl)+xlim(-1.0, 1)
save_plot(p, "prot_comparison/neutrophils.rp.sepsis1_control", 12, 8)

p = scatterNamedVectors(rna_logfc_sepsis_ctrl, prot_logfc_sepsis2_ctrl)+xlim(-1.0, 1)
save_plot(p, "prot_comparison/neutrophils.rp.sepsis2_control", 12, 8)


p = scatterNamedVectors(rna_logfc_sepsis_ctrl, prot_sepsis_combined_ctrl)+xlim(-1.0, 1)
save_plot(p, "prot_comparison/neutrophils.rpcombined.sepsis12_control", 12, 8)


p = scatterNamedVectors(prot_logfc_sepsis1_ctrl, prot_logfc_sepsis2_ctrl, name1="Block 1 Protein", name2="Block 2 Protein")
save_plot(p, "prot_comparison/neutrophils.prot.sepsis1_sepsis2", 12, 12)

p = scatterNamedVectors(prot_logfc_sepsis1_ctrl, prot_logfc_sepsis2_ctrl, name1="Block 1 Protein", name2="Block 2 Protein")
save_plot(p, "prot_comparison/neutrophils.prot.sepsis1_sepsis2", 12, 12)






#
##
###
#### Additional Plots
###
##
#

source("/mnt/input/own/projekte/pekayvaz/scripts/functions.R")

cond_expr_df = getExtendedExpressionData(obj.integrated.patient_thr, assay="RNA", group.by="condition")
cond_expr_neutro_df = getExtendedExpressionData(subset(obj.integrated.patient_thr, idents %in% c(0,1,2,3,9)), assay="RNA", group.by="condition")
write.table(cond_expr_df, "flowsets/expression_by_condition.tsv", quote=F, sep="\t", row.names=F)
write.table(cond_expr_neutro_df, "flowsets/neutro_expression_by_condition.tsv", quote=F, sep="\t", row.names=F)

cond_expr_neutro_df9 = getExtendedExpressionData(subset(obj.integrated.patient_thr, idents == 9), assay="RNA", group.by="condition")
write.table(cond_expr_neutro_df9, "flowsets/neutro9_expression_by_condition.tsv", quote=F, sep="\t", row.names=F)

p=VlnPlot(obj.integrated.patient_thr, "CD177", split.by="condition")
save_plot(p, "add_plots/vplot_cd177", 10, 3)

plotElems_cd177 = list()
plotElems_cd177[["control"]] = list(cells=cellIDForClusters(obj.integrated.patient_thr, "condition", c("control")), label="control")
plotElems_cd177[["sepsis"]] = list(cells=cellIDForClusters(obj.integrated.patient_thr, "condition", c("sepsis")), label="sepsis")
plotElems_cd177[["thrombus"]] = list(cells=cellIDForClusters(obj.integrated.patient_thr, "condition", c("thrombus")), label="thrombus")

p=makeSideBySideDotPlot(obj.integrated.patient_thr, plotElems_cd177, featureGenes=c("CD177"), group.by="idents", cols=c("red", "grey", "blue"))
save_plot(p, "add_plots/cd177", 5,12)

p=makeSideBySideDotPlot(subset(obj.integrated.patient_thr, idents %in% c(0,1,2,3,9)), plotElems_cd177, featureGenes=c("CD177"), group.by="idents", cols=c("red", "grey", "blue"))
save_plot(p, "add_plots/cd177_neutro", 5,12)

makePieTrieCounts( subset(obj.integrated.patient_thr, idents %in% c(0,1,2,3,9)), "add_plots/neutro_cluster_abundance_by_condition", "idents", "condition", size.text = 10, repel.direction="x", max.overlaps=10, fig.height=16, use.palette=NULL)
makePieTrieCounts( obj.integrated.patient_thr, "add_plots/cluster_abundance_by_condition", "cellnames_quick", "condition", size.text = 10, repel.direction="x", max.overlaps=10, fig.height=16, use.palette=NULL)

library(ComplexUpset)
allEntries = list()

baseDF = data.frame(gene=unique(rownames(obj.integrated.patient_thr)))
rownames(baseDF) = baseDF$gene

unionDF = data.frame(gene=unique(rownames(obj.integrated.patient_thr)))
rownames(unionDF) = unionDF$gene

for (cluster in names(de.sepsis.control))
{


  if (!(cluster %in% c(0,1,2,3,9)))
  {
    next()
  }

  if (cluster %in% names(de.thrombus.control))
  {

    cluster = as.character(cluster)

    thrombus.genes = de.thrombus.control[[cluster]]
    thrombus.genes.up = thrombus.genes[thrombus.genes$avg_log2FC > 0 & thrombus.genes$p_val_adj < 0.05, ]
    thrombus.genes.down = thrombus.genes[thrombus.genes$avg_log2FC < 0 & thrombus.genes$p_val_adj < 0.05, ]


    sepsis.genes = de.sepsis.control[[cluster]]
    sepsis.genes.up = sepsis.genes[sepsis.genes$avg_log2FC > 0 & sepsis.genes$p_val_adj < 0.05, ]
    sepsis.genes.down = sepsis.genes[sepsis.genes$avg_log2FC < 0 & sepsis.genes$p_val_adj < 0.05, ]

    baseDF[paste("Sepsis", "up", cluster)] = 0
    baseDF[paste("Sepsis", "down", cluster)] = 0
    baseDF[paste("Stroke", "up", cluster)] = 0
    baseDF[paste("Stroke", "down", cluster)] = 0

    baseDF[ sepsis.genes.up$gene,  paste("Sepsis", "up", cluster)] = 1
    baseDF[ sepsis.genes.down$gene,  paste("Sepsis", "down", cluster)] = 1
    baseDF[ thrombus.genes.up$gene,  paste("Stroke", "up", cluster)] = 1
    baseDF[ thrombus.genes.down$gene,  paste("Stroke", "down", cluster)] = 1


    if (!"Sepsis up" %in% names(unionDF))
    {
      unionDF[paste("Sepsis", "up")] = 0
      unionDF[paste("Sepsis", "down")] = 0
      unionDF[paste("Stroke", "up")] = 0
      unionDF[paste("Stroke", "down")] = 0
    }
    
    unionDF[ sepsis.genes.up$gene,  paste("Sepsis", "up")] = 1
    unionDF[ sepsis.genes.down$gene,  paste("Sepsis", "down")] = 1
    unionDF[ thrombus.genes.up$gene,  paste("Stroke", "up")] = 1
    unionDF[ thrombus.genes.down$gene,  paste("Stroke", "down")] = 1


  }
}
baseDF$gene=NULL
unionDF$gene = NULL

baseDF = baseDF[rowSums(baseDF) > 0,]
unionDF = unionDF[rowSums(unionDF) > 0,]

p=ComplexUpset::upset(baseDF, setdiff(colnames(baseDF), c("gene")), name='Comparison', width_ratio=0.05, min_size=5, sort_sets='ascending')
save_plot(p, "add_plots/upset_01239", 16,8)

p=ComplexUpset::upset(unionDF, setdiff(colnames(unionDF), c("gene")), name='Comparison', width_ratio=0.05, min_size=5, sort_sets='ascending')
save_plot(p, "add_plots/upset_union_01239", 16,8)







markers.use.tt=subset(exprDFOuts[["patint_thr"]]$exprdf ,avg_log2FC>0&p_val_adj<0.01&!startsWith(gene, "MT-")&!startsWith(gene, "RP"))
finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:20)
print(finalMarkers.use.tt)

p_dp_genes_idents = DotPlot(obj.integrated.patient_thr, features = unique(finalMarkers.use.tt$gene), assay="RNA", group.by="cellnames_quick")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))
save_plot(p_dp_genes_idents, paste("add_plots", "dotplot_cellnames_genes", sep="/"), 45, 8)






countByAFAbs = getCellCountDF(obj.integrated.patient_thr, prefix="", select_by = group.by, group_by="patient", relative=F, show.percent=F, outname=paste(outfolder, "ia_count_by_patient.tsv", sep="/"))


makeBarCountPlot( obj.integrated.patient_thr, "test", "idents", "patient", size.text = 10, repel.direction="x", max.overlaps=10, fig.height=16)

makeBarCountPlot( obj.integrated.patient_thr, "test", "cellnames_quick", "condition", size.text = 10, repel.direction="x", max.overlaps=10, fig.height=16)








