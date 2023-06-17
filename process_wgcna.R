## Supplementary Data for Kazer et al., 
## "An Integrated Single-Cell Analysis of Multicellular Immune Dynamics during Hyper-Acute HIV-1 Infection."
##
## Generated on 01/17/20 by SWK
## Modified 2021 by MJ for other manuscript

## Required packages for code
require(WGCNA)
require(flashClust)
require(Hmisc)
require(dplyr)

## Required packages for saving results
require(openxlsx)
require(ggplot2)
require(cowplot)
require(Seurat)
library(RColorBrewer)

#######################################################
## Functions for Module Discovery Adapted from WGCNA ##
#######################################################
combine_plot_grid = function(..., ncol=1, rel_heights=NULL)
{
  inplot = list(...)
  print(inplot)
  
  dataList = list()
  for (x in inplot)
  {
    print(x)
    dataList[[length(dataList)+1]] = x$data
  }
  
  p=cowplot::plot_grid(inplot, ncol = ncol, rel_heights = rel_heights)
  
  p$data = dataList
  
  return(p)
}


combine_plot_grid_list = function(plotlist, ...)
{
 
  dataList = list()
  for (i in 1:length(plotlist))
  {
    dataList[[i]] = plotlist[[i]]$data
  }
  
  p=cowplot::plot_grid(plotlist=plotlist, ...)
  
  p$data = dataList
  
  return(p)
}


save_plot = function(plotobj,outname, fig.width, fig.height)
{
  print(paste(outname, fig.width, fig.height))

  gc()
  
  fname=paste(outname, "png", sep=".")
  print(paste("Saving to file", fname))
  png(filename=fname, width = fig.width, height = fig.height, units = 'in', res = 300)#width = fig.width*100, height=fig.height*100)
  plot(plotobj)
  dev.off()

  gc()
  
  fname=paste(outname, "pdf", sep=".")
  print(paste("Saving to file", fname))
  pdf(file=fname, width = fig.width, height=fig.height)
  plot(plotobj)
  dev.off()
  
  gc()

  fname=paste(outname, "svg", sep=".")
  print(paste("Saving to file", fname))
  svglite::svglite(file = fname, width = fig.width, height = fig.height)
  plot(plotobj)
  dev.off()

  gc()
  
  if (class(plotobj$data) %in% c("list"))
  {
    print("list case")
    for (i in 1:length(plotobj$data))
    {
      fname = paste(outname,i, "data", sep=".")
      print(paste("Saving to file", fname))
      
      if (class(plotobj$data[[i]]) %in% c("list"))
      {
          print("multi list case")
          for (j in 1:length(plotobj$data[[i]]))
          {
              fname = paste(outname,i, j, "data", sep=".")
              print(paste("Saving to file", fname, class(plotobj$data[[i]][[j]])))

              if (class(plotobj$data[[i]][[j]]) %in% c("list", "waiver"))
              {
                next()
              }
              write.table(plotobj$data[[i]][[j]], fname, row.names = TRUE, sep="\t")    

          }
      } else {
          
          tryCatch(write.table(plotobj$data[[i]], fname, row.names = TRUE, sep="\t"), error = function(e) NULL)
      }
      
    }
  } else {
    
      fname = paste(outname,"data", sep=".")
      print(paste("Saving to file", fname))

      write.table(plotobj$data, paste(outname, "data", sep="."), row.names = TRUE, sep="\t")
  }
  
  return(plotobj)
}


# Choosing the appropriate power for generating the adjacency matrix.
FindPower <- function(datExpr){
  #choose soft-threshold power
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,
                        corOptions = list(use = 'p', method = "pearson"),networkType = "signed")
  
  # Plot the results
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
  
  # Red line corresponds to using an R^2 cut-off
  abline(h=0.80,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

# Generating the adjacency matrix and performing clustering
ClusterTOM <- function(datExpr, softPower, outname){
  #dev.off()
  #Calclute the adjacency matrix
  adj= adjacency(datExpr,type = "signed", power = softPower);
  
  #Turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations.
  TOM=TOMsimilarity(adj, TOMType = "signed");
  
  colnames(TOM) = rownames(TOM) = colnames(datExpr)
  dissTOM=1-TOM
  
  #Hierarchical clustering of the genes based on the TOM dissimilarity measure
  geneTree = flashClust(as.dist(dissTOM),method="complete");
  
  #Plot the resulting clustering tree (dendrogram)
  fname=paste(outname, "png", sep=".")
  print(paste("Saving to file", fname))
  png(filename=fname, width = 8, height = 6, units = 'in', res = 300)
  plot(geneTree, xlab="", sub="",cex=0.3)
  dev.off()
  
  return(list(adj=adj,dissTOM = dissTOM, geneTree = geneTree)) #returns list with dissimilarity TOM, and the clustered gene tree.
}

# Cut the resulting clustering dendrogram using the "tree" method for cutreeDynamic. Minimum module size can be specified.
CutTOMTree <- function(datExpr, dissTOM, geneTree, minModuleSize = 10){
  #dev.off()
  # Module identification using dynamic tree cut, you can also choose the hybrid method
  dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
  #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
  
  #Get the module labels and the size of each module. Lable 0 is reserved for unassigned genes
  print(table(dynamicMods))
  
  #Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  
  #Set the diagonal of the dissimilarity to NA 
  diag(dissTOM) = NA;
  
  #extract modules
  module_colors= setdiff(unique(dynamicColors), "grey")
  modules = lapply(module_colors, function(x){colnames(datExpr)[which(dynamicColors==x)]})
  names(modules) = module_colors
  return(list(dyanmicColors = dynamicColors, modules = modules)) #returns list with module colors, and the modules themselves
}

# Merge modules with low dissimilarity. Cutoff for dissimilarity merge can be specified
MergeSimilarModules <- function(datExpr, dynamicColors, geneTree, MEDissThres = 0.5){
  #cacluate eigengenes
  MEList = moduleEigengenes(datExpr, colors=dynamicColors, excludeGrey=TRUE)
  MEs = MEList$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs, use="pairwise.complete.obs");

  print(is.na(MEDiss) %>% table())
  
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h = MEDissThres, lwd=2, col="red")
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  #plot showing how merged modules exist
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  #extract merged modules
  merged_module_colors= setdiff(unique(mergedColors), "grey")
  merged_modules = lapply(merged_module_colors, function(x){colnames(datExpr)[which(mergedColors==x)]})
  names(merged_modules) = merged_module_colors
  
  return(list(mergedColors = mergedColors, merged_modules = merged_modules)) #returns list with merged colors, and the merged modules themselves
}

############################################################################################################
## Functions for assaying module significance against background and temporal variability in module score ##
############################################################################################################

## Test to determine if the genes within the module are truly the least dissimilar compared to randomly generated modules of the same size.
TestModuleSignificance <- function(mod, dissTOM, expr.data, n_perm = 10000, pval = 0.05, n.bin = 10){
  #vectorize the actual distribution of (dis)similarities, and remove zeros!
  true.diss = as.vector(dissTOM[mod,mod])
  true.diss = true.diss[-which(true.diss == 0)]
  
  #size of module for permutations
  mod.size = length(mod)
  
  #bin all genes by expression
  expr.avg = rowMeans(expr.data)
  expr.avg = expr.avg[order(expr.avg)]
  expr.avg.cut = as.numeric(x = cut2(x = expr.avg, m=round(length(expr.avg)/n.bin)))
  names(expr.avg.cut) = names(expr.avg)
  
  #create a table of binnings of all genes and of our module genes
  all.bin.table = table(expr.avg.cut)
  mod.bin.table = table(expr.avg.cut[mod])
  
  #randomly generate module with same expression binning structure and run t.test and record results
  test.results = data.frame(statistic = rep(NA, n_perm), pval = rep(NA, n_perm)) #container for results
print("Starting permutations")
  for (i in 1:n_perm){ #by permutation
    random.mod = list() #create an empty list we will fill with gene names (each element will be gene names by bin)
    
    #for each bin of the mod bin table, randomly select that number of genes from the full set with that expression
    for (j in names(mod.bin.table)){ 
      bin.genes = sample(names(expr.avg.cut)[which(expr.avg.cut == as.numeric(j))], mod.bin.table[j], replace = FALSE)
      random.mod[[as.numeric(j)]] = bin.genes #stick those genes into the random.mod list
    }
    #unlist and vectorize the distribution of (dis)similarities (remember to remove zeros)
    random.mod = unlist(random.mod)
    random.diss = as.vector(dissTOM[random.mod,random.mod])
    random.diss = random.diss[-which(random.diss == 0)]
    
    #perform a one-sided wilcox.test and record the statistic and p-value.
    #Note, IMPORTANT: here we perform the test asking if the true diss is LESS THAN the random diss, as we are trying to minimize dissimilarity
    test = wilcox.test(x = true.diss, y = random.diss, alternative = "less")
    test.results[i,] = c(test$statistic, test$p.value)
  }
  
  #correct for multiple hypothesis testing, and then report the proportion of bad tests
  test.results$FDR = p.adjust(test.results$pval, method = "fdr")
  num.failed.tests = sum(test.results$FDR > pval)
  print(paste(paste(num.failed.tests, n_perm, sep="/"), "permutations failed the Mann-Whitney test.", sep=" "))
  
  #is the percentage of failed tests less than or equal to the p-val specified?
  return(num.failed.tests/n_perm <= pval) #returns a vector of booleans indicating if each module was significant based on the specific p.val
}

## Test for variation in module score as a function of time. Compares many samplings of scores between pre-infection and each time point.
# Here, our time points are labeled as specified in tps; pre = pre-infection, 2001 = 0 Weeks, 2002 = 1 Week, ..., 2024 = 6 Months, 2024 = 1 Year
# meta.data here is the meta.data data.frame that exists within a Seurat object
TestModuleTemporalVariation <- function(tps = c("pre","2001","2002","2003","2004","2005","2024","2048"), #which time points to run the function on
                                        meta.data, sample.size = 50, ntest = 1000, name.of.feature = "M"){ #name.of.feature refers to prefix for columns with module scores

  #apply to run over one module at a time
  mod.min.pvals <- apply(meta.data[,grep(name.of.feature, colnames(meta.data), value=TRUE)], 2, function(mod){
    #get the scores for this cluster separated by time points in a list
    scores.by.tp = lapply(tps, function(time) mod[meta.data$TimePoint == time])
    print(as.vector(unlist(lapply(scores.by.tp, length))))
    minAllGroups = min(as.vector(unlist(lapply(scores.by.tp, length))))
    
    #run the wilcox test ntest times with sample.size between each timepoint and pre-infection, report the average p-val
    wilcox.pval.by.tp = mapply(scores = scores.by.tp[-1], tp = tps[-1], SIMPLIFY=FALSE, function(scores,tp){
      
      subsample.size = min(c(minAllGroups-1, length(scores)-1, sample.size))
      #print(paste(tp, length(scores), subsample.size))


      test.pval = replicate(ntest, expr={
        wilcox.test(x=sample(scores.by.tp[[1]], subsample.size),y=sample(scores, subsample.size))$p.value}) # CHANGED
      return(mean(test.pval))
    })
    #return the min p-value for all comparison within that cluster, FDR corrected
    min.pval = min(p.adjust(unlist(wilcox.pval.by.tp)))
    return(min.pval)
  })
  return(mod.min.pvals) #returns the vector of minimum average p-value across all time point comparisons for each module
}




na.pad <- function(x,len){
    x[1:len]
}

makePaddedDataFrame <- function(l,...){
    maxlen <- max(sapply(l,length))
    data.frame(lapply(l,na.pad,len=maxlen),...)
}

process_with_groups = function( inobj, inPCs, outfolder, quitAfterJackstraw=F, max.per.pc=250, assay="integrated_gex")
{

  outfolder = paste(outfolder, max.per.pc, sep="_")
  print(outfolder)

  dir.create(file.path(outfolder))

  p=DimPlot(object = inobj, reduction="umap")
  save_plot(p, paste(outfolder, "/", "dimplot", sep=""), fig.width=15, fig.height=6)


    print(paste("Jackstraw on assay", assay))
    inobj = JackStraw(inobj, dims=20, assay=assay)
    inobj = ScoreJackStraw(inobj, dims=1:20)

    p=JackStrawPlot(object = inobj, dims = 1:20)
    save_plot(p, paste(outfolder, "/", "jackstraw", sep=""), fig.width=15, fig.height=6)


  S.inobj = c()
  pcList = list()

  print("Fetching PC genes")
  for (i in 1:inPCs)
  {
    print(i)

    selPCGenes = PCASigGenes(inobj, pcs.use = i, max.per.pc=max.per.pc)
    print(paste(i, length(selPCGenes)))
    pcList[[paste("PC_", i)]] = selPCGenes

    S.inobj = c(S.inobj, selPCGenes)
  }

  genesOfInterest = unique(c("MT2A", "ISG15", "LY6E", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "IFI44L", "IFI6", "MX1", "IFI27",  "IFI44L", "RSAD2", "SIGLEC1", "IFIT1", "ISG15"))

  goiIntersect = intersect(genesOfInterest, unique(S.inobj))
  print(paste("Overlap genesOfInterest", length(genesOfInterest), "of", length(goiIntersect)))
  print(goiIntersect)

  padded.pcList = makePaddedDataFrame(pcList)
  write.table(padded.pcList, paste(outfolder, "/", "pcgenes.tsv", sep=""), sep="\t", row.names=F, quote = F)
  outxlsx = paste(outfolder, "/", "pcgenes.xlsx", sep="")
  if (file.exists(outxlsx))
  {
    file.remove(outxlsx)
  }
  write.xlsx(pcList, outxlsx) #save to file


  S.inobj = unique(S.inobj)
  print(paste("Total Genes for Analysis", length(S.inobj)))

  #S.CD4 = ProjectPCA(inobj, genes.print = 4, do.print = FALSE)
  #S.inobj = PCASigGenes(inobj, pcs.use = 1:inPCs, pval.cut = 0.1, max.per.pc=50)
  #print(S.inobj)

  if (quitAfterJackstraw)
  {
    return(NULL)
  }


  inobj.top.genes = S.inobj#unique(PCTopGenes(S.CD4, pc.use=1:mono.nPCs, num.genes=50, do.balanced=TRUE, use.full=TRUE)) # <- does not exist in any modern Seurat
  inobj.mat.top.genes = as.matrix(GetAssayData(inobj,assay="RNA")[inobj.top.genes,])


  FindPower(datExpr = t(inobj.mat.top.genes))

  # Columns correspond to genes and rows to samples.
  inobj.ClusterTOM = ClusterTOM(datExpr = t(inobj.mat.top.genes), softPower = 7, paste(outfolder, "/", "tom_clustering", sep=""))
  inobj.ColorsModules = CutTOMTree(datExpr = t(inobj.mat.top.genes), dissTOM = inobj.ClusterTOM$dissTOM,
                                geneTree = inobj.ClusterTOM$geneTree, minModuleSize = 3)


  inobj.MergedModules = MergeSimilarModules(datExpr = t(inobj.mat.top.genes), 
                                          dynamicColors = inobj.ColorsModules$dyanmicColors,
                                          geneTree = inobj.ClusterTOM$geneTree, MEDissThres = 0.5)

  inobj.MergedModules$merged_modules


  inobj.Modules.isSig = sapply(inobj.MergedModules$merged_modules, function(module){
    TestModuleSignificance(mod = module, dissTOM = inobj.ClusterTOM$dissTOM, expr.data = inobj.mat.top.genes,
                          n_perm = 100, pval = 0.05, n.bin = 10)
  })

  inobj.Sig.Modules = inobj.MergedModules$merged_modules[inobj.Modules.isSig] # keep those modules that are significant


  inobj = AddModuleScore(inobj, features = inobj.Sig.Modules, assay="RNA", ctrl=5, name = "M")
  ctrlTPs = inobj$mpoint
  inobj$TimePoint = ctrlTPs
  unique(ctrlTPs)

  p=DotPlot(inobj, unique(genesOfInterest), group.by="TimePoint", assay="RNA")
  save_plot(p, paste(outfolder, "/", "isg_dotplot", sep=""), fig.width=15, fig.height=6)

  Timepoint.palette = brewer.pal(4, "Dark2")

  inobj.gene.clusters.names = grep("M", colnames(inobj@meta.data), value = TRUE)
  inobj.geneclusts.time = lapply(inobj.gene.clusters.names, function(c_name){
    p = ggplot(inobj@meta.data, aes_string(x="TimePoint", y=c_name, fill="TimePoint")) +
      geom_boxplot() + scale_fill_manual(values = Timepoint.palette) + theme(legend.position = "none")
    return(p)
  })
  p=combine_plot_grid_list(plotlist = inobj.geneclusts.time)
  save_plot(p, paste(outfolder, "/", "geneclust", sep=""), fig.width=12, fig.height=12)

  print(table(inobj@meta.data$TimePoint))


  inobj.clust.wilcox.tests.minp = TestModuleTemporalVariation(meta.data = inobj@meta.data, sample.size = 150, ntest=1000, tps = c("0", "1", "2"))
  inobj.clust.wilcox.tests.minp[is.na(inobj.clust.wilcox.tests.minp)] = 1
  print(inobj.clust.wilcox.tests.minp) #print the results

  inobj.Sig.Var.Modules = inobj.Sig.Modules[inobj.clust.wilcox.tests.minp <= 0.05]
  inobj.Sig.Var.Modules = inobj.Sig.Modules[!is.na(names(inobj.Sig.Modules))]

  onobj.sig.var.lengths = lapply(inobj.Sig.Var.Modules, function(x) {length(x)})

  inobj.Sig.Var.Modules = inobj.Sig.Var.Modules[ onobj.sig.var.lengths >= 5 ]

  print( inobj.Sig.Var.Modules )

  print("Before returning results")

  
  if (length(inobj.Sig.Var.Modules))
  {

    print("Got results to report!")
    inobj@meta.data = inobj@meta.data[,-grep("M", colnames(inobj@meta.data), ignore.case = FALSE)]

    print("Adding module score")
    inobj = AddModuleScore(inobj, features = inobj.Sig.Var.Modules, ctrl=5, name = "M", assay="RNA")

    print("Getting final module scores")
    inobj.gene.clusters.names.final = grep("M", colnames(inobj@meta.data), ignore.case = FALSE, value = TRUE)

    # rename the modules with the appropraite naming scheme (e.g. M1.CD4)
    names(inobj.Sig.Var.Modules) = inobj.gene.clusters.names.final

    print(inobj.Sig.Var.Modules)
    
    outxlsx = paste(outfolder, "/", "final_modules.xlsx", sep="")

    if (file.exists(outxlsx))
    {
      file.remove(outxlsx)
    }

    write.xlsx(inobj.Sig.Var.Modules, outxlsx) #save to file
    padded.inobj.Sig.Var.Modules = makePaddedDataFrame(inobj.Sig.Var.Modules)
    write.table(padded.inobj.Sig.Var.Modules, paste(outfolder, "/", "final_modules.tsv", sep=""), sep="\t", row.names=F, quote = F)


    # plot box & whiskers of final modules
    inobj.geneclusts.time.final = lapply(names(inobj.Sig.Var.Modules), function(c_name){
      p = ggplot(inobj@meta.data, aes_string(x="TimePoint", y=c_name, fill="TimePoint")) +
        geom_boxplot() + scale_fill_manual(values = Timepoint.palette) + theme(legend.position = "none")
      return(p)
    })
    p=combine_plot_grid_list(plotlist = inobj.geneclusts.time.final)
    save_plot(p, paste(outfolder, "/", "geneclust_time", sep=""), fig.width=12, fig.height=12)

    inobj$TimePointStr = paste("TP", inobj$TimePoint, sep="")

    for (mname in names(inobj.Sig.Var.Modules))
    {
    p=DotPlot(inobj, features=inobj.Sig.Var.Modules[[mname]], group.by="TimePointStr", assay="RNA")+theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=0.5))
    save_plot(p, paste(outfolder, "/", "dotplot_", mname, sep=""), fig.width=30, fig.height=6)


    p=DoHeatmap(inobj, features=inobj.Sig.Var.Modules[[mname]], group.by="TimePointStr", slot="data")+ scale_fill_gradientn(colors = c("white", "red"))
    save_plot(p, paste(outfolder, "/", "heatmap_", mname, sep=""), fig.width=20, fig.height=15)


    }

  }


  return(list("sobj"=inobj, "sigvarmods"=inobj.Sig.Var.Modules))
}


outfolder="wgcna"
if (!dir.exists(outfolder))
{
  print(paste("Creating folder", outfolder))
  dir.create(outfolder, recursive = TRUE)
}


obj.integrated = readRDS("../obj.integrated.patient_thr.final.Rds")
unique(obj.integrated$idents)

tplist = obj.integrated$condition

tplist[tplist == "control"] = "0"
tplist[tplist == "thrombus"] = "1"
tplist[tplist == "sepsis"] = "2"

obj.integrated$mpoint = tplist
obj.integrated$mpoint = factor(obj.integrated$mpoint, levels=c("0", "1", "2"))

unique(obj.integrated$mpoint)

reslist_500 = list()
for (cname in sort((unique(obj.integrated$idents)), decreasing=TRUE))
{

    gc()
    scobj = subset(obj.integrated, idents == cname)
    name = paste(paste(outfolder,"idents", sep="/"), cname, sep="_")
    
    
    maxPerPC = 250
    res = process_with_groups(scobj, 8, name, quitAfterJackstraw=F, max.per.pc=maxPerPC, assay="integrated_gex")
    reslist_500[[cname]] = res
}
saveRDS(reslist_500, "../reslist_500_wgcna.Rds")

give.n <- function(x){
   return(c(y = mean(x), label = length(x)))
}

makeGenePlot = function(res, outname)
{
    inobj.Sig.Var.Modules = res$sigvarmods
    inobj = res$sobj
  Timepoint.palette = brewer.pal(4, "Dark2")

    inobj.geneclusts.time.final = lapply(names(inobj.Sig.Var.Modules), function(c_name){
    
    plotData = inobj@meta.data[, c("TimePoint", c_name, names(inobj.Sig.Var.Modules))]

    p = ggplot(plotData, aes_string(x="TimePoint", y=c_name, fill="TimePoint")) +
        geom_boxplot()+ stat_summary(fun.data = give.n, geom = "text") + scale_fill_manual(values = Timepoint.palette) + theme(legend.position = "none")
      return(p)
    })

    p2=combine_plot_grid_list(plotlist = inobj.geneclusts.time.final)
    save_plot(p2, paste(outname, "geneclust_time", sep=""), fig.width=12, fig.height=12)
    return(inobj.geneclusts.time.final)
}

p=makeGenePlot(reslist_500, "sourcedata_n/all_500_new_")

save.image("../covid_rev2_analyses.Rdata", compress=F)


get_gm_genes = function( resobj )
{

  objList = data.frame(matrix(nrow = 0, ncol = 2)) 
  colnames(objList) = c("module", "gene")

    for (mname in names(resobj[["sigvarmods"]]))
    {
      for (mgene in resobj[["sigvarmods"]][[mname]])
      {
        objList[nrow(objList) + 1, ] = c(mname, mgene)
      }
    }

    return(objList)
}


get_gm_genes_reslist = function(reslist)
{

  outDF = NULL

    for (mname in names(reslist))
    {
        mDF = get_gm_genes(reslist[[mname]])
        mDF$analysis = mname

        if (is.null(outDF))
        {
          outDF = mDF
        } else {
          outDF = rbind(outDF, mDF)
        }

    }

    return(outDF)
}

allModuleGenes_pooled = get_gm_genes_reslist(reslist_500)
allModuleGenes = allModuleGenes_pooled

write.table(allModuleGenes, paste("allModuleGenes.tsv", sep=""), sep="\t", row.names=F, quote = F)

get_gm_medscore = function( resobj )
{

  outDF = NULL

    for (mname in names(resobj[["sigvarmods"]]))
    {

      moduleScoresTP = resobj[["sobj"]]@meta.data[, c("TimePoint", mname)]
      colnames(moduleScoresTP) = c("TimePoint", "GM")

      summedDF = moduleScoresTP %>%
      group_by(TimePoint) %>%
      summarise(vmean = mean(GM), vsd = sd(GM), vmedian = median(GM))
      
      summedDF$module = mname
      
      if (is.null(outDF))
      {
        outDF = summedDF
      } else {
        outDF = rbind(outDF, summedDF)
      }
        
    }

    return(outDF)
}

get_gm_medscore_reslist = function(reslist)
{

  outDF = NULL

    for (mname in names(reslist))
    {
        mDF = get_gm_medscore(reslist[[mname]])
        mDF$analysis = mname

        if (is.null(outDF))
        {
          outDF = mDF
        } else {
          outDF = rbind(outDF, mDF)
        }

    }

    return(outDF)
}

allModuleScores_pooled = get_gm_medscore_reslist(reslist_500)
allModuleScores = allModuleScores_pooled
write.table(allModuleScores, paste("allModuleScores.tsv", sep=""), sep="\t", row.names=F, quote = F)

relExpsList = list(c("bcells.asympt_500", "M12"),
c("bcells.sympt_500", "M11"),
c("cd4.asympt_500", "M9"),
c("cd4.sympt_500", "M10"),
c("mono.asympt_500", "M2"),
c("mono.sympt_500", "M3"),
c("nk.asympt_500", "M7"),
c("nk.sympt_500", "M5"))

relExpsDF = NULL
for (expPair in relExpsList)
{

  pDF = allModuleScores[allModuleScores$analysis == expPair[1] &  allModuleScores$module == expPair[2], ]
  print(paste(expPair, nrow(pDF)))

  if (is.null(relExpsDF))
  {
    relExpsDF = pDF
  } else {
    relExpsDF = rbind(relExpsDF, pDF)
  }

}

library(data.table)
relExpsDF.asympt = relExpsDF[relExpsDF$analysis %like% ".asympt",]
relExpsDF.sympt = relExpsDF[relExpsDF$analysis %like% "\\.sympt",]


scale.colors <- c("#4575b4",
          "darkgrey",
          "#d73027")
p = ggplot(relExpsDF.asympt, aes(TimePoint, analysis, colour = vmedian, size = vmedian)) +
  geom_point(shape=18) +
  scale_colour_gradientn(colours = scale.colors, limits=c(-0.7, 0.5)) +
  theme(legend.position="bottom") +
  scale_size(range = c(1,15), limits=c(-0.7, 0.5)) +
  guides(color= guide_legend(title="Median Score"), size=guide_legend(title="Median Score")) + ggtitle("Module Score Plot Asympt+Ctrl")+
  theme(axis.line = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text=element_text(size=20))+ xlab("Time Point") + ylab("Subset")
save_plot(p, "module_score_plot_asympt", fig.width=10,fig.height=6)

p = ggplot(relExpsDF.sympt, aes(TimePoint, analysis, colour = vmedian, size = vmedian)) +
  geom_point(shape=18) +
  scale_colour_gradientn(colours = scale.colors, limits=c(-0.7, 0.5)) +
  theme(legend.position="bottom") +
  scale_size(range = c(1,15), limits=c(-0.7, 0.5)) +
  guides(color= guide_legend(title="Median Score"), size=guide_legend(title="Median Score")) + ggtitle("Module Score Plot Sympt+Ctrl")+
  theme(axis.line = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text=element_text(size=20)) + xlab("Time Point") + ylab("Subset")
save_plot(p, "module_score_plot_sympt", fig.width=10,fig.height=6)


save.image("../covid_rev2_analyses2.Rdata", compress=F)

