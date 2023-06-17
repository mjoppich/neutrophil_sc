
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library("stringr")
library(SeuratDisk)
library(SeuratWrappers)

# nohup ~/Rscript.sh combined_velocity_analysis/combined_library ../seurat_rds <all looms>

args = commandArgs(trailingOnly=TRUE)

outDirPrefix = args[1]
seuratObjFile = args[2]
loomFolder = args[3]

outDirPrefix = "patint_thr"
seuratObjFile="../../obj.integrated.patient_thr.final.Rds"
loomFolder="../../loom_folder/"

print(outDirPrefix)
print(seuratObjFile)
print(loomFolder)

#
## reading seurat
#
print("Reading Seurat")
obj.seurat = readRDS(seuratObjFile)

DefaultAssay(obj.seurat) = "RNA"
obj.seurat


seuratizeCellnames = function( velocytoCellnames, samplename ) {

    
    newCellnames = unlist(lapply(str_split(velocytoCellnames, ":"), function(x){
        x[1] = samplename
        x[2] = str_replace_all(x[2], "x", "-1")
        return(paste(x, collapse="_"))
    }))
    

    return(newCellnames)
}

#
## reading looms
#
print("Reading Velocities")

loomfiles = Sys.glob(paste(loomFolder, "*.loom", sep="/"))
print(loomfiles)

foundCellsTotal = 0

inLooms = list()
for (infile in loomfiles)
{
    print(infile)
    indat = ReadVelocity(file=infile)

    samplename = basename(infile)
    samplename = str_replace(samplename, ".loom", "")

    print(samplename)
    print(dim(indat$spliced))
    
    colnames(indat$spliced) = seuratizeCellnames(colnames(indat$spliced), samplename)
    colnames(indat$unspliced) = seuratizeCellnames(colnames(indat$unspliced), samplename)

    inLooms[[samplename]] = indat

    cellsNotFound = length(setdiff(colnames(indat$spliced), colnames(obj.seurat)))
    cellsFound = length(intersect(colnames(indat$spliced), colnames(obj.seurat)))

    
    print(paste("Cells not found: ", cellsNotFound))
    print(paste("Cells found: ", cellsFound))

    foundCellsTotal = foundCellsTotal + cellsFound
}

print(paste("Found cells", foundCellsTotal))
print(obj.seurat)

if (length(inLooms) == 1)
{

    unspliced.matrix = inLooms[[names(inLooms)[1]]]$unspliced

} else {

    unspliced.matrix = cbind(inLooms[[ names(inLooms)[1] ]]$unspliced, inLooms[[ names(inLooms)[2] ]]$unspliced)

    if (length(inLooms) > 2)
    {

        for (i in 3:length(inLooms))
        {
            print(names(inLooms)[i])
            unspliced.matrix = cbind(unspliced.matrix, inLooms[[names(inLooms)[i] ]]$unspliced)
        }

    }

}

dim(unspliced.matrix)

length(intersect(rownames(indat$spliced), rownames(obj.seurat)))
length(setdiff(rownames(indat$spliced), rownames(obj.seurat)))
length(setdiff(rownames(obj.seurat), rownames(indat$spliced)))

missingCells = setdiff(colnames(obj.seurat), colnames(unspliced.matrix))

table(unlist(lapply(str_split(missingCells, "_"), function(x){return(paste(x[1], x[2], sep="_"))})))


print(paste("Missing cells:", length(missingCells)))

commonCells = intersect(colnames(unspliced.matrix), colnames(obj.seurat))
commonGenes = intersect(rownames(indat$spliced), rownames(obj.seurat))

obj.common=obj.seurat[commonGenes, commonCells]
obj.common[["unspliced"]] <- CreateAssayObject(counts = unspliced.matrix[commonGenes,commonCells])

obj.common$identsstr = as.character(obj.common$idents)

SaveLoom(obj.common, paste(outDirPrefix, ".loom", sep=""), overwrite = TRUE, verbose = TRUE)
print("unspliced")
write.csv(x = t(as.matrix(obj.common@assays$unspliced@counts)), file = paste(outDirPrefix, "_unspliced.csv", sep=""))

print("RNA")
write.csv(x = t(as.matrix(obj.common@assays$RNA@counts)), file = paste(outDirPrefix, "_RNA.csv", sep=""))

print("UMAP")
write.csv(Embeddings(obj.common, reduction = "umap"), file = paste(outDirPrefix, "_umap_embeddings.csv", sep=""))

print("META")
write.csv(x = obj.common@meta.data, file = paste(outDirPrefix, "_meta.csv", sep=""))

