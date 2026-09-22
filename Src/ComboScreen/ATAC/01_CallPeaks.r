#!/usr/bin/env Rscript


source("Main.R")
source("Utilities.R")
source("Conf.R")
library(ArchR)

setwd("/data/Osaka_Omics/scATAC/")


addArchRThreads(threads = 90) 

addArchRGenome("hg38")

# selectedCells = read.csv("/home/eraslab1/Projects/EpigeneticInhibitors/TextFiles/RNA_cells.csv", row.names=1)


# selectedCells$sampleName = sapply(rownames(selectedCells), 
#                                   function(x){strsplit(x,"_")[[1]][2]})
# selectedCells$cellBarcode = sapply(rownames(selectedCells), 
#                                    function(x){strsplit(x,"-")[[1]][1]})
# selectedCells$selCells = paste0("atac_fragments_",selectedCells$sampleName,"#", selectedCells$cellBarcode,"-1")
# rownames(selectedCells) = selectedCells$selCells

# selectedCells[selectedCells$Compound == "STIM", "dosage"] = "10uM"
# selectedCells[selectedCells$Compound == "NOSTIM", "dosage"] = "10uM"

# ArrowFiles = paste0(atacFilesDir ,"ArchRArrowFiles/atac_fragments_", 1:28,".arrow")

# projOsakaMomics <- ArchRProject(
#   ArrowFiles = ArrowFiles, 
#   outputDirectory = paste0(atacFilesDir ,"OsakaMomicsATAC"),
#   copyArrows = FALSE 
# )


# proj = projOsakaMomics[selectedCells$selCells[selectedCells$selCells %in% projOsakaMomics$cellNames], ]

# rownames(selectedCells) = selectedCells$selCells
# selectedCells = selectedCells[proj$cellNames,]

# proj$dosage = selectedCells$dosage
# proj$Compound = selectedCells$Compound
# proj$Target = selectedCells$Target
# proj$sampleName = selectedCells$sampleName
# proj$cellBarcode = selectedCells$cellBarcode


# selectedCells$TSSEnrichment =  getCellColData(proj, select = "TSSEnrichment")$TSSEnrichment

# nCells = data.frame(table(selectedCells[,c("dosage","Compound")]))

# nCells$compoundDosage = paste0(nCells$Compound, "_", nCells$dosage)

# depletedCells = nCells[nCells$Freq < 100,]

# proj$compoundDosage = paste0(proj$Compound, "_", proj$dosage)
# proj = proj[proj$compoundDosage %ni% depletedCells$compoundDosage,]

# projGrouped <- addGroupCoverages(ArchRProj = proj, 
#                                  minCells = 95,
#                                  groupBy = "compoundDosage",
#                                  threads=90,
#                                  minReplicates=5)

# saveArchRProject(ArchRProj = projGrouped, 
#                  outputDirectory = paste0(atacFilesDir ,"OsakaMomicsATAC"), 
#                  load = FALSE)

projGrouped <- loadArchRProject(path = paste0(atacFilesDir ,"OsakaMomicsATAC"), force = FALSE, showLogo = TRUE)

pathToMacs2 <- findMacs2()

projPeaks <- addReproduciblePeakSet(
    ArchRProj = projGrouped, 
    groupBy = "compoundDosage",
    reproducibility = "5",
    maxPeaks=300000,
    threads=90,
    minCells=90,
    cutOff=0.01,
    method="q",
    pathToMacs2 = pathToMacs2
)

print("aloooooo")

projPeaks2 <- addPeakMatrix(projPeaks)

# print("aloooooo 2222222")
# # # k <- data.frame(getPeakSet(projPeaks2))

# # write.csv(k, paste0(atacFilesDir ,"OsakaMomicsATAC/AllPeakSet_v2.csv"), row.names=FALSE)


saveArchRProject(ArchRProj = projPeaks2, outputDirectory = paste0(atacFilesDir ,"OsakaMomicsATAC"), load = FALSE)



