#!/usr/bin/env Rscript


source("Main.R")
source("Utilities.R")
source("Conf.R")
library(ArchR)
library(parallel)
# library(foreach)
# library(doParallel)


setwd(projectDir)


generateDiffPeaks <- function(forGroup, bckGroup, dosageFor, dosageBack, filePath){
    print(paste0(forGroup,"_",dosageFor))
    print(paste0(bckGroup,"_",dosageBack))
    markerTest <- getMarkerFeatures(
          ArchRProj = myProj, 
          useMatrix = "PeakMatrix",
          groupBy = "compoundDosage",
          testMethod = "wilcoxon",
          bias = c("TSSEnrichment", "log10(nFrags)"),
          useGroups = paste0(forGroup,"_",dosageFor),
          bgdGroups = paste0(bckGroup,"_",dosageBack)
        )

    resPeaks = rowData(markerTest)
    resPeaks$Log2FC = unlist(assays(markerTest)$Log2FC)
    resPeaks$Pval = unlist(assays(markerTest)$Pval)
    resPeaks$FDR = unlist(assays(markerTest)$FDR)
    resPeaks = resPeaks[!is.na(resPeaks$Pval),]
    resPeaks = resPeaks[!is.na(resPeaks$FDR),]
    resPeaks$Compound = forGroup
    resPeaks$Background = paste0(bckGroup,"_",dosageBack)
    resPeaks$dosage = dosageFor
    
    saveRDS(resPeaks, paste0(filePath,forGroup, "_", dosageFor,"_",bckGroup,"_",dosageBack ,".rds"))
    
    #saveRDS(markerTest, paste0(filePath, "SummarizedExperiments/",forGroup, "_", dosage, "_SumExp.rds"))

}




myProj = loadArchRProject(path = "/data/Osaka_Omics/scATAC/OsakaMomicsATAC", force = FALSE, showLogo = TRUE)
deneme <- ArchR::getCellNames(myProj)
print(length(deneme))
myMetaData = getCellColData(myProj)
dosage="10nM"

r <- mclapply(1:length(unique(myMetaData$Compound)), function(i) {
         if(paste0(unique(myMetaData$Compound)[i], "_", dosage) %in% unique(myMetaData$compoundDosage)){
             generateDiffPeaks(unique(myMetaData$Compound)[i],
                           "NOSTIM",
                           dosageFor=dosage,
                           dosageBack="10uM",
                           filePath = "/data/Osaka_Omics/scATAC/OsakaMomicsATAC/DiffPeaks/NOSTIM/")  
         }     
 }, mc.cores = 80)      



# myDat = DataFrame()

# for(i in 1:length(unique(myMetaData$Compound))){

#     if(paste0(unique(myMetaData$Compound)[i], "_", "10nM") %in% unique(myMetaData$compoundDosage)){
    
#           tmp <- readRDS(paste0("/data/Osaka_Omics/scATAC/OsakaMomicsATAC/DiffPeaks/10nM/",
#                             unique(myMetaData$Compound)[i],
#                                 "_", 
#                                 "10nM",
#                                 ".rds"))
#           print(paste0(unique(myMetaData$Compound)[i], "_", "10nM"))
#           myDat <- rbind(myDat, tmp)
        
#           print(dim(myDat))
    
#     }

# }


# saveRDS(myDat, paste0("/data/Osaka_Omics/scATAC/OsakaMomicsATAC/DiffPeaks/AllDiffPeaks_10nM.rds"))