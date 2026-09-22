.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(source) & is.null(installInfo)){
      if(tolower(source) == "cran"){
        installInfo <- paste0('install.packages("',x,'")')
      }else if(tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      }else{
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}


.safelapply <- function(..., threads = 1, preschedule = FALSE){

  if(tolower(.Platform$OS.type) == "windows"){
    threads <- 1
  }

  if(threads > 1){
    .requirePackage("parallel", source = "cran")
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)

    errorMsg <- list()

    for(i in seq_along(o)){ #Make Sure this doesnt explode!
      if(inherits(o[[i]], "try-error")){
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
      }
    }

    if(length(errorMsg) != 0){

      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)

    }

  }else{

    o <- lapply(...)

  }

  o

}
                                
                                
availableSeqnames <- function(ArrowFiles = NULL, subGroup = "Fragments", threads = getArchRThreads()){
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  seqList <- .safelapply(seq_along(ArrowFiles), function(x){
    seqnames <- h5ls(ArrowFiles[x]) %>% {.[.$group==paste0("/",subGroup),]$name}
    seqnames <- seqnames[!grepl("Info", seqnames)]
    seqnames
  }, threads = threads)
  if(!all(unlist(lapply(seq_along(seqList), function(x) identical(seqList[[x]],seqList[[1]]))))){
    stop("Not All Seqnames Identical!")
  }
  o <- h5closeAll()
  return(paste0(seqList[[1]]))
}    
                        
                        
.h5read <- function(
  file = NULL,
  name = NULL,
  method = "fast",
  index = NULL,
  start = NULL,
  block = NULL,
  count = NULL
  ){

  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- H5Fopen(file, "H5F_ACC_RDONLY")
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}
                        
                        
                        
.validArrow <- function(ArrowFile = NULL){
  o <- h5closeAll()
  if(h5read(ArrowFile,"Class")!="Arrow"){
    warning(
      "This file is not a valid ArrowFile, this most likely is a bug with previous function where the class was not added.\n",
      "To fix your ArrowFiles :\n",
      "\tlapply(getArrowFiles(ArchRProj), function(x) h5write(obj = 'Arrow', file = x, name = 'Class'))",
      "\nThis will be an error in future versions."
    )
  }
  o <- h5closeAll()
  return(ArrowFile)
}

                        
.getFragsFromArrow <- function(
  ArrowFile = NULL, 
  chr = NULL, 
  out = "GRanges", 
  cellNames = NULL, 
  method = "fast"
  ){


  o <- h5closeAll()
  ArrowFile <- .validArrow(ArrowFile)
  

  #Get Sample Name
  sampleName <- .h5read(ArrowFile, paste0("Metadata/Sample"), method = method)

  o <- h5closeAll()
  nFrags <- sum(.h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method))

  if(nFrags==0){
    if(tolower(out)=="granges"){
      output <- GRanges(seqnames = chr, IRanges(start = 1, end = 1), RG = "tmp")
      output <- output[-1,]
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- output[-1,]
    }
    return(output)
  }

  if(is.null(cellNames) | tolower(method) == "fast"){
    
    output <- .h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), method = method) %>% 
      {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$RG <- Rle(
      values = paste0(sampleName, "#", .h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues"), method = method)), 
      lengths = .h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method)
    )
    if(!is.null(cellNames)){
      output <- output[BiocGenerics::which(mcols(output)$RG %bcin% cellNames)]
    }

  }else{
    
    if(!any(cellNames %in% .availableCells(ArrowFile))){

      stop("None of input cellNames are in ArrowFile availableCells!")

    }else{

      barRle <- Rle(h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues")), h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths")))
      barRle@values <- paste0(sampleName, "#", barRle@values)
      idx <- BiocGenerics::which(barRle %bcin% cellNames)
      if(length(idx) > 0){
        output <- h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), index = list(idx, 1:2)) %>% 
          {IRanges(start = .[,1], width = .[,2])}
        mcols(output)$RG <- barRle[idx]
      }else{
        output <- IRanges(start = 1, end = 1)
        mcols(output)$RG <- c("tmp")
        output <- output[-1,]
      }
    }

  }
  
  o <- h5closeAll()

  if(tolower(out)=="granges"){
    if(length(output) > 0){
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)    
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)
      output <- output[-1,]
    }
  }

  return(output)
}

                        
                        
createMyGroupBW <- function(
  i = NULL, 
  cellGroups = NULL,
  ArrowFiles = NULL, 
  cellsInArrow = NULL, 
  availableChr = NULL,
  chromLengths = NULL, 
  tiles = NULL,
  ceiling = 10,
  tileSize = 1,
  normMethod = "none",
  bwDir = bwDir2,
  threads = 1
  ){


  #Cells
  cellGroupi <- cellGroups[[i]]
  #print(sum(normBy[cellGroupi, 1]))

  #Bigwig File!
  covFile <- file.path(bwDir, paste0(make.names(names(cellGroups)[i]),
                                     "-TileSize-",tileSize,
                                     "-normMethod-",normMethod,"-ArchR.bw"))
  rmf <- .suppressAll(file.remove(covFile))

  covList <- .safelapply(seq_along(availableChr), function(k){

    it <- 0

    for(j in seq_along(ArrowFiles)){
      cellsInI <- sum(cellsInArrow[[names(ArrowFiles)[j]]] %in% cellGroupi)
      if(cellsInI > 0){
        it <- it + 1
        if(it == 1){
          fragik <- .getFragsFromArrow(ArrowFiles[j], 
                                        chr = availableChr[k], 
                                        out = "GRanges", 
                                        cellNames = cellGroupi)
        }else{
          fragik <- c(fragik, .getFragsFromArrow(ArrowFiles[j], 
                                                  chr = availableChr[k], 
                                                  out = "GRanges", 
                                                  cellNames = cellGroupi))
        }
      }
    }

    tilesk <- tiles[[k]]

    if(length(fragik) == 0){

      tilesk$reads <- 0

    }else{

      #N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      if (nTiles%%1 != 0) {
          nTiles <- trunc(nTiles) + 1
      }

      #Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$RG, cellGroupi)
      
      mat <- Matrix::sparseMatrix(
          i = c(trunc(start(fragik) / tileSize), trunc(end(fragik) / tileSize)) + 1,
          j = as.vector(c(matchID, matchID)),
          x = rep(1,  2*length(fragik)),
          dims = c(nTiles, length(cellGroupi))
        )

      if(!is.null(ceiling)){
        mat@x[mat@x > ceiling] <- ceiling
      }
      
      mat <- Matrix::rowSums(mat)
      
      rm(fragik, matchID)
         
      tilesk$reads <- mat

      if(tolower(normMethod) %in% c("ncells")){
        tilesk$reads <- tilesk$reads / length(cellGroupi)
      }else if(tolower(normMethod) %in% c("none")){
      }else{
        stop("NormMethod not recognized!")
      }

    }

    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]

    tilesk

  }, threads = threads)

  names(covList) <- availableChr

  covList <- as(covList, "RleList")

  rtracklayer::export.bw(object = covList, con = covFile)

  return(covFile)

}
                        
.suppressAll <- function(expr = NULL){
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}                   
                        
getConditionBigWig <- function(ArchRProj, myCondition, maxCells=2000){

    ArrowFiles <- getArrowFiles(ArchRProj)
    
    Groups <- getCellColData(ArchRProj = ArchRProj, 
                         select = "myCond", 
                         drop = TRUE)
    
    Cells <- ArchRProj$cellNames
    
    cellGroups <- split(Cells, Groups)
    
    
    if(!is.null(maxCells)){

        gnames <- names(cellGroups)

        cellGroups <- lapply(seq_along(cellGroups), function(x){
          if(length(cellGroups[[x]]) > maxCells){
            sample(cellGroups[[x]], maxCells)
          }else{
            cellGroups[[x]]
          }
        })

        names(cellGroups) <- gnames

      }


    bwDir1 <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs")
    bwDir2 <- file.path(getOutputDirectory(ArchRProj), "GroupBigWigs", myCondition)
    
    
    dir.create(bwDir2, showWarnings = FALSE)

    
    o <- suppressWarnings(file.remove(list.files(bwDir2, full.names = TRUE)))

    cellsInArrow <- split(
    rownames(getCellColData(ArchRProj)), 
    stringr::str_split(rownames(getCellColData(ArchRProj)), 
                       pattern="\\#", 
                       simplify=TRUE)[,1])
    
    availableChr <- availableSeqnames(head(getArrowFiles(ArchRProj)))
    
    chromLengths <- getChromLengths(ArchRProj)
    chromSizes <- getChromSizes(ArchRProj)
    
    
    allChromTiles <- readRDS("/data/Osaka_Omics/scATAC/OsakaMomicsATAC/myAllTiles.rds")

    createMyGroupBW(
      i = 1, 
      cellGroups,
      ArrowFiles, 
      cellsInArrow, 
      availableChr,
      chromLengths, 
      tiles = allChromTiles,
      ceiling = 40,
      tileSize = 1,
      normMethod = "none",
      bwDir = bwDir2,
      threads = 90
      )

}
