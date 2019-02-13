#----------------------------------------------------------------------------
#transform metabolite information to quantitative module information
setGeneric(name = "metTransform",
           def = function(annotation.result = "annotation.result.csv",
                          group,
                          control.group = "W30",
                          threads = 3,
                          trans.to = c("module", "pathway"),
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          path = "."){
             #
trans.to <- match.arg(trans.to)
species <- match.arg(species)
# file.check
if(missing(group)) {
stop("You must give the names of group which you want to process.")}
# suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
file <- dir(path)
if(all(file != annotation.result)){
  stop("annotation.result should be in the work directory")
}

if(all(file != "tags.result2")){
  stop("tags.result2 from metABM should be in the work directory")
}

if(all(file != "tags2.after.annotation")){
  stop("tags2.after.annotation from metABM should be in the work directory")
}

if(all(file != "sample.info.csv")){
  stop("sample.info (column 1: sample.name and columnd 2: group) should be in the work directory")
}

if(trans.to == "module"){
if(all(file != "module.result")){
  stop("module.result from metModule should be in the work directory")
}
}

sample.info <-
  readr::read_csv(
    file.path(path, "sample.info.csv"),
    progress = FALSE,
    col_types = readr::cols()
  )

sample.info <- as.data.frame(sample.info)

if(any(!is.element(group, unique(sample.info[,2])))){
  idx <- which(!is.element(group, unique(sample.info[,2])))
  stop(paste(group[idx], collapse = " and "), ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
}

sample.info <- sample.info[sample.info[,2] %in% group,]

annotation.result <-
  readr::read_csv(file.path(path, annotation.result),
                  progress = FALSE,
                  col_types = readr::cols())

sample <- annotation.result[,match(sample.info[,1], colnames(annotation.result))]

sample <- as.data.frame(sample)

load(file.path(path, "tags.result2"))
load(file.path(path, "tags2.after.annotation"))

  sample <- sample[, order(colnames(sample))]
  sample.info <- sample.info[order(sample.info[, 1]), ]

  as.normal <- rep(FALSE, ncol(sample))
  as.normal[grep(control.group, sample.info[, 2])] <- TRUE

  if(trans.to == "module"){
  load(file.path(path, "module.result"))

  gs <- module.result$All.metabolite.id
  gs <-
    lapply(gs, function(x)
      strsplit(x, split = ";")[[1]])
  pathwaynames <- as.list(module.result$Module.name)
  pathway.data <- list(gs, pathwaynames)
  names(pathway.data) <- c("gs", "pathwaynames")
  }

  if(trans.to == "pathway"){

    switch(species,
           "dme" = {data("dme.kegg.pathway", envir = environment())
             met.pathway <- dme.kegg.pathway
             rm(dme.kegg.pathway)},
           "hsa" = {data("hsa.kegg.pathway", envir = environment())
             met.pathway <- hsa.kegg.pathway
             rm(hsa.kegg.pathway)},
           "mmu" = {data("mmu.kegg.pathway", envir = environment())
             met.pathway <- mmu.kegg.pathway
             rm(mmu.kegg.pathway)},
           "rat" = {data("rat.kegg.pathway", envir = environment())
             met.pathway <- rat.kegg.pathway
             rm(rat.kegg.pathway)},
           "bta" = {data("bta.kegg.pathway", envir = environment())
             met.pathway <- bta.kegg.pathway
             rm(bta.kegg.pathway)},
           "gga" = {data("gga.kegg.pathway", envir = environment())
             met.pathway <- gga.kegg.pathway
             rm(gga.kegg.pathway)},
           "dre" = {data("dre.kegg.pathway", envir = environment())
             met.pathway <- dre.kegg.pathway
             rm(dre.kegg.pathway)},
           "cel" = {data("cel.kegg.pathway", envir = environment())
             met.pathway <- cel.kegg.pathway
             rm(cel.kegg.pathway)},
           "sce" = {data("sce.kegg.pathway", envir = environment())
             met.pathway <- sce.kegg.pathway
             rm(sce.kegg.pathway)},
           "ath" = {data("ath.kegg.pathway", envir = environment())
             met.pathway <- ath.kegg.pathway
             rm(ath.kegg.pathway)},
           "smm" = {data("smm.kegg.pathway", envir = environment())
             met.pathway <- smm.kegg.pathway
             rm(smm.kegg.pathway)},
           "pfa" = {data("pfa.kegg.pathway", envir = environment())
             met.pathway <- pfa.kegg.pathway
             rm(pfa.kegg.pathway)},
           "tbr" = {data("tbr.kegg.pathway", envir = environment())
             met.pathway <- tbr.kegg.pathway
             rm(tbr.kegg.pathway)},
           "eco" = {data("eco.kegg.pathway", envir = environment())
             met.pathway <- eco.kegg.pathway
             rm(eco.kegg.pathway)},
           "ppu" = {data("ppu.kegg.pathway", envir = environment())
             met.pathway <- ppu.kegg.pathway
             rm(ppu.kegg.pathway)},
           "syf" = {data("syf.kegg.pathway", envir = environment())
             met.pathway <- syf.kegg.pathway
             rm(syf.kegg.pathway)})









    gs <- met.pathway
    gs <- lapply(gs, function(x){
      as.matrix(x, ncol = 1)
    })
    pathwaynames <- as.list(names(met.pathway))
    pathway.data <- list(gs, pathwaynames)
    names(pathway.data) <- c("gs", "pathwaynames")
  }

  peak.name <- showTags2(tags2.after.annotation, slot = "name")
  rownames(sample) <- peak.name
  cat("\n")
  if(trans.to == "module"){
    cat("Transform metabolite information to quantitative subnetwork information\n")
  }else{
    cat("Transform metabolite information to quantitative pathway information\n")
  }

  module.score <- metabolite2module(
    sample = sample,
    tags = tags.result2,
    pathway.data = pathway.data,
    as.normal = as.normal,
    path = path,
    threads = threads
  )

  rm(module.score)
  cat("\n")
  cat("metTransform is done.\n")
})













