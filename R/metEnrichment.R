#' @title metEnrichment
#' @description Metabolic enrichment analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param annotation.table The metabolite annotation result which used for
#' pathway enrichment analysis. The data frame should contain two column,
#' the column names are peak.name and KEGG.ID.
#' @param sample.pos Positive annotation.result from ms2Annotation.
#' @param sample.neg Negative annotation.result from ms2Annotation.
#' @param group The group you want to use.
#' @param column Column, "hilic" or "rp".
#' @param candidate.num Top candidate.num metabolites for each peak.
#' @param sample.info Sample information.
#' @param polarity The polarity mode of data. "positive", "negative" or "both".
#' @param pos.path Directory of positive results.
#' @param neg.path Directory of negative results.
#' @param output.path The directory to output results.
#' @param species The species. "hsa" is Homo sapiens (human),
#' "dme" is Drosophlia melanogaster (fruit fly),
#' "mmu" is Mus musculus (mouse), "rat" is Rattus norvegicus (rat),
#' "bta" is Bos taurus (cow), "gga" is Gallus domesticus (chicken),
#' "dre" is Danio rerio (zebrafish), "cel" is Caenorharomyces elegans (nematode),
#'  "sce" is Saccharomyces cerevisaiae (yeast), "ath" is Arabidopsis thaliana (thale cress),
#'  "smm" is Schistosoma mansoni, "pfa" is Plasmodum falciparum 3D7,
#'  "tbr" is Trypanosoma brucei, "eco" is Escherichia coli K-12 MG1655,
#' "ppu" is Pseudomonas putida KT2440, "syf" is Synechococcus elongatus.
#' @export


# load("marker")
# load("sample.pos")
# load("sample.neg")
# load("group")
# load("sample.info")
# load("species")
# load("pos.path")
# load("neg.path")
# load("output.path")
#
# metEnrichment1(annotation.table = marker,
#                sample.pos = sample.pos,
#                sample.neg = sample.neg,
#                group = group,
#                column = column,
#                sample.info = sample.info,
#                polarity = "both",
#                species = species,
#                pos.path = pos.path,
#                neg.path = neg.path,
#                output.path = output.path)


setGeneric(
  name = "metEnrichment",
  def = function(annotation.table,
                 sample.pos,
                 sample.neg,
                 group,
                 column = c("hilic", "rp"),
                 candidate.num = 5,
                 sample.info = sample.info,
                 polarity = c("positive", "negative", "both"),
                 species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                             "dre", "cel", "sce", "ath", "smm", "pfa",
                             "tbr", "eco", "ppu", "syf"),
                 pos.path = pos.path,
                 neg.path = neg.path,
                 output.path = output.path) {



    polarity <- match.arg(polarity)
    column <- match.arg(column)
    species <- match.arg(species)


    ###########################################################################
    #creat folder
    dir.create(file.path(output.path, "intermediate_data"))
    dir.create(file.path(output.path, "Quantitative_information"))

    ##pathway enrichment analysis
    id <- as.character(annotation.table[,2])
    id <- id[!is.na(id)]

    if(length(id) == 0){
      stop("There is no markers.")
    }

    id <- unique(unlist(strsplit(id, split = ";")))

    enrichment.result <- mseAnalysis(metabolite.id = id, species = species)

    if(is.null(enrichment.result)) {return(NULL)}

    write.csv(enrichment.result, file.path(output.path, "Pathway.enrichment.analysis.csv"))

    ####overview
    msea <- new(Class = 'moduleInfo',
                Module.name = "Pathway enrichment analysis",
                Module.size = length(id),
                Module.impact = 0,
                p.value = 0,
                Detected.number = 0,
                Hidden.number = 0,
                All.metabolite.id = id,
                Detected.ID = id,
                Hidden.ID = id,
                All.metabolite.name = "no",
                Detected.metabolite.name = "no",
                Hidden.metabolite.name = "no",
                msea = enrichment.result
    )

    save(msea, file = file.path(output.path, "intermediate_data", "msea"), compress = "xz")

    # save(dn.msea, file = file.path(output.path, "intermediate_data", "dn.msea"), compress = "xz")

    pdf(file.path(output.path, "Pathway.enrichment.MSEA.pdf"),
        width = 7,
        height = 7)
    moduleplot(object = msea, cex.label = 1, n = 20)
    dev.off()

    ##
    pdf(file.path(output.path, "Pathway.enrichment.overview.pdf"),
        width = 7, height = 7)
    par(mar = c(5, 5, 4, 2))
    pathwayPlot(mse.data = msea@msea)
    dev.off()

    ##quantitative analysis of all pathways
    temp.path1 <-
      file.path(output.path, "Boxplot")
    temp.path2 <-
      file.path(output.path, "Heatmap")

    dir.create(path = temp.path1)
    dir.create(path = temp.path2)

    ###quantitative analysis
    if(polarity == "positive" | polarity == "both"){
      load(file.path(pos.path, "MRN_annotation_result",
                     "intermediate_data", "tags.result"))
      tags.result.pos <- tags.result
      if(polarity == "both"){
        tags.result.pos$name <- paste(tags.result.pos$name, "POS", sep = "_")
      }

      ###only select top candidate.num metabolite for peaks
      tags.result.pos <- lapply(unique(tags.result.pos$name), function(x){
        temp.idx <- which(tags.result.pos$name == x)
        temp.data <- tags.result.pos[temp.idx,,drop = FALSE]
        if(nrow(temp.data) > candidate.num) temp.data <- temp.data[1:candidate.num,]
        temp.data
      })

      tags.result.pos <- do.call(rbind, tags.result.pos)

    }



    if(polarity == "negative" | polarity == "both"){
      load(file.path(neg.path, "MRN_annotation_result",
                     "intermediate_data", "tags.result"))
      tags.result.neg <- tags.result
      if(polarity == "both"){
        tags.result.neg$name <- paste(tags.result.neg$name, "NEG", sep = "_")
      }

      ###only select top candidate.num metabolite for peaks
      tags.result.neg <- lapply(unique(tags.result.neg$name), function(x){
        temp.idx <- which(tags.result.neg$name == x)
        temp.data <- tags.result.neg[temp.idx,,drop = FALSE]
        if(nrow(temp.data) > candidate.num) temp.data <- temp.data[1:candidate.num,]
        temp.data
      })

      tags.result.neg <- do.call(rbind, tags.result.neg)

    }


    switch(polarity,
           "positive" = {tags.result <- tags.result.pos},
           "negative" = {tags.result <- tags.result.neg},
           "both" = {
             tags.result <- rbind(tags.result.pos,
                                  tags.result.neg)})


    anno <- tags.result
    anno <- anno[anno$isotope == "[M]", , drop = FALSE]

    anno <- anno[which(anno$name %in% as.character(annotation.table[,1])),]
    ##if one metabolite matches multiple peaks, then the dysregulated peak is used.


    ###get group.data
    group.data <- getPathway(species = species, type = "metabolite")

    ## remove the pathway which has less than 3 metabolite
    remove.idx <- unname(which(unlist(lapply(group.data, function(x){
      sum(x %in% unique(anno$to)) < 3
    }))))

    if(length(remove.idx) > 0){
      group.data <- group.data[-remove.idx]
    }

    data("kegg.rpair2", envir = environment())
    if(length(group.data) != 0){

      ##transform peak to pathway quantitative information

      singleTransform(annotation.result.pos = sample.pos,
                      annotation.result.neg = sample.neg,
                      anno = anno,
                      polarity = polarity,
                      group.data = group.data,
                      data.type = "metabolomics",
                      sample.info = sample.info,
                      group = group,
                      trans.to = "pathway",
                      scale = TRUE,
                      scale.method = "pareto",
                      species = species,
                      which.peak = "intensity",
                      column = column,
                      method = "mean",
                      output.path = file.path(output.path, "Quantitative_information"),
                      metabolic.network = kegg.rpair2)


      file.rename(from = file.path(output.path,
                                   "Quantitative_information/DNA.pathway.quantitative.result.csv"),
                  to = file.path(output.path,
                                 "Quantitative_information/pathway.quantitative.result.csv"))

      #box plot and heatmap for module
      pathway.quantitative.result <-
        readr::read_csv(file.path(output.path, "Quantitative_information",
                                  "pathway.quantitative.result.csv"),
                        progress = FALSE, col_types = readr::cols())

      pathway.quantitative.result <- as.data.frame(pathway.quantitative.result)

      pathway.node.quantitative.result <-
        readr::read_csv(file.path(output.path, "Quantitative_information",
                                  "pathway.node.quantitative.result.csv"),
                        progress = FALSE, col_types = readr::cols())

      pathway.node.quantitative.result <- as.data.frame(pathway.node.quantitative.result)

      pathway.sample <- pathway.quantitative.result[,-1]
      pathway.name.id <- pathway.quantitative.result[,1]
      pathway.name <- unlist(lapply(strsplit(pathway.name.id, split = ";"), function(x) x[1]))
      pathway.id <- unlist(lapply(strsplit(pathway.name.id, split = ";"), function(x) x[2]))
      pathway.name <- gsub(pattern = "[,|(|)|/]", replacement = "", pathway.name)
      rownames(pathway.sample) <- pathway.name



      cat("\n")
      cat("Calculate p values of pathways.\n")

      pathway.p.value <- uniTest(sample = pathway.sample,
                                 sample.info = sample.info,
                                 uni.test = "t",
                                 correct = TRUE)

      cat("\n")

      cat("Calculate fold changes of pathways.\n")
      pathway.fc <- foldChange(sample = pathway.sample,
                               sample.info = sample.info,
                               group = group)


      samplePlot(sample = pathway.sample,
                 sample.info = sample.info,
                 group = group,
                 output.path = temp.path1,
                 p.value = pathway.p.value,
                 fc = pathway.fc,
                 beeswarm = TRUE)

      pathway.sample1 <- t(apply(pathway.sample, 1, as.numeric))
      colnames(pathway.sample1) <- colnames(pathway.sample)
      rownames(pathway.sample1) <- rownames(pathway.sample)


      if(nrow(pathway.sample) > 1){
        pdf(file.path(output.path, "pathway.heatmap.pdf"),
            width = 7, height = 7)
        par(mar = c(5, 5 ,4, 2))
        heatMap(sample = pathway.sample1, sample.info = sample.info, group = group)
        dev.off()
      }

      ##heatmap for each pahtway
      group.data <- group.data[match(pathway.name.id, names(group.data))]

      for(i in 1:length(group.data)){
        # cat(i); cat(" ")
        temp.pathway <- group.data[[i]]
        temp.pathway.name <- pathway.name[i]
        temp.idx <- match(temp.pathway, pathway.node.quantitative.result$ID)
        temp.idx <- temp.idx[!is.na(temp.idx)]
        if(length(temp.idx) < 3) next()
        temp.sample <- pathway.node.quantitative.result[temp.idx,,drop = FALSE]

        rownames(temp.sample) <- temp.sample$compound.name

        pdf(file.path(temp.path2, paste(temp.pathway.name, "pdf",sep = ".")),
            width = 7, height = 7)
        par(mar = c(5, 5 ,4, 2))
        heatMap(sample = temp.sample,
                sample.info = sample.info,
                group = group)
        dev.off()

      }

    }

    cat("\n")
    cat("metEnrichment is done\n")
  })




