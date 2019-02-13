#' @title metModule
#' @description Metabolic module analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param annotation.result.pos Positive annotation.result from ms2Annotation.
#' @param annotation.result.neg Negative annotation.result from ms2Annotation.
#' @param tags2.pos Tags2 result of positive data.
#' @param tags2.neg Tags2 result of negative data.
#' @param p.value.pos Positive p values.
#' @param p.value.neg Negative p values.
#' @param foldchange.pos Positive fold change.
#' @param foldchange.neg Negative fold change.
#' @param rt.result Predicted RT results from rtPrediction.
#' @param pos.path Directory of positive results.
#' @param neg.path Directory of negative results.
#' @param sample.info Sample information.
#' @param polarity "positive", "negative" or "both".
#' @param group The group you want to use.
#' @param max.isotope The max number of isotope.
#' @param uni.test The method of univariate test.
#' @param column hilic or rp.
#' @param correct Correct p value or not.
#' @param p.cutoff The cutoff of p value.
#' @param mz.tol The mz tolerance of KEGG matching.
#' @param rt.tol1 The RT tolerance of isotope annotation (second).
#' @param rt.tol2 The RT tolerance of KEGG matching (\%).
#' @param cor.tol The tolerance of correlation.
#' @param int.tol The tolerance of intensity ratio (\%).
#' @param threads The number of threads.
#' @param output.path The directory to output results.
#' @param use.old.null Use old null distribution of not.
#' @param species The species. "hsa" is Homo sapiens (human),
#' "dme" is Drosophlia melanogaster (fruit fly),
#' "mmu" is Mus musculus (mouse), "rat" is Rattus norvegicus (rat),
#' "bta" is Bos taurus (cow), "gga" is Gallus domesticus (chicken),
#' "dre" is Danio rerio (zebrafish), "cel" is Caenorharomyces elegans (nematode),
#'  "sce" is Saccharomyces cerevisaiae (yeast), "ath" is Arabidopsis thaliana (thale cress),
#'  "smm" is Schistosoma mansoni, "pfa" is Plasmodum falciparum 3D7,
#'  "tbr" is Trypanosoma brucei, "eco" is Escherichia coli K-12 MG1655,
#' "ppu" is Pseudomonas putida KT2440, "syf" is Synechococcus elongatus.
#' @param use.all.kegg.id Use all annotations from KEGG database.
#' Default is FALSE.
#' @param only.mrn.annotation Only use MEN annotations or not.
#' @param output.module.information Output module information ot not.
#' @export

# system.time(metModule2(group = c("Normal", "AD"), column = "rp",
#                       correct = FALSE, p.cutoff = 0.05, species = "hsa",
#                       use.old.null = TRUE, only.mrn.annotation = FALSE))



# annotation.result.pos = sample.pos
# annotation.result.neg = sample.neg
# tags2.pos = tags2.pos
# tags2.neg = tags2.neg
# p.value.pos = p.value.pos
# p.value.neg = p.value.neg
# foldchange.pos = foldchange.pos
# foldchange.neg = foldchange.neg
# rt.result = rt.result
# pos.path = pos.path
# neg.path = neg.path
# sample.info = sample.info
# polarity = polarity
# group = group
# max.isotope = max.isotope
# uni.test = uni.test
# column = column
# output.path = output.path
# correct = correct
# p.cutoff = p.cutoff
# mz.tol = mz.tol
# rt.tol1 = rt.tol1# absolute second
# rt.tol2 = rt.tol2#relative %
# cor.tol = cor.tol
# int.tol = int.tol
# threads = threads
# use.old.null = use.old.null1
# species = species
# use.all.kegg.id = use.all.kegg.id
# only.mrn.annotation = only.mrn.annotation
# output.module.information = output.module.information



setGeneric(
  name = "metModule",
  def = function(annotation.result.pos,
                 annotation.result.neg,
                 tags2.pos = tags2.pos,
                 tags2.neg = tags2.neg,
                 p.value.pos = p.value.pos,
                 p.value.neg = p.value.neg,
                 foldchange.pos = foldchange.pos,
                 foldchange.neg = foldchange.neg,
                 rt.result = rt.result,
                 pos.path = pos.path,
                 neg.path = neg.path,
                 sample.info = sample.info,
                 polarity = c("positive", "negative", "both"),
                 group,
                 max.isotope = 4,
                 uni.test = c("t", "wilcox", "anova"),
                 column = c("hilic", "rp"),
                 output.path = ".",
                 correct = TRUE,
                 p.cutoff = 0.01,
                 mz.tol = 25,
                 rt.tol1 = 3,# absolute, second
                 rt.tol2 = 30,#relative, %
                 cor.tol = 0,
                 int.tol = 500,
                 threads = 3,
                 use.old.null = FALSE,
                 species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                             "dre", "cel", "sce", "ath", "smm", "pfa",
                             "tbr", "eco", "ppu", "syf"),
                 use.all.kegg.id = FALSE,
                 only.mrn.annotation = FALSE,
                 output.module.information = TRUE) {

    uni.test <- match.arg(uni.test)
    polarity <- match.arg(polarity)
    column <- match.arg(column)
    species <- match.arg(species)

    ###########################################################################
    #creat folder
    dir.create(file.path(output.path, "intermediate_data"))
    if(output.module.information) dir.create(file.path(output.path, "DNA_module_information"))
    dir.create(file.path(output.path, "DNA_functional_annotation"))
    dir.create(file.path(output.path, "DNA_functional_annotation", "Quantitative_information"))


    ###########################################################################
    #check sample names of positive and negative
    if(polarity == "both"){
      ###check sample name in sample.pos and sample.neg
      if(any(sort(colnames(annotation.result.pos)) != sort(colnames(annotation.result.neg)))){
        cat("Sample names in annotation.result.pos and annotation.result.neg are not same.\n")
      }
    }


    ###########################################################################
    #get sample.pos and sample.neg
    sample.pos <- annotation.result.pos
    sample.neg <- annotation.result.neg

    ###########################################################################
    ###check for ref.as
    if(use.old.null) {
      file <- dir(file.path(output.path, "intermediate_data"))
      need.file <- "ref.as"
      file.check <- which(!need.file %in% file)
      if (length(file.check) > 0) {
        stop(paste(need.file[file.check], collapse = " & "),
             " is not in the output directory.")
      }
    }

    ###########################################################################
    #add predicted RTs to KEGG database
    kegg.rt <- rt.result[[2]]
    rm(rt.result)
    gc()

    ##load adduct table informatio
    data("kegg.adduct.pos", envir = environment())
    data("kegg.adduct.neg", envir = environment())

    if (column == "hilic") {
      kegg.adduct.pos <- kegg.adduct.pos[which(kegg.adduct.pos$Column.hilic), ]
      kegg.adduct.neg <- kegg.adduct.neg[which(kegg.adduct.neg$Column.hilic), ]
    } else{
      kegg.adduct.pos <- kegg.adduct.pos[which(kegg.adduct.pos$Column.rp), ]
      kegg.adduct.neg <- kegg.adduct.neg[which(kegg.adduct.neg$Column.rp), ]
    }

    ##Add RT information to kegg.adduct
    idx1 <- match(kegg.adduct.pos$ID, row.names(kegg.rt))
    kegg.adduct.pos <- data.frame(kegg.adduct.pos, kegg.rt[idx1, ],
                                  stringsAsFactors = FALSE)
    kegg.adduct.pos$RT[is.na(kegg.adduct.pos$RT)] <-
      median(kegg.adduct.pos$RT[!is.na(kegg.adduct.pos$RT)])

    idx2 <- match(kegg.adduct.neg$ID, row.names(kegg.rt))
    kegg.adduct.neg <- data.frame(kegg.adduct.neg, kegg.rt[idx2, ],
                                  stringsAsFactors = FALSE)
    kegg.adduct.neg$RT[is.na(kegg.adduct.neg$RT)] <-
      median(kegg.adduct.neg$RT[!is.na(kegg.adduct.neg$RT)])


    ###########################################################################
    ##identify initial modules
    ##prepare data for module analysis
    data("kegg.rpair2", envir = environment())
    data("kegg.compound", envir = environment())

    if (column == "hilic") {
      data("adduct.table.hilic", envir = environment())
      adduct.table <- adduct.table.hilic
      rm(adduct.table.hilic)
      gc()
    }

    if (column == "rp") {
      data("adduct.table.rp", envir = environment())
      adduct.table <- adduct.table.rp
      rm(adduct.table.rp)
      gc()
    }

    adduct.table.pos <- adduct.table[adduct.table$mode == "positive", ]
    adduct.table.neg <- adduct.table[adduct.table$mode == "negative", ]

    rm(list = "adduct.table")
    gc()


    return.result <-
      getSDM(
        tags2.pos = tags2.pos,
        tags2.neg = tags2.neg,
        polarity = polarity,
        metabolic.network = kegg.rpair2,
        metabolite.pos = kegg.adduct.pos,
        metabolite.neg = kegg.adduct.neg,
        kegg.compound = kegg.compound,
        p.value.pos = p.value.pos,
        p.value.neg = p.value.neg,
        p.cutoff = p.cutoff,
        adduct.table.pos = adduct.table.pos,
        adduct.table.neg = adduct.table.neg,
        mz.tol = mz.tol,
        rt.tol = rt.tol2,
        threads = threads,
        output.path = output.path,
        use.old.null = use.old.null,
        species = species,
        use.all.kegg.id = use.all.kegg.id,
        only.mrn.annotation = only.mrn.annotation)

    ##get the result of moudle analysis
    if(polarity == "positive" | polarity == "both"){
      return.result.pos <- return.result[[1]]
      new.tags.pos <- return.result.pos[[1]]
      anno.pos <- return.result.pos[[2]]
      anno.pos.from.module <- anno.pos
      save(anno.pos.from.module,
           file = file.path(output.path, "intermediate_data", "anno.pos.from.module"), compress = "xz")

      rm(list = c("return.result.pos",
                  "anno.pos",
                  "anno.pos.from.module",
                  "kegg.adduct.pos",
                  "kegg.compound",
                  "adduct.table.pos"))

      ##process new.tags.pos and new.tags.neg, for the annotation from KEGG
      if(is.null(new.tags.pos)){
        tags2.kegg.matching.pos.from.module <- tags2.pos
        rm(list = c("tags2.pos"))
        gc()
        load(file.path(pos.path,
                       "MRN_annotation_result",
                       "intermediate_data","tags.result"))
        tags.result.pos2.from.module <- tags.result
        if(polarity == "both"){
          tags.result.pos2.from.module$name <-
            paste(tags.result.pos2.from.module$name, "POS", sep = "_")
        }

        rm(list = c("tags.result"))
        gc()
        save(tags2.kegg.matching.pos.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags2.kegg.matching.pos.from.module"),
             compress = "xz")
        save(tags.result.pos2.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags.result.pos2.from.module"), compress = "xz")
      }else{
        cat('\n')
        cat("Process positive new.tags.\n")
        temp.result.pos <- processTags2(
          new.tags2 = new.tags.pos,
          sample = sample.pos,
          mz.tol = mz.tol,
          rt.tol = rt.tol1,
          cor.tol = cor.tol,
          int.tol = int.tol,
          max.isotope = max.isotope,
          polarity = "positive",
          path = output.path
        )

        new.tags.pos <- temp.result.pos[[1]]
        tags.result.pos2.from.module <- temp.result.pos[[2]]

        tags2.kegg.matching.pos.from.module <- new.tags.pos
        save(tags2.kegg.matching.pos.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags2.kegg.matching.pos.from.module"),
             compress = "xz")
        save(tags.result.pos2.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags.result.pos2.from.module"), compress = "xz")
      }

      ##change tags2.after.kegg.matching.pos to csv like ms2Annotation
      cat("\n")
      data("kegg.compound", envir = environment())

      temp.pos <- getAnnotationResult(tags2 = tags2.kegg.matching.pos.from.module,
                                      p.value = p.value.pos,
                                      correct = correct,
                                      foldchange = foldchange.pos,
                                      tags.result2 = tags.result.pos2.from.module,
                                      kegg.compound = kegg.compound)

      readr::write_csv(temp.pos,
                       file.path(output.path,
                                 "DNA.module.annotation.result.pos.csv"))
      rm(list = "temp.pos")
      gc()


      rm(list = c("tags2.kegg.matching.pos.from.module",
                  "temp.result.pos", "new.tags.pos",
                  "tags.result.pos2.from.module"))
      gc()
    }


    if(polarity == "negative" | polarity == "both"){
      return.result.neg <- return.result[[2]]
      new.tags.neg <- return.result.neg[[1]]
      anno.neg <- return.result.neg[[2]]
      anno.neg.from.module <- anno.neg
      save(anno.neg.from.module,
           file = file.path(output.path,
                            "intermediate_data",
                            "anno.neg.from.module"), compress = "xz")

      rm(list = c("return.result.neg",
                  "anno.neg",
                  "anno.neg.from.module",
                  "kegg.adduct.neg",
                  "adduct.table.neg"))

      ##negative
      if(is.null(new.tags.neg)){
        tags2.kegg.matching.neg.from.module <- tags2.neg
        rm(list = c("tags2.neg"))
        gc()
        load(file.path(neg.path,
                       "MRN_annotation_result",
                       "intermediate_data","tags.result"))
        tags.result.neg2.from.module <- tags.result
        if(polarity == "both"){
          tags.result.neg2.from.module$name <-
            paste(tags.result.neg2.from.module$name, "NEG", sep = "_")
        }

        rm(list = c("tags.result"))
        gc()
        save(tags2.kegg.matching.neg.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags2.kegg.matching.neg.from.module"),
             compress = "xz")
        save(tags.result.neg2.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags.result.neg2.from.module"), compress = "xz")
      }else{
        cat('\n')
        cat("Process negative new.tags.\n")
        temp.result.neg <- processTags2(
          new.tags2 = new.tags.neg,
          sample = sample.neg,
          mz.tol = mz.tol,
          rt.tol = rt.tol1,
          cor.tol = cor.tol,
          int.tol = int.tol,
          max.isotope = max.isotope,
          polarity = "negative",
          path = output.path
        )

        new.tags.neg <- temp.result.neg[[1]]
        tags.result.neg2.from.module <- temp.result.neg[[2]]

        tags2.kegg.matching.neg.from.module <- new.tags.neg
        save(tags2.kegg.matching.neg.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags2.kegg.matching.neg.from.module"),
             compress = "xz")
        save(tags.result.neg2.from.module,
             file = file.path(output.path,
                              "intermediate_data",
                              "tags.result.neg2.from.module"), compress = "xz")
      }


      ##change tags2.after.kegg.matching.neg to csv like ms2Annotation
      cat("\n")
      data("kegg.compound", envir = environment())

      temp.neg <- getAnnotationResult(tags2 = tags2.kegg.matching.neg.from.module,
                                      p.value = p.value.neg,
                                      correct = correct,
                                      foldchange = foldchange.neg,
                                      tags.result2 = tags.result.neg2.from.module,
                                      kegg.compound = kegg.compound)

      cat("\n")
      readr::write_csv(temp.neg, file.path(output.path,
                                           "DNA.module.annotation.result.neg.csv"))
      rm(list = "temp.neg")
      gc()


      rm(list = c("tags2.after.kegg.matching.neg.from.module",
                  "temp.result.neg", "new.tags.neg",
                  "tags.result.neg2.from.module"))
      gc()

    }



    ##--------------------------------------------------------------------------
    #output the metabolite sets enrichment results and quantitative information
    load(file.path(output.path, "intermediate_data", "module.result2"))
    load(file.path(output.path, "intermediate_data", "module.result"))

    module.p <- unlist(lapply(module.result2, function(x)x@p.value))
    temp.index <- which(module.p < 0.05)

    if (length(temp.index) > 0) {
      cat("\n")
      cat("There are", length(temp.index), "dysregulated modules out of", length(module.p), "modules.\n")

      if(output.module.information){
        dir.create(file.path(output.path, "DNA_module_information"))
        dir.create(file.path(output.path, "DNA_module_information", "Quantitative_information"))
        ##construct group.data for singleTranform
        group.data <- module.result$All.metabolite.id[temp.index]
        group.data <- lapply(group.data, function(x){
          strsplit(x, split = ";")[[1]]
        })
        names(group.data) <- module.result$Module.name[temp.index]


        ###transform metabolite information to pathway/module information
        if(polarity == "positive" | polarity == "both"){
          load(file.path(output.path, "intermediate_data", "anno.pos.from.module"))
          load(file.path(output.path, "intermediate_data", "tags.result.pos2.from.module"))
        }

        if(polarity == "negative" | polarity == "both"){
          load(file.path(output.path, "intermediate_data", "anno.neg.from.module"))
          load(file.path(output.path, "intermediate_data", "tags.result.neg2.from.module"))
        }

        switch(polarity,
               "positive" = {anno <- anno.pos.from.module
               tags.result2 <- tags.result.pos2.from.module
               rm(list = c("anno.pos.from.module", "tags.result.pos2.from.module"))},
               "negative" = {anno <- anno.neg.from.module
               tags.result2 <- tags.result.neg2.from.module
               rm(list = c("anno.neg.from.module", "tags.result.neg2.from.module"))},
               "both" = {anno <- c(anno.pos.from.module, anno.neg.from.module)
               tags.result2 <- rbind(tags.result.pos2.from.module,
                                     tags.result.neg2.from.module)
               rm(list = c("anno.pos.from.module", "tags.result.pos2.from.module",
                           "anno.neg.from.module", "tags.result.neg2.from.module"))
               })

        tags.result2 <- tags.result2[which(tags.result2$name %in% names(anno)),]



        peak.id1 <- unname(unlist(mapply(function(x, y){
          paste(x, y, sep = ";")
        },
        x = names(anno),
        y = anno)))

        peak.id2 <- paste(tags.result2$name, tags.result2$to, sep = ";")


        anno <- tags.result2[which(peak.id2 %in% peak.id1),]


        singleTransform(annotation.result.pos = sample.pos,
                        annotation.result.neg = sample.neg,
                        anno = anno,
                        polarity = polarity,
                        group.data = group.data,
                        data.type = "metabolomics",
                        sample.info = sample.info,
                        group = group,
                        trans.to = "module",
                        scale = TRUE,
                        scale.method = "pareto",
                        species = species,
                        which.peak = "intensity",
                        column = column,
                        method = "mean",
                        output.path = file.path(output.path, "DNA_module_information", "Quantitative_information"),
                        metabolic.network = kegg.rpair2)



        ####output MSEA results of all modules
        temp.path <- file.path(output.path, "DNA_module_information")
        outputModuleInformation(module.result = module.result,
                                module.result2 = module.result2,
                                output.path = temp.path,
                                temp.index = temp.index)
      }

      ##get the dysregulated network
      all.id <- module.result$All.metabolite.id[temp.index]
      all.id <- unlist(lapply(all.id, function(x) {strsplit(x, ";")}))
      detected.id <- module.result$Detected.ID[temp.index]
      detected.id <- unlist(lapply(detected.id, function(x) {strsplit(x, ";")}))
      hidden.id <- module.result$Hidden.ID[temp.index]
      hidden.id <- unlist(lapply(hidden.id, function(x) {strsplit(x, ";")}))

      dn <- igraph::subgraph(graph = kegg.rpair2, v = all.id)

      ##pathway enrichment analysis
      # mse.result <- mseAnalysis(metabolite.id = all.id, species = species)
      mse.result <- mseAnalysis(metabolite.id = detected.id, species = species)
      write.csv(mse.result, file.path(output.path,
                                      "DNA_functional_annotation",
                                      "DNA.pathway.enrichment.result.csv"),
                row.names = TRUE)

      dn.msea <- new(Class = 'moduleInfo',
                     Module.name = "Dysregulate network",
                     Module.size = length(all.id),
                     Module.impact = 0,
                     p.value = 0,
                     Detected.number = 0,
                     Hidden.number = 0,
                     All.metabolite.id = all.id,
                     Detected.ID = detected.id,
                     Hidden.ID = hidden.id,
                     All.metabolite.name = "no",
                     Detected.metabolite.name = "no",
                     Hidden.metabolite.name = "no",
                     msea = mse.result
      )

      save(dn.msea, file = file.path(output.path, "intermediate_data", "dn.msea"), compress = "xz")

      pdf(file.path(output.path, "DNA_functional_annotation", "dysregulated.network.MSEA.pdf"),
          width = 7,
          height = 7)
      moduleplot(object = dn.msea, cex.label = 1, n = 20)
      dev.off()

      ##
      pdf(file.path(output.path, "DNA_functional_annotation", "dysregulate.network.overview.pdf"),
          width = 7, height = 7)
      par(mar = c(5, 5, 4, 2))
      pathwayPlot(mse.data = dn.msea@msea)
      dev.off()

      ##export file for cytoscape
      class <- rep(NA, length(all.id))
      i.idx <- match(detected.id, all.id)
      h.idx <- match(hidden.id, all.id)
      class[i.idx] <- "Detected"
      class[h.idx] <- "Hidden"

      DNA.node.quantitative.result <- readr::read_csv(file.path(output.path,
                                                                "DNA_module_information/Quantitative_information",
                                                                "DNA.node.quantitative.result.csv"),
                                                      progress = FALSE, col_types = readr::cols())

      idx1 <- match(all.id, DNA.node.quantitative.result$ID)

      idx <- idx1

      node.peak.name <- DNA.node.quantitative.result$peak.name[idx]

      if(polarity == "both") {
        p.value <- c(p.value.pos, p.value.neg)
        foldchange <- c(foldchange.pos, foldchange.neg)
        }

      if(polarity == "positive") {
        p.value <- p.value.pos
        foldchange <- foldchange.pos
      }

      if(polarity == "negative") {
        p.value <- p.value.neg
        foldchange <- foldchange.neg
      }

      node.p <- p.value[match(node.peak.name, names(p.value))]
      node.fc <- foldchange[match(node.peak.name, names(foldchange))]
      node.fc1 <- node.fc
      node.fc1[node.fc > 1] <- "Increase"
      node.fc1[node.fc < 1] <- "Decrease"
      node.fc1[node.fc == 1] <- "No change"
      node.fc <- node.fc1

      node.p[is.na(node.p)] <-  1
      node.fc[is.na(node.fc)] <- "ND"

      ##module information
      module <- lapply(module.result2, function(x){x@All.metabolite.id})
      names(module) <- unlist(lapply(module.result2, function(x) {x@Module.name}))

      module.group <- unlist(lapply(all.id, function(x){
        names(which(unlist(lapply(module, function(y) {is.element(x, y)}))))[1]
      }))


      ##pathway information
      pathway <- getPathway(species = species, type = "metabolite")

      pathway.group <- unlist(lapply(all.id, function(x){
        strsplit(names(which(unlist(lapply(pathway,
                                           function(y) {is.element(x, y)}))))[1], split = ";")[[1]][1]
      }))


      attribute <- data.frame(all.id, class, node.p, node.fc,
                              module.group, pathway.group,
                              stringsAsFactors = FALSE)
      colnames(attribute) <- c("ID", "Class", "P", "FC", "Module", "Pathway")

      log10P <- -log(attribute$P, 10)
      attribute <- data.frame(attribute, log10P, stringsAsFactors = FALSE)

      for.cyto <- igraph::as_edgelist(dn)
      for.cyto <- data.frame(for.cyto, stringsAsFactors = FALSE)
      colnames(for.cyto) <- c("From", "To")

      ###change ID to metabolite name
      data("kegg.compound", envir = environment())
      for.cyto[,1] <-
        unlist(lapply(strsplit(kegg.compound$Name[match(for.cyto[,1], kegg.compound$ID)], split = ";"), function(x){
          x[1]
        }))

      for.cyto[,2] <-
        unlist(lapply(strsplit(kegg.compound$Name[match(for.cyto[,2], kegg.compound$ID)], split = ";"), function(x){
          x[1]
        }))

      attribute[,1] <-
        unlist(lapply(strsplit(kegg.compound$Name[match(attribute[,1], kegg.compound$ID)], split = ";"), function(x){
          x[1]
        }))


      temp.path <- file.path(output.path, "DNA_functional_annotation", "Cytoscape_data")
      dir.create(temp.path)
      write.table(for.cyto, file.path(temp.path,
                                      "dysregulated.networks.cytoscape.txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)

      write.table(attribute, file.path(temp.path,
                                       "dysregulated.networks.attribute.txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)



      rm(list = c("all.id", "detected.id", "hidden.id", "dn", "mse.result",
                  "class", "i.idx", "h.idx", "attribute", "temp.graph",
                  "edge", "for.cyto"))
      rm(list = c("module", "pathway"))
      gc()


      #box plot and heatmap for module
      if(output.module.information) {
      temp.path1 <- file.path(output.path, "DNA_module_information", "Dysregulated_module_boxplot")
      temp.path2 <- file.path(output.path, "DNA_module_information", "Dysregulated_module_heatmap")

      dir.create(temp.path1)
      dir.create(temp.path2)

      DNA.module.quantitative.result <-
        readr::read_csv(file.path(output.path, "DNA_module_information/Quantitative_information",
                                  "DNA.module.quantitative.result.csv"),
                        progress = FALSE, col_types = readr::cols())

      DNA.module.quantitative.result <- as.data.frame(DNA.module.quantitative.result)

      DNA.node.quantitative.result <-
        readr::read_csv(file.path(output.path, "DNA_module_information/Quantitative_information",
                                  "DNA.node.quantitative.result.csv"),
                        progress = FALSE, col_types = readr::cols())

      DNA.node.quantitative.result <- as.data.frame(DNA.node.quantitative.result)

      module.name <- DNA.module.quantitative.result[,1]
      module.sample <- DNA.module.quantitative.result[,-1]
      rownames(module.sample) <- module.name

      cat("\n")
      cat("Calculate p values of modules.\n")
      module.p.value <- uniTest(sample = module.sample,
                                sample.info = sample.info,
                                uni.test = "t",
                                correct = TRUE)

      cat("\n")
      cat("Calculate fold changes of modules.\n")
      module.fc <- foldChange(sample = module.sample,
                              sample.info = sample.info,
                              group = group)


      samplePlot(sample = module.sample,
                 sample.info = sample.info,
                 group = group,
                 output.path = temp.path1,
                 p.value = module.p.value,
                 fc = module.fc,
                 beeswarm = TRUE)

      module.sample1 <- t(apply(module.sample, 1, as.numeric))
      colnames(module.sample1) <- colnames(module.sample)
      rownames(module.sample1) <- rownames(module.sample)


      pdf(file.path(temp.path2, "Dysregulated.modules.pdf"),
          width = 7, height = 7)
      par(mar = c(5, 5 ,4, 2))
      heatMap(sample = module.sample1, sample.info = sample.info, group = group)
      dev.off()

      ##heatmap for each module
      for(i in 1:length(group.data)){
        temp.module <- group.data[[i]]
        temp.module.name <- names(group.data)[i]
        temp.idx <- match(temp.module, DNA.node.quantitative.result$ID)
        temp.idx <- temp.idx[!is.na(temp.idx)]
        temp.sample <- DNA.node.quantitative.result[temp.idx,,drop = FALSE]

        rownames(temp.sample) <- temp.sample$compound.name

        pdf(file.path(temp.path2, paste("Dysregulated",temp.module.name, "pdf",sep = ".")),
            width = 7, height = 7)
        par(mar = c(5, 5 ,4, 2))
        heatMap(sample = temp.sample,
                sample.info = sample.info,
                group = group)
        dev.off()

      }
      }
    }else{
      cat("\n")
      cat("There are no dysregulated modules with p value less than 0.05.\n")
    }

      ##quantitative analysis of all pathways
      temp.path3 <- file.path(output.path, "DNA_functional_annotation", "Dysregulated_network_boxplot")
      temp.path4 <- file.path(output.path, "DNA_functional_annotation", "Dysregulated_network_heatmap")
      temp.path5 <- file.path(output.path, "DNA_functional_annotation", "Quantitative_information")

      dir.create(temp.path3)
      dir.create(temp.path4)
      dir.create(temp.path5)

      ###quantitative analysis
      if(length(temp.index) > 0){
        if(polarity == "positive" | polarity == "both"){
          load(file.path(output.path, "intermediate_data", "anno.pos.from.module"))
          load(file.path(output.path, "intermediate_data", "tags.result.pos2.from.module"))
        }

        if(polarity == "negative" | polarity == "both"){
          load(file.path(output.path, "intermediate_data", "anno.neg.from.module"))
          load(file.path(output.path, "intermediate_data", "tags.result.neg2.from.module"))
        }

        switch(polarity,
               "positive" = {anno <- anno.pos.from.module
               tags.result2 <- tags.result.pos2.from.module
               rm(list = c("anno.pos.from.module", "tags.result.pos2.from.module"))},
               "negative" = {anno <- anno.neg.from.module
               tags.result2 <- tags.result.neg2.from.module
               rm(list = c("anno.neg.from.module", "tags.result.neg2.from.module"))},
               "both" = {anno <- c(anno.pos.from.module, anno.neg.from.module)
               tags.result2 <- rbind(tags.result.pos2.from.module,
                                     tags.result.neg2.from.module)
               rm(list = c("anno.pos.from.module", "tags.result.pos2.from.module",
                           "anno.neg.from.module", "tags.result.neg2.from.module"))
               })


        # tags.result2 <- tags.result2[which(tags.result2$name %in% names(anno)),]
        peak.id1 <- unname(unlist(mapply(function(x, y){
          paste(x, y, sep = ";")
        },
        x = names(anno),
        y = anno)))

        peak.id2 <- paste(tags.result2$name, tags.result2$to, sep = ";")

        tags.result2_1 <- tags.result2[match(peak.id1, peak.id2),, drop = FALSE]
        tags.result2_2 <- tags.result2[-which(tags.result2$to %in% unique(unlist(anno))),, drop = FALSE]
        anno <- rbind(tags.result2_1, tags.result2_2)
        anno <- anno[anno$isotope == "[M]",, drop = FALSE]
      }else{
        if(polarity == "positive" | polarity == "both"){
          load(file.path(output.path, "intermediate_data", "tags.result.pos2.from.module"))
        }

        if(polarity == "negative" | polarity == "both"){
          load(file.path(output.path, "intermediate_data", "tags.result.neg2.from.module"))
        }

        switch(polarity,
               "positive" = {tags.result2 <- tags.result.pos2.from.module
               rm(list = c("tags.result.pos2.from.module"))},
               "negative" = {tags.result2 <- tags.result.neg2.from.module
               rm(list = c("tags.result.neg2.from.module"))},
               "both" = {
               tags.result2 <- rbind(tags.result.pos2.from.module,
                                     tags.result.neg2.from.module)
               rm(list = c("tags.result.pos2.from.module",
                           "tags.result.neg2.from.module"))
               })
        anno <- tags.result2
        anno <- anno[anno$isotope == "[M]",, drop = FALSE]
      }


###get group.data
     group.data <- getPathway(species = species, type = "metabolite")

     ## remove the pathway which has less than 3 metabolite
     remove.idx <- unname(which(unlist(lapply(group.data, function(x){
      sum(x %in% unique(anno$to)) < 3
     }))))

     if(length(remove.idx) > 0){
      group.data <- group.data[-remove.idx]
     }

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
                       output.path = temp.path5,
                       metabolic.network = kegg.rpair2)


     }

     #box plot and heatmap for module
     DNA.pathway.quantitative.result <-
       readr::read_csv(file.path(temp.path5,
                                 "DNA.pathway.quantitative.result.csv"),
                       progress = FALSE, col_types = readr::cols())

     DNA.pathway.quantitative.result <- as.data.frame(DNA.pathway.quantitative.result)

     pathway.node.quantitative.result <-
       readr::read_csv(file.path(temp.path5,
                                 "pathway.node.quantitative.result.csv"),
                       progress = FALSE, col_types = readr::cols())

     pathway.node.quantitative.result <- as.data.frame(pathway.node.quantitative.result)

     pathway.sample <- DNA.pathway.quantitative.result[,-1]
     pathway.name.id <- DNA.pathway.quantitative.result[,1]
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
                output.path = temp.path3,
                p.value = pathway.p.value,
                fc = pathway.fc,
                beeswarm = TRUE)

     pathway.sample1 <- t(apply(pathway.sample, 1, as.numeric))
     colnames(pathway.sample1) <- colnames(pathway.sample)
     rownames(pathway.sample1) <- rownames(pathway.sample)


     pdf(file.path(temp.path4, "Dysregulated.pathways.pdf"),
         width = 7, height = 7)
     par(mar = c(5, 5 ,4, 2))
     heatMap(sample = pathway.sample1, sample.info = sample.info, group = group)
     dev.off()

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

       pdf(file.path(temp.path4, paste(temp.pathway.name, "pdf",sep = ".")),
           width = 7, height = 7)
       par(mar = c(5, 5 ,4, 2))
       heatMap(sample = temp.sample,
               sample.info = sample.info,
               group = group)
       dev.off()

     }
    cat("\n")
    cat("metModule is done\n")
  })




##----------------------------------------------------------------------------
# title getSDM
# description Get dysregulated modules.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param tags2.pos Positive data from readAnnotation.
# param tags2.neg Negative data from readAnnotation.
# param metabolic.network kegg.rpair2.
# param metabolite.pos Positive kegg.adduct.
# param metabolite.neg Negative kegg.adduct.
# param kegg.compound kegg.compound
# param p.value.pos The p value of positive peaks.
# param p.value.neg The p value of negative peaks.
# param p.cutoff p value cutoff.
# param adduct.table.pos The positive adduct table.
# param adduct.table.neg The negative adduct table.
# param mz.tol The mz tolerance of KEGG matching.
# param rt.tol THe RT tolerance of KEGG matching (\%).
# param threads How many threads do you want to use? Default is the number of
# your PC threads - 3.
# param path The directory.
# param output.path The directory of outputing results.
# param use.old.null Use old null distribution of not.
# param species The species.

# tags2.pos = tags2.pos
# tags2.neg = tags2.neg
# polarity = "positive"
# metabolic.network = kegg.rpair2
# metabolite.pos = kegg.adduct.pos
# metabolite.neg = kegg.adduct.neg
# kegg.compound = kegg.compound
# p.value.pos = p.value.pos
# p.value.neg = p.value.neg
# p.cutoff = p.cutoff
# adduct.table.pos = adduct.table.pos
# adduct.table.neg = adduct.table.neg
# mz.tol = mz.tol
# rt.tol = rt.tol2
# threads = threads
# path = path
# output.path = output.path
# use.old.null = use.old.null
# species = species
# use.all.kegg.id = use.all.kegg.id
# only.mrn.annotation = only.mrn.annotation


setGeneric(name = "getSDM",
           def = function(tags2.pos,
                          tags2.neg,
                          polarity,
                          metabolic.network = kegg.rpair2,
                          metabolite.pos = kegg.adduct.pos,
                          metabolite.neg = kegg.adduct.neg,
                          kegg.compound = kegg.compound,
                          p.value.pos,
                          p.value.neg,
                          p.cutoff = 0.01,
                          adduct.table.pos,
                          adduct.table.neg,
                          mz.tol = 25,
                          rt.tol = 30,
                          threads = 3,
                          output.path = ".",
                          use.old.null = TRUE,
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          use.all.kegg.id = FALSE,
                          only.mrn.annotation = FALSE){

             if(sum(p.value.pos < p.cutoff) == 0 & sum(p.value.neg < p.cutoff) == 0) {
               stop("No dysregulated peaks.\n")
             }
             species <- match.arg(species)

             if(polarity == "positive" | polarity == "both"){
               cat("\n")
               cat("Positive.\n")
               result.pos <- tags2Result(tags2 = tags2.pos, score.cutoff = 0)
             }


             if(polarity == "negative" | polarity == "both"){
               cat("\n")
               cat("Negative.\n")
               result.neg <- tags2Result(tags2 = tags2.neg, score.cutoff = 0)
             }

             if(polarity == "positive") result.neg <- NULL
             if(polarity == "negative") result.pos <- NULL


             ##find the significant changed peaks
             if(only.mrn.annotation){
               if(polarity == "positive" | polarity == "both") index.pos <- which(p.value.pos < p.cutoff & !is.na(result.pos$to))
               if(polarity == "negative" | polarity == "both") index.neg <- which(p.value.neg < p.cutoff & !is.na(result.neg$to))

               switch(polarity,
                     "positive" = {index <- index.pos},
                     "negative" = {index <- index.neg},
                     "both" = {index <- c(index.pos, index.neg)})

               if(length(index) == 0){
                 warning("The peaks with p value less than 0.05 have no MRN based annotations.\nSo only.mrn.annotation has been set as FALSE.")
                 if(polarity == "positive" | polarity == "both") index.pos <- which(p.value.pos < p.cutoff)
                 if(polarity == "negative" | polarity == "both") index.neg <- which(p.value.neg < p.cutoff)
                 only.mrn.annotation <- FALSE
               }
             }else{
               if(polarity == "positive" | polarity == "both") index.pos <- which(p.value.pos < p.cutoff)
               if(polarity == "negative" | polarity == "both") index.neg <- which(p.value.neg < p.cutoff)
             }

             if(polarity == "positive" | polarity == "both") size.pos <- length(index.pos)
             if(polarity == "negative" | polarity == "both") size.neg <- length(index.neg)

             if(polarity == "positive" | polarity == "both"){
               result.pos1 <- result.pos[index.pos,]
               ##remove isotope annotation form result
               result.pos1 <- removeIsotopeFromResult(result = result.pos1)
             }

             if(polarity == "negative" | polarity == "both"){
               result.neg1 <- result.neg[index.neg,]
               result.neg1 <- removeIsotopeFromResult(result = result.neg1)
             }


             if(polarity == "positive" | polarity == "both"){
               #-----------------------------------------------------------------
               ####for positive mode
               to.pos1 <- result.pos1$to
               names(to.pos1) <- result.pos$name[index.pos]
               rm(index.pos)
               gc()
               to.pos1 <- lapply(to.pos1, function(x) {
                 if(is.na(x)) {return(x)
                 }else{
                   x <- strsplit(x, split = ";")[[1]]
                   x <- x[!duplicated(x)]
                   x <- paste(x, collapse = ";")
                   x
                 }
               })

               names(to.pos1) <- result.pos1$name

               anno.pos1 <- to.pos1
               ###kegg mz and rt matching
               mz.pos <- as.numeric(result.pos1$mz)
               rt.pos <- as.numeric(result.pos1$rt)
               cat("\n")
               cat("Annotate positive peaks utilizing KEGG database.\n")
               match.result.pos <- keggMatch(mz = mz.pos, rt = rt.pos,
                                             metabolite = metabolite.pos,
                                             mz.tol = mz.tol, rt.tol = rt.tol,
                                             polarity = "positive")

               ##anno.pos2 is the annotation results from KEGG according to mz and rt
               anno.pos2 <- lapply(match.result.pos, function(x) {
                 if(is.null(x)) return(NA)
                 temp <- unique(x$ID)
                 temp <- paste(temp, collapse = ";")
                 temp
               })

               names(anno.pos2) <- result.pos1$name

               if(use.all.kegg.id){
                 anno.pos <- mapply(function(x,y){
                   if(is.na(x)) {
                     if(is.na(y)) {
                       return(NA)
                     }else{
                       y <- strsplit(y, split = ";")[[1]]
                       y <- y[!is.na(y)]
                       return(y)
                     }
                   }else{
                     x <- strsplit(x, split = ";")[[1]]
                     x <- x[!is.na(x)]
                     y <- ifelse(is.na(y), y, strsplit(y, split = ";")[[1]])
                     # y <- strsplit(y, split = ";")[[1]]
                     y <- y[!is.na(y)]
                     x <- unique(c(x, y))
                     return(x)
                   }
                 },
                 x = anno.pos1,
                 y = anno.pos2)
               }else{
                 anno.pos <- mapply(function(x,y){
                   if(is.na(x)) {
                     if(is.na(y)) {
                       return(NA)
                     }else{
                       y <- strsplit(y, split = ";")[[1]]
                       y <- y[!is.na(y)]
                       return(y)
                     }
                   }else{
                     x <- strsplit(x, split = ";")[[1]]
                     x <- x[!is.na(x)]
                     return(x)
                   }
                 },
                 x = anno.pos1,
                 y = anno.pos2)
               }



               names(anno.pos) <- result.pos1$name

               index.pos <- unname(which(unlist(lapply(anno.pos, function(x) {!all(is.na(x))}))))
               result.pos2 <- result.pos1[index.pos,]
               rm(result.pos1)
               gc()
               anno.pos <- anno.pos[index.pos]
               anno.pos1 <- anno.pos1[index.pos]
               anno.pos2 <- anno.pos2[index.pos]

               match.result.pos2 <- match.result.pos[index.pos]
               rm(match.result.pos)
               gc()
             }


             ####for negative mode
             if(polarity == "negative" | polarity == "both"){
               to.neg1 <- result.neg1$to
               names(to.neg1) <- result.neg$name[index.neg]
               rm(index.neg)
               gc()
               to.neg1 <- lapply(to.neg1, function(x) {
                 if(is.na(x)) {return(x)
                 }else{
                   x <- strsplit(x, split = ";")[[1]]
                   x <- x[!duplicated(x)]
                   x <- paste(x, collapse = ";")
                   x
                 }
               })
               names(to.neg1) <- result.neg1$name
               anno.neg1 <- to.neg1
               ###kegg mz and rt matching
               mz.neg <- as.numeric(result.neg1$mz)
               rt.neg <- as.numeric(result.neg1$rt)
               cat("\n")
               cat("Annotate negative peaks utilizing KEGG database.\n")
               match.result.neg <- keggMatch(mz = mz.neg, rt = rt.neg,
                                             metabolite = metabolite.neg,
                                             mz.tol = mz.tol,rt.tol = rt.tol,
                                             polarity = "negative")

               ##anno.neg2 is the annotation results from KEGG according to mz and rt
               anno.neg2 <- lapply(match.result.neg, function(x) {
                 if(is.null(x)) return(NA)
                 temp <- unique(x$ID)
                 temp <- paste(temp, collapse = ";")
                 temp
               })

               names(anno.neg2) <- result.neg1$name

               if(use.all.kegg.id){
                 anno.neg <- mapply(function(x,y){
                   if(is.na(x)) {
                     if(is.na(y)) {
                       return(NA)
                     }else{
                       y <- strsplit(y, split = ";")[[1]]
                       y <- y[!is.na(y)]
                       return(y)
                     }
                   }else{
                     x <- strsplit(x, split = ";")[[1]]
                     x <- x[!is.na(x)]
                     y <- ifelse(is.na(y), y, strsplit(y, split = ";")[[1]])
                     # y <- strsplit(y, split = ";")[[1]]
                     y <- y[!is.na(y)]
                     x <- unique(c(x, y))
                     return(x)
                   }
                 },
                 x = anno.neg1,
                 y = anno.neg2)
               }else{
                 anno.neg <- mapply(function(x,y){
                   if(is.na(x)) {
                     if(is.na(y)) {
                       return(NA)
                     }else{
                       y <- strsplit(y, split = ";")[[1]]
                       y <- y[!is.na(y)]
                       return(y)
                     }
                   }else{
                     x <- strsplit(x, split = ";")[[1]]
                     x <- x[!is.na(x)]
                     return(x)
                   }
                 },
                 x = anno.neg1,
                 y = anno.neg2)
               }

               names(anno.neg) <- result.neg1$name

               index.neg <- unname(which(unlist(lapply(anno.neg, function(x) {!all(is.na(x))}))))
               result.neg2 <- result.neg1[index.neg,]
               rm(result.neg1)
               gc()
               anno.neg <- anno.neg[index.neg]
               anno.neg1 <- anno.neg1[index.neg]
               anno.neg2 <- anno.neg2[index.neg]

               match.result.neg2 <- match.result.neg[index.neg]
               rm(match.result.neg)
               gc()
             }



             #unique.id is the annotations of significant peaks
             switch(polarity,
                    "positive" = {anno <- anno.pos},
                    "negative" = {anno <- anno.neg},
                    "both" = {anno <- c(anno.pos, anno.neg)})
             unique.id <- unique(unname(unlist(anno)))
             unique.id <- unique.id[which(unique.id %in% igraph::V(metabolic.network)$name)]

             cat("\n")
             cat("Identify initial modules from dysregulated network.\n")
             cat("\n")
             module.result <- getModule(node.id = unique.id,
                                        total.id.number = length(unique.id),
                                        metabolic.network = metabolic.network,
                                        threads = threads,
                                        max.reaction = 3)
             ##module is the node name for each module, hidden.id is the hidden metabolites
             module <- module.result[[1]]
             activity.score <- module.result[[2]]
             hidden.id <- module.result[[3]]
             impact <- module.result[[4]]
             rm(module.result)
             gc()

             ##get NULL distribution
             cat("\n")
             cat("Get pseudo modules.\n")

             if(!use.old.null){
               ref.activity.score <- getNullDistribution(result.pos = result.pos,
                                                          result.neg = result.neg,
                                                          polarity,
                                                          mz.tol = mz.tol,
                                                          rt.tol = rt.tol,
                                                          metabolite.pos = metabolite.pos,
                                                          metabolite.neg = metabolite.neg,
                                                          metabolic.network = metabolic.network,
                                                          size.pos = size.pos,
                                                          size.neg = size.neg,
                                                          times = 20,
                                                          threads = threads,
                                                          use.all.kegg.id = use.all.kegg.id,
                                                          only.mrn.annotation = only.mrn.annotation,
                                                          progressbar = FALSE)

               ref.as <- unlist(ref.activity.score) + 0.00000001
               save(ref.as, file = file.path(output.path, "intermediate_data", "ref.as"), compress = "xz")
               rm(ref.activity.score)
               gc()
             }else{
               cat("\n")
               cat("Use old pseudo scores.\n")
               load(file.path(output.path, "intermediate_data", "ref.as"))
             }

             activity.score <- activity.score + 0.00000001

             ###
             para <- MASS::fitdistr(x = ref.as, densfun = "gamma")[[1]]
             # standard.distribution <- rgamma(n = length(ref.as), shape = para[1], rate = para[2])
             standard.distribution <- rgamma(n = 100000, shape = para[1], rate = para[2])
             # ks.test(x = ref.as, y = standard.distribution)


             cdf <- ecdf(x = standard.distribution)
             rm(list = "standard.distribution")
             gc()

             module.p <- 1 - cdf(activity.score)
             rm(list = c("standard.distribution", "cdf", "ref.as", "para"))
             gc()

             ###output result
             #order module accorind to p value
             module <- module[order(module.p)]
             impact <- impact[order(module.p)]
             module.p <- module.p[order(module.p)]


             ## get information of each module
             module.result <- lapply(module, function(x){
               temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)

               i.n <- sum(x %in% unique.id)
               h.n <- sum(x %in% hidden.id)
               t.n <- length(x)

               i.id <- x[which(x %in% unique.id)]
               h.id <- x[which(x %in% hidden.id)]

               i.name <- kegg.compound$Name[match(i.id, kegg.compound$ID)]
               i.name <- unlist(lapply(i.name, function(x) strsplit(x, split = ";")[[1]][1]))
               i.name <- paste(i.name, collapse = ";")

               h.name <- kegg.compound$Name[match(h.id, kegg.compound$ID)]
               h.name <- unlist(lapply(h.name, function(x) strsplit(x, split = ";")[[1]][1]))
               h.name <- paste(h.name, collapse = ";")

               t.name <- kegg.compound$Name[match(x, kegg.compound$ID)]
               t.name <- unlist(lapply(t.name, function(x) strsplit(x, split = ";")[[1]][1]))
               t.name <- paste(t.name, collapse = ";")

               i.id <- paste(i.id, collapse = ";")
               h.id <- paste(h.id, collapse = ";")
               t.id <- paste(x, collapse = ";")

               temp.result <- c(t.n, i.n, h.n, t.id, i.id, h.id, t.name, i.name, h.name)
               rm(list = c('t.n', 'i.n', 'h.n', 't.id', 'i.id', 'h.id', 't.name', 'i.name', 'h.name'))
               gc()
               temp.result
             }
             )


             module.result <- do.call(rbind, module.result)
             module.name <- paste("module", 1:length(module), sep = "")
             module.result <- data.frame(module.name, impact, module.p, module.result,
                                         stringsAsFactors = FALSE)
             colnames(module.result) <- c("Module.name", "Module.impact","p.value", "Module.size",
                                          "Detected.number", "Hidden.number",
                                          "All.metabolite.id",
                                          "Detected.ID", "Hidden.ID",
                                          "All.metabolite.name", "Detected.name",
                                          "Hidden.name")

             ###metabolite sets enrichment analysis for modules
             pbapply::pboptions(style = 1)
             cat("\n")
             cat("Functional annotation for dysregulated modules.\n")

             mse.result <- pbapply::pblapply(module, function(x){
               mseAnalysis(metabolite.id = x, species = species)
             })

             #add the mse.result to module.result
             pathway.name.id <- lapply(mse.result, function(x){
               if(is.null(x)) return(c(NA, NA))
               name.id <- rownames(x)
               name <- unlist(lapply(name.id, function(x) {
                 strsplit(x, ";")[[1]][1]
               }))
               id <- unlist(lapply(name.id, function(x) {
                 strsplit(x, ";")[[1]][2]
               }))
               name <- paste(name, collapse = ";")
               id <- paste(id, collapse = ";")
               return(c(name, id))
             })

             pathway.name.id <- do.call(rbind, pathway.name.id)
             colnames(pathway.name.id) <- c("pathway.name", "pathway.id")
             module.result <-data.frame(module.result, pathway.name.id,
                                        stringsAsFactors = FALSE)
             rm(list = "pathway.name.id")
             gc()

             save(module.result, file = file.path(output.path, "intermediate_data", "module.result"), compress = "xz")
             # write.csv(module.result, file.path(output.path, "module_information", "module.result.csv"), row.names = FALSE)


             ##transform module.result to moduleInfo
             module.result1 <- apply(module.result, 1, list)

             module.result2 <- lapply(module.result1, function(module){
               module <- module[[1]]
               new(Class = "moduleInfo",
                   Module.name = unname(module["Module.name"]),
                   Module.size = as.numeric(unname(module["Module.size"])),
                   Module.impact = as.numeric(unname(module["Module.impact"])),
                   p.value = as.numeric(unname(module["p.value"])),
                   Detected.number = as.numeric(unname(module["Detected.number"])),
                   Hidden.number = as.numeric(unname(module["Hidden.number"])),
                   All.metabolite.id = strsplit(unname(module["All.metabolite.id"]), ";")[[1]],
                   Detected.ID = strsplit(unname(module["Detected.ID"]), ";")[[1]],
                   Hidden.ID = strsplit(unname(module["Hidden.ID"]), ";")[[1]],
                   All.metabolite.name = strsplit(unname(module["All.metabolite.name"]), ";")[[1]],
                   Detected.metabolite.name = strsplit(unname(module["Detected.metabolite.name"]), ";")[[1]],
                   Hidden.metabolite.name = strsplit(unname(module["Hidden.metabolite.name"]), ";")[[1]],
                   msea = data.frame()
               )
             })


             mse.result <- lapply(mse.result, function(x){
               if(is.null(x)) return(data.frame())
               x
             })

             module.result2 <- mapply(function(x, y){
               x@msea <- y
               x
             },
             x = module.result2,
             y = mse.result)

             save(module.result2, file = file.path(output.path, "intermediate_data","module.result2"), compress = "xz")

             ##using module to filter MRN annotation and KEGG annotation result
             index <- which(module.p < 0.05)
             if(length(index) == 0) {
               cat("\n")
               cat("There are no initial modules with p valuess less than 0.05.\n")
             }
             if(length(index) > 0){
               ##for positive peaks
               cat("\n")
               cat("Filter annotations according to dysregulated modules.\n")

               ##positive
               if(polarity == "positive" | polarity == "both"){
                 return.result.pos <- filterAnnotationByEnrichment(module = module,
                                                                   index = index,
                                                                   anno1 = anno.pos1,
                                                                   anno2 = anno.pos2,
                                                                   result2 = result.pos2, tags2 = tags2.pos,
                                                                   match.result2 = match.result.pos2,
                                                                   mz.tol = mz.tol,
                                                                   rt.tol = rt.tol)
                 annotation1 <- return.result.pos[[5]]
                 annotation2 <- return.result.pos[[6]]
                 return.result.pos <- return.result.pos[-c(5,6)]

                 MRN.annotation.of.dysregulated.peak.filtered.by.module.pos <- annotation1
                 KEGG.annotation.of.dysregulated.peak.filtered.by.module.pos <- annotation2
                 save(MRN.annotation.of.dysregulated.peak.filtered.by.module.pos,
                      file = file.path(output.path, "intermediate_data",
                                                    "MRN.annotation.of.dysregulated.peak.filtered.by.module.pos"))

                 save(KEGG.annotation.of.dysregulated.peak.filtered.by.module.pos,
                      file = file.path(output.path, "intermediate_data",
                                       "KEGG.annotation.of.dysregulated.peak.filtered.by.module.pos"))

                 rm(list = c("annotation1", "annotation2"))

               }

               ##negative

               if(polarity == "negative" | polarity == "both"){
                 return.result.neg <- filterAnnotationByEnrichment(module = module,
                                                                   index = index,
                                                                   anno1 = anno.neg1,
                                                                   anno2 = anno.neg2,
                                                                   result2 = result.neg2,
                                                                   tags2 = tags2.neg,
                                                                   match.result2 = match.result.neg2,
                                                                   mz.tol = mz.tol,
                                                                   rt.tol = rt.tol)
                 annotation1 <- return.result.neg[[5]]
                 annotation2 <- return.result.neg[[6]]
                 return.result.neg <- return.result.neg[-c(5,6)]

                 MRN.annotation.of.dysregulated.peak.filtered.by.module.neg <- annotation1
                 KEGG.annotation.of.dysregulated.peak.filtered.by.module.neg <- annotation2
                 save(MRN.annotation.of.dysregulated.peak.filtered.by.module.neg,
                      file = file.path(output.path, "intermediate_data",
                                       "MRN.annotation.of.dysregulated.peak.filtered.by.module.neg"))

                 save(KEGG.annotation.of.dysregulated.peak.filtered.by.module.neg,
                      file = file.path(output.path, "intermediate_data",
                                       "KEGG.annotation.of.dysregulated.peak.filtered.by.module.neg"))

                 rm(list = c("annotation1", "annotation2"))

               }

             }


             if(polarity == "both"){
               names(return.result.pos) <- c("tags2.pos", "anno.pos", "anno.pos1", "anno.pos2")
               names(return.result.neg) <- c("tags2.neg", "anno.neg", "anno.neg1", "anno.neg2")
               return.result <- list(return.result.pos, return.result.neg)
               names(return.result) <- c("positive", "negative")
               rm(list = c("return.result.pos", "return.result.neg"))
             }

             if(polarity == "positive"){
               return.result <- list(return.result.pos, NULL)
               rm(list = c("return.result.pos"))
             }

             if(polarity == "negative"){
               return.result <- list(NULL, return.result.neg)
               rm(list = c("return.result.neg"))
             }
             return.result <- return.result
           })


#------------------------------------------------------------------------------
# title getNullDistribution
# description Get NULL distribution from a network.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param result.pos The positive result from tags2.
# param result.neg The negative result from tags2.
# param metabolite.pos Positive kegg.adduct.
# param metabolite.neg Negative kegg.adduct.
# param metabolic.network kegg.rpair2.
# param polarity The polarity.
# param size.pos The size of positive peaks shoube be selected.
# param size.neg The size of negative peaks shoube be selected.
# param times repeat times
# param threads How many threads do you want to use?
# param mz.tol The tolerance of mz.
# param rt.tol The tolerance of RT (\%).
# param use.all.kegg.id Use all the KEGG annotations or not.
# return The NULL activity scores.


setGeneric(name = "getNullDistribution",
           def = function(result.pos,
                          result.neg,
                          polarity,
                          metabolite.pos,
                          metabolite.neg,
                          metabolic.network = kegg.rpair2,
                          size.pos,
                          size.neg,
                          times = 20,
                          threads = 5,
                          mz.tol = 25,
                          rt.tol = 30,
                          use.all.kegg.id = FALSE,
                          only.mrn.annotation = FALSE,
                          progressbar = FALSE){

             if(only.mrn.annotation){
               if(polarity == "positive" | polarity == "both") result.pos <- result.pos[!is.na(result.pos$to),]
               if(polarity == "negative" | polarity == "both") result.neg <- result.neg[!is.na(result.neg$to),]
             }
             # pbapply::pboptions(type="timer", style = 1)
             acti.score <- vector(mode = "list", length = times)

             for(i in seq_len(times)){
               # if(i %% 10 == 0 | i == 1) cat(i); cat(" ")
               cat(i,"/",times); cat("\n")
               if(polarity == "positive" | polarity == "both") {
                 idx.pos <- sample(1:nrow(result.pos), size.pos)
                 result.pos1 <- result.pos[idx.pos,]
                 result.pos1 <- removeIsotopeFromResult(result = result.pos1)
                 temp.mz.pos <- as.numeric(result.pos$mz[idx.pos])
                 temp.rt.pos <- as.numeric(result.pos$rt[idx.pos])
                 match.result.pos <- keggMatch(mz = temp.mz.pos,
                                               rt = temp.rt.pos,
                                               metabolite = metabolite.pos,
                                               mz.tol = mz.tol,rt.tol = rt.tol,
                                               polarity = "positive")
                 anno1.pos <- result.pos1$to

                 anno2.pos <- lapply(match.result.pos, function(x) {
                   if(is.null(x)) return(NA)
                   # unique(x$ID)
                   paste(unique(x$ID), collapse = ";")
                 })

                 if(use.all.kegg.id){
                   anno.pos <- mapply(function(x,y){
                     if(is.na(x)) {
                       if(is.na(y)) {
                         return(NA)
                       }else{
                         y <- strsplit(y, split = ";")[[1]]
                         y <- y[!is.na(y)]
                         return(y)
                       }
                     }else{
                       x <- strsplit(x, split = ";")[[1]]
                       x <- x[!is.na(x)]
                       y <- ifelse(is.na(y), y, strsplit(y, split = ";")[[1]])
                       # y <- strsplit(y, split = ";")[[1]]
                       y <- y[!is.na(y)]
                       x <- unique(c(x, y))
                       return(x)
                     }
                   },
                   x = anno1.pos,
                   y = anno2.pos)
                 }else{
                   anno.pos <- mapply(function(x,y){
                     if(is.na(x)) {
                       if(is.na(y)) {
                         return(NA)
                       }else{
                         y <- strsplit(y, split = ";")[[1]]
                         y <- y[!is.na(y)]
                         return(y)
                       }
                     }else{
                       x <- strsplit(x, split = ";")[[1]]
                       x <- x[!is.na(x)]
                       return(x)
                     }
                   },
                   x = anno1.pos,
                   y = anno2.pos)
                 }


                 anno.pos <- unlist(anno.pos)
                 anno.pos <- anno.pos[!is.na(anno.pos)]
                 anno.pos <- unique(anno.pos)
                 rm(list = c("anno1.pos", "anno2.pos"))
                 gc()
               }

               if(polarity == "negative" | polarity == "both"){
                 idx.neg <- sample(1:nrow(result.neg), size.neg)
                 result.neg1 <- result.neg[idx.neg,]
                 result.neg1 <- removeIsotopeFromResult(result = result.neg1)
                 temp.mz.neg <- as.numeric(result.neg$mz[idx.neg])
                 temp.rt.neg <- as.numeric(result.neg$rt[idx.neg])
                 match.result.neg <- keggMatch(mz = temp.mz.neg, rt = temp.rt.neg,
                                               metabolite = metabolite.neg,
                                               mz.tol = mz.tol,rt.tol = rt.tol,
                                               polarity = "negative")


                 anno1.neg <- result.neg1$to
                 # anno1.neg <- lapply(result.neg1$to, function(x){
                 #   if(is.na(x)) {return(NA)}else{strsplit(x, split = ";")[[1]]}
                 # })

                 anno2.neg <- lapply(match.result.neg, function(x) {
                   if(is.null(x)) return(NA)
                   unique(x$ID)
                 })

                 if(use.all.kegg.id){
                   anno.neg <- mapply(function(x,y){
                     if(is.na(x)) {
                       if(is.na(y)) {
                         return(NA)
                       }else{
                         y <- strsplit(y, split = ";")[[1]]
                         y <- y[!is.na(y)]
                         return(y)
                       }
                     }else{
                       x <- strsplit(x, split = ";")[[1]]
                       x <- x[!is.na(x)]
                       y <- ifelse(is.na(y), y, strsplit(y, split = ";")[[1]])
                       # y <- strsplit(y, split = ";")[[1]]
                       y <- y[!is.na(y)]
                       x <- unique(c(x, y))
                       return(x)
                     }
                   },
                   x = anno1.neg,
                   y = anno2.neg)
                 }else{
                   anno.neg <- mapply(function(x,y){
                     if(is.na(x)) {
                       if(is.na(y)) {
                         return(NA)
                       }else{
                         y <- strsplit(y, split = ";")[[1]]
                         y <- y[!is.na(y)]
                         return(y)
                       }
                     }else{
                       x <- strsplit(x, split = ";")[[1]]
                       x <- x[!is.na(x)]
                       return(x)
                     }
                   },
                   x = anno1.neg,
                   y = anno2.neg)
                 }

                 anno.neg <- unlist(anno.neg)
                 anno.neg <- anno.neg[!is.na(anno.neg)]
                 anno.neg <- unique(anno.neg)
                 rm(list = c("anno1.neg", "anno2.neg"))
                 gc()
               }

               #---------------------------------------------------------------

               switch(polarity,
                      "positive" = {anno <- anno.pos},
                      "negative" = {anno <- anno.neg},
                      "both" = {anno <- c(anno.pos, anno.neg)})

               unique.id <- anno[which(anno %in% igraph::V(metabolic.network)$name)]

               ac.score <- getModule(node.id = unique.id,
                                     metabolic.network = metabolic.network,
                                     total.id.number = length(unique.id),
                                     threads = threads,
                                     progressbar = FALSE,
                                     max.reaction = 3,
                                     calculate.impact = FALSE)[[2]]

               acti.score[[i]] <- unlist(ac.score)
             }
             acti.score <- acti.score
           })






setGeneric(name = "outputModuleInformation",
           def = function(module.result,
                          module.result2,
                          output.path,
                          temp.index){
             temp.path <- file.path(output.path, "Module_functional_information")
             dir.create(temp.path)

             for (i in temp.index) {
               temp.module <- module.result2[[i]]
               temp.module.msea <- temp.module@msea
               if (nrow(temp.module.msea) == 0) next()
               temp.module.name <- temp.module@Module.name
               pdf(file.path(temp.path, paste(temp.module.name, "MSEA.pdf", sep = ".")),
                   width = 7,
                   height = 7)
               moduleplot(object = temp.module, cex.label = 1, n = 20)
               dev.off()
               write.csv(temp.module@msea,
                         file.path(temp.path,
                                   paste(temp.module.name, "MSEA.csv", sep = ".")))
               rm(list = c("temp.module", "temp.module.msea", "temp.module.name"))
               gc()

             }


             pdf(file.path(output.path, "dysregulated.module.overview.pdf"),
                 width = 7, height = 7)
             par(mar = c(5, 5, 4, 2))
             moduleScatterPlot(module.result = module.result, p.cutoff = 0.05)
             dev.off()

             write.csv(module.result, file.path(output.path, "dysregulated.module.information.csv"))

           })




#--------------------------------------------------

# param mz The vector of mz.
# param rt The vector of RT.
# param metabolite kegg.adduct.
# param mz.tol The mz tolerance of matching.
# param rt.tol The RT tolerance of matching.
# param polarity The polarity.
# param threads How many threads do you want to use? Default is the number of
# your PC threads - 3.
# param progressbar Show progressbar or not.
# return KEGG matching result

setGeneric(name = "keggMatch",
           def = function(mz,
                          rt,
                          metabolite,
                          mz.tol = 25,
                          rt.tol = 30,
                          polarity = c("positive", "negative")){
             options(warn = -1)
             polarity <- match.arg(polarity)
             metabolite <- metabolite[metabolite$Mode == polarity,]
             metabolite <- metabolite[metabolite$Accurate.mass<=max(mz),]

             #----------------------------------------------------------------
             metabolite.mz <- sapply(metabolite$Accurate.mass,function(mass){ifelse(mass>=400,mass,400)})
             result <- mapply(function(x, y){
               # mz.error <- abs(x - metabolite$Accurate.mass)*10^6/metabolite$Accurate.mass
               mz.error <- abs(x - metabolite$Accurate.mass)*10^6/metabolite.mz
               rt.error <- abs(y - metabolite$RT)*100/metabolite$RT
               index <- which(mz.error <= mz.tol & rt.error <= rt.tol)
               if(length(index) == 0) return(NULL)
               temp <- data.frame(metabolite[index,,drop = FALSE],
                                  mz.error[index], rt.error[index],
                                  stringsAsFactors = FALSE)
               colnames(temp)[c(14,15)] <- c("mz.error", "rt.error")
               list(temp)
             },
             x = mz,
             y = rt)

             result <- lapply(result, function(x) {
               if(is.null(x)) return(NULL)
               if(length(x) == 1) return(x[[1]])
               return(x)
             })

             result <- result
           })



# setGeneric(name = "keggMatch",
#            def = function(mz,
#                           rt,
#                           metabolite,
#                           mz.tol = 25,
#                           rt.tol = 30,
#                           polarity = c("positive", "negative"),
#                           threads = 3,
#                           progressbar = TRUE){
#              options(warn = -1)
#              polarity <- match.arg(polarity)
#              metabolite <- metabolite[metabolite$Mode == polarity,]
#
#              info <- data.frame(mz,rt, stringsAsFactors = FALSE)
#              info <- apply(info, 1, list)
#
#              temp.fun <- function(x, metabolite, mz.tol, rt.tol){
#                temp.mz <- x[[1]][[1]]
#                temp.rt <- x[[1]][[2]]
#                mz.error <- abs(metabolite$Accurate.mass - temp.mz)*10^6/temp.mz
#                rt.error <- abs(metabolite$RT - temp.rt)*100/temp.rt
#                temp.index <- which(mz.error <= mz.tol & rt.error <= rt.tol)
#                if(length(temp.index) == 0) return(NULL)
#                temp.result <- data.frame(metabolite, mz.error,
#                                          rt.error,
#                                          stringsAsFactors = FALSE)[temp.index,,drop = FALSE]
#                rm(list = c("temp.mz", "temp.rt", "mz.error", "rt.error", "temp.index"))
#                temp.result
#              }
#
#              # if(progressbar) cat("KEGG matching\n")
#
#              result <-
#                BiocParallel::bplapply(info, temp.fun,
#                                       BPPARAM = BiocParallel::SnowParam(workers = threads,
#                                                                         progressbar = progressbar),
#                                       metabolite = metabolite,
#                                       mz.tol = mz.tol,
#                                       rt.tol = rt.tol)
#              result <- result
#            })
#
#
# setGeneric(name = "keggMatch1",
#            def = function(mz,
#                           rt,
#                           metabolite,
#                           mz.tol = 25,
#                           rt.tol = 30,
#                           polarity = c("positive", "negative"),
#                           threads = 3,
#                           progressbar = TRUE){
#              options(warn = -1)
#              polarity <- match.arg(polarity)
#              metabolite <- metabolite[metabolite$Mode == polarity,]
#
#              #----------------------------------------------------------------
#
#              system.time(mz.error <- lapply(mz, function(x){
#                abs(x - metabolite$Accurate.mass)*10^6/x
#              }))
#
#              system.time(rt.error <- lapply(rt, function(x){
#                abs(x - metabolite$RT)*100/x
#              }))
#
#              system.time(index <- mapply(function(x, y){
#              which(x <= mz.tol & y <= rt.tol)
#              },
#              x = mz.error,
#              y = rt.error))
#
#              system.time(result <- mapply(function(x, y, z){
#              if(length(z) == 0) return(NULL)
#                temp <- data.frame(metabolite[z,,drop = FALSE], x[z], y[z],
#                           stringsAsFactors = FALSE)
#                colnames(temp)[c(14,15)] <- c("mz.error", "rt.error")
#                temp
#
#              },
#              x = mz.error,
#              y = rt.error,
#              z = index))
#
#              result <- result
#            })






# kegg.distance <- distances(graph = kegg.rpair2, v = V(kegg.rpair2)$name, to = V(kegg.rpair2)$name)
# kegg.distance[which(is.infinite(kegg.distance), arr.ind = TRUE)] <- 30
# save(kegg.distance, file = "kegg.distance")

#------------------------------------------------------------------------------
# title getModule
# description Get modules from a network.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param node.id The vector of node ID.
# param total.id.number The number of node.
# param metabolic.network kegg.rpair2.
# param threads How many threads do you want to use? Default is the number of
#  your PC threads - 3.
# para max.reaction The max raction number to find hidden nodes.
# param progressbar Show progressbar or not.
# param calculate.impact Calculate impact or not.
# return The module information.

setGeneric(name = "getModule",
           function(node.id,
                    total.id.number,
                    metabolic.network = kegg.rpair2,
                    threads = 3,
                    max.reaction = 3,
                    progressbar = TRUE,
                    calculate.impact = TRUE){

             unique.id <- unique(node.id)
             ##find the hidden metabolite which can connect two input metabolites in up to
             ## three reaction
             temp.fun <- function(name, metabolic.network,
                                  unique.id, max.reaction){
               options(warn = -1)
               library(igraph, quietly = TRUE,
                       logical.return = FALSE, warn.conflicts = FALSE)
               idx <- match(name, unique.id)
               to <- unique.id[-c(1:idx)]

               distance <- igraph::distances(graph = metabolic.network, v = name, to = to)[1,]
               distance[is.infinite(distance)] <- 5
               to <- to[which(distance <= max.reaction)]
               if(length(to) == 0) return(NULL)
               result <-
                 igraph::shortest_paths(graph = metabolic.network,
                                        from = name, to = to)[[1]]
               # len <- unlist(lapply(result, function(x) length(x)))
               # result <- result[len<=4]
               result <- lapply(result, names)
               result <- unique(unlist(result))
               result <- setdiff(result, unique.id)
               rm(list = c("idx", "to", "distance"))
               gc()
               return(result)
             }

             if(progressbar) cat("Find hidden metabolites.\n")
             hidden.id <-
               BiocParallel::bplapply(unique.id[-length(unique.id)],
                                      FUN = temp.fun,
                                      BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                        progressbar = progressbar),
                                      metabolic.network = metabolic.network,
                                      unique.id = unique.id,
                                      max.reaction = max.reaction)

             hidden.id <- unlist(hidden.id)
             hidden.id <- unique(hidden.id)

             subnetwork <- igraph::subgraph(graph = metabolic.network,
                                            v = c(unique.id, hidden.id))

             # fc <- igraph::cluster_fast_greedy(subnetwork)
             wt <- igraph::cluster_walktrap(subnetwork)

             node.name <- wt$names
             membership <- wt$membership
             membership1 <- table(membership)

             group <- lapply(names(membership1), function(x) {
               node.name[which(membership == as.numeric(x))]
             })

             ##remove group which less than 3 nodes
             group <- group[unlist(lapply(group, length)) > 3]

             ###the group whose node more than 100 should be re-subgraph
             idx <- which(unlist(lapply(group, length)) > 100)

             while(length(idx) > 0){
               large.group <- group[idx]
               group <- group[-idx]
               new.group <- lapply(large.group, function(x) {
                 temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)
                 temp.wt <- igraph::cluster_walktrap(temp.graph)
                 temp.node.name <- temp.wt$names
                 temp.membership <- temp.wt$membership
                 temp.membership1 <- table(temp.membership)

                 temp.group <- lapply(names(temp.membership1), function(x) {
                   temp.node.name[which(temp.membership == as.numeric(x))]
                 })

                 ##remove group which less than 3 nodes
                 temp.group <- temp.group[unlist(lapply(temp.group, length)) > 3]
                 temp.group
               })

               if(length(new.group) == 1){
                 new.group <- new.group[[1]]
               }else{
                 new.group1 <- new.group[[1]]
                 for(i in 2:length(new.group)){
                   new.group1 <- c(new.group1, new.group[[i]])
                 }
                 new.group <- new.group1
               }
               group <- c(group, new.group)
               idx <- which(unlist(lapply(group, length)) > 100)
             }
             # save(group, file = "group")
             ##remove group which edge number is less than node number
             remove.idx <- which(unlist(lapply(group, function(x){
               temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)
               node.number <- length(igraph::V(temp.graph))
               edge.number <- length(igraph::E(temp.graph))
               diff.number <- node.number - edge.number
               degree <- igraph::degree(graph = temp.graph, v = x)
               if(diff.number > 0 & all(degree < 3)){
                 return(TRUE)
               }else{
                 FALSE
               }
             })))

             group <- group[-remove.idx]

             ##remove the hidden nodes whose degree is less than 2
             if(progressbar) cat("Remove the hidden metabolites whose degrees are less than two.\n")
             group <- lapply(group, function(x){
               temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)
               degree <- igraph::degree(graph = temp.graph, v = x)

               class <- sapply(names(degree), function(x){
                 ifelse(is.element(x, unique.id), "Detetcted", "Hidden")
               })

               remove.idx <- which(degree == 1 & class == "Hidden")

               while(length(remove.idx) > 0){
                 x <- x[-remove.idx]
                 temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)

                 degree <- igraph::degree(graph = temp.graph, v = x)

                 class <- sapply(names(degree), function(x){
                   ifelse(is.element(x, unique.id), "Detetcted", "Hidden")
                 })
                 remove.idx <- which(degree == 1 & class == "Hidden")

               }
               x
             })

             group <- group[unlist(lapply(group, length)) > 3]

             # pbapply::pboptions(type = "timer", style = 1)
             if(progressbar) {cat("\n")
               cat("Calculate activity scores of dysregulated modules.\n")}

             temp.fun <- function(x, unique.id, hidden.id, total.id.number,
                                  metabolic.network = metabolic.network){
               library(igraph, quietly = TRUE,
                       logical.return = FALSE, warn.conflicts = FALSE)
               options(warn = -1)

               temp.graph <- igraph::subgraph(graph = metabolic.network, v = x)
               input.number <- sum(x %in% unique.id)
               hidden.number <- sum(x %in% hidden.id)

               m <- length(igraph::E(metabolic.network))
               e <- length(igraph::E(temp.graph))
               node.name <- names(igraph::V(temp.graph))

               value <- sapply(node.name, function(node){
                 temp.value <- sapply(node.name, function(x){
                   temp1 <- igraph::degree(graph = metabolic.network, v = x)
                   temp2 <- igraph::degree(graph = metabolic.network, v = node)
                   temp1*temp2/(4*m^2)
                 })
                 sum(temp.value)
               })

               modularity <- e/m - sum(value)
               # modularity
               q <- modularity * sqrt(total.id.number/length(x))
               a <- input.number * q/length(x)
               a
             }


             system.time(activity.score <-
                           BiocParallel::bplapply(group, temp.fun,
                                                  BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                                    progressbar = progressbar),
                                                  unique.id = unique.id,
                                                  hidden.id = hidden.id,
                                                  total.id.number = total.id.number,
                                                  metabolic.network = metabolic.network))

             activity.score <- unlist(activity.score)

             if(calculate.impact){
               ##Module impact
               if(progressbar) cat("Calculate impacts of dysregulated modules.\n")
               impact <- lapply(group, function(temp.group){
                 temp.graph <- igraph::subgraph(graph = metabolic.network,
                                                v = temp.group)
                 centrality <- getCentrality(graph = temp.graph, type = "d")
                 temp.idx <- which(temp.group %in% unique.id)
                 impact <- sum(centrality[temp.idx])/sum(centrality)
                 impact
               })
               impact <- unlist(impact)
               module <- group
               result <- list(module, activity.score, hidden.id, impact)
               names(result) <- c("module", "activity.score", "hidden.id", "impact")
             }else{
               module <- group
               result <- list(module, activity.score, hidden.id)
               names(result) <- c("module", "activity.score", "hidden.id")
             }
             result <- result
           })
