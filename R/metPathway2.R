#' @title metPathway2
#' @description Metabolic pathway analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param annotation.result.pos Positive annotation.result from ms2Annotation.
#' @param annotation.result.neg Negative annotation.result from ms2Annotation.
#' @param group The group you want to use.
#' @param max.isotope The max number of isotope.
#' @param uni.test The method of univariate test.
#' @param column hilic or rp.
#' @param pos.path The positive directory.
#' @param neg.path The negative directory.
#' @param correct Correct p value or not.
#' @param p.cutoff The cutoff of p value.
#' @param mz.tol The mz tolerance of KEGG matching.
#' @param rt.tol1 The RT tolerance of isotope annotation (second).
#' @param rt.tol2 The RT tolerance of KEGG matching (\%).
#' @param cor.tol THe tolerance of correlation.
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
#' @export
#

# annotation.result.pos = "ms2.match.annotation.result.csv"
# annotation.result.neg = "ms2.match.annotation.result.csv"
# group = c("Normal", "AD")
# max.isotope = 4
# uni.test = "t"
# column = "hilic"
# pos.path = file.path(".", "POS")
# neg.path = file.path(".", "NEG")
# output.path = file.path(".", "POS and NEG", "Dysregulated_network_analysis_result")
# correct = FALSE
# p.cutoff = 0.05
# mz.tol = 25
# rt.tol1 = 3# absolute, second
# rt.tol2 = 30#relative, %
# cor.tol = 0
# int.tol = 500
# threads = 3
# use.old.null = FALSE
# species = "hsa"
# use.all.kegg.id = FALSE
# only.mrn.annotation = FALSE
#


# metPathway2(annotation.result.pos = "ms2.match.annotation.result.csv",
#          annotation.result.neg = "ms2.match.annotation.result.csv",
#          group = c("Normal", "AD"),
#          max.isotope = 4,
#          uni.test = "t",
#          column = "hilic",
#          pos.path = file.path(".", "POS"),
#          neg.path = file.path(".", "NEG"),
#          output.path = file.path(".", "POS and NEG", "Dysregulated_network_analysis_result"),
#          correct = FALSE,
#          p.cutoff = 0.05,
#          mz.tol = 25,
#          rt.tol1 = 3,# absolute, second
#          rt.tol2 = 30,#relative, %
#          cor.tol = 0,
#          int.tol = 500,
#          threads = 3,
#          use.old.null = TRUE,
#          species = "hsa",
#          use.all.kegg.id = FALSE,
#          only.mrn.annotation = FALSE)


setGeneric(
  name = "metPathway2",
  def = function(annotation.result.pos = "ms2.match.annotation.result.csv",
                 annotation.result.neg = "ms2.match.annotation.result.csv",
                 group,
                 max.isotope = 4,
                 uni.test = c("t", "wilcox", "anova"),
                 column = c("hilic", "rp"),
                 pos.path = file.path(".", "POS"),
                 neg.path = file.path(".", "NEG"),
                 output.path = file.path(".", "POS and NEG", "Dysregulated_network_analysis_result"),
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
                 only.mrn.annotation = FALSE) {

    now.path <- getwd()
    now.path1 <- strsplit(now.path, split = "/")[[1]]
    if(tail(now.path1, 1) == "POS" | tail(now.path1, 1) == "NEG"){
      new.path <- paste(now.path1[-length(now.path1)], collapse = "/")
      pos.path <- file.path(new.path, "POS")
      neg.path = file.path(new.path, "NEG")
      output.path = file.path(new.path, "POS and NEG", "Dysregulated_network_analysis_result")
    }else{
      cat("Your work directory is: ", now.path, ". \nPlease make sure the work directory is 'POS' or 'NEG', or above them.\n", sep = "")
    }
    dir.create(file.path(new.path, "POS and NEG"))
    dir.create(output.path)
    dir.create(file.path(output.path, "intermediate_data"))
    # dir.create(file.path(output.path, "Pathway_information"))
    dir.create(file.path(output.path, "pathway_information"))

    if(missing(group)) {
      stop("You must give the names of group which you want to process.")}

    uni.test <- match.arg(uni.test)
    column <- match.arg(column)
    species <- match.arg(species)

    metPathway.parameters <- c(paste(group, collapse = ";"), max.isotope,
                              uni.test, column,
                              correct,
                              p.cutoff, mz.tol,
                              rt.tol1, rt.tol2,
                              threads, use.old.null, species)
    metPathway.parameters <- data.frame(c("group", "max.isotope",
                                         "uni.test", "column",
                                         "correct",
                                         "p.cutoff", "mz.tol",
                                         "rt.tol1", "rt.tol2",
                                         "threads", "use.old.null",
                                         "species"),
                                       metPathway.parameters, stringsAsFactors = FALSE)
    colnames(metPathway.parameters) <- c('parameter', "value")
    write.csv(metPathway.parameters, file.path(output.path, "metPathway2.parameters.csv"),
              row.names = FALSE)

    # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))

    #check for need files
    ####positive
    file <- c(dir(file.path(pos.path, "MS2_match_result_POS","intermediate_data")),
              dir(file.path(pos.path, "MRN_annotation_result_POS","intermediate_data")))
    need.file <- c(annotation.result.pos,
                   "tags2.after.redundancy.remove",
                   "sample.info.csv", "rt.result")

    file.check <- which(!need.file %in% file)
    if (length(file.check) > 0) {
      stop(paste(need.file[file.check], collapse = " & "),
           " are not in the positive directory.")
    }


    ####negative
    file <- c(dir(file.path(neg.path, "MS2_match_result_NEG","intermediate_data")),
              dir(file.path(neg.path, "MRN_annotation_result_NEG","intermediate_data")))
    need.file <- c(annotation.result.neg,
                   "tags2.after.redundancy.remove",
                   "sample.info.csv", "rt.result")

    file.check <- which(!need.file %in% file)
    if (length(file.check) > 0) {
      stop(paste(need.file[file.check], collapse = " & "),
           " are not in the negative directory.")
    }

    ###check for temp.mse.result
    if(use.old.null) {
      file <- dir(file.path(output.path, "intermediate_data"))
      need.file <- "temp.mse.result"
      file.check <- which(!need.file %in% file)
      if (length(file.check) > 0) {
        stop(paste(need.file[file.check], collapse = " & "),
             " is not in the output directory.")
      }
    }

    load(file.path(pos.path, "MRN_annotation_result_POS","intermediate_data",
                   "tags2.after.redundancy.remove"))
    tags2.after.redundancy.remove.pos <- tags2.after.redundancy.remove
    rm(list="tags2.after.redundancy.remove")
    gc()

    #change peak name
    tags2.after.redundancy.remove.pos <-
      lapply(tags2.after.redundancy.remove.pos, function(x){
        x@name <- paste(x@name, "POS", sep = "_")
        x
      })

    load(file.path(neg.path, "MRN_annotation_result_NEG", "intermediate_data",
                   "tags2.after.redundancy.remove"))
    tags2.after.redundancy.remove.neg <- tags2.after.redundancy.remove
    rm(list="tags2.after.redundancy.remove")
    gc()

    #change paek name
    tags2.after.redundancy.remove.neg <-
      lapply(tags2.after.redundancy.remove.neg, function(x){
        x@name <- paste(x@name, "NEG", sep = "_")
        x
      })

    load(file.path(pos.path, "MRN_annotation_result_POS","intermediate_data", "rt.result"))
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


    sample.info <-
      readr::read_csv(
        file.path(pos.path, "MS2_match_result_POS","intermediate_data", "sample.info.csv"),
        progress = FALSE,
        col_types = readr::cols()
      )

    sample.info <- as.data.frame(sample.info)

    file.copy(file.path(pos.path, "MS2_match_result_POS","intermediate_data", "sample.info.csv"),
              to = file.path(output.path, "intermediate_data"))


    if(any(!is.element(group, unique(sample.info[,2])))){
      idx <- which(!is.element(group, unique(sample.info[,2])))
      stop(paste(group[idx], collapse = " and "),
           ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
    }

    sample.info <- sample.info[sample.info[,2] %in% group,]

    annotation.result.pos <-
      readr::read_csv(file.path(pos.path, "MS2_match_result_POS", annotation.result.pos),
                      progress = FALSE,
                      col_types = readr::cols())

    readr::write_csv(annotation.result.pos,
                     file.path(output.path, "intermediate_data", "ms2.match.annotation.result.pos.csv"))

    sample.pos <- annotation.result.pos[,match(sample.info[,1],
                                               colnames(annotation.result.pos))]

    sample.pos <- as.data.frame(sample.pos)
    rownames(sample.pos) <- paste(annotation.result.pos$name, "POS", sep = "_")

    rm(list="annotation.result.pos")
    gc()

    annotation.result.neg <-
      readr::read_csv(file.path(neg.path, "MS2_match_result_NEG",annotation.result.neg),
                      progress = FALSE,
                      col_types = readr::cols())

    readr::write_csv(annotation.result.neg,
                     file.path(output.path, "intermediate_data", "ms2.match.annotation.result.neg.csv"))

    sample.neg <- annotation.result.neg[,match(sample.info[,1],
                                               colnames(annotation.result.neg))]
    sample.neg <- as.data.frame(sample.neg)
    rownames(sample.neg) <- paste(annotation.result.neg$name, "NEG", sep = "_")
    rm(list = "annotation.result.neg")
    gc()

    ###check sample name in sample.pos and sample.neg
    if(any(sort(colnames(sample.pos)) != sort(colnames(sample.neg)))){
      stop("Sample names in annotation.result.pos and annotation.result.neg are not same.\n")
    }


    sample <- rbind(sample.pos, sample.neg)
    # rm(list = c("sample.pos", "sample.neg"))
    if(any(dir(file.path(output.path, "intermediate_data")) == "p.value")){
      load(file.path(output.path, "intermediate_data", "p.value"))
    }else{
      cat("Calculate p values of peaks.\n")
      p.value <- uniTest(
        sample = as.matrix(sample),
        sample.info = sample.info,
        uni.test = uni.test,
        correct = correct
      )
      names(p.value) <- annotation.result$name
    }


    num.peak <- sum(p.value <= p.cutoff)
    cat("There are ", num.peak, " out of ", nrow(sample), " peaks with p values less than ",
        p.cutoff, ".\n", sep = "")

    while(num.peak <= 30){
      p.cutoff <- readline("Do you want to change the p.cutoff or stop this analysis?
                           (type a new p.cutoff or type b to stop analysis)")
      if(p.cutoff == "b") stop()
      p.cutoff <- as.numeric(p.cutoff)
      num.peak <- sum(p.value <= p.cutoff)
      cat("There are ", num.peak, " out of ", nrow(sample), " peaks with p values less than ",
          p.cutoff, ".\n", sep = "")
    }

    p.value.pos <- p.value[grep("POS", names(p.value))]
    p.value.neg <- p.value[grep("NEG", names(p.value))]


    if(any(dir(file.path(output.path, "intermediate_data")) == "foldchange")){
      load(file.path(output.path, "intermediate_data", "foldchange"))
    }else{
      cat("\n")
      cat("Calculate fold change of peaks.\n")
      foldchange <- foldChange(sample = sample,
                               sample.info = sample.info,
                               by.what = "median", group = group)
      names(foldchange) <- annotation.result$name
    }

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

    adduct.table.pos <-
      adduct.table[adduct.table$mode == "positive", ]
    adduct.table.neg <-
      adduct.table[adduct.table$mode == "negative", ]

    rm(list = "adduct.table")
    gc()



    #-----------------------------------------------------------------------
    #begin find the pathway using pathway enrichment method.
    return.result <-
    findPathways2(tags2.pos = tags2.after.redundancy.remove.pos,
                  tags2.neg = tags2.after.redundancy.remove.neg,
                  metabolic.network = kegg.rpair2,
                  metabolite.pos = kegg.adduct.pos,
                  metabolite.neg = kegg.adduct.neg,
                  kegg.compound = kegg.compound,
                  p.value.pos = p.value.pos,
                  p.value.neg = p.value.neg,
                  p.cutoff = p.cutoff,
                  adduct.table.pos = adduct.table.pos,
                  adduct.table.neg = adduct.table.neg,
                  mz.tol = mz.tol, rt.tol = rt.tol2,
                  threads = threads,
                  path = path, output.path = output.path,
                  use.old.null = use.old.null,
                  species = species,
                  use.all.kegg.id = use.all.kegg.id,
                  only.mrn.annotation = only.mrn.annotation)

    return.result.pos <- return.result[[1]]
    return.result.neg <- return.result[[2]]

    new.tags.pos <- return.result.pos[[1]]
    anno.pos <- return.result.pos[[2]]
    anno.pos.from.pathway <- anno.pos
    save(anno.pos.from.pathway,
         file = file.path(output.path,
                          "intermediate_data",
                          "anno.pos.from.pathway"), compress = "xz")

    new.tags.neg <- return.result.neg[[1]]
    anno.neg <- return.result.neg[[2]]
    anno.neg.from.pathway <- anno.neg
    save(anno.neg.from.pathway,
         file = file.path(output.path,
                          "intermediate_data",
                          "anno.neg.from.pathway"), compress = "xz")

    rm(list = c("return.result","return.result.pos", "return.result.neg",
                "anno.pos", "anno.neg",
                "anno.pos.from.pathway", "anno.neg.from.pathway",
                "kegg.adduct.pos", "kegg.adduct.neg",
                "kegg.compound",
                "adduct.table.pos", "adduct.table.neg"))
    gc()

    #----------------------------------------------------------------------------
    ##process new.tags.pos and new.tags.neg, for the annotation from KEGG
    if(is.null(new.tags.pos)){
      tags2.after.kegg.matching.pos.from.pathway <- tags2.after.redundancy.remove.pos
      rm(list = c("tags2.after.redundancy.remove.pos"))
      gc()
      load(file.path(pos.path,
                     "MRN_annotation_result_POS",
                     "intermediate_data","tags.result"))
      tags.result.pos2.from.pathway <- tags.result
      tags.result.pos2.from.pathway$name <-
        paste(tags.result.pos2.from.pathway$name, "POS", sep = "_")
      rm(list = c("tags.result"))
      gc()
      save(tags2.after.kegg.matching.pos.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags2.after.kegg.matching.pos.from.pathway"),
           compress = "xz")
      save(tags.result.pos2.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags.result.pos2.from.pathway"), compress = "xz")
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
      tags.result.pos2.from.pathway <- temp.result.pos[[2]]

      tags2.after.kegg.matching.pos.from.pathway <- new.tags.pos
      save(tags2.after.kegg.matching.pos.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags2.after.kegg.matching.pos.from.pathway"),
           compress = "xz")
      save(tags.result.pos2.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags.result.pos2.from.pathway"), compress = "xz")
    }


    ##change tags2.after.kegg.matching.pos to csv like ms2Annotation
    cat("\n")
    load(file.path(pos.path, "Dysregulated_network_analysis_result_POS", "intermediate_data", "p.value"))
    load(file.path(pos.path, "Dysregulated_network_analysis_result_POS", "intermediate_data", "foldchange"))
    data("kegg.compound", envir = environment())

    temp.pos <- getAnnotationResult(tags2 = tags2.after.kegg.matching.pos.from.pathway,
                                    p.value = p.value,
                                    correct = correct,
                                    foldchange = foldchange,
                                    tags.result2 = tags.result.pos2.from.pathway,
                                    kegg.compound = kegg.compound)

    cat("\n")
    readr::write_csv(temp.pos,
                     file.path(output.path,
                               "DNA.annotation.result.pos.from.pahtway.csv"))
    rm(list = "temp.pos")
    gc()


    rm(list = c("tags2.after.kegg.matching.pos.from.pathway",
                "temp.result.pos", "new.tags.pos",
                "tags.result.pos2.from.pathway"))
    gc()


    ##negative
    if(is.null(new.tags.neg)){
      tags2.after.kegg.matching.neg.from.pathway <- tags2.after.redundancy.remove.neg
      rm(list = c("tags2.after.redundancy.remove.neg"))
      gc()
      load(file.path(neg.path,
                     "MRN_annotation_result_NEG",
                     "intermediate_data","tags.result"))
      tags.result.neg2.from.pathway <- tags.result
      tags.result.neg2.from.pathway$name <-
        paste(tags.result.neg2.from.pathway$name, "NEG", sep = "_")
      rm(list = c("tags.result"))
      gc()
      save(tags2.after.kegg.matching.neg.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags2.after.kegg.matching.neg.from.pathway"),
           compress = "xz")
      save(tags.result.neg2.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags.result.neg2.from.pathway"), compress = "xz")
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
      tags.result.neg2.from.pathway <- temp.result.neg[[2]]

      tags2.after.kegg.matching.neg.from.pathway <- new.tags.neg
      save(tags2.after.kegg.matching.neg.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags2.after.kegg.matching.neg.from.pathway"),
           compress = "xz")
      save(tags.result.neg2.from.pathway,
           file = file.path(output.path,
                            "intermediate_data",
                            "tags.result.neg2.from.pathway"), compress = "xz")
    }

    ##change tags2.after.kegg.matching.neg to csv like ms2Annotation
    cat("\n")
    load(file.path(neg.path, "Dysregulated_network_analysis_result_NEG", "intermediate_data", "p.value"))
    load(file.path(neg.path, "Dysregulated_network_analysis_result_NEG", "intermediate_data", "foldchange"))
    data("kegg.compound", envir = environment())

    temp.neg <- getAnnotationResult(tags2 = tags2.after.kegg.matching.neg.from.pathway,
                                    p.value = p.value,
                                    correct = correct,
                                    foldchange = foldchange,
                                    tags.result2 = tags.result.neg2.from.pathway,
                                    kegg.compound = kegg.compound)

    cat("\n")
    readr::write_csv(temp.neg, file.path(output.path,
                                         "DNA.annotation.result.neg.from.pathway.csv"))
    rm(list = "temp.neg")
    gc()


    rm(list = c("tags2.after.kegg.matching.neg.from.pathway",
                "temp.result.neg", "new.tags.neg",
                "tags.result.neg2.from.pathway"))
    gc()



    ###transform metabolite information to pathway/Pathway information
    singleTransform2(annotation.result.pos = "ms2.match.annotation.result.pos.csv",
                     annotation.result.neg = "ms2.match.annotation.result.neg.csv",
                     data.type = "metabolomics",
                     sample.info = "sample.info.csv",
                     group = group, trans.to = "pathway",
                     scale = TRUE,
                     scale.method = "pareto",
                     species = species,
                     which.peak = "intensity",
                     column = column, method = "mean",
                     path = file.path(output.path, "intermediate_data"),
                     output.path = file.path(output.path, "quantitative_information"),
                     metabolic.network = kegg.rpair2)

    load(file.path(output.path, "quantitative_information", "node.quantitative.data"))
    node.quantitative.data2 <- node.quantitative.data
    save(node.quantitative.data2, file = file.path(output.path,
                                                   "quantitative_information",
                                                   "node.quantitative.data2"), compress = "xz")
    readr::write_csv(node.quantitative.data2, file.path(output.path,
                                                        "quantitative_information",
                                                        "node.quantitative.data2.csv"))

    unlink(file.path(output.path, "quantitative_information",
                     "node.quantitative.data"), recursive = TRUE)
    unlink(file.path(output.path,  "quantitative_information",
                     "node.quantitative.data.csv"), recursive = TRUE)


    ###use node quantitative data2 to change pathway.result
    ##add the detailed information to mse.result
    switch(species,
           "hsa" = {data("hsa.kegg.pathway", envir = environment())
             M = hsa.kegg.pathway},
           "dme" = {data("dme.kegg.pathway", envir = environment())
             M = dme.kegg.pathway},
           "mmu" = {data("mmu.kegg.pathway", envir = environment())
             M = mmu.kegg.pathway},
           "rat" = {data("rat.kegg.pathway", envir = environment())
             M = rat.kegg.pathway},
           "bta" = {data("bta.kegg.pathway", envir = environment())
             M = bta.kegg.pathway},
           "gga" = {data("gga.kegg.pathway", envir = environment())
             M = gga.kegg.pathway},
           "dre" = {data("dre.kegg.pathway", envir = environment())
             M = dre.kegg.pathway},
           "cel" = {data("cel.kegg.pathway", envir = environment())
             M = cel.kegg.pathway},
           "sce" = {data("sce.kegg.pathway", envir = environment())
             M = sce.kegg.pathway},
           "ath" = {data("ath.kegg.pathway", envir = environment())
             M = ath.kegg.pathway},
           "smm" = {data("smm.kegg.pathway", envir = environment())
             M = smm.kegg.pathway},
           "pfa" = {data("pfa.kegg.pathway", envir = environment())
             M = pfa.kegg.pathway},
           "tbr" = {data("tbr.kegg.pathway", envir = environment())
             M = tbr.kegg.pathway},
           "eco" = {data("eco.kegg.pathway", envir = environment())
             M = eco.kegg.pathway},
           "ppu" = {data("ppu.kegg.pathway", envir = environment())
             M = ppu.kegg.pathway},
           "syf" = {data("syf.kegg.pathway", envir = environment())
             M = syf.kegg.pathway}
    )


    ##--------------------------------------------------------------------------
    #----------------------------------------------------------------------------

    pathway.data <- M
    rm(list = "M")

    overlap <- unname(unlist(lapply(pathway.data, function(x){
      length(intersect(x, node.quantitative.data$to))
    })))

    pathway.length <- unname(unlist(lapply(pathway.data, length)))

    detected.metabolite.id <- unname(unlist(lapply(pathway.data, function(x){
      paste(intersect(x, node.quantitative.data$to), collapse = ";")
    })))

    pathway.name <- unlist(lapply(strsplit(names(pathway.data), split = ";"), function(x) x[1]))
    pathway.id <- unlist(lapply(strsplit(names(pathway.data), split = ";"), function(x) x[2]))

    load(file.path(output.path, "intermediate_data", "pathway.result"))
    pathway.result2 <- data.frame(pathway.name, pathway.id, pathway.length, overlap, detected.metabolite.id,
                                  stringsAsFactors = FALSE)

    pathway.result2 <- pathway.result2[match(pathway.result$pathway.id, pathway.result2$pathway.id),]

    save(pathway.result2, file = file.path(output.path, "intermediate_data", "pathway.result2"),
         compress = "xz")

    write.csv(pathway.result2, file.path(output.path, "pathway_information", "pathway.result2.csv"), row.names = FALSE)



    temp.path2 <- file.path(output.path, "pathway_information", "boxplot")
    dir.create(temp.path2)

    load(file.path(output.path,  "quantitative_information", "pathway.quantitative.data"))

    ###only select the pathway which at least have more than 3 detected metabolites
    temp.idx <- match(paste(pathway.result2$pathway.name[which(pathway.result2$overlap >= 3)],
                            pathway.result2$pathway.id[which(pathway.result2$overlap >= 3)], sep = ";"),
                      pathway.quantitative.data[,1])

    pathway.quantitative.data <- pathway.quantitative.data[temp.idx,]

    pathway.name.id <- pathway.quantitative.data[,1]
    pathway.name <- unname(unlist(lapply(pathway.name.id, function(x) strsplit(x, split = ";")[[1]][1])))
    pathway.id <- unname(unlist(lapply(pathway.name.id, function(x) strsplit(x, split = ";")[[1]][2])))
    pathway.sample <- pathway.quantitative.data[,-1]
    rownames(pathway.sample) <- pathway.name

    cat("\n")
    cat("Calculate p values of pathways.\n")
    pathway.p.value <- uniTest(sample = pathway.sample,
                               sample.info = sample.info, uni.test = "t",
                               correct = TRUE)

    cat("\n")
    cat("Calculate fold changes of pathways.\n")
    pathway.fc <- foldChange(sample = pathway.sample,
                             sample.info = sample.info,
                             group = group)


    samplePlot(sample = pathway.sample, sample.info = sample.info, file.name = pathway.id,
               group = group, output.path = temp.path2, p.value = pathway.p.value, fc = pathway.fc,
               beeswarm = TRUE)

    pathway.sample1 <- t(apply(pathway.sample, 1, as.numeric))
    colnames(pathway.sample1) <- colnames(pathway.sample)
    rownames(pathway.sample1) <- rownames(pathway.sample)

    if(nrow(pathway.sample1) >= 2){
      pdf(file.path(output.path, "pathway_information", "pathway.heatmap.pdf"),
          width = 7, height = 7)
      par(mar = c(5, 5 ,4, 2))
      heatMap(sample = pathway.sample1, sample.info = sample.info, group = group)
      dev.off()
    }

    cat("\n")
    cat("metPathway2 is done\n")
  })






####--------------------------------------------------------------------------

# title findPathways
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


setGeneric(name = "findPathways2",
           def = function(tags2.pos,
                          tags2.neg,
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
                          path = ".",
                          output.path = ".",
                          use.old.null = FALSE,
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          use.all.kegg.id = FALSE,
                          only.mrn.annotation = FALSE){
             if(sum(p.value.pos < p.cutoff) == 0 & sum(p.value.neg < p.cutoff) == 0) {
               stop("No significantly dysregulated peaks.\n")
             }

             species <- match.arg(species)

             cat("\n")
             cat("Positive.\n")
             result.pos <- tags2Result(tags2 = tags2.pos, score.cutoff = 0)

             cat("\n")
             cat("Negative.\n")
             result.neg <- tags2Result(tags2 = tags2.neg, score.cutoff = 0)

             ##find the significant changed peaks
             if(only.mrn.annotation){
               index.pos <- which(p.value.pos < p.cutoff & !is.na(result.pos$to))
               index.neg <- which(p.value.neg < p.cutoff & !is.na(result.neg$to))
               if(length(index.pos) ==0 & length(index.neg) == 0){
                 warning("The peaks with p value less than 0.05 have no MRN based annotations.\nSo only.mrn.annotation has been set as FALSE.")
                 index.pos <- which(p.value.pos < p.cutoff)
                 index.neg <- which(p.value.neg < p.cutoff)
                 only.mrn.annotation <- FALSE
               }
             }else{
               index.pos <- which(p.value.pos < p.cutoff)
               index.neg <- which(p.value.neg < p.cutoff)
             }

             ##size is the number of dysregulated peaks
             size.pos <- length(index.pos)
             size.neg <- length(index.neg)

             result.pos1 <- result.pos[index.pos,]
             result.neg1 <- result.neg[index.neg,]
             ##remove isotope annotation form result
             result.pos1 <- removeIsotopeFromResult(result = result.pos1)
             result.neg1 <- removeIsotopeFromResult(result = result.neg1)

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


             ####for negative mode
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

             #unique.id is the annotations of dysregulated peaks
             anno <- c(anno.pos, anno.neg)
             unique.id <- unique(unname(unlist(anno)))
             unique.id <- unique.id[which(unique.id %in% igraph::V(metabolic.network)$name)]

             cat("\n")
             cat("Identify initial pathways.\n")
             mse.result <- mseAnalysis(metabolite.id = unique.id, species = species)

             ##get NULL distribution
             cat("\n")
             cat("Get pseudo pathways.\n")

             if(!use.old.null){
               ##temp.fun is function for multiple threads calculation
               temp.fun <- function(idx, result.pos, result.neg,
                                    size.pos,  size.neg,
                                    only.mrn.annotation,
                                    metabolite.pos, metabolite.neg,
                                    mz.tol, rt.tol,
                                    use.all.kegg.id, species,
                                    keggMatch,
                                    mseAnalysis1,
                                    M = M,
                                    removeIsotopeFromResult){
                 temp.idx.pos <- sample(1:nrow(result.pos), size.pos)
                 temp.idx.neg <- sample(1:nrow(result.neg), size.neg)
                 if(only.mrn.annotation){
                   temp.idx.pos <- temp.idx.pos[which(!is.na(result.pos$to[temp.idx.pos]))]
                   temp.idx.neg <- temp.idx.neg[which(!is.na(result.neg$to[temp.idx.neg]))]
                 }

                 temp.result.pos1 <- result.pos[temp.idx.pos,]
                 temp.result.neg1 <- result.neg[temp.idx.neg,]

                 temp.result.pos1 <- removeIsotopeFromResult(result = temp.result.pos1)
                 temp.result.neg1 <- removeIsotopeFromResult(result = temp.result.neg1)

                 temp.to.pos1 <- temp.result.pos1$to
                 temp.to.neg1 <- temp.result.neg1$to
                 names(temp.to.pos1) <- result.pos$name[temp.idx.pos]
                 names(temp.to.neg1) <- result.neg$name[temp.idx.neg]

                 temp.to.pos1 <- lapply(temp.to.pos1, function(x) {
                   if(is.na(x)) {return(x)
                   }else{
                     x <- strsplit(x, split = ";")[[1]]
                     x <- x[!duplicated(x)]
                     x <- paste(x, collapse = ";")
                     x
                   }
                 })
                 names(temp.to.pos1) <- temp.result.pos1$name
                 temp.anno.pos1 <- temp.to.pos1

                 temp.to.neg1 <- lapply(temp.to.neg1, function(x) {
                   if(is.na(x)) {return(x)
                   }else{
                     x <- strsplit(x, split = ";")[[1]]
                     x <- x[!duplicated(x)]
                     x <- paste(x, collapse = ";")
                     x
                   }
                 })
                 names(temp.to.neg1) <- temp.result.neg1$name
                 temp.anno.neg1 <- temp.to.neg1

                 ###kegg mz and rt matching
                 temp.mz.pos <- as.numeric(temp.result.pos1$mz)
                 temp.rt.pos <- as.numeric(temp.result.pos1$rt)

                 temp.mz.neg <- as.numeric(temp.result.neg1$mz)
                 temp.rt.neg <- as.numeric(temp.result.neg1$rt)


                 temp.match.result.pos <- keggMatch(mz = temp.mz.pos,
                                                    rt = temp.rt.pos,
                                                    metabolite = metabolite.pos,
                                                    mz.tol = mz.tol,rt.tol = rt.tol,
                                                    polarity = "positive")

                 temp.match.result.neg <- keggMatch(mz = temp.mz.neg, rt = temp.rt.neg,
                                                    metabolite = metabolite.neg,
                                                    mz.tol = mz.tol,rt.tol = rt.tol,
                                                    polarity = "negative")


                 ##temp.anno.pos2 is the annotation results from KEGG according to mz and rt
                 temp.anno.pos2 <- lapply(temp.match.result.pos, function(x) {
                   if(is.null(x)) return(NA)
                   temp <- unique(x$ID)
                   temp <- paste(temp, collapse = ";")
                   temp
                 })

                 names(temp.anno.pos2) <- temp.result.pos1$name

                 if(use.all.kegg.id){
                   temp.anno.pos <- mapply(function(x,y){
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
                   x = temp.anno.pos1,
                   y = temp.anno.pos2)
                 }else{
                   temp.anno.pos <- mapply(function(x,y){
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
                   x = temp.anno.pos1,
                   y = temp.anno.pos2)
                 }

                 names(temp.anno.pos) <- temp.result.pos1$name

                 ##temp.anno.neg2 is the annotation results from KEGG according to mz and rt
                 temp.anno.neg2 <- lapply(temp.match.result.neg, function(x) {
                   if(is.null(x)) return(NA)
                   temp <- unique(x$ID)
                   temp <- paste(temp, collapse = ";")
                   temp
                 })

                 names(temp.anno.neg2) <- temp.result.neg1$name

                 if(use.all.kegg.id){
                   temp.anno.neg <- mapply(function(x,y){
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
                   x = temp.anno.neg1,
                   y = temp.anno.neg2)
                 }else{
                   temp.anno.neg <- mapply(function(x,y){
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
                   x = temp.anno.neg1,
                   y = temp.anno.neg2)
                 }

                 names(temp.anno.neg) <- temp.result.neg1$name



                 temp.index.pos <- unname(which(unlist(lapply(temp.anno.pos, function(x) {!all(is.na(x))}))))
                 temp.result.pos2 <- temp.result.pos1[temp.index.pos,]
                 rm(temp.result.pos1)
                 gc()
                 temp.anno.pos <- temp.anno.pos[temp.index.pos]
                 temp.anno.pos1 <- temp.anno.pos1[temp.index.pos]
                 temp.anno.pos2 <- temp.anno.pos2[temp.index.pos]

                 temp.match.result.pos2 <- temp.match.result.pos[temp.index.pos]
                 rm(temp.match.result.pos)
                 gc()


                 temp.index.neg <- unname(which(unlist(lapply(temp.anno.neg, function(x) {!all(is.na(x))}))))
                 temp.result.neg2 <- temp.result.neg1[temp.index.neg,]
                 rm(temp.result.neg1)
                 gc()
                 temp.anno.neg <- temp.anno.neg[temp.index.neg]
                 temp.anno.neg1 <- temp.anno.neg1[temp.index.neg]
                 temp.anno.neg2 <- temp.anno.neg2[temp.index.neg]

                 temp.match.result.neg2 <- temp.match.result.neg[temp.index.neg]
                 rm(temp.match.result.neg)
                 gc()

                 #unique.id is the annotations of significant peaks
                 temp.unique.id <- unique(unname(unlist(c(temp.anno.pos, temp.anno.neg))))
                 temp.mse.result <- mseAnalysis1(metabolite.id = temp.unique.id, M = M)
                 temp.mse.result

               }

               #######load pathway data
               switch(species,
                      "hsa" = {data("hsa.kegg.pathway", envir = environment())
                        M = hsa.kegg.pathway},
                      "dme" = {data("dme.kegg.pathway", envir = environment())
                        M = dme.kegg.pathway},
                      "mmu" = {data("mmu.kegg.pathway", envir = environment())
                        M = mmu.kegg.pathway},
                      "rat" = {data("rat.kegg.pathway", envir = environment())
                        M = rat.kegg.pathway},
                      "bta" = {data("bta.kegg.pathway", envir = environment())
                        M = bta.kegg.pathway},
                      "gga" = {data("gga.kegg.pathway", envir = environment())
                        M = gga.kegg.pathway},
                      "dre" = {data("dre.kegg.pathway", envir = environment())
                        M = dre.kegg.pathway},
                      "cel" = {data("cel.kegg.pathway", envir = environment())
                        M = cel.kegg.pathway},
                      "sce" = {data("sce.kegg.pathway", envir = environment())
                        M = sce.kegg.pathway},
                      "ath" = {data("ath.kegg.pathway", envir = environment())
                        M = ath.kegg.pathway},
                      "smm" = {data("smm.kegg.pathway", envir = environment())
                        M = smm.kegg.pathway},
                      "pfa" = {data("pfa.kegg.pathway", envir = environment())
                        M = pfa.kegg.pathway},
                      "tbr" = {data("tbr.kegg.pathway", envir = environment())
                        M = tbr.kegg.pathway},
                      "eco" = {data("eco.kegg.pathway", envir = environment())
                        M = eco.kegg.pathway},
                      "ppu" = {data("ppu.kegg.pathway", envir = environment())
                        M = ppu.kegg.pathway},
                      "syf" = {data("syf.kegg.pathway", envir = environment())
                        M = syf.kegg.pathway}
               )

               temp.mse.result <- BiocParallel::bplapply(1:100, temp.fun,
                                                         BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                                           progressbar = TRUE),
                                                         result.pos = result.pos,
                                                         result.neg = result.neg,
                                                         size.pos = size.pos,
                                                         size.neg = size.neg,
                                                         only.mrn.annotation = only.mrn.annotation,
                                                         metabolite.pos = metabolite.pos,
                                                         metabolite.neg = metabolite.neg,
                                                         mz.tol = mz.tol,
                                                         rt.tol = rt.tol,
                                                         use.all.kegg.id = use.all.kegg.id,
                                                         species = species,
                                                         keggMatch = keggMatch,
                                                         mseAnalysis1 = mseAnalysis1,
                                                         M = M,
                                                         removeIsotopeFromResult)
               save(temp.mse.result, file = file.path(output.path, "intermediate_data","temp.mse.result"), compress = "xz")
             }else{
               cat("\n")
               cat("Use old pseudo scores.\n")
               cat("\n")
               load(file.path(output.path, "intermediate_data","temp.mse.result"))
             }

             ##mse.result is from mse analysis using annotation of dysregulated peaks
             mse.result <- mse.result[order(rownames(mse.result)),]

             #temp.mse.result is from mse analyssis using annotation of
             #reference peaks
             temp.mse.result <- lapply(temp.mse.result, function(x){
               x <- x[order(rownames(x)),]
               x
             })

             p.overlap <- lapply(1:nrow(mse.result), function(idx){
               temp <- lapply(temp.mse.result, function(x){
                 as.numeric(x[idx,c(1,4,5)])
               })
               temp <- do.call(rbind, temp)
               temp.p <- temp[,1]
               temp.overlap <- temp[,3]/temp[,2]
               return.result <- list(temp.p, temp.overlap)
               names(return.result) <- c("p", "overlap")
               return.result
             })

             refer.p <- lapply(p.overlap, function(x) x[[1]])
             refer.log.p <- lapply(refer.p, function(x) -log(x,10))

             refer.overlap <- lapply(p.overlap, function(x) x[[2]])

             names(refer.log.p) <- names(refer.overlap) <- rownames(mse.result)


             p.adjust.p <- mapply(function(ref.p, p){
               random.data <- runif(n = 101, min = 0, max = sd(ref.p)/2)
               check <- tryCatch({
                 para <- MASS::fitdistr(x = ref.p + random.data[1:100], densfun = "gamma")[[1]]
               }, error = function(e){
                 NA
               })
               if(is.na(check)) return(1)
               para <- MASS::fitdistr(x = ref.p + random.data[1:100], densfun = "gamma")[[1]]
               standard.distribution <- rgamma(n = 100000, shape = para[1], rate = para[2])
               cdf <- ecdf(x = standard.distribution)
               adjust.p <- 1 - cdf(p+random.data[101])
             },
             ref.p = refer.log.p,
             p = -log(mse.result[,1],10))



             p.adjust.overlap <- mapply(function(ref.overlap, overlap){
               random.data <- runif(n = 101, min = 0, max = sd(ref.overlap)/2)
               check <- tryCatch({
                 para <- MASS::fitdistr(x = ref.overlap + random.data[1:100], densfun = "gamma")[[1]]
               }, error = function(e){
                 NA
               })
               if(is.na(check)) return(1)
               para <- MASS::fitdistr(x = ref.overlap + random.data[1:100], densfun = "gamma")[[1]]
               standard.distribution <- rgamma(n = 100000, shape = para[1], rate = para[2])
               cdf <- ecdf(x = standard.distribution)
               adjust.p <- 1 - cdf(overlap+random.data[101])
             },
             ref.overlap = refer.overlap,
             overlap = mse.result[,5]/mse.result[,4])



             mse.result <- data.frame(mse.result, p.adjust.p, p.adjust.overlap,
                                      stringsAsFactors = FALSE)

             colnames(mse.result)[c(6,7)] <- c("adjusted.p.value.by.pvalue", "adjusted.p.value.by.overlap")
             mse.result <- mse.result[order(mse.result$adjusted.p.value.by.pvalue),]

             ##add the detailed information to mse.result
             switch(species,
                    "hsa" = {data("hsa.kegg.pathway", envir = environment())
                      M = hsa.kegg.pathway},
                    "dme" = {data("dme.kegg.pathway", envir = environment())
                      M = dme.kegg.pathway},
                    "mmu" = {data("mmu.kegg.pathway", envir = environment())
                      M = mmu.kegg.pathway},
                    "rat" = {data("rat.kegg.pathway", envir = environment())
                      M = rat.kegg.pathway},
                    "bta" = {data("bta.kegg.pathway", envir = environment())
                      M = bta.kegg.pathway},
                    "gga" = {data("gga.kegg.pathway", envir = environment())
                      M = gga.kegg.pathway},
                    "dre" = {data("dre.kegg.pathway", envir = environment())
                      M = dre.kegg.pathway},
                    "cel" = {data("cel.kegg.pathway", envir = environment())
                      M = cel.kegg.pathway},
                    "sce" = {data("sce.kegg.pathway", envir = environment())
                      M = sce.kegg.pathway},
                    "ath" = {data("ath.kegg.pathway", envir = environment())
                      M = ath.kegg.pathway},
                    "smm" = {data("smm.kegg.pathway", envir = environment())
                      M = smm.kegg.pathway},
                    "pfa" = {data("pfa.kegg.pathway", envir = environment())
                      M = pfa.kegg.pathway},
                    "tbr" = {data("tbr.kegg.pathway", envir = environment())
                      M = tbr.kegg.pathway},
                    "eco" = {data("eco.kegg.pathway", envir = environment())
                      M = eco.kegg.pathway},
                    "ppu" = {data("ppu.kegg.pathway", envir = environment())
                      M = ppu.kegg.pathway},
                    "syf" = {data("syf.kegg.pathway", envir = environment())
                      M = syf.kegg.pathway}
             )

             data("kegg.compound", envir = environment())

             pathway.data <- M
             rm(M)
             gc()

             pathway.name.id <- rownames(mse.result)
             pathway.name <- unlist(lapply(strsplit(pathway.name.id, split = ";"), function(x) x[1]))
             pathway.id <- unlist(lapply(strsplit(pathway.name.id, split = ";"), function(x) x[2]))

             temp.info <- lapply(pathway.name.id, function(x){
               temp.idx <- match(x, names(pathway.data))
               detected.id <- intersect(pathway.data[[temp.idx]], unique.id)
               detected.name <- kegg.compound$Name[match(detected.id, kegg.compound$ID)]
               detected.name <- unlist(lapply(strsplit(detected.name, split = ";"), function(x) x[1]))
               hidden.id <- setdiff(pathway.data[[temp.idx]], detected.id)
               hidden.name <- kegg.compound$Name[match(hidden.id, kegg.compound$ID)]
               hidden.name <- unlist(lapply(strsplit(hidden.name, split = ";"), function(x) x[1]))
               detected.id <- paste(detected.id, collapse = ";")
               detected.name <- paste(detected.name, collapse = ";")
               hidden.id <- paste(hidden.id, collapse = ";")
               hidden.name <- paste(detected.name, collapse = ";")
               return(c(detected.id, detected.name, hidden.id, hidden.name))
             })

             temp.info <- do.call(rbind, temp.info)
             colnames(temp.info) <- c("detected.metabolite.id", "detected.metabolite.name", "hidden.metabolite.id", "hidden.metabolite.name")
             mse.result <- data.frame(pathway.name, pathway.id, mse.result, temp.info,
                                      stringsAsFactors = FALSE)
             write.csv(mse.result, file.path(output.path, "pathway_information","pathway.result.csv"), row.names = FALSE)
             write.csv(mse.result, file.path(output.path, "pathway_information","pathway.result.csv"), row.names = FALSE)
             pathway.result <- mse.result
             save(pathway.result, file = file.path(output.path, "intermediate_data","pathway.result"), compress = "xz")

             rm(list = c("pathway.data", "pathway.result"))
             gc()

             ##filter annotation accoring to enrichment pathways
             index <- which(mse.result$adjusted.p.value.by.overlap < 0.05 | mse.result$adjusted.p.value.by.pvalue < 0.05)

             if(length(index) == 0){
               cat("There are no intial pathways with p values less than 0.05.\n")
             }else{
               cat("Filter annotations according to dysregulated pathways.\n")

               annotation.pos1 <- vector(mode = "list", length = length(index))
               annotation.pos2 <- vector(mode = "list", length = length(index))

               annotation.neg1 <- vector(mode = "list", length = length(index))
               annotation.neg2 <- vector(mode = "list", length = length(index))

               for(i in index){
                 # positive
                 temp.pathway <- mse.result$detected.metabolite.id[i]
                 temp.pathway <- strsplit(temp.pathway, split = ";")[[1]]

                 anno1.result.pos <- lapply(anno.pos1, function(x){
                   if(is.na(x)){return(NA)}
                   x <- strsplit(x, split = ";")[[1]]
                   idx <- which(x %in% temp.pathway)
                   if(length(idx) == 0) return(NA)
                   return(paste(x[idx], collapse = ";"))
                 })
                 anno1.result.pos <- unlist(anno1.result.pos)

                 anno2.result.pos <- lapply(anno.pos2, function(x){
                   if(is.na(x)){return(NA)}
                   x <- strsplit(x, split = ";")[[1]]
                   idx <- which(x %in% temp.pathway)
                   if(length(idx) == 0) return(NA)
                   return(paste(x[idx], collapse = ";"))
                 })

                 anno2.result.pos <- unlist(anno2.result.pos)

                 annotation.pos1[[i]] <- anno1.result.pos
                 annotation.pos2[[i]] <- anno2.result.pos

                 # negative
                 temp.pathway <- mse.result$detected.metabolite.id[i]
                 temp.pathway <- strsplit(temp.pathway, split = ";")[[1]]

                 anno1.result.neg <- lapply(anno.neg1, function(x){
                   if(is.na(x)){return(NA)}
                   x <- strsplit(x, split = ";")[[1]]
                   idx <- which(x %in% temp.pathway)
                   if(length(idx) == 0) return(NA)
                   return(paste(x[idx], collapse = ";"))
                 })
                 anno1.result.neg <- unlist(anno1.result.neg)

                 anno2.result.neg <- lapply(anno.neg2, function(x){
                   if(is.na(x)){return(NA)}
                   x <- strsplit(x, split = ";")[[1]]
                   idx <- which(x %in% temp.pathway)
                   if(length(idx) == 0) return(NA)
                   return(paste(x[idx], collapse = ";"))
                 })

                 anno2.result.neg <- unlist(anno2.result.neg)

                 annotation.neg1[[i]] <- anno1.result.neg
                 annotation.neg2[[i]] <- anno2.result.neg
               }


               ###positive
               annotation.pos1 <- do.call(cbind, annotation.pos1)
               annotation.pos2 <- do.call(cbind, annotation.pos2)

               peak.name.pos <- names(anno.pos1)

               annotation.pos1 <- data.frame(peak.name.pos,
                                             unlist(anno.pos1), annotation.pos1,
                                             stringsAsFactors = FALSE)
               annotation.pos2 <- data.frame(peak.name.pos,
                                             unlist(anno.pos2), annotation.pos2,
                                             stringsAsFactors = FALSE)

               perfect.annotation.pos1 <- apply(annotation.pos1, 1, function(x){
                 x <- as.character(x)
                 if(is.na(x[2])) return(NA)
                 temp.idx <- which(!is.na(x[-c(1:2)]))
                 if(length(temp.idx) == 0) return(x[2])
                 if(length(temp.idx) == 1) return(x[-c(1:2)][temp.idx])
                 return(paste(x[-c(1:2)][temp.idx], collapse = ";"))
               })

               perfect.annotation.pos2 <- apply(annotation.pos2, 1, function(x){
                 x <- as.character(x)
                 if(is.na(x[2])) return(NA)
                 temp.idx <- which(!is.na(x[-c(1:2)]))
                 if(length(temp.idx) == 0) return(NA)
                 return(paste(x[-c(1:2)][temp.idx], collapse = ";"))
               })


               colnames(annotation.pos1) <- c("Peak.name","MRN.annotation",
                                              mse.result$pathway.name[index])
               colnames(annotation.pos2) <- c("Peak.name","KEGG.annotation",
                                              mse.result$pathway.name[index])

               annotation.pos1 <- data.frame(annotation.pos1[,1:2], perfect.annotation.pos1,
                                             annotation.pos1[,-c(1:2)], stringsAsFactors = FALSE)
               annotation.pos2 <- data.frame(annotation.pos2[,1:2], perfect.annotation.pos2,
                                             annotation.pos2[,-c(1:2)], stringsAsFactors = FALSE)

               colnames(annotation.pos1)[3] <- colnames(annotation.pos2)[3] <- "Annotation"

               write.csv(annotation.pos1,
                         file.path(output.path, "intermediate_data",
                                   'Positive.significant.peaks.annotation.MRN.from.pathway.csv'),
                         row.names = FALSE)
               write.csv(annotation.pos2,
                         file.path(output.path, "intermediate_data",
                                   'Positive.significant.peaks.annotation.KEGG.from.pathway.csv'),
                         row.names = FALSE)


               ###negative
               annotation.neg1 <- do.call(cbind, annotation.neg1)
               annotation.neg2 <- do.call(cbind, annotation.neg2)

               peak.name.neg <- names(anno.neg1)

               annotation.neg1 <- data.frame(peak.name.neg,
                                             unlist(anno.neg1), annotation.neg1,
                                             stringsAsFactors = FALSE)
               annotation.neg2 <- data.frame(peak.name.neg,
                                             unlist(anno.neg2), annotation.neg2,
                                             stringsAsFactors = FALSE)

               perfect.annotation.neg1 <- apply(annotation.neg1, 1, function(x){
                 x <- as.character(x)
                 if(is.na(x[2])) return(NA)
                 temp.idx <- which(!is.na(x[-c(1:2)]))
                 if(length(temp.idx) == 0) return(x[2])
                 if(length(temp.idx) == 1) return(x[-c(1:2)][temp.idx])
                 return(paste(x[-c(1:2)][temp.idx], collapse = ";"))
               })


               perfect.annotation.neg2 <- apply(annotation.neg2, 1, function(x){
                 x <- as.character(x)
                 if(is.na(x[2])) return(NA)
                 temp.idx <- which(!is.na(x[-c(1:2)]))
                 if(length(temp.idx) == 0) return(NA)
                 return(paste(x[-c(1:2)][temp.idx], collapse = ";"))
               })


               colnames(annotation.neg1) <- c("Peak.name","MRN.annotation",
                                              mse.result$pathway.name[index])
               colnames(annotation.neg2) <- c("Peak.name","KEGG.annotation",
                                              mse.result$pathway.name[index])

               annotation.neg1 <- data.frame(annotation.neg1[,1:2], perfect.annotation.neg1,
                                             annotation.neg1[,-c(1:2)], stringsAsFactors = FALSE)
               annotation.neg2 <- data.frame(annotation.neg2[,1:2], perfect.annotation.neg2,
                                             annotation.neg2[,-c(1:2)], stringsAsFactors = FALSE)

               colnames(annotation.neg1)[3] <- colnames(annotation.neg2)[3] <- "Annotation"

               write.csv(annotation.neg1,
                         file.path(output.path, "intermediate_data",
                                   'Negative.significant.peaks.annotation.MRN.from.pathway.csv'),
                         row.names = FALSE)
               write.csv(annotation.neg2,
                         file.path(output.path, "intermediate_data",
                                   'Negative.significant.peaks.annotation.KEGG.from.pathway.csv'),
                         row.names = FALSE)


               ##new anno should be from prefer.annotation1 and prefer.annotation2
               ##new anno should be from prefer.annotation1 and prefer.annotation2
               anno.pos <- mapply(function(x,y){
                 if(!is.na(x)) return(list(x))
                 return(list(y))
               },
               x = perfect.annotation.pos1,
               y = perfect.annotation.pos2)

               anno.pos <- lapply(anno.pos, function(x){
                 if(is.na(x)) return(NA)
                 strsplit(x, split = ";")[[1]]
               })

               anno.neg <- mapply(function(x,y){
                 if(!is.na(x)) return(list(x))
                 return(list(y))
               },
               x = perfect.annotation.neg1,
               y = perfect.annotation.neg2)

               anno.neg <- lapply(anno.neg, function(x){
                 if(is.na(x)) return(NA)
                 strsplit(x, split = ";")[[1]]
               })


               ##add those result to tags2 positive
               peak.name.pos <- showTags2(tags2.pos, slot = "name")
               index.pos <- match(annotation.pos1$Peak.name, peak.name.pos)

               for(i in 1:length(index.pos)){
                 temp.anno <- annotation.pos1[i,3]
                 if(is.na(temp.anno)) next
                 temp.anno <- strsplit(temp.anno, split = ";")[[1]]
                 temp.tags <- tags2.pos[[index.pos[i]]]
                 temp.tags.annotation <- temp.tags@annotation
                 temp.tags.to <- unlist(lapply(temp.tags.annotation, function(x) x$to))
                 temp.idx <- which(temp.tags.to %in% temp.anno)
                 temp.tags.annotation <- temp.tags.annotation[temp.idx]
                 tags2.pos[[index.pos[i]]]@annotation <- temp.tags.annotation
                 rm(list = c("temp.anno", "temp.tags", "temp.tags.annotation", "temp.tags.to",
                             "temp.idx", "temp.tags.annotation"))
                 gc()
               }

               ##add those result to tags2 negative
               peak.name.neg <- showTags2(tags2.neg, slot = "name")
               index.neg <- match(annotation.neg1$Peak.name, peak.name.neg)

               for(i in 1:length(index.neg)){
                 temp.anno <- annotation.neg1[i,3]
                 if(is.na(temp.anno)) next
                 temp.anno <- strsplit(temp.anno, split = ";")[[1]]
                 temp.tags <- tags2.neg[[index.neg[i]]]
                 temp.tags.annotation <- temp.tags@annotation
                 temp.tags.to <- unlist(lapply(temp.tags.annotation, function(x) x$to))
                 temp.idx <- which(temp.tags.to %in% temp.anno)
                 temp.tags.annotation <- temp.tags.annotation[temp.idx]
                 tags2.neg[[index.neg[i]]]@annotation <- temp.tags.annotation
                 rm(list = c("temp.anno", "temp.tags", "temp.tags.annotation", "temp.tags.to",
                             "temp.idx", "temp.tags.annotation"))
                 gc()
               }



               ##construct kegg.result for tags2 positive
               index.pos <- which(!is.na(annotation.pos2$Annotation) & is.na(annotation.pos1$Annotation))
               if(length(index.pos) > 0){
                 kegg.result.pos <- vector(mode = "list", length = length(index.pos))
                 for(i in index.pos){
                   temp.name <- annotation.pos2$Peak.name[i]
                   temp.anno <- strsplit(annotation.pos2$Annotation[i], split = ";")[[1]]
                   temp.match.result <- match.result.pos2[[i]]
                   temp.id <- temp.match.result$ID
                   temp.idx <- which(temp.id %in% temp.anno)
                   temp.match.result <- temp.match.result[temp.idx,]
                   peakName <- temp.name
                   peakIndex <- match(temp.name, peak.name.pos)
                   peakMz <- showTags2(tags2.pos[peakIndex], slot = "mz")
                   peakRT <- showTags2(tags2.pos[peakIndex], slot = "rt")
                   # mzError.ppm <- abs(temp.match.result$Accurate.mass - peakMz)*10^6/peakMz
                   mzError.ppm <- abs(temp.match.result$Accurate.mass - peakMz)*10^6/ifelse(peakMz >= 400, peakMz, 400)
                   rtError <- abs(temp.match.result$RT - peakRT)*100/peakRT
                   peakID <- temp.match.result$ID
                   adduct <- temp.match.result$Adduct
                   charge <- temp.match.result$Charge
                   theoreticalMZ <- temp.match.result$Accurate.mass
                   theoreticalRT <- temp.match.result$RT
                   temp.kegg.result <- data.frame(peakName, peakIndex, peakMz, peakRT,
                                                  mzError.ppm,rtError,
                                                  peakID, adduct, charge,
                                                  theoreticalMZ, theoreticalRT,
                                                  stringsAsFactors = FALSE)
                   kegg.result.pos[[match(i, index.pos)]] <- temp.kegg.result
                 }


                 for(i in 1:length(kegg.result.pos)){
                   tags2.pos <- kegg2peakInfo(kegg.result = kegg.result.pos[[i]],
                                              mz.tol = mz.tol,
                                              rt.tol = rt.tol,
                                              weight.mz = 0.5, weight.rt = 0.5,weight.dp = 0,
                                              peak.info = tags2.pos,
                                              kegg.compound = kegg.compound)
                 }


                 ##construct kegg.result for tags2 negative
                 index.neg <- which(!is.na(annotation.neg2$Annotation) & is.na(annotation.neg1$Annotation))
                 if(length(index.neg) > 0){
                   kegg.result.neg <- vector(mode = "list", length = length(index.neg))
                   for(i in index.neg){
                     temp.name <- annotation.neg2$Peak.name[i]
                     temp.anno <- strsplit(annotation.neg2$Annotation[i], split = ";")[[1]]
                     temp.match.result <- match.result.neg2[[i]]
                     temp.id <- temp.match.result$ID
                     temp.idx <- which(temp.id %in% temp.anno)
                     temp.match.result <- temp.match.result[temp.idx,]
                     peakName <- temp.name
                     peakIndex <- match(temp.name, peak.name.neg)
                     peakMz <- showTags2(tags2.neg[peakIndex], slot = "mz")
                     peakRT <- showTags2(tags2.neg[peakIndex], slot = "rt")
                     # mzError.ppm <- abs(temp.match.result$Accurate.mass - peakMz)*10^6/peakMz
                     mzError.ppm <- abs(temp.match.result$Accurate.mass - peakMz)*10^6/ifelse(peakMz>=400,peakMz, 400)
                     rtError <- abs(temp.match.result$RT - peakRT)*100/peakRT
                     peakID <- temp.match.result$ID
                     adduct <- temp.match.result$Adduct
                     charge <- temp.match.result$Charge
                     theoreticalMZ <- temp.match.result$Accurate.mass
                     theoreticalRT <- temp.match.result$RT
                     temp.kegg.result <- data.frame(peakName, peakIndex, peakMz, peakRT,
                                                    mzError.ppm,rtError,
                                                    peakID, adduct, charge,
                                                    theoreticalMZ, theoreticalRT,
                                                    stringsAsFactors = FALSE)
                     kegg.result.neg[[match(i, index.neg)]] <- temp.kegg.result
                   }


                   for(i in 1:length(kegg.result.neg)){
                     tags2.neg <- kegg2peakInfo(kegg.result = kegg.result.neg[[i]],
                                                mz.tol = mz.tol,
                                                rt.tol = rt.tol,
                                                weight.mz = 0.5, weight.rt = 0.5,weight.dp = 0,
                                                peak.info = tags2.neg,
                                                kegg.compound = kegg.compound)
                   }
                 }

               }else{
                 tags2.pos <- NULL
                 tags2.neg <- NULL
               }
             }


             return.result.pos <- list(tags2.pos, anno.pos, anno.pos1, anno.pos2)
             names(return.result.pos) <- c("tags2.pos", "anno.pos", "anno.pos1", "anno.pos2")

             return.result.neg <- list(tags2.neg, anno.neg, anno.neg1, anno.neg2)
             names(return.result.neg) <- c("tags2.neg", "anno.neg", "anno.neg1", "anno.neg2")
             return.result <- list(return.result.pos, return.result.neg)
             names(return.result) <- c("positive", "negative")
             return.result <- return.result
           })



