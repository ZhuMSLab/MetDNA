#' @title metPathway
#' @description Metabolic pathway analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param annotation.result Annotation.result from ms2Annotation.
#' @param group The group you want to use.
#' @param max.isotope The max number of isotope.
#' @param uni.test The method of univariate test.
#' @param column hilic or rp.
#' @param polarity positive or negative.
#' @param path The work directory.
#' @param output.path The directory for outputting results.
#' @param correct Correct p value or not.
#' @param p.cutoff The cutoff of p value. Default is 0.01.
#' @param mz.tol The mz tolerance of KEGG matching. Defalut is 25 ppm.
#' @param rt.tol1 The RT tolerance of isotope annotation (second). Default is
#' 3 s.
#' @param rt.tol2 The RT tolerance of KEGG matching (\%). Default is 30\%.
#' @param cor.tol THe tolerance of correlation.
#' @param int.tol The tolerance of intensity ratio for
#' isotope annotation(\%).
#' @param threads The number of threads.
#' @param use.old.null Use old null distribution of not. Default is FALSE.
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




# annotation.result = "ms2.match.annotation.result.csv"
# group = c("Normal", "AD")
# max.isotope = 4
# uni.test = "t"
# column = "rp"
# polarity = "positive"
# path =
#   paste(file.path(".", "MRN_annotation_result"),
#         ifelse(polarity == "positive", "POS", 'NEG'), sep = "_")
# output.path = paste(file.path(".", "Dysregulated_network_analysis_result"),
#                     ifelse(polarity == "positive", "POS", "NEG"), sep = '_')
# correct = FALSE
# p.cutoff = 0.05
# mz.tol = 25
# rt.tol1 = 3# absolute, second
# rt.tol2 = 30#relative, %
# cor.tol = 0
# int.tol = 500
# threads = 3
# # use.old.null = FALSE,
# species = "hsa"
# use.all.kegg.id = FALSE
# only.mrn.annotation = FALSE
# use.old.null = TRUE
#
# system.time(metPathway(annotation.result = annotation.result, group = group, max.isotope = max.isotope, column = "rp",
#                        polarity = "positive", correct = FALSE, p.cutoff = 0.05, species = "hsa", use.all.kegg.id = FALSE,
#                        only.mrn.annotation = FALSE))



setGeneric(name = "metPathway",
           def = function(annotation.result = "ms2.match.annotation.result.csv",
                          group,
                          max.isotope = 4,
                          uni.test = c("t", "wilcox", "anova"),
                          column = c("hilic", "rp"),
                          polarity = c("positive", "negative"),
                          path =
                            paste(file.path(".", "MRN_annotation_result"),
                                  ifelse(polarity == "positive", "POS", 'NEG'), sep = "_"),
                          output.path = paste(file.path(".", "Dysregulated_network_analysis_result"),
                                              ifelse(polarity == "positive", "POS", "NEG"), sep = '_'),
                          correct = TRUE,
                          p.cutoff = 0.01,
                          mz.tol = 25,
                          rt.tol1 = 3,# absolute, second
                          rt.tol2 = 30,#relative, %
                          cor.tol = 0,
                          int.tol = 500,
                          threads = 3,
                          # use.old.null = FALSE,
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          use.all.kegg.id = FALSE,
                          only.mrn.annotation = FALSE,
                          use.old.null = FALSE) {

             path <- path[1]
             output.path <- output.path[1]
             dir.create(output.path)
             dir.create(file.path(output.path, "intermediate_data"))
             # dir.create(file.path(output.path, "module_information"))
             dir.create(file.path(output.path, "pathway_information"))
             if(missing(group)) {
               stop("You must give the names of group which you want to process.")}

             uni.test <- match.arg(uni.test)
             column <- match.arg(column)
             polarity <- match.arg(polarity)
             species <- match.arg(species)

             metPathway.parameters <- c(paste(group, collapse = ";"), max.isotope,
                                       uni.test, column,
                                       polarity, correct,
                                       p.cutoff, mz.tol, rt.tol1, rt.tol2,
                                       threads, species)
             metPathway.parameters <- data.frame(c("group", "max.isotope",
                                                  "uni.test", "column",
                                                  'polarity', "correct",
                                                  "p.cutoff", "mz.tol", "rt.tol1", "rt.tol2",
                                                  "threads",
                                                  "species"),
                                                metPathway.parameters,
                                                stringsAsFactors = FALSE)
             colnames(metPathway.parameters) <- c('parameter', "value")
             write.csv(metPathway.parameters, file.path(output.path, "metPathway.parameters.csv"),
                       row.names = FALSE)
             rm(list = "metPathway.parameters")
             gc()
             # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))

             #check for need files
             file1 <- dir(file.path(path, "intermediate_data"))
               need.file <- c(annotation.result, "tags2.after.redundancy.remove",
                              "sample.info.csv", "rt.result")

               file.check <- which(!need.file %in% file1)
               if (length(file.check) > 0) {
                 stop(paste(need.file[file.check], collapse = " & "),
                      " are not in the directory.")
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


             rm(list = c("file1"))
             gc()

             load(file.path(path, "intermediate_data","tags2.after.redundancy.remove"))
             load(file.path(path, "intermediate_data","rt.result"))
             kegg.rt <- rt.result[[2]]

             rm(rt.result)
             gc()

             if (polarity == "positive") {
               data("kegg.adduct.pos", envir = environment())
               kegg.adduct <- kegg.adduct.pos
               rm(list = "kegg.adduct.pos")
               gc()
             } else{
               data("kegg.adduct.neg", envir = environment())
               kegg.adduct <- kegg.adduct.neg
               rm(list = "kegg.adduct.neg")
               gc()
             }

             if (column == "hilic") {
               kegg.adduct <- kegg.adduct[which(kegg.adduct$Column.hilic), ]
             } else{
               kegg.adduct <- kegg.adduct[which(kegg.adduct$Column.rp), ]
             }


             ##Add RT information to kegg.adduct
             idx <- match(kegg.adduct$ID, row.names(kegg.rt))
             kegg.adduct <- data.frame(kegg.adduct, kegg.rt[idx, ],
                                       stringsAsFactors = FALSE)

             kegg.adduct$RT[is.na(kegg.adduct$RT)] <-
               median(kegg.adduct$RT[!is.na(kegg.adduct$RT)])

             rm(list = c("idx","kegg.rt"))
             gc()

             ##sample.info
             sample.info <-
               readr::read_csv(
                 file.path(path, "intermediate_data","sample.info.csv"),
                 progress = FALSE,
                 col_types = readr::cols()
               )

             sample.info <- as.data.frame(sample.info)

             file.copy(file.path(path, "intermediate_data", "sample.info.csv"),
                       to = file.path(output.path, "intermediate_data"))

             if(any(!is.element(group, unique(sample.info[,2])))){
               idx <- which(!is.element(group, unique(sample.info[,2])))
               stop(paste(group[idx], collapse = " and "),
                    ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
             }

             sample.info <- sample.info[sample.info[,2] %in% group,]

             annotation.result <-
               readr::read_csv(file.path(path, "intermediate_data",annotation.result),
                               progress = FALSE,
                               col_types = readr::cols())

             file.copy(file.path(path, "intermediate_data", "ms2.match.annotation.result.csv"),
                       to = file.path(output.path, "intermediate_data"))

             sample <-
               annotation.result[,match(sample.info[,1], colnames(annotation.result))]

             sample <- as.data.frame(sample)
             rownames(sample) <- annotation.result$name


             #--------------------------------------------------------------------------
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

             rm(annotation.result)
             gc()

             num.peak <- sum(p.value <= p.cutoff)
             cat("\n")
             cat("There are ", num.peak, " out of ", nrow(sample), " peaks with p values less than ",
                 p.cutoff, ".\n", sep = "")
             cat("\n")

             while(num.peak <= 30){
               p.cutoff <- readline("Do you want to change the p.cutoff or stop this analysis?
                                    (type a new p.cutoff or type b to stop analysis)")
               if(p.cutoff == "b") stop()
               p.cutoff <- as.numeric(p.cutoff)
               num.peak <- sum(p.value <= p.cutoff)
               cat("There are ", num.peak, " out of ", nrow(sample), " peaks with p values less than ",
                   p.cutoff, ".\n", sep = "")
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

             adduct.table <-
               adduct.table[adduct.table$mode == polarity, ]

             #-----------------------------------------------------------------------
             #begin find the pathway using pathway enrichment method.
             return.result <-
               findPathway(tags2 = tags2.after.redundancy.remove,
                 metabolic.network = kegg.rpair2,
                 metabolite = kegg.adduct,
                 kegg.compound = kegg.compound,
                 p.value = p.value,
                 p.cutoff = p.cutoff,
                 adduct.table = adduct.table,
                 polarity = polarity,
                 mz.tol = mz.tol,
                 rt.tol = rt.tol2,
                 threads = threads,
                 path = path,
                 output.path = output.path,
                 use.old.null = use.old.null,
                 species = species,
                 use.all.kegg.id = use.all.kegg.id,
                 only.mrn.annotation = only.mrn.annotation)

             new.tags <- return.result[[1]]
             anno <- return.result[[2]]
             anno.from.pathway <- anno
             save(anno.from.pathway,
                  file = file.path(output.path,
                                   "intermediate_data","anno.from.pathway"),
                  compress = "xz")

             rm(list = c("return.result", "anno", "anno.from.pathway",
                         "kegg.adduct", "kegg.compound", "adduct.table"))
             gc()


             #process new.tags, for the annotation form KEGG.
             if(is.null(new.tags)){
               tags2.after.kegg.matching.from.pathway <-
                 tags2.after.redundancy.remove
               rm(list = c("tags2.after.redundancy.remove"))
               gc()
               load(file.path(path, "intermediate_data","tags.result"))
               tags.result2.from.pathway <- tags.result
               rm(list = c("tags.result"))
               gc()
               save(tags2.after.kegg.matching.from.pathway,
                    file = file.path(output.path,
                                     "intermediate_data",
                                     "tags2.after.kegg.matching.from.pathway"),
                    compress = "xz")
               save(tags.result2.from.pathway,
                    file = file.path(output.path,
                                     "intermediate_data",
                                     "tags.result2.from.pathway"), compress = "xz")
             }else{
               cat('\n')
               cat("Process new.tags.\n")
               temp.result <- processTags2(
                 new.tags2 = new.tags,
                 sample = sample,
                 mz.tol = mz.tol,
                 rt.tol = rt.tol1,
                 cor.tol = cor.tol,
                 int.tol = int.tol,
                 max.isotope = max.isotope,
                 polarity = polarity,
                 path = output.path
               )

               tags2.after.kegg.matching.from.pathway <- temp.result[[1]]
               tags.result2.from.pathway <- temp.result[[2]]

               rm(list = c("temp.result"))

               save(tags2.after.kegg.matching.from.pathway,
                    file = file.path(output.path,
                                     "intermediate_data",
                                     "tags2.after.kegg.matching.from.pathway"),
                    compress = "xz")
               # if(!is.na(tags.result2.from.pathway[[1]])){
               #   # file.copy(file.path(output.path, "redun2"), file.path(output.path, "intermediate_data"))
               #   # unlink(x = file.path(output.path, "redun2"), recursive = TRUE)
               # }else{
               #   load(file.path(path, "intermediate_data","tags.result"))
               #   tags.result2.from.pathway <- tags.result
               #   rm(list = c("tags.result"))
               #   gc()
               # }
               save(tags.result2.from.pathway,
                    file = file.path(output.path,
                                     "intermediate_data",
                                     "tags.result2.from.pathway"), compress = "xz")
             }

             ##change tags2.after.kegg.matching to csv like ms2Annotation

             cat("\n")
             # load(file.path(output.path, "intermediate_data", "p.value"))
             # load(file.path(output.path, "intermediate_data", "foldchange"))
             data("kegg.compound", envir = environment())

             temp <- getAnnotationResult(tags2 = tags2.after.kegg.matching.from.pathway,
                                         p.value = p.value,
                                         correct = correct,
                                         foldchange = foldchange,
                                         tags.result2 = tags.result2.from.pathway,
                                         kegg.compound = kegg.compound)

             readr::write_csv(temp,
                              file.path(output.path,
                                        "DNA.pathway.annotation.result.csv"))

             rm(list = c("temp", "kegg.compound"))
             gc()

             rm(list = c("tags2.after.kegg.matching.from.pathway",
                         "temp.result", "new.tags",
                         "tags.result2.from.pathway"))
             gc()
             #
             singleTransform(annotation.result = "ms2.match.annotation.result.csv",
                             data.type = "metabolomics",
                             sample.info = "sample.info.csv",
                             group = group, trans.to = "pathway",
                             scale = TRUE,
                             scale.method = "pareto",
                             species = species,
                             which.peak = "intensity",
                             polarity = polarity,
                             column = column, method = "mean",
                             path = file.path(output.path, "intermediate_data"),
                             output.path = file.path(output.path, "quantitative_information"),
                             metabolic.network = kegg.rpair2)

             load(file.path(output.path,
                            "quantitative_information",
                            "node.quantitative.data"))
             node.quantitative.data2 <- node.quantitative.data
             save(node.quantitative.data2, file = file.path(output.path,
                                                            "quantitative_information",
                                                            "node.quantitative.data2"), compress = "xz")
             readr::write_csv(node.quantitative.data2, file.path(output.path,
                                                                 "quantitative_information",
                                                                 "node.quantitative.data2.csv"))

             unlink(file.path(output.path, "quantitative_information",
                              "node.quantitative.data"), recursive = TRUE)
             unlink(file.path(output.path, "quantitative_information",
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
             cat("metPathway is done\n")
           })


######-----------------------------------------------------------------------
# tags2 = tags2.after.redundancy.remove
# metabolic.network = kegg.rpair2
# metabolite = kegg.adduct
# kegg.compound = kegg.compound
# p.value = p.value
# p.cutoff = p.cutoff
# adduct.table = adduct.table
# polarity = polarity
# mz.tol = mz.tol
# rt.tol = rt.tol2
# threads = threads
# path = path
# output.path = output.path
# use.old.null = use.old.null
# species = species
# use.all.kegg.id = use.all.kegg.id
# only.mrn.annotation = only.mrn.annotation

# system.time(test <- findPathway(tags2 = tags2.after.redundancy.remove, metabolic.network = kegg.rpair2, metabolite = kegg.adduct.neg,
#                     kegg.compound = kegg.compound, p.value = p.value,p.cutoff = p.cutoff, adduct.table = adduct.table,
#                     polarity = polarity, mz.tol = mz.tol, rt.tol = rt.tol, threads = threads, path = path,
#                     output.path = output.path, use.old.null = use.old.null, species = species, use.all.kegg.id = use.all.kegg.id,
#                     only.mrn.annotation = FALSE))


setGeneric(name = "findPathway", def = function(tags2,
                                                metabolic.network,
                                                metabolite,
                                                kegg.compound,
                                                p.value,
                                                p.cutoff = 0.01,
                                                adduct.table,
                                                polarity = c("positive", "negative"),
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
  if(sum(p.value < p.cutoff) == 0) stop("No dysregulated peaks.\n")

  polarity <- match.arg(polarity)
  species <- match.arg(species)
  adduct.table <- adduct.table[adduct.table$mode == polarity,]
  metabolite <- metabolite[metabolite$Mode == polarity,]

  result <- tags2Result(tags2 = tags2, score.cutoff = 0)
  # rm(tags2)
  ##find the significant changed peaks
  if(only.mrn.annotation){
    index <-which(p.value < p.cutoff & !is.na(result$to))
    if(length(index) ==0){
      warning("The peaks with p value less than 0.05 have no MRN based annotations.\nSo only.mrn.annotation has been set as FALSE.")
      index <- which(p.value < p.cutoff)
      only.mrn.annotation <- FALSE
    }
  }else{
    index <- which(p.value < p.cutoff)
  }

  ##size is the number of dysregulated peaks
  size <- length(index)

  result1 <- result[index,]
  ##remove isotope annotation form result, this function is in tools
  result1 <- removeIsotopeFromResult(result = result1)

  to1 <- result1$to
  names(to1) <- result$name[index]
  rm(index)
  gc()
  to1 <- lapply(to1, function(x) {
    if(is.na(x)) {return(x)
    }else{
      x <- strsplit(x, split = ";")[[1]]
      x <- x[!duplicated(x)]
      x <- paste(x, collapse = ";")
      x
    }
  })
  names(to1) <- result1$name
  anno1 <- to1
  ###kegg mz and rt matching
  mz <- as.numeric(result1$mz)
  rt <- as.numeric(result1$rt)
  cat("\n")
  cat("Annotate peaks utilizing KEGG database.\n")
  match.result <- keggMatch(mz = mz, rt = rt,
                            metabolite = metabolite,
                            mz.tol = mz.tol,rt.tol = rt.tol,
                            polarity = polarity)

  ##anno2 is the annotation results from KEGG according to mz and rt
  anno2 <- lapply(match.result, function(x) {
    if(is.null(x)) return(NA)
    temp <- unique(x$ID)
    temp <- paste(temp, collapse = ";")
    temp
  })

  names(anno2) <- result1$name

  if(use.all.kegg.id){
    anno <- mapply(function(x,y){
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
    x = anno1,
    y = anno2)
  }else{
    anno <- mapply(function(x,y){
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
    x = anno1,
    y = anno2)
  }

  names(anno) <- result1$name

  index <- unname(which(unlist(lapply(anno, function(x) {!all(is.na(x))}))))
  result2 <- result1[index,]
  rm(result1)
  gc()
  anno <- anno[index]
  anno1 <- anno1[index]
  anno2 <- anno2[index]

  match.result2 <- match.result[index]
  rm(match.result)
  gc()


  #unique.id is the annotations of dysregulated peaks
  unique.id <- unique(unname(unlist(anno)))
  # unique.id <- unique.id[which(unique.id %in% igraph::V(metabolic.network)$name)]
  cat("\n")
  cat("Identify initial pathways.\n")

  mse.result <- mseAnalysis(metabolite.id = unique.id, species = species)

  ##get NULL distribution
  cat("\n")
  cat("Get pseudo pathways.\n")
  if(!use.old.null){
    ##temp.fun is function for multiple threads calculation
    temp.fun <- function(idx, result, size, only.mrn.annotation,
                         metabolite, mz.tol, rt.tol, polarity,
                         use.all.kegg.id,
                         keggMatch,
                         mseAnalysis1,
                         M = M,
                         removeIsotopeFromResult){
      temp.idx <- sample(1:nrow(result), size)
      if(only.mrn.annotation){
        temp.idx <- temp.idx[which(!is.na(result$to[temp.idx]))]
      }

      temp.result1 <- result[temp.idx,]
      ##remove isotope annotation form result
      temp.result1 <- removeIsotopeFromResult(result = temp.result1)

      temp.to1 <- temp.result1$to
      names(temp.to1) <- result$name[temp.idx]

      temp.to1 <- lapply(temp.to1, function(x) {
        if(is.na(x)) {return(x)
        }else{
          x <- strsplit(x, split = ";")[[1]]
          x <- x[!duplicated(x)]
          x <- paste(x, collapse = ";")
          x
        }
      })
      names(temp.to1) <- temp.result1$name
      temp.anno1 <- temp.to1

      ###kegg mz and rt matching
      temp.mz <- as.numeric(temp.result1$mz)
      temp.rt <- as.numeric(temp.result1$rt)


      temp.match.result <- keggMatch(mz = temp.mz, rt = temp.rt,
                                     metabolite = metabolite,
                                     mz.tol = mz.tol,rt.tol = rt.tol,
                                     polarity = polarity)

      ##anno2 is the annotation results from KEGG according to mz and rt
      temp.anno2 <- lapply(temp.match.result, function(x) {
        if(is.null(x)) return(NA)
        temp <- unique(x$ID)
        temp <- paste(temp, collapse = ";")
        temp
      })

      names(temp.anno2) <- temp.result1$name

      if(use.all.kegg.id){
        temp.anno <- mapply(function(x,y){
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
        x = temp.anno1,
        y = temp.anno2)
      }else{
        temp.anno <- mapply(function(x,y){
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
        x = temp.anno1,
        y = temp.anno2)
      }

      names(temp.anno) <- temp.result1$name

      temp.index <- unname(which(unlist(lapply(temp.anno, function(x) {!all(is.na(x))}))))
      temp.result2 <- temp.result1[temp.index,]
      rm(temp.result1)
      gc()
      temp.anno <- temp.anno[temp.index]
      temp.anno1 <- temp.anno1[temp.index]
      temp.anno2 <- temp.anno2[temp.index]

      temp.match.result2 <- temp.match.result[temp.index]
      rm(temp.match.result)
      gc()

      #unique.id is the annotations of significant peaks
      temp.unique.id <- unique(unname(unlist(temp.anno)))
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
                                              result = result,
                                              size = size,
                                              only.mrn.annotation = only.mrn.annotation,
                                              metabolite = metabolite,
                                              mz.tol = mz.tol,
                                              rt.tol = rt.tol,
                                              polarity = polarity,
                                              use.all.kegg.id = use.all.kegg.id,
                                              keggMatch = keggMatch,
                                              mseAnalysis1 = mseAnalysis1,
                                              M = M,
                                              removeIsotopeFromResult = removeIsotopeFromResult)
    save(temp.mse.result, file = file.path(output.path, "intermediate_data","temp.mse.result"), compress = "xz")
  }else{
    cat("\n")
    cat("Use old pseudo scores.\n")
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
  pathway.result <- mse.result
  save(pathway.result, file = file.path(output.path, "intermediate_data","pathway.result"), compress = "xz")

  # rm(list = c("pathway.result"))
  gc()

  ##filter annotation accoring to enrichment pathways
  index <- which(mse.result$adjusted.p.value.by.overlap < 0.05 | mse.result$adjusted.p.value.by.pvalue < 0.05)

  if(length(index) == 0){
    cat("There are no initial pathways with p values less than 0.05.\n")
  }else{
    cat("Filter annotations according to dysregulated pathways.\n")

    annotation1 <- vector(mode = "list", length = length(index))
    annotation2 <- vector(mode = "list", length = length(index))

    for(i in index){
      # cat(i); cat(" ")
      temp.pathway <- mse.result$detected.metabolite.id[i]
      temp.pathway <- strsplit(temp.pathway, split = ";")[[1]]
      anno1.result <- lapply(anno1, function(x){
        if(is.na(x)){return(NA)}
        x <- strsplit(x, split = ";")[[1]]
        idx <- which(x %in% temp.pathway)
        if(length(idx) == 0) return(NA)
        return(paste(x[idx], collapse = ";"))
      })
      anno1.result <- unlist(anno1.result)

      anno2.result <- lapply(anno2, function(x){
        if(is.na(x)){return(NA)}
        x <- strsplit(x, split = ";")[[1]]
        idx <- which(x %in% temp.pathway)
        if(length(idx) == 0) return(NA)
        return(paste(x[idx], collapse = ";"))
      })

      anno2.result <- unlist(anno2.result)

      annotation1[[i]] <- anno1.result
      annotation2[[i]] <- anno2.result
    }

    annotation1 <- do.call(cbind, annotation1)
    annotation2 <- do.call(cbind, annotation2)

    peak.name <- names(anno1)

    annotation1 <- data.frame(peak.name,
                              unlist(anno1), annotation1, stringsAsFactors = FALSE)
    annotation2 <- data.frame(peak.name,
                              unlist(anno2), annotation2, stringsAsFactors = FALSE)

    perfect.annotation1 <- apply(annotation1, 1, function(x){
      x <- as.character(x)
      if(is.na(x[2])) return(NA)
      temp.idx <- which(!is.na(x[-c(1:2)]))
      if(length(temp.idx) == 0) return(x[2])
      if(length(temp.idx) == 1) return(x[-c(1:2)][temp.idx])
      return(paste(x[-c(1:2)][temp.idx], collapse = ";"))
    })


    perfect.annotation2 <- apply(annotation2, 1, function(x){
      x <- as.character(x)
      if(is.na(x[2])) return(NA)
      temp.idx <- which(!is.na(x[-c(1:2)]))
      if(length(temp.idx) == 0) return(NA)
      return(paste(x[-c(1:2)][temp.idx], collapse = ";"))
    })

    ##new anno should be from prefer.annotation1 and prefer.annotation2
    anno <- mapply(function(x,y){
      if(!is.na(x)) return(list(x))
      return(list(y))
    },
    x = perfect.annotation1,
    y = perfect.annotation2)

    anno <- lapply(anno, function(x){
      if(is.na(x)) return(NA)
      strsplit(x, split = ";")[[1]]
    })

    names(anno) <- names(anno1)
    anno <- anno[-which(unlist(lapply(anno, function(x) all(is.na(x)))))]




    colnames(annotation1) <- c("Peak.name","MRN.annotation",
                               mse.result$pathway.name[index])
    colnames(annotation2) <- c("Peak.name","KEGG.annotation",
                               mse.result$pathway.name[index])

    annotation1 <- data.frame(annotation1[,1:2], perfect.annotation1,
                              annotation1[,-c(1:2)], stringsAsFactors = FALSE)
    annotation2 <- data.frame(annotation2[,1:2], perfect.annotation2,
                              annotation2[,-c(1:2)], stringsAsFactors = FALSE)

    colnames(annotation1)[3] <- colnames(annotation2)[3] <- "Annotation"

    write.csv(annotation1,
              file.path(output.path, "intermediate_data",
                        'Significant.peaks.annotation.MRN.from.pathway.csv'),
              row.names = FALSE)
    write.csv(annotation2,
              file.path(output.path, "intermediate_data",
                        'Significant.peaks.annotation.KEGG.from.pathway.csv'),
              row.names = FALSE)


    ##add those result to tags2
    peak.name <- showTags2(tags2, slot = "name")
    index <- match(annotation1$Peak.name, peak.name)

    for(i in 1:length(index)){
      temp.anno <- annotation1[i,3]
      if(is.na(temp.anno)) next
      temp.anno <- strsplit(temp.anno, split = ";")[[1]]
      temp.tags <- tags2[[index[i]]]
      temp.tags.annotation <- temp.tags@annotation
      temp.tags.to <- unlist(lapply(temp.tags.annotation, function(x) x$to))
      temp.idx <- which(temp.tags.to %in% temp.anno)
      temp.tags.annotation <- temp.tags.annotation[temp.idx]
      tags2[[index[i]]]@annotation <- temp.tags.annotation
      rm(list = c("temp.anno", "temp.tags", "temp.tags.annotation", "temp.tags.to",
                  "temp.idx", "temp.tags.annotation"))
      gc()
    }

    ##construct kegg.result for tags2
    index <- which(!is.na(annotation2$Annotation) & is.na(annotation1$Annotation))
    if(length(index) > 0){
      kegg.result <- vector(mode = "list", length = length(index))
      for(i in index){
        temp.name <- annotation2$Peak.name[i]
        temp.anno <- strsplit(annotation2$Annotation[i], split = ";")[[1]]
        temp.match.result <- match.result2[[i]]
        temp.id <- temp.match.result$ID
        temp.idx <- which(temp.id %in% temp.anno)
        temp.match.result <- temp.match.result[temp.idx,]
        peakName <- temp.name
        peakIndex <- match(temp.name, peak.name)
        peakMz <- showTags2(tags2[peakIndex], slot = "mz")
        peakRT <- showTags2(tags2[peakIndex], slot = "rt")
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
        kegg.result[[match(i, index)]] <- temp.kegg.result
      }


      for(i in 1:length(kegg.result)){
        tags2 <- kegg2peakInfo(kegg.result = kegg.result[[i]], mz.tol = mz.tol,
                               rt.tol = rt.tol,
                               weight.mz = 0.5, weight.rt = 0.5,weight.dp = 0,
                               peak.info = tags2,
                               kegg.compound = kegg.compound)
      }
    }else{
      tags2 <- NULL
    }
  }

  return.result <- list(tags2, anno, anno1, anno2)
  names(return.result) <- c("tags2", "anno", "anno1", "annot2")
  return.result <- return.result

})





