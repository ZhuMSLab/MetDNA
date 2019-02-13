
#
# system.time(metdnaForOneMode(polarity = "positive",
#                 column = "rp",
#                 group = c("Normal", "AD"),
#                 uni.test = "t",
#                 path = "/home/jasper/Downloads/7f68899f/POS",
#                 correct = FALSE,
#                 p.cutoff = 0.05,
#                 species = "hsa",
#                 dn.analysis = TRUE,
#                 use.old.null = TRUE,
#                 pathway.enrichment = TRUE))


setGeneric(name = "metdnaForOneMode",
           def = function(ms1.data = "data.csv",
                          sample.info = "sample.info.csv",
                          mz.tol = 25,
                          rt.tol.for.ms1.ms2.match = 10,#second
                          path = ".",
                          log.path = ".",
                          instrument = c("SciexTripleTOF", "AgilentQTOF",
                                         "OtherQTOF", "ThermoOrbitrap"),
                          ms2.type = c("mgf", "msp", "mzXML"),
                          polarity = c("positive", "negative", "both"),
                          column = c("hilic", "rp"),
                          ce = c("30", "10", "20", "35,15", "40", "50"),
                          ms2.match.plot = TRUE,
                          candidate.num = 5,
                          ###
                          prefer.adduct = ifelse(polarity == "positive", "M+H", "M-H"),
                          use.default.md = TRUE,
                          threads = 3,
                          max.isotope = 4,
                          rt.tol1 = 3,#second
                          rt.tol2 = 30,#%
                          cor.tol = 0,
                          int.tol = 500,
                          dp.tol = 0.5,
                          max.step = 3,
                          score.cutoff = 0,
                          remain = FALSE,
                          remain.per = 0.5,
                          seed.neighbor.match.plot = FALSE,
                          ###
                          group,
                          uni.test = c("t", "wilcox"),
                          correct = TRUE,
                          p.cutoff = 0.01,
                          use.old.null = FALSE,
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          use.all.kegg.id = FALSE,
                          only.mrn.annotation = FALSE,
                          check.data = TRUE,
                          ms2.annotation = TRUE,
                          mrn.annotation = TRUE,
                          dn.analysis = FALSE,
                          pathway.enrichment = FALSE,
                          analysis.report = TRUE
           ){



             #parameter decision
             instrument <- match.arg(instrument)
             ms2.type <- match.arg(ms2.type)
             ##change the ms2.type to mgf if there are mgf files, and change
             ## ms2.type to msp if there are msp files.
             temp.file <- dir(path)
             postfix <- getPostfix(x = temp.file)
             postfix <- unique(postfix[!is.na(postfix)])

             if(any(postfix == "mgf")){
               ms2.type <- "mgf"
             }

             if(any(postfix == "msp")){
               ms2.type <- "msp"
             }

             polarity <- match.arg(polarity)
             column <- match.arg(column)
             ce <- match.arg(ce)
             uni.test <- match.arg(uni.test)
             species <- match.arg(species)


             #--------------------------------------------------------------------
             #check group
             old.sample.info.name <- sample.info
             old.data.name <- ms1.data
             temp.error <- errorDisplay(
               sample.info <- read.csv(file.path(path, sample.info), stringsAsFactors = FALSE),
               error.info = paste("Error: There is no", sample.info, "in", path, ". Please check it.\n"))
             if(class(temp.error) == "try-error") stop("Error")

             # sample.info <- read.csv(file.path(path, sample.info), stringsAsFactors = FALSE)
             if(missing(group)){
               cat("You don't provide the group names.\n")
               cat("You don't provide the group names.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               group <- sort(unique(as.character(sample.info[,2])))
               if(length(group) == 2){
                 cat("The group:", group, ", will be used.\n")
                 cat("The group:", group, ", will be used.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               }else{
                 cat(paste("Error: There are", length(group), "groups (", group, ") in your sample.info. Please provide the group names you want to process.\n"), file = file.path(log.path, "run.log.txt"), append = TRUE)
                 stop(paste("There are", length(group), "groups (", group, ") in your sample.info. Please provide the group names you want to process.\n"))
               }
             }else{
               temp.idx <- which(!(group %in% as.character(sample.info[,2])))
               if(length(temp.idx) > 0){
                 cat(paste('Error: ',group[temp.idx], "is/are not in your sample.info. (", unique(sample.info[,2]), ")\n"), file = file.path(log.path, "run.log.txt"), append = TRUE)
                 stop(paste(group[temp.idx], "is/are not in your sample.info. (", unique(sample.info[,2]), ")\n"))
               }
             }

             ##check ms2.annotation.result and MRN.annotation.result exist or not.
               if(any(dir(file.path(path, ifelse(polarity == "positive",
                                                 "MS2_match_result",
                                                 "MS2_match_result"))) == "ms2.match.annotation.result.csv")){
                 check.data <- FALSE
                 ms2.annotation <- FALSE
               }

               if(any(dir(file.path(path, ifelse(polarity == "positive",
                                                 "MRN_annotation_result",
                                                 "MRN_annotation_result"))) == "MRN.annotation.result.csv")){
                 check.data <- FALSE
                 ms2.annotation <- FALSE
                 mrn.annotation <- FALSE
               }

             ###check data
             cat("1. Check data.\n")
             cat("1. Check data.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
             cat("------------------------------------------------------------\n")
             if(check.data){
               switch(ms2.type,
                      "mgf" = {ms2.file <- grep("mgf", dir(path), value = TRUE)},
                      "msp" = {ms2.file <- grep("msp", dir(path), value = TRUE)},
                      "mzXML" = {ms2.file <- grep("mzXML", dir(path), value = TRUE)})


               check.result <- try({checkData(data = ms1.data,
                                         sample.info = old.sample.info.name,
                                         ms2.type = ms2.type,
                                         path = path)}, silent = TRUE)

               if(class(check.result) == "try-error"){
                 cat(check.result[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
                 stop(check.result[[1]])
               }


               if(any(as.numeric(check.result[,4]) > 0)){
                 cat('Error: Please check your data to make sure that they are valid.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
                 stop('Please check your data to make sure that they are valid.\n')
               }
             }else{
               cat('Skip data check.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Skip data check.\n")
             }


             ###ms2Annotation
             cat("\n")
             cat('2. MS/MS match annotation.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
             cat("2. MS/MS match annotation.\n")
             cat("------------------------------------------------------------\n")
             if(ms2.annotation){

               temp.error <- try({ms2Annotation(ms1.file = ms1.data,
                                                     sample.info = old.sample.info.name,
                                                     mz.tol = mz.tol,
                                                     rt.tol = rt.tol.for.ms1.ms2.match,
                                                     path = path,
                                                     output.path = file.path(path, "MS2_match_result"),
                                                     instrument = instrument,
                                                     ms2.type = ms2.type,
                                                     polarity = polarity,
                                                     column = column,
                                                     ce = ce,
                                                     ms2.match.plot = ms2.match.plot)}, silent = TRUE)

             if(class(temp.error) == "try-error") {
               cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
               stop(temp.error[[1]])
             }

             }else{
               cat('Skip MS/MS match annotation.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Skip MS/MS match annotation.\n")
             }

             ###metABM
             if(ms2.annotation) cat("\n")
             cat("\n")
             cat('3. Metabolic reaction network based metabolite annotation.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
             cat("3. Metabolic reaction network based metabolite annotation.\n")
             cat("------------------------------------------------------------\n")
             if(mrn.annotation){


               temp.error <- try({metABM(annotation.result = "ms2.match.annotation.result.csv",
                      ms2.data = "ms2",
                      prefer.adduct = prefer.adduct,
                      use.default.md = use.default.md,
                      column = column,
                      polarity = polarity,
                      threads = threads,
                      path = file.path(path, "MS2_match_result"),
                      output.path = file.path(path, "MRN_annotation_result"),
                      max.isotope = max.isotope,
                      mz.tol = mz.tol,
                      rt.tol1 = rt.tol1,
                      rt.tol2 = rt.tol2,
                      cor.tol = cor.tol,
                      int.tol = int.tol,
                      dp.tol = dp.tol,
                      max.step = max.step,
                      score.cutoff = score.cutoff,
                      remain = remain,
                      remain.per = remain.per,
                      seed.neighbor.match.plot = seed.neighbor.match.plot,
                      candidate.num = candidate.num)}, silent = TRUE)

               if(class(temp.error) == "try-error") {
                 cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
                 stop(temp.error[[1]])
               }
             }else{
               cat('Skip metabolic reaction network based annotation.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Skip metabolic reaction network based annotation.\n")
             }


             ###metModule
             cat("\n")
             cat('4. Dysregulated network analysis.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
             cat("4. Dysregulated network analysis.\n")

               output.path <- file.path(path, "Dysregulated_network_analysis_result")
               dir.create(output.path)
               dir.create(file.path(output.path, "intermediate_data"))

             cat("------------------------------------------------------------\n")
             if(dn.analysis){
               cat('Module analysis.\n', file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Module analysis.\n")
               ##get the directory
               ##read data
               annotation.result <- readr::read_csv(file.path(path, ms1.data),
                                                        progress = FALSE,
                                                        col_types = readr::cols())
               annotation.result <- as.data.frame(annotation.result)
               sample.info <- readr::read_csv(file.path(path, old.sample.info.name),
                                                  progress = FALSE,
                                                  col_types = readr::cols())
               sample.info <- as.data.frame(sample.info)


               ###########################################################################
               #load tags2.data and rt.result
               #positive
               load(file.path(path, "MRN_annotation_result","intermediate_data","tags2.after.redundancy.remove"))
               load(file.path(path, "MRN_annotation_result","intermediate_data", "rt.result"))
               tags2 <- tags2.after.redundancy.remove
               rm(list="tags2.after.redundancy.remove")
               gc()


               ##############################################################################
               #Calculate p value and fold change
               if(any(!is.element(group, unique(sample.info[,2])))){
                 idx <- which(!is.element(group, unique(sample.info[,2])))
                 cat("Error: ",paste(group[idx], collapse = " and "),
                     ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 stop(paste(group[idx], collapse = " and "),
                      ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
               }

               sample.info <- sample.info[sample.info[,2] %in% group,]

               sample <- annotation.result[,match(sample.info[,1],
                                                          colnames(annotation.result))]

               sample <- as.data.frame(sample)
                 rownames(sample) <- annotation.result$name
                 cat("\n")
                 cat("Calculate p values of peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of peaks.\n")


                 if(ncol(sample) <= 3){
                   cat("At least one group has only on sample. So the p values are assigned 1.\n")
                   p.value <- rep(1, nrow(sample))
                 }else{
                   p.value <- uniTest(
                     sample = as.matrix(sample),
                     sample.info = sample.info,
                     uni.test = uni.test,
                     correct = correct)
                 }

                 names(p.value) <- rownames(sample)
                 p.value[is.na(p.value)] <- 1


               num.peak <- sum(p.value <= p.cutoff)
               cat("\n")
               cat("There are ", num.peak, " out of ",
                   length(p.value), " peaks with p values ",
                   ifelse(correct, '(FDR correct)', "(No correct)"),
                   " less than ",
                   p.cutoff, ".\n", sep = "", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("There are ", num.peak, " out of ",
                   length(p.value), " peaks with p values ",
                   ifelse(correct, '(FDR correct)', "(No correct)"),
                   " less than ",
                   p.cutoff, ".\n", sep = "")


               ##if peak number with p value less than p.cutoff is less than 10
               count <- 1
               # while(num.peak <= 10 & count <= 5){
               #   count <- count + 1
               #   cat("\n")
               #   cat("Error: Peaks with p values < p.cutoff is less than 10, please change the correct method and the p.cutoff.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               #   stop("Peaks with p values < p.cutoff is less than 10, please change the correct method and the p.cutoff.\n")
               #
               #   p.cutoff <- readline("p.cutoff: ")
               #   p.cutoff <- as.numeric(p.cutoff)
               #
               #   if(correct){
               #     cat("\n")
               #     correct <- readline("correct (Type 1 for TRUE or 0 for FALSE): ")
               #     while(correct != "1" & correct != "0"){
               #       cat("\n")
               #       correct <- readline("correct (Type 1 for TRUE or 0 for FALSE): ")
               #     }
               #     correct <- ifelse(correct == "1", TRUE, FALSE)
               #   }
               #
               #
               #   ##re-calculate p value
               #     cat("\n")
               #     cat("Re-calculate p values of peaks.\n")
               #     p.value <- uniTest(
               #       sample = as.matrix(sample),
               #       sample.info = sample.info,
               #       uni.test = uni.test,
               #       correct = correct)
               #     names(p.value) <- rownames(sample)
               #     p.value[is.na(p.value)] <- 1
               #
               #
               #   num.peak <- sum(p.value <= p.cutoff)
               #   cat("\n")
               #   cat("There are ", num.peak, " out of ",
               #       length(p.value), " peaks with p values ",
               #       ifelse(correct, '(FDR correct)', "(No correct)"),
               #       " less than ",
               #       p.cutoff, ".\n", sep = "")
               # }


                 cat("\n")
                 cat("Calculate fold change of peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of peaks.\n")
                 foldchange <- foldChange(sample = sample,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange) <- rownames(sample)

                 save(p.value, file = file.path(output.path, "intermediate_data", "p.value"))
                 save(foldchange, file = file.path(output.path, "intermediate_data", "foldchange"))


               ###volcano plot
               pdf(file = file.path(output.path, "volcano.plot.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               volcanoPlot(p.value = p.value, fc = foldchange, correct = correct,
                           p.cutoff = p.cutoff, fc.cutoff = 1, cex = 0.3)
               dev.off()

               if(polarity == "positive"){
                 sample.pos <- sample
                 sample.neg <- NULL
                 tags2.pos <- tags2
                 tags2.neg <- NULL
                 p.value.pos <- p.value
                 p.value.neg <- NULL
                 foldchange.pos <- foldchange
                 foldchange.neg <- NULL
                 pos.path <- path
                 neg.path <- NULL
               }else{
                 sample.neg <- sample
                 sample.pos <- NULL
                 tags2.neg <- tags2
                 tags2.pos <- NULL
                 p.value.neg <- p.value
                 p.value.pos <- NULL
                 foldchange.neg <- foldchange
                 foldchange.pos <- NULL
                 neg.path <- path
                 pos.path <- NULL
               }


               if(sum(p.value < p.cutoff) >= 10){

                 temp.error <- try({metModule(annotation.result.pos = sample.pos,
                           annotation.result.neg = sample.neg,
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
                           polarity = polarity,
                           group = group,
                           max.isotope = max.isotope,
                           uni.test = uni.test,
                           column = column,
                           output.path = output.path,
                           correct = correct,
                           p.cutoff = p.cutoff,
                           mz.tol = mz.tol,
                           rt.tol1 = rt.tol1,# absolute, second
                           rt.tol2 = rt.tol2,#relative, %
                           cor.tol = cor.tol,
                           int.tol = int.tol,
                           threads = threads,
                           use.old.null = use.old.null,
                           species = species,
                           use.all.kegg.id = use.all.kegg.id,
                           only.mrn.annotation = only.mrn.annotation,
                           output.module.information = TRUE)}, silent = TRUE)
                 if(class(temp.error) == "try-error") {
                   cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
                   stop(temp.error[[1]])
                 }
               }

             }else{
               cat("Skip dysregulated network analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Skip dysregulated network analysis.\n")

               ##read data
                 annotation.result <- readr::read_csv(file.path(path, ms1.data),
                                                          progress = FALSE,
                                                          col_types = readr::cols())
                 annotation.result <- as.data.frame(annotation.result)
                 sample.info <- readr::read_csv(file.path(path, old.sample.info.name),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
                 sample.info <- as.data.frame(sample.info)


               #Calculate p value and fold change
               if(any(!is.element(group, unique(sample.info[,2])))){
                 idx <- which(!is.element(group, unique(sample.info[,2])))
                 stop(paste(group[idx], collapse = " and "),
                      ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
               }

               sample.info <- sample.info[sample.info[,2] %in% group,]


                 sample <- annotation.result[,match(sample.info[,1],
                                                            colnames(annotation.result))]

                 sample <- as.data.frame(sample)
                   rownames(sample) <- annotation.result$name
                   cat("\n")
                   cat("Calculate p values of peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of peaks.\n")

                 if(ncol(sample) <= 3){
                   cat("At least one group has only on sample. So the p values are assigned 1.\n")
                   p.value <- rep(1, nrow(sample))
                 }else{
                   p.value <- uniTest(
                     sample = as.matrix(sample),
                     sample.info = sample.info,
                     uni.test = uni.test,
                     correct = correct)
                 }
                 names(p.value) <- rownames(sample)
                 p.value[is.na(p.value)] <- 1

                 cat("\n")
                 cat("Calculate fold change of peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of peaks.\n")
                 foldchange <- foldChange(sample = sample,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange) <- rownames(sample)

                 save(p.value, file = file.path(output.path, "intermediate_data", "p.value"))
                 save(foldchange, file = file.path(output.path, "intermediate_data", "foldchange"))
             }




             ###pathway enrichment analysis
             cat("\n")
             cat("5. Pathway enrichment analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
             cat("5. Pathway enrichment analysis.\n")

               output.path <- file.path(path,"Pathway_enrichment_analysis_result")
             dir.create(output.path)
             dir.create(file.path(output.path, "intermediate_data"))

             cat("------------------------------------------------------------\n")
             if(pathway.enrichment){
               cat("Pathway enrichment analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Pathway enrichment analysis.\n")
               cat("\n")

               ##read data
               annotation.result <- readr::read_csv(file.path(path, ms1.data),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
               annotation.result <- as.data.frame(annotation.result)
               sample.info <- readr::read_csv(file.path(path, old.sample.info.name),
                                              progress = FALSE,
                                              col_types = readr::cols())
               sample.info <- as.data.frame(sample.info)

               ###########################################################################
               #load tags2.data and rt.result
               #positive
               load(file.path(path, "MRN_annotation_result","intermediate_data","tags2.after.redundancy.remove"))
               load(file.path(path, "MRN_annotation_result","intermediate_data", "rt.result"))
               tags2 <- tags2.after.redundancy.remove
               rm(list="tags2.after.redundancy.remove")
               gc()


               ##############################################################################
               #Calculate p value and fold change
               if(any(!is.element(group, unique(sample.info[,2])))){
                 idx <- which(!is.element(group, unique(sample.info[,2])))
                 stop(paste(group[idx], collapse = " and "),
                      ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
               }

               sample.info <- sample.info[sample.info[,2] %in% group,]

               sample <- annotation.result[,match(sample.info[,1],
                                                  colnames(annotation.result))]

               sample <- as.data.frame(sample)
               rownames(sample) <- annotation.result$name

               cat("Calculate p values of peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Calculate p values of peaks.\n")

               if(ncol(sample) <= 3){
                 cat("At least one group has only one sample. So the p values are assigned 1.\n")
                 p.value <- rep(1, nrow(sample))
               }else{
                 p.value <- uniTest(
                   sample = as.matrix(sample),
                   sample.info = sample.info,
                   uni.test = uni.test,
                   correct = correct)
               }

               names(p.value) <- rownames(sample)
               p.value[is.na(p.value)] <- 1


               num.peak <- sum(p.value <= p.cutoff)
               cat("There are ", num.peak, " out of ",
                   length(p.value), " peaks with p values ",
                   ifelse(correct, '(FDR correct)', "(No correct)"),
                   " less than ",
                   p.cutoff, ".\n", sep = "", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("There are ", num.peak, " out of ",
                   length(p.value), " peaks with p values ",
                   ifelse(correct, '(FDR correct)', "(No correct)"),
                   " less than ",
                   p.cutoff, ".\n", sep = "")


               ##if peak number with p value less than p.cutoff is less than 10

               # count <- 1
               # while(num.peak <= 10 & count <= 5){
               #   count <- count + 1
               #   cat("\n")
               #   # cat("Warning: Peaks with p values < p.cutoff is less than 10, please change the correct method and the p.cutoff.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               #   # warning("Peaks with p values < p.cutoff is less than 10, please change the correct method and the p.cutoff.\n")
               #
               #   p.cutoff <- readline("p.cutoff: ")
               #   p.cutoff <- as.numeric(p.cutoff)
               #
               #   if(correct){
               #     cat("\n")
               #     correct <- readline("correct (Type 1 for TRUE or 0 for FALSE): ")
               #     while(correct != "1" & correct != "0"){
               #       cat("\n")
               #       correct <- readline("correct (Type 1 for TRUE or 0 for FALSE): ")
               #     }
               #     correct <- ifelse(correct == "1", TRUE, FALSE)
               #   }
               #
               #
               #   ##re-calculate p value
               #   cat("\n")
               #   cat("Re-calculate p values of peaks.\n")
               #   p.value <- uniTest(
               #     sample = as.matrix(sample),
               #     sample.info = sample.info,
               #     uni.test = uni.test,
               #     correct = correct)
               #   names(p.value) <- rownames(sample)
               #   p.value[is.na(p.value)] <- 1
               #
               #
               #   num.peak <- sum(p.value <= p.cutoff)
               #   cat("\n")
               #   cat("There are ", num.peak, " out of ",
               #       length(p.value), " peaks with p values ",
               #       ifelse(correct, '(FDR correct)', "(No correct)"),
               #       " less than ",
               #       p.cutoff, ".\n", sep = "")
               # }


               cat("\n")
               cat("Calculate fold change of peaks.\n")
               foldchange <- foldChange(sample = sample,
                                        sample.info = sample.info,
                                        by.what = "median",
                                        group = group)
               names(foldchange) <- rownames(sample)

               save(p.value, file = file.path(output.path, "intermediate_data", "p.value"))
               save(foldchange, file = file.path(output.path, "intermediate_data", "foldchange"))


               ###volcano plot
               pdf(file = file.path(output.path, "volcano.plot.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               volcanoPlot(p.value = p.value, fc = foldchange, correct = correct,
                           p.cutoff = p.cutoff, fc.cutoff = 1, cex = 0.3)
               dev.off()

               if(polarity == "positive"){
                 sample.pos <- sample
                 sample.neg <- NULL
                 pos.path <- path
                 neg.path <- NULL
               }else{
                 sample.neg <- sample
                 sample.pos <- NULL
                 neg.path <- path
                 pos.path <- NULL
               }


               if(sum(p.value < p.cutoff) >= 10){
                 ##Do the user provide the annotation or not.
                   if(any(dir(path) == "marker.csv")) {
                     marker <- read.csv(file.path(path, "marker.csv"), stringsAsFactors = FALSE)
                   }else{
                     marker <- NULL
                   }


                 if(is.null(marker)){
                   cat("You don't provide the markers.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                   cat("You don't provide the markers.\n")

                   index <- which(p.value < p.cutoff)
                   temp.tags2 <- tags2[index]

                   marker <- tags2Result(tags2 = temp.tags2,
                                          candidate.num = candidate.num)
                   marker <- as.data.frame(marker)
                   marker <- marker[which(!is.na(marker$to)),]
                   marker <- marker[,c("name", "to")]
                   colnames(marker) <- c("peak.name", "KEGG.ID")
                 }

                 temp.error <- try({metEnrichment(annotation.table = marker,
                               sample.pos = sample.pos,
                               sample.neg = sample.neg,
                               group = group,
                               column = column,
                               candidate.num = candidate.num,
                               sample.info = sample.info,
                               polarity = polarity,
                               species = species,
                               pos.path = pos.path,
                               neg.path = neg.path,
                               output.path = output.path)}, silent = TRUE)

                 if(class(temp.error) == "try-error") {
                   cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
                   # stop(temp.error[[1]])
                 }
               }
             }else{
               cat("Skip pathway enrichment analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Skip pathway enrichment analysis.\n")
             }

             ##analysis report
             if(analysis.report){
               if(any(dir(file.path(output.path, "intermediate_data")) == "module.result")){
                 cat("\n")
                 cat("6. Analysis report generation.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("6. Analysis report generation.\n")
                 cat("------------------------------------------------------------\n")

                 if(polarity == "positive"){
                   pos.path <- path
                   neg.path <- NULL
                 }else{
                   neg.path <- path
                   pos.path <- NULL
                 }

                 #
                 temp.error <- try({analysisGeneration(polarity = polarity,
                                                       pos.path = pos.path,
                                                       neg.path = neg.path,
                                                       sample.info = old.sample.info.name,
                                                       data.pos = old.data.name,
                                                       data.neg = old.data.name,
                                                       correct = correct,
                                                       p.cutoff = p.cutoff,
                                                       candidate.num = candidate.num,
                                                       type = "html",
                                                       output.path = ifelse(polarity == "positive", pos.path, neg.path))},
                                   silent = TRUE)

                 if(class(temp.error) == "try-error") {
                   cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
                 }
               }

               # ##analysis report

               cat("\n")
               cat("6. Analysis report generation.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("6. Analysis report generation.\n")
               cat("------------------------------------------------------------\n")

               if(polarity == "positive"){
                 pos.path <- path
                 neg.path <- NULL
               }else{
                 neg.path <- path
                 pos.path <- NULL
               }

               temp.error <- try({analysisGeneration(polarity = polarity,
                                                     pos.path = pos.path,
                                                     neg.path = neg.path,
                                                     sample.info = old.sample.info.name,
                                                     data.pos = old.data.name,
                                                     data.neg = old.data.name,
                                                     correct = correct,
                                                     p.cutoff = p.cutoff,
                                                     candidate.num = candidate.num,
                                                     type = "html",
                                                     output.path = ifelse(polarity == "positive", pos.path, neg.path))}, silent = TRUE)

               if(class(temp.error) == "try-error") {
                 cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
               }
             }







           })









