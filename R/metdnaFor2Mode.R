# system.time(metdnaFor2Mode(ms1.data.pos = "data.csv",
#                            ms1.data.neg = "data.csv",
#                            sample.info.pos = "sample.info.csv",
#                            sample.info.neg = "sample.info.csv",
#                            mz.tol = 25,
#                            rt.tol.for.ms1.ms2.match = 10,#second
#                            path = "/home/jasper/Downloads/7f68899f/POS and NEG",
#                            pos.path = "/home/jasper/Downloads/7f68899f/POS",
#                            neg.path = "/home/jasper/Downloads/7f68899f/NEG",
#                            ms2.type = "mgf",
#                            column = "rp",
#                            ce = "30",
#                            ###
#                            group = c("Normal", "AD"),
#                            uni.test = "t",
#                            correct = FALSE,
#                            p.cutoff = 0.1,
#                            use.old.null = FALSE,
#                            species = "hsa",
#                            dn.analysis = FALSE,
#                            pathway.enrichment = FALSE))


setGeneric(name = "metdnaFor2Mode",
           def = function(ms1.data.pos = "data.csv",
                          ms1.data.neg = "data.csv",
                          sample.info.pos = "sample.info.csv",
                          sample.info.neg = "sample.info.csv",
                          mz.tol = 25,
                          rt.tol.for.ms1.ms2.match = 10,#second
                          path = ".",
                          log.path = ".",
                          pos.path = ".",
                          neg.path = ".",
                          ms2.type = c("mgf", "msp", "mzXML"),
                          column = c("hilic", "rp"),
                          ce = c("30", "10", "15", "20", "25", "35",
                                 "35,15", "40", "45", "50", "55",
                                 "60", "65", "70"),
                          ###
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
                          ###
                          group,
                          uni.test = c("t", "wilcox"),
                          correct = TRUE,
                          p.cutoff = 0.01,
                          use.old.null = FALSE,
                          candidate.num = 5,
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          use.all.kegg.id = FALSE,
                          only.mrn.annotation = FALSE,
                          dn.analysis = FALSE,
                          pathway.enrichment = TRUE
           ){


             #--------------------------------------------------------------------
             #parameter decision
             ms2.type <- match.arg(ms2.type)
             column <- match.arg(column)
             ce <- match.arg(ce)
             uni.test <- match.arg(uni.test)
             species <- match.arg(species)

             old.sample.info.pos.name <- sample.info.pos
             old.sample.info.neg.name <- sample.info.neg

             ###metModule
             cat("\n")
             cat("Dysregulated network analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
             cat("Dysregulated network analysis.\n")

             output.path = file.path(path, "Dysregulated_network_analysis_result")
             dir.create(output.path)
             dir.create(file.path(output.path, "intermediate_data"))

             cat("------------------------------------------------------------\n")
             if(dn.analysis){
               cat("Module analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Module analysis.\n")
               cat("\n")
               ##get the directory
               ##read data
                 annotation.result.pos <- readr::read_csv(file.path(pos.path, ms1.data.pos),
                                                          progress = FALSE,
                                                          col_types = readr::cols())
                 annotation.result.pos <- as.data.frame(annotation.result.pos)
                 sample.info.pos <- readr::read_csv(file.path(pos.path, old.sample.info.pos.name),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
                 sample.info.pos <- as.data.frame(sample.info.pos)


                 annotation.result.neg <- readr::read_csv(file.path(neg.path, ms1.data.neg),
                                                          progress = FALSE,
                                                          col_types = readr::cols())
                 annotation.result.neg <- as.data.frame(annotation.result.neg)
                 sample.info.neg <- readr::read_csv(file.path(neg.path, old.sample.info.neg.name),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
                 sample.info.neg <- as.data.frame(sample.info.neg)


                 if(any(sort(sample.info.pos$sample.name) != sort(sample.info.neg$sample.name))){
                   cat("Error: Sample names in positive and negative mode are not same.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                   stop("Sample names in positive and negative mode are not same.\n")
                 }
                 sample.info <- sample.info.pos



               ###########################################################################
               #load tags2.data and rt.result
               #positive
                 load(file.path(pos.path, "MRN_annotation_result","intermediate_data","tags2.after.redundancy.remove"))
                 load(file.path(pos.path, "MRN_annotation_result","intermediate_data", "rt.result"))
                 tags2.pos <- tags2.after.redundancy.remove
                 rm(list="tags2.after.redundancy.remove")
                 gc()
                 #change peak name
                   tags2.pos <-
                     lapply(tags2.pos, function(x){
                       x@name <- paste(x@name, "POS", sep = "_")
                       x
                     })

               ##negative
                 load(file.path(neg.path, "MRN_annotation_result", "intermediate_data", "tags2.after.redundancy.remove"))
                 load(file.path(neg.path, "MRN_annotation_result","intermediate_data", "rt.result"))
                 tags2.neg <- tags2.after.redundancy.remove
                 rm(list="tags2.after.redundancy.remove")
                 gc()

                   #change peak name
                   tags2.neg <-
                     lapply(tags2.neg, function(x){
                       x@name <- paste(x@name, "NEG", sep = "_")
                       x
                     })

               ##############################################################################
               #Calculate p value and fold change
               if(any(!is.element(group, unique(sample.info[,2])))){
                 idx <- which(!is.element(group, unique(sample.info[,2])))
                 stop(paste(group[idx], collapse = " and "),
                      ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
               }

               sample.info <- sample.info[sample.info[,2] %in% group,]

                 sample.pos <- annotation.result.pos[,match(sample.info[,1],
                                                            colnames(annotation.result.pos))]

                 sample.pos <- as.data.frame(sample.pos)
                 rownames(sample.pos) <- paste(annotation.result.pos$name, "POS", sep = "_")


                 sample.neg <- annotation.result.neg[,match(sample.info[,1],
                                                            colnames(annotation.result.neg))]
                 sample.neg <- as.data.frame(sample.neg)
                 rownames(sample.neg) <- paste(annotation.result.neg$name, "NEG", sep = "_")

                cat("\n")
                cat("Calculate p values of positive peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of positive peaks.\n")
                 p.value.pos <- uniTest(
                   sample = as.matrix(sample.pos),
                   sample.info = sample.info,
                   uni.test = uni.test,
                   correct = correct)
                 names(p.value.pos) <- rownames(sample.pos)
                 p.value.pos[is.na(p.value.pos)] <- 1

                 cat("\n")
                 cat("Calculate p values of negative peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of negative peaks.\n")
                 p.value.neg<- uniTest(
                   sample = as.matrix(sample.neg),
                   sample.info = sample.info,
                   uni.test = uni.test,
                   correct = correct)
                 names(p.value.neg) <- rownames(sample.neg)
                 p.value.neg[is.na(p.value.neg)] <- 1


               p.value <- c(p.value.pos, p.value.neg)

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


               ##if peak number with p value less than p.cutoff is less than 30
               count <- 1
               while(num.peak <= 10 & count <= 5){
                 count <- count + 1
                 cat("\n")
                 cat("Error: Peaks with p values < p.cutoff is less than 10, please change the correct method and the p.cutoff.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 stop("Peaks with p values < p.cutoff is less than 10, please change the correct method and the p.cutoff.\n")

                 p.cutoff <- readline("p.cutoff: ")
                 p.cutoff <- as.numeric(p.cutoff)

                 if(correct){
                   cat("\n")
                   correct <- readline("correct (Type 1 for TRUE or 0 for FALSE): ")
                   while(correct != "1" & correct != "0"){
                     cat("\n")
                     correct <- readline("correct (Type 1 for TRUE or 0 for FALSE): ")
                   }
                   correct <- ifelse(correct == "1", TRUE, FALSE)
                 }


                 ##re-calculate p value
                   cat("\n")
                   cat("Re-calculate p values of positive peaks.\n")
                   p.value.pos <- uniTest(
                     sample = as.matrix(sample.pos),
                     sample.info = sample.info,
                     uni.test = uni.test,
                     correct = correct)
                   names(p.value.pos) <- rownames(sample.pos)
                   p.value.pos[is.na(p.value.pos)] <- 1


                   cat("\n")
                   cat("Re-calculate p values of negative peaks.\n")
                   p.value.neg <- uniTest(
                     sample = as.matrix(sample.neg),
                     sample.info = sample.info,
                     uni.test = uni.test,
                     correct = correct)
                   names(p.value.neg) <- rownames(sample.neg)
                   p.value.neg[is.na(p.value.neg)] <- 1

                   p.value <- c(p.value.pos, p.value.neg)


                 num.peak <- sum(p.value <= p.cutoff)
                 cat("\n")
                 cat("There are ", num.peak, " out of ",
                     length(p.value), " peaks with p values ",
                     ifelse(correct, '(FDR correct)', "(No correct)"),
                     " less than ",
                     p.cutoff, ".\n", sep = "")
               }

                 cat("\n")
                 cat("Calculate fold change of positive peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of positive peaks.\n")
                 foldchange.pos <- foldChange(sample = sample.pos,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange.pos) <- rownames(sample.pos)


                 cat("\n")
                 cat("Calculate fold change of negative peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of negative peaks.\n")
                 foldchange.neg <- foldChange(sample = sample.neg,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange.neg) <- rownames(sample.neg)

                 foldchange <- c(foldchange.pos, foldchange.neg)


                 save(p.value.pos, file = file.path(output.path, "intermediate_data", "p.value.pos"))
                 save(foldchange.pos, file = file.path(output.path, "intermediate_data", "foldchange.pos"))
                 save(p.value.neg, file = file.path(output.path, "intermediate_data", "p.value.neg"))
                 save(foldchange.neg, file = file.path(output.path, "intermediate_data", "foldchange.neg"))


               ###volcano plot
               pdf(file = file.path(output.path, "volcano.plot.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               volcanoPlot(p.value = p.value, fc = foldchange, correct = correct,
                           p.cutoff = p.cutoff, fc.cutoff = 1, cex = 0.3)
               dev.off()

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
                           polarity = "both",
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
                 annotation.result.pos <- readr::read_csv(file.path(pos.path, ms1.data.pos),
                                                          progress = FALSE,
                                                          col_types = readr::cols())
                 annotation.result.pos <- as.data.frame(annotation.result.pos)
                 sample.info.pos <- readr::read_csv(file.path(pos.path, old.sample.info.pos.name),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
                 sample.info.pos <- as.data.frame(sample.info.pos)

                 annotation.result.neg <- readr::read_csv(file.path(neg.path, ms1.data.neg),
                                                          progress = FALSE,
                                                          col_types = readr::cols())
                 annotation.result.neg <- as.data.frame(annotation.result.neg)
                 sample.info.neg <- readr::read_csv(file.path(neg.path, old.sample.info.neg.name),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
                 sample.info.neg <- as.data.frame(sample.info.neg)


                 if(any(sort(sample.info.pos$sample.name) != sort(sample.info.neg$sample.name))){
                   cat("Error: Sample names in positive and negative mode are not same.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                   stop("Sample names in positive and negative mode are not same.\n")
                 }
                 sample.info <- sample.info.pos


               #Calculate p value and fold change
               if(any(!is.element(group, unique(sample.info[,2])))){
                 idx <- which(!is.element(group, unique(sample.info[,2])))

                 cat(paste("Error: ",group[idx], collapse = " and "),
                     ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.", file = file.path(log.path, "run.log.txt"), append = TRUE)

                 stop(paste(group[idx], collapse = " and "),
                      ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
               }

               sample.info <- sample.info[sample.info[,2] %in% group,]

                 sample.pos <- annotation.result.pos[,match(sample.info[,1],
                                                            colnames(annotation.result.pos))]

                 sample.pos <- as.data.frame(sample.pos)
                   rownames(sample.pos) <- paste(annotation.result.pos$name, "POS", sep = "_")



                 sample.neg <- annotation.result.neg[,match(sample.info[,1],
                                                            colnames(annotation.result.neg))]
                 sample.neg <- as.data.frame(sample.neg)
                 rownames(sample.neg) <- paste(annotation.result.neg$name, "NEG", sep = "_")

                 cat("\n")
                 cat("Calculate p values of positive peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of positive peaks.\n")
                 p.value.pos <- uniTest(
                   sample = as.matrix(sample.pos),
                   sample.info = sample.info,
                   uni.test = uni.test,
                   correct = correct)
                 names(p.value.pos) <- rownames(sample.pos)
                 p.value.pos[is.na(p.value.pos)] <- 1

                 cat("\n")
                 cat("Calculate p values of negative peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of negative peaks.\n")
                 p.value.neg<- uniTest(
                   sample = as.matrix(sample.neg),
                   sample.info = sample.info,
                   uni.test = uni.test,
                   correct = correct)
                 names(p.value.neg) <- rownames(sample.neg)
                 p.value.neg[is.na(p.value.neg)] <- 1

                 p.value <- c(p.value.pos, p.value.neg)

                 cat("\n")
                 cat("Calculate fold change of positive peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of positive peaks.\n")
                 foldchange.pos <- foldChange(sample = sample.pos,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange.pos) <- rownames(sample.pos)

                 cat("\n")
                 cat("Calculate fold change of negative peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of negative peaks.\n")
                 foldchange.neg <- foldChange(sample = sample.neg,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange.neg) <- rownames(sample.neg)

                 foldchange <- c(foldchange.pos, foldchange.neg)

                 save(p.value.pos, file = file.path(output.path, "intermediate_data", "p.value.pos"))
                 save(foldchange.pos, file = file.path(output.path, "intermediate_data", "foldchange.pos"))
                 save(p.value.neg, file = file.path(output.path, "intermediate_data", "p.value.neg"))
                 save(foldchange.neg, file = file.path(output.path, "intermediate_data", "foldchange.neg"))
             }

             ###pathway enrichment analysis
             cat("\n")
             cat("Pathway enrichment analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
             cat("Pathway enrichment analysis.\n")

             output.path <- file.path(path,"Pathway_enrichment_analysis_result")
             dir.create(output.path)
             dir.create(file.path(output.path, "intermediate_data"))

             cat("------------------------------------------------------------\n")
             if(pathway.enrichment){
               cat("\n")
               cat("Pathway enrichment analysis.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Pathway enrichment analysis.\n")
               cat("\n")
               ##get the directory
               ##read data
                 annotation.result.pos <- readr::read_csv(file.path(pos.path, ms1.data.pos),
                                                          progress = FALSE,
                                                          col_types = readr::cols())
                 annotation.result.pos <- as.data.frame(annotation.result.pos)
                 annotation.result.pos$name <- paste(annotation.result.pos$name, "POS", sep = "_")

                 sample.info.pos <- readr::read_csv(file.path(pos.path, old.sample.info.pos.name),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
                 sample.info.pos <- as.data.frame(sample.info.pos)


                 annotation.result.neg <- readr::read_csv(file.path(neg.path, ms1.data.neg),
                                                          progress = FALSE,
                                                          col_types = readr::cols())
                 annotation.result.neg <- as.data.frame(annotation.result.neg)
                 annotation.result.neg$name <- paste(annotation.result.neg$name, "NEG", sep = "_")

                 sample.info.neg <- readr::read_csv(file.path(neg.path, old.sample.info.neg.name),
                                                    progress = FALSE,
                                                    col_types = readr::cols())
                 sample.info.neg <- as.data.frame(sample.info.neg)


                 if(any(sort(sample.info.pos$sample.name) != sort(sample.info.neg$sample.name))){
                   cat("Error: Sample names in positive and negative mode are not same.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                   stop("Sample names in positive and negative mode are not same.\n")
                 }
                 sample.info <- sample.info.pos



               ###########################################################################
               #load tags2.data
               #positive
                 load(file.path(pos.path, "MRN_annotation_result","intermediate_data","tags2.after.redundancy.remove"))
                 tags2.pos <- tags2.after.redundancy.remove
                 rm(list="tags2.after.redundancy.remove")
                 gc()
                 #change peak name
                   tags2.pos <-
                     lapply(tags2.pos, function(x){
                       x@name <- paste(x@name, "POS", sep = "_")
                       x
                     })


                 ##remove the peaks which have no annotations
                 remove.index.pos <- which(showTags2(tags2 = tags2.pos, slot = "annotation.len") == 0)
                 remove.name.pos <- showTags2(tags2 = tags2.pos, slot = "name")[remove.index.pos]
                 tags2.pos <- tags2.pos[-remove.index.pos]
                 annotation.result.pos <- annotation.result.pos[-which(annotation.result.pos$name %in% remove.name.pos),]

               ##negative
                 load(file.path(neg.path, "MRN_annotation_result", "intermediate_data", "tags2.after.redundancy.remove"))
                 tags2.neg <- tags2.after.redundancy.remove
                 rm(list="tags2.after.redundancy.remove")
                 gc()

                   #change peak name
                   tags2.neg <-
                     lapply(tags2.neg, function(x){
                       x@name <- paste(x@name, "NEG", sep = "_")
                       x
                     })


                 ##remove the peaks which have no annotations
                   remove.index.neg <- which(showTags2(tags2 = tags2.neg, slot = "annotation.len") == 0)
                   remove.name.neg <- showTags2(tags2 = tags2.neg, slot = "name")[remove.index.neg]
                   tags2.neg <- tags2.neg[-remove.index.neg]
                   annotation.result.neg <- annotation.result.neg[-which(annotation.result.neg$name %in% remove.name.neg),]

               ##############################################################################
               #Calculate p value and fold change
               if(any(!is.element(group, unique(sample.info[,2])))){
                 idx <- which(!is.element(group, unique(sample.info[,2])))

                 cat(paste("Error: ",group[idx], collapse = " and "),
                     ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.", file = file.path(log.path, "run.log.txt"), append = TRUE)

                 stop(paste(group[idx], collapse = " and "),
                      ifelse(length(idx) > 1, " are ", " is "), "not in your sample.info.")
               }

               sample.info <- sample.info[sample.info[,2] %in% group,]

                 sample.pos <- annotation.result.pos[,match(sample.info[,1],
                                                            colnames(annotation.result.pos))]

                 sample.pos <- as.data.frame(sample.pos)
                 rownames(sample.pos) <- annotation.result.pos$name

                 sample.neg <- annotation.result.neg[,match(sample.info[,1],
                                                            colnames(annotation.result.neg))]
                 sample.neg <- as.data.frame(sample.neg)
                 rownames(sample.neg) <- annotation.result.neg$name
                 cat("\n")
                 cat("Calculate p values of positive peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of positive peaks.\n")
                 p.value.pos <- uniTest(
                   sample = as.matrix(sample.pos),
                   sample.info = sample.info,
                   uni.test = uni.test,
                   correct = correct)
                 names(p.value.pos) <- rownames(sample.pos)
                 p.value.pos[is.na(p.value.pos)] <- 1
                 cat("\n")
                 cat("Calculate p values of negative peaks.\n", file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate p values of negative peaks.\n")
                 p.value.neg<- uniTest(
                   sample = as.matrix(sample.neg),
                   sample.info = sample.info,
                   uni.test = uni.test,
                   correct = correct)
                 names(p.value.neg) <- rownames(sample.neg)
                 p.value.neg[is.na(p.value.neg)] <- 1

                 p.value <- c(p.value.pos, p.value.neg)

                 tags2 <- c(tags2.pos, tags2.neg)

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


               ##if peak number with p value less than p.cutoff is less than 30
               # count <- 1
               # while(num.peak <= 10 & count <= 5){
               #   count <- count + 1
               #   cat("\n")
               #   cat("Error: Peaks with p values < p.cutoff is less than 10, please change the correct method and the p.cutoff.\n",
               #       file = file.path(log.path, "run.log.txt"), append = TRUE)
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
               #     cat("Re-calculate p values of positive peaks.\n")
               #     p.value.pos <- uniTest(
               #       sample = as.matrix(sample.pos),
               #       sample.info = sample.info,
               #       uni.test = uni.test,
               #       correct = correct)
               #     names(p.value.pos) <- rownames(sample.pos)
               #     p.value.pos[is.na(p.value.pos)] <- 1
               #
               #     cat("\n")
               #     cat("Re-calculate p values of negative peaks.\n")
               #     p.value.neg <- uniTest(
               #       sample = as.matrix(sample.neg),
               #       sample.info = sample.info,
               #       uni.test = uni.test,
               #       correct = correct)
               #     names(p.value.neg) <- rownames(sample.neg)
               #     p.value.neg[is.na(p.value.neg)] <- 1
               #
               #     p.value <- c(p.value.pos, p.value.neg)
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
                 cat("Calculate fold change of positive peaks.\n",
                     file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of positive peaks.\n")
                 foldchange.pos <- foldChange(sample = sample.pos,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange.pos) <- rownames(sample.pos)


                 cat("\n")
                 cat("Calculate fold change of negative peaks.\n",
                     file = file.path(log.path, "run.log.txt"), append = TRUE)
                 cat("Calculate fold change of negative peaks.\n")
                 foldchange.neg <- foldChange(sample = sample.neg,
                                              sample.info = sample.info,
                                              by.what = "median",
                                              group = group)
                 names(foldchange.neg) <- rownames(sample.neg)


               foldchange <- c(foldchange.pos, foldchange.neg)

                 save(p.value.pos, file = file.path(output.path, "intermediate_data", "p.value.pos"))
                 save(foldchange.pos, file = file.path(output.path, "intermediate_data", "foldchange.pos"))
                 save(p.value.neg, file = file.path(output.path, "intermediate_data", "p.value.neg"))
                 save(foldchange.neg, file = file.path(output.path, "intermediate_data", "foldchange.neg"))


               ###volcano plot
               pdf(file = file.path(output.path, "volcano.plot.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               volcanoPlot(p.value = p.value, fc = foldchange, correct = correct,
                           p.cutoff = p.cutoff, fc.cutoff = 1, cex = 0.3)
               dev.off()

               if(sum(p.value < p.cutoff) >= 10){
                 ##Do the user provide the annotation or not.
                   if(any(dir(pos.path) == "marker.csv")) {
                     marker.pos <- read.csv(file.path(pos.path, "marker.csv"), stringsAsFactors = FALSE)
                     marker.pos[,1] <- paste(marker.pos[,1], "POS", sep = "_")
                   }else{
                     marker.pos <- NULL
                   }

                   if(any(dir(neg.path) == "marker.csv")) {
                     marker.neg <- read.csv(file.path(neg.path, "marker.csv"), stringsAsFactors = FALSE)
                     marker.neg[,1] <- paste(marker.neg[,1], "NEG", sep = "_")
                   }else{
                     marker.neg <- NULL
                   }


                 marker <- rbind(marker.pos, marker.neg)

                 if(is.null(marker)){
                   cat("\n")
                   cat("You don't provide the markers.\n",
                       file = file.path(log.path, "run.log.txt"), append = TRUE)
                   cat("You don't provide the markers.\n")
                   index <- which(p.value < p.cutoff)
                   temp.tags2 <- tags2[index]
                   marker <- tags2Result(tags2 = temp.tags2,
                                         candidate.num = candidate.num)
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
                               polarity = "both",
                               species = species,
                               pos.path = pos.path,
                               neg.path = neg.path,
                               output.path = output.path)}, silent = TRUE)
                    if(class(temp.error) == "try-error") {
                      cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
                      stop(temp.error[[1]])
                    }

               }
             }else{
               cat("Skip pathway enrichment analysis.\n",
                   file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Skip pathway enrichment analysis.\n")
             }



             ##analysis report


             if(any(dir(file.path(output.path, "intermediate_data")) == "module.result")){
               cat("\n")
               cat("Analysis report generation.\n",
                   file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Analysis report generation.\n")
               cat("------------------------------------------------------------\n")
               #

               temp.error <- try({analysisGeneration(polarity = "both",
                              pos.path = pos.path,
                              neg.path = neg.path,
                              output.path = path,
                              sample.info = old.sample.info.pos.name,
                              data.pos = ms1.data.pos,
                              data.neg = ms1.data.neg,
                              correct = correct,
                              p.cutoff = p.cutoff,
                              candidate.num = candidate.num,
                              type = "html")}, silent = TRUE)

               if(class(temp.error) == "try-error") {
                 cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
               }

             }


             # ##analysis report
             # if(any(dir(output.path) == "Pathway.enrichment.analysis.csv")){
               cat("\n")
               cat("Analysis report generation.\n",
                   file = file.path(log.path, "run.log.txt"), append = TRUE)
               cat("Analysis report generation.\n")
               cat("------------------------------------------------------------\n")

               temp.error <- try({analysisGeneration(polarity = "both",
                                                 pos.path = pos.path,
                                                 neg.path = neg.path,
                                                 output.path = path,
                                                 sample.info = old.sample.info.pos.name,
                                                 data.pos = ms1.data.pos,
                                                 data.neg = ms1.data.neg,
                                                 correct = correct,
                                                 p.cutoff = p.cutoff,
                                                 candidate.num = candidate.num,
                                                 type = "html")}, silent = TRUE)

               if(class(temp.error) == "try-error") {
                 cat(temp.error[[1]], file = file.path(log.path, "run.log.txt"), append = TRUE)
               }
             # }
           })









