#' @title analysisGeneration
#' @description Generate analysis report
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param polarity The polarity
#' @param pos.path The directory of positive data.
#' @param neg.path The directory of negative data.
#' @param data.pos The name of positive data.
#' @param data.neg The name of negative data.
#' @param sample.info The name of sample.info.
#' @param output.path Output directory.
#' @param correct Correct p value or not.
#' @param p.cutoff P value cutoff.
#' @param candidate.num Top candidate.num metabolites for each peak.
#' @param type Analysis report type.
#' @return Analysis report.
#' @export
#


# analysisGeneration(polarity = "both",
#                    pos.path = pos.path,
#                    neg.path = neg.path,
#                    output.path = path,
#                    sample.info = old.sample.info.pos.name,
#                    data.pos = ms1.data.pos,
#                    data.neg = ms1.data.neg,
#                    correct = correct,
#                    p.cutoff = p.cutoff,
#                    type = "html")


setGeneric(name = "analysisGeneration",
           def = function(
             polarity = c("positive", "negative", "both"),
             pos.path = ".",
             neg.path = ".",
             data.pos = "data.csv",
             data.neg = "data.csv",
             sample.info = "sample.info.csv",
             output.path = ".",
             correct = TRUE,
             p.cutoff = 0.01,
             candidate.num = candidate.num,
             type = c("html", "pdf", "all")){

             options(warn = -1)
             polarity <- match.arg(polarity)
             type <- match.arg(type)

             ###plot a null plot for volcanoplot.and.pathway.overview.png
             dir.create(file.path(output.path, "Analysis_report"))
             if(all(dir(file.path(output.path, "Analysis_report")) != "volcanoplot.and.pathway.overview.png")){
               png(filename = file.path(output.path, "Analysis_report/volcanoplot.and.pathway.overview.png"),
                   width = 14, height = 7, res = 200, units = "in")
               plot(0 ,xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                    col = "white")
               legend("topleft", legend = "No pathway enrichment analysis.", cex = 1.5, bty = "n")
               dev.off()
             }

             ##get the directory
             now.path <- ifelse(polarity == "positive", pos.path, neg.path)
             now.path1 <- strsplit(now.path, split = "/")[[1]]
             now.path1 <- paste(now.path1[-length(now.path1)], collapse = "/")
             # if(tail(now.path1, 1) == "POS" | tail(now.path1, 1) == "NEG" | tail(now.path1, 1) == "POS and NEG"){
             #   new.path <- paste(now.path1[-length(now.path1)], collapse = "/")
             # }else{
             #   new.path <- now.path
             # }
             # neg.path <- file.path(new.path, "NEG")
             path <- file.path(output.path, "Dysregulated_network_analysis_result")
             path1 <- file.path(output.path, "Pathway_enrichment_analysis_result")

             output.path <- file.path(output.path, "Analysis_report")

             dir.create(output.path)

             ####get the template of analysisGeneration
             if(polarity == "both"){
               rmarkdown::draft(file.path(output.path, "analysis_report"),
                                template = "MetDNA_POS_and_NEG", package = "analysisReport",
                                create_dir = TRUE, edit = FALSE)

               file.copy(from =
                           file.path(output.path, "analysis_report",
                                     c("MetDNA.logo.png", "MetDNA.template.Rmd")),
                         output.path, overwrite = TRUE)

               unlink(file.path(output.path, "analysis_report"), recursive = TRUE)
             }else{
               rmarkdown::draft(file.path(output.path, "analysis_report"),
                                template = "MetDNA_POS_or_NEG", package = "analysisReport",
                                create_dir = TRUE, edit = FALSE)

               file.copy(from =
                           file.path(output.path, "analysis_report",
                                     c("MetDNA.logo.png", "MetDNA.template.Rmd")),
                         output.path, overwrite = TRUE)

               unlink(file.path(output.path, "analysis_report"), recursive = TRUE)
             }

             #------------------------------------------------------------------
             ###parameters
             if(polarity == "both"){
               file.copy(from = file.path(pos.path, "MetDNA.parameters.csv"),
                         to = file.path(output.path, "MetDNA.parameters.POS.csv"), overwrite = TRUE)
               # file.copy(from = file.path(neg.path, "MetDNA.parameters.csv"),
               #           to = file.path(output.path, "MetDNA.parameters.NEG.csv"), overwrite = TRUE)

               parameter.pos <- read.csv(file.path(output.path, "MetDNA.parameters.POS.csv"), stringsAsFactors = FALSE)
               # parameter.neg <- read.csv(file.path(output.path, "MetDNA.parameters.NEG.csv"), stringsAsFactors = FALSE)
               parameter.pos <- parameter.pos[c(5:7, 13:18),-1]
               # parameter.neg <- parameter.neg[c(5:7, 13:18),-1]
               meaning <- c("Ionization Polarity",
                            "Liquid Chromatograph",
                            "Collision energy",
                            "Group of your data",
                            "Univariate Statistics",
                            "P-value Adjustment",
                            "Cutoff of P-value",
                            "Species",
                            "MS Instrument")
               parameter.pos$parameter <- meaning
               # parameter.neg$parameter <- meaning
               colnames(parameter.pos) <- c("Patameter", "Setting")
               # colnames(parameter.neg) <- c("Patameter", "Value")
               # rownames(parameter.pos) <- rownames(parameter.neg) <- NULL
               rownames(parameter.pos) <- NULL
               parameter <- parameter.pos
               parameter[1,2] <- "both"
               parameter <- parameter[c(1,2,9,3,4,5,8,7,6),]
               group <- parameter[5,2]
               group <- strsplit(x = group, split = ",")[[1]]
               ctr.group <- group[1]
               case.group <- group[2]
               group <- data.frame(c("Control Group", "Case Group"),c(ctr.group, case.group),
                                   stringsAsFactors = FALSE)
               colnames(group) <- c("Patameter", "Setting")
               parameter <- rbind(parameter[1:4,], group, parameter[6:9,])
               rownames(parameter) <- NULL

               parameter[1,2] <- parameterMap(name = parameter[1,2])
               parameter[2,2] <- parameterMap(name = parameter[2,2])
               parameter[3,2] <- parameterMap(name = parameter[3,2])
               parameter[7,2] <- parameterMap(name = parameter[7,2])
               parameter[8,2] <- parameterMap(name = parameter[8,2])
               parameter[10,2] <- parameterMap(name = parameter[10,2])


               save(parameter, file = file.path(output.path, "parameter"), compress = "xz")
               # save(parameter.neg, file = file.path(output.path, "parameter.neg"), compress = "xz")

               unlink(x = file.path(output.path, "MetDNA.parameters.POS.csv"), recursive = TRUE)
               # unlink(x = file.path(output.path, "MetDNA.parameters.NEG.csv"), recursive = TRUE)
               # parameter <- parameter.pos
             }else{
               file.copy(from = file.path(ifelse(polarity=="positive", pos.path, neg.path),"MetDNA.parameters.csv"), to = output.path)
               parameter <- read.csv(file.path(output.path, "MetDNA.parameters.csv"), stringsAsFactors = FALSE)
               parameter <- parameter[c(5:7, 13:18),-1]
               meaning <- c("Ionization Polarity",
                            "Liquid Chromatograph",
                            "Collision energy",
                            "Group of your data",
                            "Univariate Statistics",
                            "P-value Adjustment",
                            "Cutoff of P-value",
                            "Species",
                            "MS Instrument")
               parameter$parameter <- meaning
               colnames(parameter) <- c("Patameter", "Setting")
               rownames(parameter) <- NULL
               parameter <- parameter[c(1,2,9,3,4,5,8,7,6),]
               group <- parameter[5,2]
               group <- strsplit(x = group, split = ",")[[1]]
               ctr.group <- group[1]
               case.group <- group[2]
               group <- data.frame(c("Control Group", "Case Group"),c(ctr.group, case.group),
                                   stringsAsFactors = FALSE)
               colnames(group) <- c("Patameter", "Setting")
               parameter <- rbind(parameter[1:4,], group, parameter[6:9,])
               rownames(parameter) <- NULL
               parameter[1,2] <- parameterMap(name = parameter[1,2])
               parameter[2,2] <- parameterMap(name = parameter[2,2])
               parameter[3,2] <- parameterMap(name = parameter[3,2])
               parameter[7,2] <- parameterMap(name = parameter[7,2])
               parameter[8,2] <- parameterMap(name = parameter[8,2])
               parameter[10,2] <- parameterMap(name = parameter[10,2])
               save(parameter, file = file.path(output.path, "parameter"), compress = "xz")
               unlink(x = file.path(output.path, "MetDNA.parameters.csv"), recursive = TRUE)
             }




             #-----------------------------------------------------------------------
             ##data information
             dataInformation(pos.path = pos.path,
                             neg.path = neg.path,
                             sample.info = sample.info,
                             data.pos = data.pos,
                             data.neg = data.neg,
                             output.path = output.path,
                             polarity = polarity)

             #-----------------------------------------------------------------------
             ##redundancy remove
             if(polarity == "both"){
               load(file.path(pos.path, "MRN_annotation_result","intermediate_data","redun1"))
               redun1.pos <- do.call(rbind, redun1)
               load(file.path(neg.path, "MRN_annotation_result","intermediate_data","redun1"))
               redun1.neg <- do.call(rbind, redun1)

               pdf(file = file.path(output.path, "redundancy.removal.in.metABM.POS.pdf"),
                   width = 7, height = 7)
               par(mar = c(5,5,4,2))
               redundancyPlot(redun = redun1.pos)
               dev.off()

               pdf(file = file.path(output.path, "redundancy.removal.in.metABM.NEG.pdf"),
                   width = 7, height = 7)
               par(mar = c(5,5,4,2))
               redundancyPlot(redun = redun1.neg)
               dev.off()


               ######################
               load(file.path(pos.path, "MRN_annotation_result","intermediate_data","tags.result"))
               tags.result.pos <- tags.result


               tags.result.pos <- removeTagsResult(tags.result = tags.result.pos,
                                candidate.num = candidate.num)


               pdf(file.path(output.path, "annotation.information.POS.pdf"),
                   width = 7, height = 7)
               par(mar = c(5, 5, 4, 2))
               annotationPlot(tags.result = tags.result.pos)
               dev.off()

               load(file.path(neg.path, "MRN_annotation_result","intermediate_data","tags.result"))

               tags.result.neg <- tags.result


               tags.result.neg <- removeTagsResult(tags.result = tags.result.neg,
                                                   candidate.num = candidate.num)



               pdf(file.path(output.path, "annotation.information.NEG.pdf"),
                   width = 7, height = 7)
               par(mar = c(5, 5, 4, 2))
               annotationPlot(tags.result = tags.result.neg)
               dev.off()

               png(file.path(output.path, "annotation.information.and.redundancy.removal.in.metABM.POS.png"),
                    width = 14, height = 7, res = 200, units = "in")
               layout(mat = matrix(c(1,2), ncol = 2))
               par(mar = c(5, 5, 4, 2))
               annotationPlot(tags.result = tags.result.pos)
               redundancyPlot(redun = redun1.pos)
               layout(1)
               dev.off()

               png(file.path(output.path, "annotation.information.and.redundancy.removal.in.metABM.NEG.png"),
                    width = 14, height = 7, res = 200, units = "in")
               layout(mat = matrix(c(1,2), ncol = 2))
               par(mar = c(5, 5, 4, 2))
               annotationPlot(tags.result = tags.result.neg)
               redundancyPlot(redun = redun1.neg)
               layout(1)
               dev.off()


               ###recursive annotation result
               ###annotation grade result
               pdf(file.path(output.path, "annotation.grade.barplot.POS.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               gradePlot(tags.result = tags.result.pos, type = "barplot")
               dev.off()


               pdf(file.path(output.path, "annotation.grade.pie.POS.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               gradePlot(tags.result = tags.result.pos, type = "pie")
               dev.off()

               pdf(file.path(output.path, "annotation.grade.barplot.NEG.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               gradePlot(tags.result = tags.result.neg, type = "barplot")
               dev.off()


               pdf(file.path(output.path, "annotation.grade.pie.NEG.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               gradePlot(tags.result = tags.result.neg, type = "pie")
               dev.off()


               png(file.path(output.path, "annotation.grade.POS.png"),
                    width = 14, height = 7, res = 200, units = "in")
               par(mar = c(5,5,4,2))
               layout(mat = matrix(c(1,2), ncol = 2))
               gradePlot(tags.result = tags.result.pos, type = "barplot")
               gradePlot(tags.result = tags.result.pos, type = "pie")
               layout(1)
               dev.off()

               png(file.path(output.path, "annotation.grade.NEG.png"),
                    width = 14, height = 7, res = 200, units = "in")
               par(mar = c(5,5,4,2))
               layout(mat = matrix(c(1,2), ncol = 2))
               gradePlot(tags.result = tags.result.neg, type = "barplot")
               gradePlot(tags.result = tags.result.neg, type = "pie")
               layout(1)
               dev.off()

             }else{
               load(file.path(ifelse(polarity == "positive", pos.path, neg.path), "MRN_annotation_result/intermediate_data/redun1"))
               redun1 <- do.call(rbind, redun1)
               pdf(file = file.path(output.path, "redundancy.removal.in.metABM.pdf"),
                   width = 7, height = 7)
               par(mar = c(5,5,4,2))
               redundancyPlot(redun = redun1)
               dev.off()


               ######################
               load(file.path(ifelse(polarity == "positive", pos.path, neg.path), "MRN_annotation_result/intermediate_data/tags.result"))


               tags.result <- removeTagsResult(tags.result = tags.result,
                                                   candidate.num = candidate.num)

               pdf(file.path(output.path, "annotation.information.pdf"),
                   width = 7, height = 7)
               par(mar = c(5, 5, 4, 2))
               annotationPlot(tags.result = tags.result)
               dev.off()


               png(file.path(output.path, "annotation.information.and.redundancy.removal.in.metABM.png"),
                    width = 14, height = 7, res = 200, units = "in")
               layout(mat = matrix(c(1,2), ncol = 2))
               par(mar = c(5, 5, 4, 2))
               annotationPlot(tags.result = tags.result)
               redundancyPlot(redun = redun1)
               layout(1)
               dev.off()


               ###recursive annotation result
               ###annotation grade result
               pdf(file.path(output.path, "annotation.grade.barplot.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               gradePlot(tags.result = tags.result, type = "barplot")
               dev.off()


               pdf(file.path(output.path, "annotation.grade.pie.pdf"), width = 7, height = 7)
               par(mar = c(5,5,4,2))
               gradePlot(tags.result = tags.result, type = "pie")
               dev.off()


               png(file.path(output.path, "annotation.grade.png"),
                    width = 14, height = 7, res = 200, units = "in")
               par(mar = c(5,5,4,2))
               layout(mat = matrix(c(1,2), ncol = 2))
               gradePlot(tags.result = tags.result, type = "barplot")
               gradePlot(tags.result = tags.result, type = "pie")
               layout(1)
               dev.off()

             }


             ####dysregulated network analysis
             ##p.value and foldchange
             if(polarity == "positive"){
               load(file.path(pos.path,
                              "Pathway_enrichment_analysis_result/intermediate_data",
                              "p.value"))
               load(file.path(pos.path,
                              "Pathway_enrichment_analysis_result/intermediate_data",
                              "foldchange"))
             }

             if(polarity == "negative"){
               load(file.path(neg.path,
                              "Pathway_enrichment_analysis_result/intermediate_data",
                              "p.value"))
               load(file.path(neg.path,
                              "Pathway_enrichment_analysis_result/intermediate_data",
                              "foldchange"))
             }

             if(polarity == "both"){
               load(file.path(path1, "intermediate_data", "p.value.pos"))
               load(file.path(path1, "intermediate_data", "foldchange.pos"))
               load(file.path(path1, "intermediate_data", "p.value.neg"))
               load(file.path(path1, "intermediate_data", "foldchange.neg"))
             }



             if(polarity == "both"){
               p.value <- c(p.value.pos, p.value.neg)
               foldchange <- c(foldchange.pos, foldchange.neg)
             }else{
               p.value <- list(p.value)[[1]]
               foldchange <- list(foldchange)[[1]]
             }


             ####dysregulated network analysis
             if(polarity == 'both'){
               if(any(dir(file.path(path, "intermediate_data")) == "dn.msea")){
                 load(file.path(path, "intermediate_data/dn.msea"))
               }
             }

             if(polarity == 'positive'){
               if(any(dir(file.path(pos.path, "Dysregulated_network_analysis_result/intermediate_data")) == "dn.msea")){
                 load(file.path(pos.path, "Dysregulated_network_analysis_result/intermediate_data/dn.msea"))
               }
             }

             if(polarity == 'negative'){
               if(any(dir(file.path(neg.path, "Dysregulated_network_analysis_result/intermediate_data")) == "dn.msea")){
                 load(file.path(neg.path, "Dysregulated_network_analysis_result/intermediate_data/dn.msea"))
               }
             }

             ####pathway enrichment analysis
             if(polarity == 'both'){
               if(any(dir(file.path(path1, "intermediate_data")) == "msea")){
                 load(file.path(path1, "intermediate_data/msea"))
               }
             }

             if(polarity == 'positive'){
               if(any(dir(file.path(pos.path, "Pathway_enrichment_analysis_result/intermediate_data")) == "msea")){
                 load(file.path(pos.path, "Pathway_enrichment_analysis_result/intermediate_data/msea"))
               }
             }

             if(polarity == 'negative'){
               if(any(dir(file.path(neg.path, "Pathway_enrichment_analysis_result/intermediate_data")) == "msea")){
                 load(file.path(neg.path, "Pathway_enrichment_analysis_result/intermediate_data/msea"))
               }
             }


             if(exists(x = "dn.msea")){
               save(dn.msea, file = file.path(output.path, "dn.msea"), compress = "xz")

               ##volcano plot and pathway.overview
               png(file.path(output.path, "volcanoplot.and.pathway.overview.png"),
                    width = 14, height = 7, res = 200, units = "in")
               par(mar = c(5,5,4,2))
               layout(mat = matrix(c(1,2), ncol = 2))
               volcanoPlot(p.value = p.value, fc = foldchange, correct = correct, p.cutoff = p.cutoff, cex = 0.5)
               pathwayPlot(mse.data = dn.msea@msea)
               layout(1)
               dev.off()
             }else{
               png(file.path(output.path, "volcanoplot.and.pathway.overview.png"),
                    width = 14, height = 7, res = 200, units = "in")
               par(mar = c(5,5,4,2))
               layout(mat = matrix(c(1,2), ncol = 2))
               volcanoPlot(p.value = p.value, fc = foldchange, correct = correct, p.cutoff = p.cutoff, cex = 0.5)
               plot(0, col = "white", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
               legend("top", legend = "No dysregultated network", bty = "n", cex =  1.5)
               layout(1)
               dev.off()
             }


             if(exists(x = "msea")){
               save(msea, file = file.path(output.path, "msea"), compress = "xz")

               ##volcano plot and pathway.overview
               png(file.path(output.path, "volcanoplot.and.pathway.overview.png"),
                    width = 14, height = 7, res = 200, units = "in")
               par(mar = c(5,5,4,2))
               layout(mat = matrix(c(1,2), ncol = 2))
               volcanoPlot(p.value = p.value, fc = foldchange,
                           correct = correct, p.cutoff = p.cutoff, cex = 0.5)
               pathwayPlot(mse.data = msea@msea)
               layout(1)
               dev.off()
             }else{
               png(file.path(output.path, "volcanoplot.and.pathway.overview.png"),
                    width = 14, height = 7, res = 200, units = "in")
               par(mar = c(5,5,4,2))
               layout(mat = matrix(c(1,2), ncol = 2))
               volcanoPlot(p.value = p.value, fc = foldchange, correct = correct, p.cutoff = p.cutoff, cex = 0.5)
               plot(0, col = "white", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
               legend("top", legend = "No dysregultated network", bty = "n", cex =  1.5)
               layout(1)
               dev.off()
             }

             #------------------------------------------------------------------------------
             ###module information

             ##transform rmd to HTML or pdf
             if(type == "html"){
               rmarkdown::render(file.path(output.path, "MetDNA.template.Rmd"), rmarkdown::html_document())
             }

             if(type == "pdf"){
               rmarkdown::render(file.path(output.path, "MetDNA.template.Rmd"), rmarkdown::pdf_document())
             }

             if(type == "all"){
               rmarkdown::render(file.path(output.path, "MetDNA.template.Rmd"), rmarkdown::html_document())
               rmarkdown::render(file.path(output.path, "MetDNA.template.Rmd"), rmarkdown::pdf_document())
             }

             ##remove the png
             file <- dir(output.path)
             file.remove(file.path(output.path, grep("png", file, value = TRUE)))
             file.remove(file.path(output.path, "dn.msea"))
             file.remove(file.path(output.path, "msea"))
             file.remove(file.path(output.path, "parameter"))
             file.remove(file.path(output.path, "parameter.pos"))
             file.remove(file.path(output.path, "parameter.neg"))
             file.remove(file.path(output.path, "MetDNA.logo.png"))
             file.remove(file.path(output.path, "MetDNA.template.Rmd"))

             ###rename
             file.rename(from = file.path(output.path, "MetDNA.template.html"),
                         to = file.path(output.path, "MetDNA.analysis.report.html"))

             file.rename(from = file.path(output.path, "annotation.grade.barplot.pdf"),
                         to = file.path(output.path, "Identification.grade.barplot.pdf"))

             file.rename(from = file.path(output.path, "annotation.grade.pie.pdf"),
                         to = file.path(output.path, "Identification.grade.pie.pdf"))

             file.rename(from = file.path(output.path, "annotation.information.pdf"),
                         to = file.path(output.path, "Identification.information.pdf"))

             file.rename(from = file.path(output.path, "peak.distribution.plot.pos.pdf"),
                         to = file.path(output.path, "Peak.intensity.profile.pos.pdf"))
             file.rename(from = file.path(output.path, "peak.distribution.plot.neg.pdf"),
                         to = file.path(output.path, "Peak.intensity.profile.neg.pdf"))

             file.rename(from = file.path(output.path, "redundancy.removal.in.metABM.pdf"),
                         to = file.path(output.path, "Peak.and.metabolite.redundancy.pdf"))

             file.rename(from = file.path(output.path, "annotation.grade.barplot.NEG.pdf"),
                         to = file.path(output.path, "Identification.grade.barplot.neg.pdf"))

             file.rename(from = file.path(output.path, "annotation.grade.barplot.POS.pdf"),
                         to = file.path(output.path, "Identification.grade.barplot.pos.pdf"))

             file.rename(from = file.path(output.path, "annotation.grade.pie.NEG.pdf"),
                         to = file.path(output.path, "Identification.grade.pie.neg.pdf"))

             file.rename(from = file.path(output.path, "annotation.grade.pie.POS.pdf"),
                         to = file.path(output.path, "Identification.grade.pie.pos.pdf"))

             file.rename(from = file.path(output.path, "annotation.information.NEG.pdf"),
                         to = file.path(output.path, "Identification.information.neg.pdf"))

             file.rename(from = file.path(output.path, "annotation.information.POS.pdf"),
                         to = file.path(output.path, "Identification.information.pos.pdf"))

             file.rename(from = file.path(output.path, "redundancy.removal.in.metABM.NEG.pdf"),
                         to = file.path(output.path, "Peak.and.metabolite.redundancy.neg.pdf"))

             file.rename(from = file.path(output.path, "redundancy.removal.in.metABM.POS.pdf"),
                         to = file.path(output.path, "Peak.and.metabolite.redundancy.pos.pdf"))

             cat("Analysis report is done.\n")
           })


# ##module.score
setGeneric(name = "redundancyPlot",
           def = function(redun){
             plot(redun[,1], col = "#729ECE", type = "o",
                  ylim = c(range(redun)[1] * 0.8, range(redun)[2] * 1.2),
                  pch = 19, lwd = 2, cex = 1.5,
                  xlab = "Recursive times",
                  ylab = 'Redundancy', cex.lab = 1.8,
                  cex.axis = 1.5, xaxt = "n")

             par(xpd = FALSE)
             axis(side = 1, at = c(1:nrow(redun)),
                  cex.axis = 1.5, labels = c("Raw", 1:(nrow(redun) - 1)))
             points(redun[,2], col = "#FF9E4A", type = "o", pch = 17, lwd = 2, cex = 1.5)


             maptools::pointLabel(x = rep(c(1:nrow(redun)), 2),
                                  y = c(redun[,1],redun[,2]),
                                  labels = as.character(c(round(redun[,1], 2),
                                                          round(redun[,2], 2))),
                                  cex = 1.3)

             abline(v = 3, lty = 2, lwd = 1.5, col = "tomato")
             abline(v = 4, lty = 2, lwd = 1.5, col = "tomato")
             legend("topright",
                    legend = c("Peak redundancy", "Metabolite redundancy"),
                    pch = c(19, 17), col = c("#729ECE", "#FF9E4A"),
                    bty = "n", pt.cex = 1.3, cex = 1.3, lty = 2, lwd = 1.5)
           })






setGeneric(name = "annotationPlot",
           def = function(tags.result){
             level <- as.numeric(tags.result$level)

             seed.number <- unlist(lapply(sort(unique(level)), function(x) {
               temp <- tags.result[level == x,]
               temp <- temp[temp$as.seed == "TRUE",]
               length(unique(temp$to))
               # length(unique(tags.result$to[level == x]))
             }))


             temp <- lapply(sort(unique(level)), function(x){
               temp <- tags.result[level == x,]
               unique(temp$to)
             })

             cummulative.metabolite.number <- unlist(lapply(1:length(temp), function(x){
               length(unique(unlist(temp[1:x])))
             }))


             temp <- barplot(seed.number,
                             ylab = "#Metabolite",
                             cex.axis = 1.5,cex.lab = 1.8,
                             ylim = c(0,max(cummulative.metabolite.number)),
                             border = NA,
                             col = "#729ECE",, names.arg = c(1:length(seed.number)) - 1,
                             xlab = "Round", cex.names = 1.5)
             par(xpd = FALSE)
             abline(h=0)
             par(xpd = TRUE)
             points(temp[,1], y = cummulative.metabolite.number, pch = 17,
                    col = "#ED665D", cex = 1, type = "o", lwd = 2)

             maptools::pointLabel(x = temp[,1],
                                  y = cummulative.metabolite.number,
                                  labels = as.character(cummulative.metabolite.number),
                                  cex = 0.8)
             par(xpd = FALSE)

             legend("right", legend = c("Metabolite number", "Cummulative metabolite number"),
                    col = c("#729ECE", "#ED665D"), bty = "n", pch = c(15, 17), pt.cex = 1.3,
                    cex = 1.3, lty = 2, lwd = 2)
           })






setGeneric(name = "gradePlot",
           def = function(tags.result, type = c("barplot", "pie")){
             type <- match.arg(type)
             unique.id <- unique(tags.result$to)
             grade <- unname(sapply(unique.id, function(temp.id) {
               temp.idx <- grep(temp.id, tags.result$group)
               sort(tags.result$Confidence[temp.idx])[1]
             }))
             grade <- table(grade)

             if(type == "barplot"){
               temp <- barplot(grade, border = NA, cex.lab = 1.8,
                               cex.axis = 1.5, col = c("#729ECE", "#FF9E4A", "#AD8BC9", "#ED665D"),
                               ylab = 'Metabolite ID number', cex.names = 1.5,
                               main = "Identification grade", cex.main = 1.5)
               par(xpd = FALSE)
               abline(h = 0)
               text(x = temp, y = grade/2, labels = grade, cex = 1.3,
                    col = "white")
             }else{
               pie(x = grade, radius = 1, clockwise = TRUE,
                   col = c("#729ECE", "#FF9E4A", "#AD8BC9", "#ED665D"), border = NA,
                   cex = 1.5, main = 'Identification grade', cex.main = 1.5,
                   labels = "")
               par(new = TRUE)
               pie(1, col = "white", labels = "", radius = 0.6, border = NA)
               legend("center",
                      legend = c(paste("Grade 1:", round(grade[1]*100/sum(grade),1), "%"),
                                 paste("Grade 2:", round(grade[2]*100/sum(grade), 1), "%"),
                                 paste("Grade 3:", round(grade[3]*100/sum(grade),1), "%"),
                                 paste("Grade 4", round(grade[4]*100/sum(grade), 1), "%")),
                      bty = "n", pch = 19, pt.cex = 1.8, col = c("#729ECE", "#FF9E4A", "#AD8BC9", "#ED665D"),
                      cex = 1.3)
             }





           })



#-----------------------------------------------------------------------------
##data information
setGeneric(name = "dataInformation",
           def = function(pos.path,
                          neg.path,
                          sample.info = "sample.info.csv",
                          data.pos = "data.csv",
                          data.neg = "data.csv",
                          output.path,
                          polarity = c("positive", "negative", "both")){

             old.data.pos.name <- data.pos
             old.data.neg.name <- data.neg
             old.sample.info.name <- sample.info

             polarity <- match.arg(polarity)

             my.theme <- ggplot2::theme_bw()+
               ggplot2::theme(axis.title.x = ggplot2::element_text(size = 18),
                              axis.title.y = ggplot2::element_text(size = 18)) +
               ggplot2::theme(axis.text.x = ggplot2::element_text(size = 15),
                              axis.text.y = ggplot2::element_text(size = 15)) +
               ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) +
               ggplot2::theme(legend.text = ggplot2::element_text(size = 10))

             if(polarity == "positive" | polarity == "both"){
               data.pos <- readr::read_csv(file.path(pos.path, old.data.pos.name),
                                           col_types = readr::cols(), progress = FALSE)
               data.pos <- as.data.frame(data.pos)
               sample.info <- readr::read_csv(file.path(pos.path, old.sample.info.name),
                                              col_types = readr::cols())
               sample.info <- as.data.frame(sample.info)
               tags.idx <- match(c("name", "mz", "rt"), colnames(data.pos))
               tags.pos <- data.pos[,tags.idx]
               sample.pos <- data.pos[,-tags.idx]
               mz.pos <- tags.pos$mz
               rt.pos <- tags.pos$rt
               int.log.pos <- log(apply(sample.pos, 1, function(x){
                 median(x, na.rm = TRUE)}
               ) + 10, 10)
               rt.mz.int.pos <- data.frame(rt.pos, mz.pos, int.log.pos,
                                           stringsAsFactors = FALSE)

               rt.mz.int.pos <-
                 ggplot2::ggplot(data = rt.mz.int.pos,
                                 ggplot2::aes(x = rt.pos, y = mz.pos, colour = int.log.pos)) +
                 ggplot2::geom_point(alpha = 0.3) +
                 ggplot2::scale_color_gradient(low = "green", high = "red") +
                 ggplot2::labs(x = "Retention time (RT, second)",
                               y = "Mass to charge ratio (m/z)",
                               colour = "log10(intensity)") +
                 my.theme+ggplot2::ggtitle("Positive mode")+
                 ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20))

               #pdf
               ggplot2::ggsave(
                 filename = file.path(output.path, "peak.distribution.plot.pos.pdf"),
                 plot = rt.mz.int.pos,
                 width = 8,
                 height = 6
               )

               ##png
               ggplot2::ggsave(
                 filename = file.path(output.path, "peak.distribution.plot.png"),
                 plot = rt.mz.int.pos,
                 width = 8,
                 height = 6, dpi = 600)
             }

             ##negative
             if(polarity == "negative" | polarity == "both"){
               data.neg <- readr::read_csv(file.path(neg.path, old.data.neg.name),
                                           col_types = readr::cols(), progress = FALSE)
               data.neg <- as.data.frame(data.neg)
               sample.info <- readr::read_csv(file.path(neg.path, old.sample.info.name),
                                              col_types = readr::cols())
               sample.info <- as.data.frame(sample.info)
               tags.idx <- match(c("name", "mz", "rt"), colnames(data.neg))
               tags.neg <- data.neg[,tags.idx]
               sample.neg <- data.neg[,-tags.idx]
               mz.neg <- tags.neg$mz
               rt.neg <- tags.neg$rt
               int.log.neg <- log(apply(sample.neg, 1, function(x){
                 median(x, na.rm = TRUE)}
               ) + 10, 10)
               rt.mz.int.neg <- data.frame(rt.neg, mz.neg, int.log.neg,
                                           stringsAsFactors = FALSE)

               rt.mz.int.neg <-
                 ggplot2::ggplot(data = rt.mz.int.neg,
                                 ggplot2::aes(x = rt.neg, y = mz.neg, colour = int.log.neg)) +
                 ggplot2::geom_point(alpha = 0.3) +
                 ggplot2::scale_color_gradient(low = "green", high = "red") +
                 ggplot2::labs(x = "Retention time (RT, second)",
                               y = "Mass to charge ratio (m/z)",
                               colour = "log10(intensity)") +
                 my.theme+ggplot2::ggtitle("Negative mode")+
                 ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20))

               ##pdf
               ggplot2::ggsave(
                 filename = file.path(output.path, "peak.distribution.plot.neg.pdf"),
                 plot = rt.mz.int.neg,
                 width = 8,
                 height = 6
               )
               ##png
               # ggplot2::ggsave(
               #   filename = file.path(output.path, "peak.distribution.plot.png"),
               #   plot = rt.mz.int.neg,
               #   width = 8,
               #   height = 6,dpi = 600)


               png(file.path(output.path, "peak.distribution.plot.png"),
                   width = 8, height = 6, res = 200, units = "in")
               gridExtra::grid.arrange(rt.mz.int.neg, ncol = 1)
               dev.off()

             }


             if(polarity == "both"){
               png(file.path(output.path, "peak.distribution.plot.png"),
                    width = 16, height = 7, res = 200, units = "in")
               gridExtra::grid.arrange(rt.mz.int.pos, rt.mz.int.neg, ncol = 2)
               dev.off()
             }


           })




####parameter transformation
setGeneric(name = "parameterMap", def = function(name){
  species <- data.frame("r" = c("hsa", "mmu",
                                "rat", "bta",
                                "dme",
                                "gga",  "dre",
                                "cel",
                                "sce",
                                "ath",
                                "smm",
                                "pfa",
                                "tbr",
                                "eco",
                                "ppu", "syf"),
                        "web" = c("Homo sapiens (human)", "Mus musculus (mouse)",
                                  "Rattus norvegicus (rat)", "Bos taurus (cow)",
                                  "Drosophila melanogaster (fruit fly)",
                                  "Gallus gallus (chicken)", "Danio rerio (zebrafish)",
                                  "Caenorhabditis elegans (nematode)",
                                  "Saccharomyces cerevisiae (yeast)",
                                  "Arabidopsis thaliana (thale cress)",
                                  "Schistosoma mansoni",
                                  "Plasmodum falciparum 3D7 (Malaria)",
                                  "Trypanosoma brucei",
                                  "Escherichia coli K-12 MG1655",
                                  "Pseudomonas putida KT2440", "Synechococcus elongatus"),
                        stringsAsFactors = FALSE)


  STAT_METHOD = data.frame("r" = c('t', 'wilcox'),
                           "web" = c('Student t-test', 'Wilcox test'),
                           stringsAsFactors = FALSE)

  ADJUST_P = data.frame("r" = c('FALSE', 'TRUE'),
                        "web" = c('No', 'Yes'),
                        stringsAsFactors = FALSE)



  INSTT = data.frame("r" = c('SciexTripleTOF',
                             'AgilentQTOF',
                             'OtherQTOF',
                             'ThermoOrbitrap'),
                     "web" = c(  'Sciex TripleTOF',
                                 'Agilent QTOF',
                                 'Other QTOF',
                                 'Thermo Orbitrap (HCD)'),
                     stringsAsFactors = FALSE)


  OTHER = data.frame("r" = c('positive',
                             'negative',
                             'both',
                             "hilic",
                             "rp"),
                     "web" = c(  'Positive',
                                 'Negative',
                                 'Both',
                                 "HILIC",
                                 "RP"),
                     stringsAsFactors = FALSE)



  paramter.map <- rbind(species, STAT_METHOD, ADJUST_P, INSTT, OTHER)

  paramter.map$web[match(name, paramter.map$r)]

})

