# load("remain.idx")
# load("raw.tags2")
# load('new.tags')
# load("result2")
# write.csv(result2, "result2.csv")
# peak.name <- showTags2(raw.tags2[remain.idx], slot = "name")
# idx <- match(peak.name, result2$name)
# idx <- idx[!is.na(idx)]
# idx <- sort(idx)
# result2_2 <- result2[idx,]
# write.csv(result2_2, "result2_2.csv")
# remain.idx
#
# temp.name <- result2_2$name[temp.idx]
#
# temp.idx <- match(temp.name, showTags2(raw.tags2, "name"))
#
# id <- showTags2(raw.tags2[temp.idx], slot = "annotation.id")
# id <- lapply(id, function(x) {paste(x, collapse = ";")})
# id <- unlist(id)
#
#
#
# kegg.compound[c(1565, 25),]
#
# match(c("M130T687", "M130T733"),unlist(lapply(ms2, function(x) {x[[1]]["NAME",1]})))
# ms2_1 <- ms2[[174]][[2]]
# ms2_2 <- ms2[[172]][[2]]
#
# par(mar = c(5,5,4,2))
# plot(0,0, xlim = c(0, 150), ylim = c(-1,1), col = "white",
#      xlab = "Mass to charge ratio (m/z)", ylab = "Relative intensity",
#      cex.lab = 1.8, cex.axis = 1.5)
# abline(h = 0, lty = 2, lwd = 1.5)
# points(ms2_1[,1], ms2_1[,2]/max(ms2_1[,2]), type = "h", col = "lightseagreen", lwd = 1.5)
# points(ms2_2[,1], -ms2_2[,2]/max(ms2_2[,2]), type = "h", col = "salmon", lwd = 1.5)
#
# load("inHouse.compound")
#
# match(c('C01879',"C00025"), inHouse.compound$KEGG.ID)
# inHouse.compound[c(477, 1556),]
#
# match(c("S0010", "n00354"), rt.all$`in house ID`)
# rt.all[346,]
#
# a <- rt.result[[1]]
# match(c("S0010", "n00354"), rownames(a))
# a[c(474,1538),]
#
#
# match(c("M130T687", "M130T733"), showTags2(raw.tags2, "name"))
# raw.tags2[[553]]
# raw.tags2[[551]]
#
#
#
#
#
#
#
# load("kegg.rpair2")
#
#
#
# setwd("/home/jasper/work/metABM/figure")
# dir()
# load("kegg.rpair")
# load("kegg.rpair2.rda")
# load("inHouse.compound.rda")
# ms2.data <- ms2Annotation::zhuMetlib
#
#
#
#
# lab.id <- inHouse.compound$Lab.ID
# kegg.id <- inHouse.compound$KEGG.ID
#
# remove.idx <- which(is.na(kegg.id) | nchar(kegg.id) != 6)
#
# lab.id <- lab.id[-remove.idx]
# kegg.id <- kegg.id[-remove.idx]
#
#
# result <- list()
#
# for(i in seq_along(lab.id)){
#   cat(i); cat(" ")
# if(i == length(lab.id)) {result[[i]] <- NULL; next}
#   temp.id <- lab.id[i]
#   temp.id1 <- kegg.id[i]
#   neighbor <- getNeighbor(metabolite.id = temp.id1, step = 1, graph = kegg.rpair2)
#   if(is.null(neighbor)) {result[[i]] <- NULL; next}
#   id <- neighbor$Node.ID[which(neighbor$Node.ID %in% kegg.id[(i+1):length(kegg.id)])]
#   if(length(id) == 0) {result[[i]] <- NULL; next}
#   id1 <- lab.id[match(id, kegg.id)]
#   result[[i]] <- data.frame(temp.id, id1, stringsAsFactors = FALSE)
# }
#
#
# result <- do.call(rbind, result)
#
#
#
# dp2 <- list()
# unique.id <- unique(c(result[,1], result[,2]))
#
# for(i in 1:length(unique.id)){
#   cat(i); cat(" ")
#   temp.id <- unique.id[i]
#   idx1 <- which(result[,1] == temp.id)
#   idx2 <- which(result[,2] == temp.id)
#   idx <- c(idx1, idx2)
#
#   temp <- NULL
#   for(j in idx){
#     cat(j); cat(" ")
#     temp.id1 <- result[j,1]
#     temp.id2 <- result[j,2]
#     mass1 <- inHouse.compound$EM[match(temp.id1, inHouse.compound$Lab.ID)]
#     mass2 <- inHouse.compound$EM[match(temp.id2, inHouse.compound$Lab.ID)]
#
#     if(is.na(mass1)) {mass1 = 1}
#     if(is.na(mass2)) {mass2 = 1}
#     ms2_1 <- getMS2(libID = result[j,1], lib = ms2.data)
#     ms2_2 <- getMS2(libID = result[j,2], lib = ms2.data)
#
#     if(mass1 >= mass2) {
#      temp[match(j, idx)] <- as.numeric(IdentifyFeature(ms2_1, ms2_2,
#                                           direction = "reverse"))
#     }else{
#       temp[match(j, idx)] <- as.numeric(IdentifyFeature(ms2_2, ms2_1,
#                                           direction = "reverse"))
#     }
#
#   }
# dp2[[i]] <- temp
# }
#
# dp2 <- lapply(dp2, function(x) {
#   x[is.na(x)] <- 0
#   x[is.nan(x)] <- 0
#   x
# })
#
#
# temp1 <- dp2
# temp1 <- unlist(lapply(temp1, function(x) {max(x)}))
#
#
# remove.idx <- which(temp1 == 0)
#
# remove.id <- unique.id[remove.idx]
#
# neighbor <- result[-which(result[,1] %in% remove.id | result[,2] %in% remove.id),]
#
# result <- neighbor
#
#
# dp2 <- list()
# unique.id <- unique(c(result[,1], result[,2]))
#
# for(i in 1:length(unique.id)){
#   cat(i); cat(" ")
#   temp.id <- unique.id[i]
#   idx1 <- which(result[,1] == temp.id)
#   idx2 <- which(result[,2] == temp.id)
#   idx <- c(idx1, idx2)
#
#   temp <- NULL
#   for(j in idx){
#     cat(j); cat(" ")
#     temp.id1 <- result[j,1]
#     temp.id2 <- result[j,2]
#     mass1 <- inHouse.compound$EM[match(temp.id1, inHouse.compound$Lab.ID)]
#     mass2 <- inHouse.compound$EM[match(temp.id2, inHouse.compound$Lab.ID)]
#
#     if(is.na(mass1)) {mass1 = 1}
#     if(is.na(mass2)) {mass2 = 1}
#     ms2_1 <- getMS2(libID = result[j,1], lib = ms2.data)
#     ms2_2 <- getMS2(libID = result[j,2], lib = ms2.data)
#
#     if(mass1 >= mass2) {
#       temp[match(j, idx)] <- as.numeric(IdentifyFeature(ms2_1, ms2_2,
#                                                         direction = "reverse"))
#     }else{
#       temp[match(j, idx)] <- as.numeric(IdentifyFeature(ms2_2, ms2_1,
#                                                         direction = "reverse"))
#     }
#
#   }
#   dp2[[i]] <- temp
# }
#
# dp2 <- lapply(dp2, function(x) {
#   x[is.na(x)] <- 0
#   x[is.nan(x)] <- 0
#   x
# })
#
#
#
# new.result <- list()
# for(i in 1:length(unique.id)){
#   cat(i); cat(" ")
#   temp.id <- unique.id[i]
#   idx1 <- which(result[,1] == temp.id)
#   idx2 <- which(result[,2] == temp.id)
#   if(length(idx1) != 0) {
#     neighbor1 <- result[idx1,2]
#   }else{
#       neighbor1 <- NULL
#     }
#   if(length(idx2) != 0) {
#     neighbor2 <- result[idx2,1]
#   }else{
#       neighbor2 <- NULL
#   }
#
#   neighbor <- c(neighbor1, neighbor2)
#
#   nei.pool <- setdiff(unique.id, neighbor)
#   new.neighbor <- nei.pool[sample(1:length(nei.pool), length(neighbor))]
#   new.result[[i]] <- data.frame(temp.id, new.neighbor, stringsAsFactors = FALSE)
# }
#
#
# new.result <- do.call(rbind, new.result)
#
#
#
#
#
# dp2_1 <- list()
#
# for(i in 1:length(unique.id)){
#   cat(i); cat(" ")
#   temp.id <- unique.id[i]
#   idx1 <- which(result[,1] == temp.id)
#   idx2 <- which(result[,2] == temp.id)
#   idx <- c(idx1, idx2)
#
#   temp <- NULL
#   for(j in idx){
#     cat(j); cat(" ")
#     temp.id1 <- new.result[j,1]
#     temp.id2 <- new.result[j,2]
#     mass1 <- inHouse.compound$EM[match(temp.id1, inHouse.compound$Lab.ID)]
#     mass2 <- inHouse.compound$EM[match(temp.id2, inHouse.compound$Lab.ID)]
#
#     if(is.na(mass1)) {mass1 = 1}
#     if(is.na(mass2)) {mass2 = 1}
#     ms2_1 <- getMS2(libID = new.result[j,1], lib = ms2.data)
#     ms2_2 <- getMS2(libID = new.result[j,2], lib = ms2.data)
#
#     if(mass1 >= mass2) {
#       temp[match(j, idx)] <- as.numeric(IdentifyFeature(ms2_1, ms2_2,
#                                                         direction = "reverse"))
#     }else{
#       temp[match(j, idx)] <- as.numeric(IdentifyFeature(ms2_2, ms2_1,
#                                                         direction = "reverse"))
#     }
#
#   }
#   dp2_1[[i]] <- temp
# }
#
# dp2_1 <- lapply(dp2_1, function(x) {
#   x[is.na(x)] <- 0
#   x[is.nan(x)] <- 0
#   x
# })
#
#
#
#
#
# x <- unlist(dp2)
# y <- unlist(dp2_1)
# remove.idx <- which(x == 0)
# x <- x[-remove.idx]
# y <- y[-remove.idx]
# z <- c(x, y)
# class <- c(rep(1, length(x)), rep(0, length(y)))
# library(pROC)
# roc1 <- roc(class, z, ci = TRUE)
# plot(roc1)
#
# threshold <- roc1$thresholds
# threshold[1] <- 0
# threshold[792] <- 1
#
#
# par(mar = c(5,5,4,2))
# plot(roc1, xlab = "Specificity (1-False positive rate)", ylab = "Sensitivity (True positive rate)",
#      lwd = 3, cex.lab = 1.8, cex.axis = 1.5,
#      col = "lightseagreen")
# legend("bottom", legend = c('AUC: 0.935', "95% CI: 0.926-0.943"),
#        bty = "n", cex = 1.5)
#
# library(ROCR)
# pred <- prediction(predictions = z, labels = class)
# rp.perf <- performance(pred, "prec", "rec")
#
#
# roc.perf <- performance(pred, "tpr", "fpr")
#
#
# f1 <- performance(pred, "f")
# plot(threshold, roc1$specificities, type = "l",
#      col = "lightseagreen",lwd = 2, xlab = "Cutoff of dot product",
#      cex.lab = 1.8, cex.axis = 1.5, ylab = "Specificity/Sensiticity",
#      ylim = c(0,1))
#
#
# points(threshold, roc1$sensitivities, type = "l",
#        lwd = 2, col = "salmon")
#
# # points(f1@x.values[[1]], f1@y.values[[1]], type = "l",
# #        lwd = 2, col = "orchid4")
#
# abline(v = 0.5, lty = 3, lwd = 1.5)
#
# points(x = 0.5, roc1$specificities[455], type = "p", col = "lightseagreen",
#        cex = 1.5, pch = 19)
# points(x = 0.5, roc1$sensitivities[455], type = "p", col = "salmon",
#        cex = 1.5, pch = 19)
#
# text(x = 0.5, roc1$specificities[455], labels = "(0.5, 0.950)",
#      cex = 1.5, pos = 1)
#
# text(x = 0.5, roc1$sensitivities[455], labels = "(0.5, 0.518)",
#      cex = 1.5, pos = 1)
#
# legend("bottom", legend = c("Specificity/(1-False positive rate)",
#                                     "Sensitivity(True positive rate)"),
#        bty = "n", pch = 19, lty = 1, col = c("lightseagreen", "salmon"),
#        pt.cex = 1.5, cex = 1.5, lwd = 2)
# # points(x = 0.5, f1@y.values[[1]][331], type = "p", col = "orchid4",
# #        cex = 1.5, pch = 19)
#
#
# x1 <- dp2
# y1 <- dp2_1
#
#
# x1 <- unlist(lapply(x1, function(x) max(x)))
# y1 <- unlist(lapply(y1, function(x) max(x)))
#
#
# boxplot(x1, y1)
#
# library(beanplot)
# beanplot(x1, y1, ylab = "Dot product",
#          cex.lab = 1.8, cex.axis = 1.5, names = c("Reaction pair", "Non-RP"),
#           col = "salmon", border = "salmon")
# abline(h = 0.5, lty = 3, lwd = 1.5, col = "lightseagreen")
#
#
# x2 <- unlist(dp2)
# y2 <- unlist(dp2_1)
#
#
# boxplot(x2, y2)
#
# library(beanplot)
# beanplot(x2, y2, ylab = "Dot product",
#          cex.lab = 1.8, cex.axis = 1.5, names = c("Reaction pair", "Non-RP"),
#          col = "salmon", border = "salmon")
# abline(h = 0.5, lty = 3, lwd = 1.5, col = "lightseagreen")
#
#
# dim(new.result)
#
#
# neighbor <- result
#
#
# sum(x2>=0.5)
# length(x2) - sum(x2>=0.5)
#
# sum(y2>=0.5)
# length(y2) - sum(y2>=0.5)
# temp <- cbind(c(910, 149), c(86,409))
# temp1 <- barplot(temp, names.arg = c("Reaction pair", "Non-RP"),
#         ylab = "Number", cex.lab = 1.8, cex.axis = 1.5,
#         col = c('lightseagreen', "salmon"), border = NA,
#         cex.names = 1.5)
#
#
# text(x = 0.7, y = 346/2, labels = "DP>=0.5:\n70.0%", cex = 1.5, col = "white")
# text(x = 0.7, y = 346 + 149/2, labels = "DP<0.5:\n30.0%", cex = 1.5, col = "white")
#
# text(x = 1.9, y = 86/2, labels = "DP>=0.5:\n17.4%", cex = 1.5, col = "white")
# text(x = 1.9, y = 86 + 409/2, labels = "DP<0.5:\n82.6%", cex = 1.5, col = "white")
#
#
#
#
#
#
# ###W03 vs W30
# annotation.result <- readr::read_csv("annotation.result.csv")
# sample.info <- read.csv("sample.info.csv",stringsAsFactors = FALSE)
#
# sample.info <- sample.info[which(sample.info[,2] %in% c("W03", "W30")),]
# sample <- annotation.result[,match(sample.info[,1], colnames(annotation.result))]
#
# p.value <- uniTest(sample = sample, sample.info = sample.info, uni.test = "t", correct = FALSE)
# fc <- foldChange(sample = sample, sample.info = sample.info, by.what = "median", group = c("W03", 'W30'))
# fc[which(is.infinite(fc))] <- max(fc[!is.infinite(fc)])
#
# my.col <- rep("grey", length(p.value))
# my.col[p.value < 0.05 & fc > 1] <- 'salmon'
# my.col[p.value < 0.05 & fc < 1] <- 'lightseagreen'
#
# par(mar = c(5,5,4,2))
# plot(log(fc,2), -log(p.value, 10), col = my.col, pch = 19,
#      cex = 0.3, xlab = "log2(Fold change)",
#      ylab = "-log10(p value)",
#      cex.lab = 1.8,
#      cex.axis = 1.5
# )
# abline(v = 0, lty = 3, lwd = 1.5)
# abline(h = 1.3, lty = 3, lwd = 1.5)
#
#
# ###W03 vs W15
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/W03 vs W15")
# annotation.result <- readr::read_csv("annotation.result.csv")
# sample.info <- read.csv("sample.info.csv",stringsAsFactors = FALSE)
#
# sample.info <- sample.info[which(sample.info[,2] %in% c("W03", "W15")),]
# sample <- annotation.result[,match(sample.info[,1], colnames(annotation.result))]
#
# p.value <- uniTest(sample = sample, sample.info = sample.info, uni.test = "t", correct = FALSE)
# fc <- foldChange(sample = sample, sample.info = sample.info, by.what = "median", group = c("W03", 'W15'))
# fc[which(is.infinite(fc))] <- max(fc[!is.infinite(fc)])
#
# my.col <- rep("grey", length(p.value))
# my.col[p.value < 0.05 & fc > 1] <- 'salmon'
# my.col[p.value < 0.05 & fc < 1] <- 'lightseagreen'
#
# par(mar = c(5,5,4,2))
# plot(log(fc,2), -log(p.value, 10), col = my.col, pch = 19,
#      cex = 0.3, xlab = "log2(Fold change)",
#      ylab = "-log10(p value)",
#      cex.lab = 1.8,
#      cex.axis = 1.5
# )
# abline(v = 0, lty = 3, lwd = 1.5)
# abline(h = 1.3, lty = 3, lwd = 1.5)
#
# marker1 <- which(p.value < 0.05)
#
#
# ##W03 vs W30
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/W03 vs W30")
# annotation.result <- readr::read_csv("annotation.result.csv")
# sample.info <- read.csv("sample.info.csv",stringsAsFactors = FALSE)
#
# sample.info <- sample.info[which(sample.info[,2] %in% c("W03", "W30")),]
# sample <- annotation.result[,match(sample.info[,1], colnames(annotation.result))]
#
# p.value <- uniTest(sample = sample, sample.info = sample.info, uni.test = "t", correct = FALSE)
# fc <- foldChange(sample = sample, sample.info = sample.info, by.what = "median", group = c("W03", 'W30'))
# fc[which(is.infinite(fc))] <- max(fc[!is.infinite(fc)])
#
# my.col <- rep("grey", length(p.value))
# my.col[p.value < 0.05 & fc > 1] <- 'salmon'
# my.col[p.value < 0.05 & fc < 1] <- 'lightseagreen'
#
# par(mar = c(5,5,4,2))
# plot(log(fc,2), -log(p.value, 10), col = my.col, pch = 19,
#      cex = 0.3, xlab = "log2(Fold change)",
#      ylab = "-log10(p value)",
#      cex.lab = 1.8,
#      cex.axis = 1.5
# )
# abline(v = 0, lty = 3, lwd = 1.5)
# abline(h = 1.3, lty = 3, lwd = 1.5)
#
# marker2 <- which(p.value < 0.05)
#
# ###P03 vs P30
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/P03 vs P30")
# annotation.result <- readr::read_csv("annotation.result.csv")
# sample.info <- read.csv("sample.info.csv",stringsAsFactors = FALSE)
#
# sample.info <- sample.info[which(sample.info[,2] %in% c("P03", "P30")),]
# sample <- annotation.result[,match(sample.info[,1], colnames(annotation.result))]
#
# p.value <- uniTest(sample = sample, sample.info = sample.info, uni.test = "t", correct = FALSE)
# fc <- foldChange(sample = sample, sample.info = sample.info, by.what = "median", group = c("P03", 'P30'))
# fc[which(is.infinite(fc))] <- max(fc[!is.infinite(fc)])
#
# my.col <- rep("grey", length(p.value))
# my.col[p.value < 0.05 & fc > 1] <- 'salmon'
# my.col[p.value < 0.05 & fc < 1] <- 'lightseagreen'
#
# par(mar = c(5,5,4,2))
# plot(log(fc,2), -log(p.value, 10), col = my.col, pch = 19,
#      cex = 0.3, xlab = "log2(Fold change)",
#      ylab = "-log10(p value)",
#      cex.lab = 1.8,
#      cex.axis = 1.5
# )
# abline(v = 0, lty = 3, lwd = 1.5)
# abline(h = 1.3, lty = 3, lwd = 1.5)
#
# marker3 <- which(p.value < 0.05)
#
# ###W30 vs P30
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/W30 vs P30")
# annotation.result <- readr::read_csv("annotation.result.csv")
# sample.info <- read.csv("sample.info.csv",stringsAsFactors = FALSE)
#
# sample.info <- sample.info[which(sample.info[,2] %in% c("W30", "P30")),]
# sample <- annotation.result[,match(sample.info[,1], colnames(annotation.result))]
#
# p.value <- uniTest(sample = sample, sample.info = sample.info, uni.test = "t", correct = FALSE)
# fc <- foldChange(sample = sample, sample.info = sample.info, by.what = "median", group = c("W30", 'P30'))
# fc[which(is.infinite(fc))] <- max(fc[!is.infinite(fc)])
#
# my.col <- rep("grey", length(p.value))
# my.col[p.value < 0.05 & fc > 1] <- 'salmon'
# my.col[p.value < 0.05 & fc < 1] <- 'lightseagreen'
#
# par(mar = c(5,5,4,2))
# plot(log(fc,2), -log(p.value, 10), col = my.col, pch = 19,
#      cex = 0.3, xlab = "log2(Fold change)",
#      ylab = "-log10(p value)",
#      cex.lab = 1.8,
#      cex.axis = 1.5
# )
# abline(v = 0, lty = 3, lwd = 1.5)
# abline(h = 1.3, lty = 3, lwd = 1.5)
# marker4 <- which(p.value < 0.05)
#
#
#
#
#
# library(VennDiagram)
# a <- venn.diagram(x = list(a = marker2, b = marker3), filename = NULL,
#                   col = c("lightseagreen", "salmon"), lwd = 8,cex = 1.5,
#                   category.names = c("W30 vs W03", "P03 vs P30"), cat.cex = 1.5
#                   )
# grid.draw(a)
#
#
# a <- venn.diagram(x = list(a = marker1, b = marker2), filename = NULL,
#                   col = c("lightseagreen", "salmon"), lwd = 8,cex = 1.5,
#                   category.names = c("W03 vs W15", "W03 vs W30"), cat.cex = 1.5
# )
# grid.draw(a)
#
#
#
