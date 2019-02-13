# setwd("/home/jasper/work/MetDNA/fly neg")
#
# system.time(metABM(prefer.adduct = "M-H", column = "hilic",
#                    polarity = "negative", threads = 3))
# #
# #
# system.time(metModule(annotation.result = "annotation.result.csv",
#                       path = "/home/jasper/work/MetDNA/fly neg/metABM result",
#                       output.path = "/home/jasper/work/MetDNA/fly neg/metModule result2",
#                       polarity = "negative",
#                       group = c('W03', "W30"), uni.test = "t", column = "hilic",
#                       # path = "/home/jasper/work/MetDNA/fly pos/metABM",
#                       # output.path = "/home/jasper/work/MetDNA/fly pos/metModule",
#                       threads = 3, use.old.null = TRUE, species = "dme"))
#
#
#
# metModule2(annotation.result.pos = "annotation.result.csv",
#   annotation.result.neg = "annotation.result.csv",
#   group = c("W03", "W30"),
#   max.isotope = 4,
#   uni.test = "t",
#   column ="hilic",
#   pos.path = "/home/jasper/work/MetDNA/fly pos/metABM result",
#   neg.path = "/home/jasper/work/MetDNA/fly neg/metABM result",
#   output.path = "/home/jasper/work/MetDNA/fly pos and neg 2",
#   correct = TRUE,
#   p.cutoff = 0.01,
#   mz.tol = 25,
#   rt.tol1 = 3,# absolute, second
#   rt.tol2 = 30,#relative, %
#   cor.tol = 0,
#   int.tol = 500,
#   threads = 3,
#   use.old.null = TRUE,
#   species = "dme",
#   use.all.kegg.id = FALSE)

#
# annotation.result <- readr::read_csv("annotation.result.csv")
# sample.info <- read.csv("sample.info.csv", stringsAsFactors = FALSE)
# sample.info <- sample.info[sample.info$group %in% c('W03', "W30"),]
# # sample.info <- sample.info[sample.info$group %in% c("4w", "12w", "24w", "32w", "52w"),]
# p.value <- uniTest(sample = annotation.result, sample.info = sample.info,
#                    uni.test = "t", correct = TRUE)
#
#
# # annotation.result <- readr::read_csv("annotation.result.csv")
# # colnames(annotation.result)
# # idx1 <- grep("_4w_", colnames(annotation.result))
# # annotation.result[,idx1] <- annotation.result[,idx1]/median.4
# #
# # idx2 <- grep("_12w_", colnames(annotation.result))
# # annotation.result[,idx2] <- annotation.result[,idx2]/median.12
# #
# # idx3 <- grep("_24w_", colnames(annotation.result))
# # annotation.result[,idx3] <- annotation.result[,idx3]/median.24
# #
# # idx4 <- grep("_32w_", colnames(annotation.result))
# # annotation.result[,idx4] <- annotation.result[,idx4]/median.32
# #
# # idx5 <- grep("_52w_", colnames(annotation.result))
# # annotation.result[,idx5] <- annotation.result[,idx5]/median.52
# #
# # write.csv(annotation.result, "annotation.result.csv", row.names = FALSE)
#
# sum(p.value < 0.05)
# system.time(metModule(annotation.result = "annotation.result.csv",
#                       group = c("W03", "W30"),
#                       path = "/home/jasper/work/metABM/fly pos/MRN annotation",
#                       output.path = "/home/jasper/work/metABM/fly pos/DSN identification/W03 vs W30 test",
#                       polarity = "positive",
#                       uni.test = "t",
#                       column = "hilic",
#                       p.cutoff = 0.05,
#                       species = "dme",
#                       use.old.null = FALSE,
#                       use.all.kegg.id = FALSE,
#                       correct = TRUE, threads = 3))
#
# # system.time(metTransform(annotation.result = "annotation.result.csv",
# #                          group = c("Normal", "MCI"), control.group = "Normal"))
#
#
# system.time(singleTransform(annotation.result = "annotation.result.csv",
#                             data.type = "metabolomics", sample.info = "sample.info.csv",
#                             group = c("4w", "32w"),
#                             scale = TRUE, scale.method = "pareto",
#                             trans.to = "pathway", species = "dme", method = "sum"
#                             ))
#
#
#
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/W03 vs W30/Module_MSE analysis")
# dir()
# dn.msea <- read.csv("Dysregulated_Network_MSEA.csv", stringsAsFactors = FALSE)
# rownames(dn.msea) <- dn.msea$X
# dn.msea <- dn.msea[,-1]
# dn.msea <- new(Class = 'moduleInfo',
#                Module.name = "Dysregulate_Network",
#                Module.size = 0,
#                Module.impact = 0,
#                p.value = 0,
#                Detected.number = 0,
#                Hidden.number = 0,
#                All.metabolite.id = "1",
#                Detected.ID = "1",
#                Hidden.ID = "1",
#                All.metabolite.name = "1",
#                Detected.metabolite.name = "1",
#                Hidden.metabolite.name = "1",
#                msea = dn.msea
# )
#
#
# moduleplot(dn.msea, n = 30, mar = c(5,21,4,2), main = "(Metabolomics)\nW03 vs W30")
#
# dn <- read.table("Dysregulated_Networks.txt", header =TRUE, as.is = TRUE)
# attr <- read.table("Dysregulated_Networks_attr.txt", header =TRUE, as.is = TRUE)
# dn.id <- unique(attr$ID)
# setwd("..")
# score <- read.csv('module.score.csv', stringsAsFactors = FALSE)
# idx1 <- grep("W03", colnames(score))
# idx2 <- grep("W30", colnames(score))
#
#
# # fc <- rep("no", 14)
# fc <- apply(score, 1, function(x) {
#   median(as.numeric(x)[idx2])/median(as.numeric(x)[idx1])
# })
#
# names(fc) <- score[,1]
# load("module.result")
#
# idx <- which(module.result$p.value < 0.05)
# module.result1 <- module.result[idx,]
# all.id <- module.result1$All.metabolite.id
# all.id <- lapply(all.id, function(x){
#   strsplit(x, split = ";")[[1]]
# })
#
# score$X
#
# fc1 <- fc
# fc1[fc>1] <- "Increase"
# fc[fc < 1] <- "Decrease"
#
#
#
# fc.attr <- sapply(dn.id, function(x){
# temp <- which(unlist(lapply(all.id, function(y){
#     is.element(x, y)
#   })))
# paste("module", temp, sep = "")
# })
#
#
# fc.attr1 <- unlist(lapply(fc.attr, function(x){
#   unname(fc1[match(x, names(fc1))])
# }))
#
# fc.attr2 <- cbind(fc.attr1, fc.attr)
#
# rownames(fc.attr2) == attr$ID
#
# attr <- data.frame(attr, fc.attr2, stringsAsFactors = FALSE)
# attr[,3][which(attr[,2] == "Hidden")] <- "no"
# colnames(attr)[c(3,4)] <- c("fc", "module")
# write.table(attr, "attr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/P03 vs P30/Module_MSE analysis")
# dir()
# dn.msea <- read.csv("Dysregulated_Network_MSEA.csv", stringsAsFactors = FALSE)
# rownames(dn.msea) <- dn.msea$X
# dn.msea <- dn.msea[,-1]
# dn.msea <- new(Class = 'moduleInfo',
#                Module.name = "Dysregulate_Network",
#                Module.size = 0,
#                Module.impact = 0,
#                p.value = 0,
#                Detected.number = 0,
#                Hidden.number = 0,
#                All.metabolite.id = "1",
#                Detected.ID = "1",
#                Hidden.ID = "1",
#                All.metabolite.name = "1",
#                Detected.metabolite.name = "1",
#                Hidden.metabolite.name = "1",
#                msea = dn.msea
# )
#
#
# moduleplot(dn.msea, n = 30, mar = c(5,20,4,2), main = "(Metabolomics)\nP03 vs P30")
#
# dn <- read.table("Dysregulated_Networks.txt", header =TRUE, as.is = TRUE)
# attr <- read.table("Dysregulated_Networks_attr.txt", header =TRUE, as.is = TRUE)
# dn.id <- unique(attr$ID)
# setwd("..")
# score <- read.csv('module.score.csv', stringsAsFactors = FALSE)
# idx1 <- grep("P03", colnames(score))
# idx2 <- grep("P30", colnames(score))
#
#
# # fc <- rep("no", 14)
# fc <- apply(score, 1, function(x) {
#   median(as.numeric(x)[idx2])/median(as.numeric(x)[idx1])
# })
#
# names(fc) <- score[,1]
# load("module.result")
#
# idx <- which(module.result$p.value < 0.05)
# module.result1 <- module.result[idx,]
# all.id <- module.result1$All.metabolite.id
# all.id <- lapply(all.id, function(x){
#   strsplit(x, split = ";")[[1]]
# })
#
# score$X
#
# fc1 <- fc
# fc1[fc>1] <- "Increase"
# fc[fc < 1] <- "Decrease"
#
#
#
# fc.attr <- sapply(dn.id, function(x){
#   temp <- which(unlist(lapply(all.id, function(y){
#     is.element(x, y)
#   })))
#   paste("module", temp, sep = "")
# })
#
#
# fc.attr1 <- unlist(lapply(fc.attr, function(x){
#   unname(fc1[match(x, names(fc1))])
# }))
#
# fc.attr2 <- cbind(fc.attr1, fc.attr)
#
# rownames(fc.attr2) == attr$ID
#
# attr <- data.frame(attr, fc.attr2, stringsAsFactors = FALSE)
# attr[,3][which(attr[,2] == "Hidden")] <- "no"
# colnames(attr)[c(3,4)] <- c("fc", "module")
# write.table(attr, "attr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/W03 vs P03/Module_MSE analysis")
# dir()
# dn.msea <- read.csv("Dysregulated_Network_MSEA.csv", stringsAsFactors = FALSE)
# rownames(dn.msea) <- dn.msea$X
# dn.msea <- dn.msea[,-1]
# dn.msea <- new(Class = 'moduleInfo',
#                Module.name = "Dysregulate_Network",
#                Module.size = 0,
#                Module.impact = 0,
#                p.value = 0,
#                Detected.number = 0,
#                Hidden.number = 0,
#                All.metabolite.id = "1",
#                Detected.ID = "1",
#                Hidden.ID = "1",
#                All.metabolite.name = "1",
#                Detected.metabolite.name = "1",
#                Hidden.metabolite.name = "1",
#                msea = dn.msea
# )
#
#
# moduleplot(dn.msea, n = 30, mar = c(5,20,4,2), main = "Metabolomics\n(W03 vs P03)")
#
# dn <- read.table("Dysregulated_Networks.txt", header =TRUE, as.is = TRUE)
# attr <- read.table("Dysregulated_Networks_attr.txt", header =TRUE, as.is = TRUE)
# dn.id <- unique(attr$ID)
# setwd("..")
# score <- read.csv('module.score.csv', stringsAsFactors = FALSE)
# idx1 <- grep("W03", colnames(score))
# idx2 <- grep("P03", colnames(score))
#
#
# # fc <- rep("no", 14)
# fc <- apply(score, 1, function(x) {
#   median(as.numeric(x)[idx2])/median(as.numeric(x)[idx1])
# })
#
# names(fc) <- score[,1]
# load("module.result")
#
# idx <- which(module.result$p.value < 0.05)
# module.result1 <- module.result[idx,]
# all.id <- module.result1$All.metabolite.id
# all.id <- lapply(all.id, function(x){
#   strsplit(x, split = ";")[[1]]
# })
#
# score$X
#
# fc1 <- fc
# fc1[fc>1] <- "Increase"
# fc[fc < 1] <- "Decrease"
#
#
#
# fc.attr <- sapply(dn.id, function(x){
#   temp <- which(unlist(lapply(all.id, function(y){
#     is.element(x, y)
#   })))
#   paste("module", temp, sep = "")
# })
#
#
# fc.attr1 <- unlist(lapply(fc.attr, function(x){
#   unname(fc1[match(x, names(fc1))])
# }))
#
# fc.attr2 <- cbind(fc.attr1, fc.attr)
#
# rownames(fc.attr2) == attr$ID
#
# attr <- data.frame(attr, fc.attr2, stringsAsFactors = FALSE)
# attr[,3][which(attr[,2] == "Hidden")] <- "no"
# colnames(attr)[c(3,4)] <- c("fc", "module")
# write.table(attr, "attr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#
#
#
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/W30 vs P30/Module_MSE analysis")
# dir()
# dn.msea <- read.csv("Dysregulated_Network_MSEA.csv", stringsAsFactors = FALSE)
# rownames(dn.msea) <- dn.msea$X
# dn.msea <- dn.msea[,-1]
# dn.msea <- new(Class = 'moduleInfo',
#                Module.name = "Dysregulate_Network",
#                Module.size = 0,
#                Module.impact = 0,
#                p.value = 0,
#                Detected.number = 0,
#                Hidden.number = 0,
#                All.metabolite.id = "1",
#                Detected.ID = "1",
#                Hidden.ID = "1",
#                All.metabolite.name = "1",
#                Detected.metabolite.name = "1",
#                Hidden.metabolite.name = "1",
#                msea = dn.msea
# )
#
#
# moduleplot(dn.msea, n = 30, mar = c(5,20,4,2), main = "Metabolomics\n(W30 vs P30)")
#
# dn <- read.table("Dysregulated_Networks.txt", header =TRUE, as.is = TRUE)
# attr <- read.table("Dysregulated_Networks_attr.txt", header =TRUE, as.is = TRUE)
# dn.id <- unique(attr$ID)
# setwd("..")
# score <- read.csv('module.score.csv', stringsAsFactors = FALSE)
# idx1 <- grep("W30", colnames(score))
# idx2 <- grep("P30", colnames(score))
#
#
# # fc <- rep("no", 14)
# fc <- apply(score, 1, function(x) {
#   median(as.numeric(x)[idx2])/median(as.numeric(x)[idx1])
# })
#
# names(fc) <- score[,1]
# load("module.result")
#
# idx <- which(module.result$p.value < 0.05)
# module.result1 <- module.result[idx,]
# all.id <- module.result1$All.metabolite.id
# all.id <- lapply(all.id, function(x){
#   strsplit(x, split = ";")[[1]]
# })
#
# score$X
#
# fc1 <- fc
# fc1[fc>1] <- "Increase"
# fc[fc < 1] <- "Decrease"
#
#
#
# fc.attr <- sapply(dn.id, function(x){
#   temp <- which(unlist(lapply(all.id, function(y){
#     is.element(x, y)
#   })))
#   paste("module", temp, sep = "")
# })
#
#
# fc.attr1 <- unlist(lapply(fc.attr, function(x){
#   unname(fc1[match(x, names(fc1))])
# }))
#
# fc.attr2 <- cbind(fc.attr1, fc.attr)
#
# rownames(fc.attr2) == attr$ID
#
# attr <- data.frame(attr, fc.attr2, stringsAsFactors = FALSE)
# attr[,3][which(attr[,2] == "Hidden")] <- "no"
# colnames(attr)[c(3,4)] <- c("fc", "module")
# write.table(attr, "attr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#
#
#
#
# #####################AD
# setwd("/home/jasper/work/metABM/AD POS/DSN identification/Normal vs MCI/Module_MSE analysis")
# dir()
# dn.msea <- read.csv("Dysregulated_Network_MSEA.csv", stringsAsFactors = FALSE)
# rownames(dn.msea) <- dn.msea$X
# dn.msea <- dn.msea[,-1]
# dn.msea <- new(Class = 'moduleInfo',
#                Module.name = "Dysregulate_Network",
#                Module.size = 0,
#                Module.impact = 0,
#                p.value = 0,
#                Detected.number = 0,
#                Hidden.number = 0,
#                All.metabolite.id = "1",
#                Detected.ID = "1",
#                Hidden.ID = "1",
#                All.metabolite.name = "1",
#                Detected.metabolite.name = "1",
#                Hidden.metabolite.name = "1",
#                msea = dn.msea
# )
#
#
# moduleplot(dn.msea, n = 30, mar = c(5,20,4,2), main = "Metabolomics\n(Normal vs MCI)")
#
# dn <- read.table("Dysregulated_Networks.txt", header =TRUE, as.is = TRUE)
# attr <- read.table("Dysregulated_Networks_attr.txt", header =TRUE, as.is = TRUE)
# dn.id <- unique(attr$ID)
# setwd("..")
# score <- read.csv('module.score.csv', stringsAsFactors = FALSE)
# sample.info <- read.csv("sample.info.csv", stringsAsFactors = FALSE)
# sample.info[,1] == colnames(score)
#
# idx1 <- sample.info[,1][grep("Normal", sample.info[,2])]
# idx2 <- sample.info[,1][grep("MCI", sample.info[,2])]
#
# idx1 <- match(idx1, colnames(score))
# idx2 <- match(idx2, colnames(score))
#
# # fc <- rep("no", 14)
# fc <- apply(score, 1, function(x) {
#   median(as.numeric(x)[idx2])/median(as.numeric(x)[idx1])
# })
#
# names(fc) <- score[,1]
# load("module.result")
#
# idx <- which(module.result$p.value < 0.05)
# module.result1 <- module.result[idx,]
# all.id <- module.result1$All.metabolite.id
# all.id <- lapply(all.id, function(x){
#   strsplit(x, split = ";")[[1]]
# })
#
# score$X
#
# fc1 <- fc
# fc1[fc>1] <- "Increase"
# fc[fc < 1] <- "Decrease"
#
#
#
# fc.attr <- sapply(dn.id, function(x){
#   temp <- which(unlist(lapply(all.id, function(y){
#     is.element(x, y)
#   })))
#   paste("module", temp, sep = "")
# })
#
#
# fc.attr1 <- unlist(lapply(fc.attr, function(x){
#   unname(fc1[match(x, names(fc1))])
# }))
#
# fc.attr2 <- cbind(fc.attr1, fc.attr)
#
# rownames(fc.attr2) == attr$ID
#
# attr <- data.frame(attr, fc.attr2, stringsAsFactors = FALSE)
# attr[,3][which(attr[,2] == "Hidden")] <- "no"
# colnames(attr)[c(3,4)] <- c("fc", "module")
# write.table(attr, "attr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#
#
# #####################AD
# setwd("/home/jasper/work/metABM/AD POS/DSN identification/Normal vs AD/Module_MSE analysis")
# dir()
# dn.msea <- read.csv("Dysregulated_Network_MSEA.csv", stringsAsFactors = FALSE)
# rownames(dn.msea) <- dn.msea$X
# dn.msea <- dn.msea[,-1]
# dn.msea <- new(Class = 'moduleInfo',
#                Module.name = "Dysregulate_Network",
#                Module.size = 0,
#                Module.impact = 0,
#                p.value = 0,
#                Detected.number = 0,
#                Hidden.number = 0,
#                All.metabolite.id = "1",
#                Detected.ID = "1",
#                Hidden.ID = "1",
#                All.metabolite.name = "1",
#                Detected.metabolite.name = "1",
#                Hidden.metabolite.name = "1",
#                msea = dn.msea
# )
#
#
# moduleplot(dn.msea, n = 30, mar = c(5,20,4,2), main = "Metabolomics\n(Normal vs AD)")
#
# dn <- read.table("Dysregulated_Networks.txt", header =TRUE, as.is = TRUE)
# attr <- read.table("Dysregulated_Networks_attr.txt", header =TRUE, as.is = TRUE)
# dn.id <- unique(attr$ID)
# setwd("..")
# score <- read.csv('module.score.csv', stringsAsFactors = FALSE)
# sample.info <- read.csv("sample.info.csv", stringsAsFactors = FALSE)
# sample.info[,1] == colnames(score)
#
# idx1 <- sample.info[,1][grep("Normal", sample.info[,2])]
# idx2 <- sample.info[,1][grep("AD", sample.info[,2])]
#
# idx1 <- match(idx1, colnames(score))
# idx2 <- match(idx2, colnames(score))
#
# # fc <- rep("no", 14)
# fc <- apply(score, 1, function(x) {
#   median(as.numeric(x)[idx2])/median(as.numeric(x)[idx1])
# })
#
# names(fc) <- score[,1]
# load("module.result")
#
# idx <- which(module.result$p.value < 0.05)
# module.result1 <- module.result[idx,]
# all.id <- module.result1$All.metabolite.id
# all.id <- lapply(all.id, function(x){
#   strsplit(x, split = ";")[[1]]
# })
#
# score$X
#
# fc1 <- fc
# fc1[fc>1] <- "Increase"
# fc1[fc < 1] <- "Decrease"
#
#
#
# fc.attr <- sapply(dn.id, function(x){
#   temp <- which(unlist(lapply(all.id, function(y){
#     is.element(x, y)
#   })))
#   paste("module", temp, sep = "")
# })
#
#
# fc.attr1 <- unlist(lapply(fc.attr, function(x){
#   unname(fc1[match(x, names(fc1))])
# }))
#
# fc.attr2 <- cbind(fc.attr1, fc.attr)
#
# rownames(fc.attr2) == attr$ID
#
# attr <- data.frame(attr, fc.attr2, stringsAsFactors = FALSE)
# attr[,3][which(attr[,2] == "Hidden")] <- "no"
# colnames(attr)[c(3,4)] <- c("fc", "module")
# write.table(attr, "attr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#
#
#
#
#
#
# setwd("/home/jasper/work/metABM/fly pos/DSN identification/QC W03 W30 P03 P30")
# dir()
# # annotation.result <- readr::read_csv("annotation.result.csv")
# # sample.name <- grep("Sample", colnames(annotation.result), value = TRUE)
# # group <- stringr::str_extract(string = sample.name, pattern = "[WPE]{1}[0-9]+")
# # group[is.na(group)] <- 'QC'
#
# # sample.info <- data.frame(sample.name, group, stringsAsFactors = FALSE)
#
# # write.csv(sample.info, "sample.info.csv", row.names = FALSE)
# system.time(metTransform(annotation.result = "annotation.result.csv",threads = 4,
#                          group = c("QC","W03", "W30", "P03", 'P30'), control.group = "QC", trans.to = "pathway"))
#
#
#
#
# setwd("/home/jasper/work/metABM/AD POS/DSN identification/QC vs Normal vs MCI vs AD")
# dir()
# annotation.result <- readr::read_csv("annotation.result.csv")
# qc <- readr::read_csv("qc.pos.csv")
# qc1 <- qc[,sample(1:200, 100)]
# write.csv(qc1, "qc.csv", row.names = FALSE)
# annotation.result <- annotation.result[,grep("YL", colnames(annotation.result))]
# annotation.result <- data.frame(annotation.result, qc1, stringsAsFactors = FALSE)
# write.csv(annotation.result, "annotation.result.csv", row.names = FALSE)
# sample.info <- read.csv("sample.info.csv", stringsAsFactors = F)
# sample.info <- sample.info[sample.info$group != "QC",]
# sample.name <- sample.info$sample.name
# group <- sample.info$group
# sample.name <- c(grep("QC", colnames(annotation.result), value = TRUE), sample.name)
# group <- c(rep("QC", 100), group)
# sample.info <- data.frame(sample.name, group, stringsAsFactors = FALSE)
# write.csv(sample.info, "sample.info.csv", row.names = FALSE)
# system.time(metTransform(annotation.result = "annotation.result.csv",threads = 6,
#                          group = c("QC","Normal", "MCI", "AD"), control.group = "QC", trans.to = "pathway"))
#







# setwd("/home/jasper/work/MetDNA/annotation validation/POS/test5")
# dir()
# metABM(use.default.md = TRUE, remain = TRUE)
#
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test1", "MRN_annotation_result_POS","remain.idx"))
# idx1 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test2", "MRN_annotation_result_POS","remain.idx"))
# idx2 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test3", "MRN_annotation_result_POS","remain.idx"))
# idx3 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test4", "MRN_annotation_result_POS","remain.idx"))
# idx4 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test5", "MRN_annotation_result_POS","remain.idx"))
# idx5 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test6", "MRN_annotation_result_POS","remain.idx"))
# idx6 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test7", "MRN_annotation_result_POS","remain.idx"))
# idx7 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test8", "MRN_annotation_result_POS","remain.idx"))
# idx8 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test9", "MRN_annotation_result_POS","remain.idx"))
# idx9 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test10", "MRN_annotation_result_POS","remain.idx"))
# idx10 <- remain.idx
# load(file.path(".", "MRN_annotation_result_POS","remain.idx"))
#
# data1 <- readr::read_csv(file.path(".", "ms2_match_result_POS", "ms2.match.annotation.result.csv"))
# data1 <- as.data.frame(data1)
#
# data2 <- readr::read_csv(file.path(".", "MRN_annotation_result_POS", "MRN.annotation.result.csv"))
# data2 <- as.data.frame(data2)
#
# sum(!is.na(data1$hits.reverse))
# sum(!is.na(data2$ID))
#
#
#
# data("inHouse.compound")
#
# id1 <- unlist(lapply(strsplit(data1$hits.reverse[remain.idx], split = ";"), function(x){
#   temp <- stringr::str_extract(string = x, pattern = "labid\\{[0-9A-Za-z]+\\}")
#   temp <- gsub(pattern = "labid\\{", replacement = "", x = temp)
#   temp <- gsub(pattern = "\\}", replacement = "", x = temp)
#   paste(inHouse.compound$KEGG.ID[match(temp, inHouse.compound$Lab.ID)], collapse = ";")
# }))
#
#
# id2 <- data2$ID[remain.idx]
#
# cbind(id1, id2)
#
# save(id1, file = "id1", compress = "xz")
# save(id2, file = "id2", compress = "xz")
#
#
#
# setwd("/home/jasper/work/MetDNA/annotation validation/POS/all seed")
# dir()
# metABM(use.default.md = TRUE, remain = FALSE)
#
# data1 <- readr::read_csv(file.path(".", "ms2_match_result_POS", "ms2.match.annotation.result.csv"))
# data1 <- as.data.frame(data1)
#
# data2 <- readr::read_csv(file.path(".", "MRN_annotation_result_POS", "MRN.annotation.result.csv"))
# data2 <- as.data.frame(data2)
#
# sum(!is.na(data1$hits.reverse))
# sum(!is.na(data2$ID))
#
#
# data1 <- readr::read_csv(file.path(".", "ms2_match_result_POS", "ms2.match.annotation.result.csv"))
# data1 <- as.data.frame(data1)
#
#
#
# ##################
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test1", "MRN_annotation_result_POS","remain.idx"))
# idx1 <- remain.idx
#
# load("/home/jasper/work/MetDNA/annotation validation/POS/test1/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test1/id2")
#
# test1 <- data.frame(data1[idx1,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test1) <- c("peak.name", "id1", "id2")
# save(test1, file = "test1", compress = "xz")
#
#
# ##test2
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test2", "MRN_annotation_result_POS","remain.idx"))
# idx2 <- remain.idx
#
# load("/home/jasper/work/MetDNA/annotation validation/POS/test2/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test2/id2")
#
# test2 <- data.frame(data1[idx2,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test2) <- c("peak.name", "id1", "id2")
# save(test2, file = "test2", compress = "xz")
#
# ##test3
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test3", "MRN_annotation_result_POS","remain.idx"))
# idx3 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test3/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test3/id2")
#
# test3 <- data.frame(data1[idx3,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test3) <- c("peak.name", "id1", "id2")
# save(test3, file = "test3", compress = "xz")
#
#
# ##test4
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test4", "MRN_annotation_result_POS","remain.idx"))
# idx4 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test4/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test4/id2")
#
# test4 <- data.frame(data1[idx4,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test4) <- c("peak.name", "id1", "id2")
# save(test4, file = "test4", compress = "xz")
#
# ##test5
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test5", "MRN_annotation_result_POS","remain.idx"))
# idx5 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test5/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test5/id2")
#
# test5 <- data.frame(data1[idx5,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test5) <- c("peak.name", "id1", "id2")
# save(test5, file = "test5", compress = "xz")
#
#
#
# ##test6
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test6", "MRN_annotation_result_POS","remain.idx"))
# idx6 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test6/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test6/id2")
#
# test6 <- data.frame(data1[idx6,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test6) <- c("peak.name", "id1", "id2")
# save(test6, file = "test6", compress = "xz")
#
#
# ##test7
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test7", "MRN_annotation_result_POS","remain.idx"))
# idx7 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test7/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test7/id2")
#
# test7 <- data.frame(data1[idx7,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test7) <- c("peak.name", "id1", "id2")
# save(test7, file = "test7", compress = "xz")
#
#
#
# ##test8
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test8", "MRN_annotation_result_POS","remain.idx"))
# idx8 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test8/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test8/id2")
#
# test8 <- data.frame(data1[idx8,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test8) <- c("peak.name", "id1", "id2")
# save(test8, file = "test8", compress = "xz")
#
#
# ##test9
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test9", "MRN_annotation_result_POS","remain.idx"))
# idx9 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test9/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test9/id2")
#
# test9 <- data.frame(data1[idx9,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test9) <- c("peak.name", "id1", "id2")
# save(test9, file = "test9", compress = "xz")
#
#
# ##test9
# load(file.path("/home/jasper/work/MetDNA/annotation validation/POS/test9", "MRN_annotation_result_POS","remain.idx"))
# idx10 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/POS/test10/id1")
# load("/home/jasper/work/MetDNA/annotation validation/POS/test10/id2")
#
# test10 <- data.frame(data1[idx10,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test10) <- c("peak.name", "id1", "id2")
# save(test10, file = "test10", compress = "xz")
#
#
#
# test1 <- test1[!is.na(test1$id2),]
# test2 <- test2[!is.na(test2$id2),]
# test3 <- test3[!is.na(test3$id2),]
# test4 <- test4[!is.na(test4$id2),]
# test5 <- test5[!is.na(test5$id2),]
# test6 <- test6[!is.na(test6$id2),]
# test7 <- test7[!is.na(test7$id2),]
# test8 <- test8[!is.na(test8$id2),]
# test9 <- test9[!is.na(test9$id2),]
# test10 <- test10[!is.na(test10$id2),]
#
# dim(test1)
# dim(test2)
# dim(test3)
# dim(test4)
# dim(test5)
# dim(test6)
# dim(test7)
# dim(test8)
# dim(test9)
# dim(test10)
#
#
# unique(c(test1$peak.name,test2$peak.name,test3$peak.name,test4$peak.name,test5$peak.name,test6$peak.name,test7$peak.name,
#          test8$peak.name,test9$peak.name,test10$peak.name))
#
#
#
# test <- rbind(test1, test2, test3, test4,test5,test6,test7,test8,test9,test10)
# colnames(test) <- c("Peak.name","Real.identification", "MRN.annotation")
# save(test, file = "test")
#
#
#
# setwd("/home/jasper/work/MetDNA/annotation validation/POS/summary")
# load("test")
#
# note <- rep(NA, nrow(test))
# note[which(test$Real.identification == test$MRN.annotation)] <- "Right"
# note[is.na(note)] <- "Wrong"
# test <- data.frame(test, note, stringsAsFactors = FALSE)
# write.csv(test, "test.csv")
#
#
#
#
#
#
#
#
# ################negative
# setwd("/home/jasper/work/MetDNA/annotation validation/NEG/test3")
# dir()
# metABM(use.default.md = TRUE, remain = TRUE, polarity = "negative")
#
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test1", "MRN_annotation_result_NEG","remain.idx"))
# idx1 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test2", "MRN_annotation_result_NEG","remain.idx"))
# idx2 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test3", "MRN_annotation_result_NEG","remain.idx"))
# idx3 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test4", "MRN_annotation_result_NEG","remain.idx"))
# idx4 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test5", "MRN_annotation_result_NEG","remain.idx"))
# idx5 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test6", "MRN_annotation_result_NEG","remain.idx"))
# idx6 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test7", "MRN_annotation_result_NEG","remain.idx"))
# idx7 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test8", "MRN_annotation_result_NEG","remain.idx"))
# idx8 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test9", "MRN_annotation_result_NEG","remain.idx"))
# idx9 <- remain.idx
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test10", "MRN_annotation_result_NEG","remain.idx"))
# idx10 <- remain.idx
#
# setwd("/home/jasper/work/MetDNA/annotation validation/NEG/test10")
# load(file.path(".", "MRN_annotation_result_NEG","remain.idx"))
#
# data1 <- readr::read_csv(file.path(".", "ms2_match_result_NEG", "ms2.match.annotation.result.csv"))
# data1 <- as.data.frame(data1)
#
# data2 <- readr::read_csv(file.path(".", "MRN_annotation_result_NEG", "MRN.annotation.result.csv"))
# data2 <- as.data.frame(data2)
#
# sum(!is.na(data1$hits.reverse))
# sum(!is.na(data2$ID))
#
#
#
# data("inHouse.compound")
#
# id1 <- unlist(lapply(strsplit(data1$hits.reverse[remain.idx], split = ";"), function(x){
#   temp <- stringr::str_extract(string = x, pattern = "labid\\{[0-9A-Za-z]+\\}")
#   temp <- gsub(pattern = "labid\\{", replacement = "", x = temp)
#   temp <- gsub(pattern = "\\}", replacement = "", x = temp)
#   paste(inHouse.compound$KEGG.ID[match(temp, inHouse.compound$Lab.ID)], collapse = ";")
# }))
#
#
# id2 <- data2$ID[remain.idx]
#
# cbind(id1, id2)
#
# save(id1, file = "id1", compress = "xz")
# save(id2, file = "id2", compress = "xz")
#
#
#
# setwd("/home/jasper/work/MetDNA/annotation validation/NEG/all seed")
# dir()
# metABM(use.default.md = TRUE, remain = FALSE, polarity = "negative")
#
# data1 <- readr::read_csv(file.path(".", "ms2_match_result_NEG", "ms2.match.annotation.result.csv"))
# data1 <- as.data.frame(data1)
#
# data2 <- readr::read_csv(file.path(".", "MRN_annotation_result_NEG", "MRN.annotation.result.csv"))
# data2 <- as.data.frame(data2)
#
# sum(!is.na(data1$hits.reverse))
# sum(!is.na(data2$ID))
#
#
# data1 <- readr::read_csv(file.path(".", "ms2_match_result_NEG", "ms2.match.annotation.result.csv"))
# data1 <- as.data.frame(data1)
#
#
#
# ##################
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test1", "MRN_annotation_result_NEG","remain.idx"))
# idx1 <- remain.idx
#
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test1/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test1/id2")
#
# test1 <- data.frame(data1[idx1,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test1) <- c("peak.name", "id1", "id2")
# save(test1, file = "test1", compress = "xz")
#
#
# ##test2
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test2", "MRN_annotation_result_NEG","remain.idx"))
# idx2 <- remain.idx
#
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test2/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test2/id2")
#
# test2 <- data.frame(data1[idx2,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test2) <- c("peak.name", "id1", "id2")
# save(test2, file = "test2", compress = "xz")
#
# ##test3
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test3", "MRN_annotation_result_NEG","remain.idx"))
# idx3 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test3/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test3/id2")
#
# test3 <- data.frame(data1[idx3,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test3) <- c("peak.name", "id1", "id2")
# save(test3, file = "test3", compress = "xz")
#
#
# ##test4
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test4", "MRN_annotation_result_NEG","remain.idx"))
# idx4 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test4/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test4/id2")
#
# test4 <- data.frame(data1[idx4,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test4) <- c("peak.name", "id1", "id2")
# save(test4, file = "test4", compress = "xz")
#
# ##test5
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test5", "MRN_annotation_result_NEG","remain.idx"))
# idx5 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test5/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test5/id2")
#
# test5 <- data.frame(data1[idx5,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test5) <- c("peak.name", "id1", "id2")
# save(test5, file = "test5", compress = "xz")
#
#
#
# ##test6
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test6", "MRN_annotation_result_NEG","remain.idx"))
# idx6 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test6/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test6/id2")
#
# test6 <- data.frame(data1[idx6,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test6) <- c("peak.name", "id1", "id2")
# save(test6, file = "test6", compress = "xz")
#
#
# ##test7
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test7", "MRN_annotation_result_NEG","remain.idx"))
# idx7 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test7/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test7/id2")
#
# test7 <- data.frame(data1[idx7,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test7) <- c("peak.name", "id1", "id2")
# save(test7, file = "test7", compress = "xz")
#
#
#
# ##test8
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test8", "MRN_annotation_result_NEG","remain.idx"))
# idx8 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test8/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test8/id2")
#
# test8 <- data.frame(data1[idx8,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test8) <- c("peak.name", "id1", "id2")
# save(test8, file = "test8", compress = "xz")
#
#
# ##test9
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test9", "MRN_annotation_result_NEG","remain.idx"))
# idx9 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test9/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test9/id2")
#
# test9 <- data.frame(data1[idx9,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test9) <- c("peak.name", "id1", "id2")
# save(test9, file = "test9", compress = "xz")
#
#
# ##test9
# load(file.path("/home/jasper/work/MetDNA/annotation validation/NEG/test9", "MRN_annotation_result_NEG","remain.idx"))
# idx10 <- remain.idx
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test10/id1")
# load("/home/jasper/work/MetDNA/annotation validation/NEG/test10/id2")
#
# test10 <- data.frame(data1[idx10,c(1)], id1, id2, stringsAsFactors = FALSE)
# colnames(test10) <- c("peak.name", "id1", "id2")
# save(test10, file = "test10", compress = "xz")
#
#
#
# test1 <- test1[!is.na(test1$id2),]
# test2 <- test2[!is.na(test2$id2),]
# test3 <- test3[!is.na(test3$id2),]
# test4 <- test4[!is.na(test4$id2),]
# test5 <- test5[!is.na(test5$id2),]
# test6 <- test6[!is.na(test6$id2),]
# test7 <- test7[!is.na(test7$id2),]
# test8 <- test8[!is.na(test8$id2),]
# test9 <- test9[!is.na(test9$id2),]
# test10 <- test10[!is.na(test10$id2),]
#
# dim(test1)
# dim(test2)
# dim(test3)
# dim(test4)
# dim(test5)
# dim(test6)
# dim(test7)
# dim(test8)
# dim(test9)
# dim(test10)
#
#
# unique(c(test1$peak.name,test2$peak.name,test3$peak.name,test4$peak.name,test5$peak.name,test6$peak.name,test7$peak.name,
#          test8$peak.name,test9$peak.name,test10$peak.name))
#
#
#
# test <- rbind(test1, test2, test3, test4,test5,test6,test7,test8,test9,test10)
# colnames(test) <- c("Peak.name","Real.identification", "MRN.annotation")
# save(test, file = "test")
#
#
#
# setwd("/home/jasper/work/MetDNA/annotation validation/NEG/summary")
# load("test")
#
# note <- rep(NA, nrow(test))
# note[which(test$Real.identification == test$MRN.annotation)] <- "Right"
# note[is.na(note)] <- "Wrong"
# test <- data.frame(test, note, stringsAsFactors = FALSE)
# write.csv(test, "test.csv")





##find right seeds for annotation
# setwd("/home/jasper/work/MetDNA/demo/fly/POS/dp0.8")
#
# ms2Annotation(ms1.file = "data.csv",
#               mz.tol = 25,
#               rt.tol = 10,
#               dp.cutoff = 0.8,
#               polarity = "positive",
#               path = ".",
#               output.path = file.path(".", "MS2_match_result"),
#               ms2.type = "mgf",
#               column = "hilic",
#               ce = "30")
#
# load("/home/jasper/work/MetDNA/demo/fly/POS/dp0.8/MS2_match_result/intermediate_data/ms2")
#
# data("kegg.rpair2", envir = environment())
# data("kegg.compound", envir = environment())
# data("inHouse.compound", envir = environment())
#
#
# ##add prediction RT to inhouse compound and KEGG compound
# load("rt.result")
# inhouse.rt <- rt.result[[1]]
# kegg.rt <- rt.result[[2]]
#
# idx <- match(inHouse.compound$Lab.ID, row.names(inhouse.rt))
#
# inHouse.compound <- data.frame(inHouse.compound, inhouse.rt[idx,],
#                                stringsAsFactors = FALSE)
#
# ##Add the predicted RT to inHouse compound and KEGG compound
#
# inHouse.compound$RT[is.na(inHouse.compound$RT)] <-
#   median(inHouse.compound$RT[!is.na(inHouse.compound$RT)])
#
# idx <- match(kegg.compound$ID, row.names(kegg.rt))
#
# kegg.compound <- data.frame(kegg.compound, kegg.rt[idx,],
#                             stringsAsFactors = FALSE)
#
# kegg.compound$RT[is.na(kegg.compound$RT)] <-
#   median(kegg.compound$RT[!is.na(kegg.compound$RT)])
#
#
#
# column <- "hilic"
# polarity <- "positive"
#
# if(column == "hilic") {
#   data("adduct.table.hilic", envir = environment())
#   adduct.table <- adduct.table.hilic
#   rm(list = c("adduct.table.hilic"))
#   gc()
# }
#
# if(column == "rp") {
#   data("adduct.table.rp", envir = environment())
#   adduct.table <- adduct.table.rp
#   rm(list = c("adduct.table.rp"))
#   gc()
# }
#
# adduct.table <- adduct.table[adduct.table$mode == polarity,]
#
# path = file.path(".", "MS2_match_result")
# output.path = file.path(".", "MRN_annotation_result")
#
#
# data <- readr::read_csv(file.path(path, "ms2.match.annotation.result.csv"))
# data <- as.data.frame(data)
#
# load("rt.all")
#
#
# data1 <- removeByRT(data = data, rt.data = rt.all, rt.tol = 10, direction = "less")
#
#
# sum(!is.na(data1$hits.reverse))
#
# write.csv(data1, "data1.csv", row.names = FALSE)
#
#
# metA(data = "data1.csv",
#      ms2 = ms2,
#      adduct.table = adduct.table,
#      max.isotope = 4,
#      polarity = "positive",
#      metabolite = kegg.compound,
#      metabolic.network = kegg.rpair2,
#      inHouse.compound = inHouse.compound,
#      threads = 3,
#      score.cutoff = 0,
#      path = ".",
#      output.path = file.path(".", "right annotation"),
#      rt.filter = FALSE)
#
#
# #find wrong seeds for annotation
# dir.create("dp0.1")
# setwd("dp0.1")
# ms2Annotation(ms1.file = "data.csv",
#               mz.tol = 25,
#               rt.tol = 10,
#               dp.cutoff = 0.1,
#               polarity = "positive",
#               path = ".",
#               output.path = file.path(".", "MS2_match_result"),
#               ms2.type = c("mgf", "mzXML", "msp"),
#               column = "hilic",
#               ce = "30")
#
#
#
#
# load("/home/jasper/work/MetDNA/demo/fly/POS/dp0.1/MS2_match_result/intermediate_data/ms2")
#
# data("kegg.rpair2", envir = environment())
# data("kegg.compound", envir = environment())
# data("inHouse.compound", envir = environment())
#
# ##add prediction RT to inhouse compound and KEGG compound
# load("rt.result")
# inhouse.rt <- rt.result[[1]]
# kegg.rt <- rt.result[[2]]
#
# idx <- match(inHouse.compound$Lab.ID, row.names(inhouse.rt))
#
# inHouse.compound <- data.frame(inHouse.compound, inhouse.rt[idx,],
#                                stringsAsFactors = FALSE)
#
# ##Add the predicted RT to inHouse compound and KEGG compound
#
# inHouse.compound$RT[is.na(inHouse.compound$RT)] <-
#   median(inHouse.compound$RT[!is.na(inHouse.compound$RT)])
#
# idx <- match(kegg.compound$ID, row.names(kegg.rt))
#
# kegg.compound <- data.frame(kegg.compound, kegg.rt[idx,],
#                             stringsAsFactors = FALSE)
#
# kegg.compound$RT[is.na(kegg.compound$RT)] <-
#   median(kegg.compound$RT[!is.na(kegg.compound$RT)])
#
#
#
#
#
#
#
#
# column <- "hilic"
# polarity <- "positive"
#
# if(column == "hilic") {
#   data("adduct.table.hilic", envir = environment())
#   adduct.table <- adduct.table.hilic
#   rm(list = c("adduct.table.hilic"))
#   gc()
# }
#
# if(column == "rp") {
#   data("adduct.table.rp", envir = environment())
#   adduct.table <- adduct.table.rp
#   rm(list = c("adduct.table.rp"))
#   gc()
# }
#
# adduct.table <- adduct.table[adduct.table$mode == polarity,]
#
# path = file.path(".", "MS2_match_result")
# output.path = file.path(".", "MRN_annotation_result")
#
#
# data <- readr::read_csv(file.path(path, "ms2.match.annotation.result.csv"))
# data <- as.data.frame(data)
#
# load("rt.all")
#
#
# data2 <- removeByDP(data = data, dp.cutoff = 0.2, direction = "less")
# data2 <- removeByRT(data = data2, rt.data = rt.all, rt.tol = 240, direction = "bigger")
#
#
#
# sum(!is.na(data2$hits.reverse))
# sum(!is.na(data2$hits.forward))
#
# # write.csv(data1, "data1.csv", row.names = FALSE)
# write.csv(data2, "data2.csv", row.names = FALSE)
#
#
# metA(data = "data2.csv",
#      ms2 = ms2,
#      adduct.table = adduct.table,
#      max.isotope = 4,
#      polarity = "positive",
#      metabolite = kegg.compound,
#      metabolic.network = kegg.rpair2,
#      inHouse.compound = inHouse.compound,
#      threads = 3,
#      score.cutoff = 0,
#      path = ".",
#      output.path = file.path(".", "wrong annotation"),
#      rt.filter = FALSE)
#
#
#
# temp <- unlist(lapply(tags2.after.redundancy.remove, function(x){
#   annotation <- x@annotation
#   if(length(annotation) == 0) return(NA)
#   unlist(lapply(annotation, function(y) y$ms2.sim))
# }))
#
# load("tags.result")
# ###############negative
# setwd("/home/jasper/work/MetDNA/demo/fly/NEG")
#
#
# ##ms2 annotation first
# ms2Annotation(ms1.file = "data.csv",
#               mz.tol = 25,
#               rt.tol = 10,
#               polarity = "negative",
#               path = ".",
#               output.path = file.path(".", "MS2_match_result"),
#               ms2.type = c("mgf", "mzXML", "msp"),
#               column = "hilic",
#               ce = "30")
#
#
#
#
# load("/home/jasper/work/MetDNA/demo/fly/NEG/MS2_match_result/intermediate_data/ms2")
#
# data("kegg.rpair2", envir = environment())
# data("kegg.compound", envir = environment())
# data("inHouse.compound", envir = environment())
#
#
#
#
# ##add prediction RT to inhouse compound and KEGG compound
# load("rt.result")
# inhouse.rt <- rt.result[[1]]
# kegg.rt <- rt.result[[2]]
#
# idx <- match(inHouse.compound$Lab.ID, row.names(inhouse.rt))
#
# inHouse.compound <- data.frame(inHouse.compound, inhouse.rt[idx,],
#                                stringsAsFactors = FALSE)
#
# ##Add the predicted RT to inHouse compound and KEGG compound
#
# inHouse.compound$RT[is.na(inHouse.compound$RT)] <-
#   median(inHouse.compound$RT[!is.na(inHouse.compound$RT)])
#
# idx <- match(kegg.compound$ID, row.names(kegg.rt))
#
# kegg.compound <- data.frame(kegg.compound, kegg.rt[idx,],
#                             stringsAsFactors = FALSE)
#
# kegg.compound$RT[is.na(kegg.compound$RT)] <-
#   median(kegg.compound$RT[!is.na(kegg.compound$RT)])
#
#
#
# column <- "hilic"
# polarity <- "negative"
#
# if(column == "hilic") {
#   data("adduct.table.hilic", envir = environment())
#   adduct.table <- adduct.table.hilic
#   rm(list = c("adduct.table.hilic"))
#   gc()
# }
#
# if(column == "rp") {
#   data("adduct.table.rp", envir = environment())
#   adduct.table <- adduct.table.rp
#   rm(list = c("adduct.table.rp"))
#   gc()
# }
#
# adduct.table <- adduct.table[adduct.table$mode == polarity,]
#
# path = file.path(".", "MS2_match_result")
# output.path = file.path(".", "MRN_annotation_result")
#
#
# data <- readr::read_csv(file.path(path, "ms2.match.annotation.result.csv"))
# data <- as.data.frame(data)
#
# load("rt.all")
#
#
# data1 <- removeByRT(data = data, rt.data = rt.all, rt.tol = 10, direction = "less")
# data2 <- removeByRT(data = data, rt.data = rt.all, rt.tol = 130, direction = "bigger")
#
# sum(!is.na(data1$hits.reverse))
# sum(!is.na(data2$hits.reverse))
#
# write.csv(data1, "data1.csv", row.names = FALSE)
# write.csv(data2, "data2.csv", row.names = FALSE)
#
# metA(data = "data1.csv",
#      ms2 = ms2,
#      adduct.table = adduct.table,
#      max.isotope = 4,
#      polarity = "negative",
#      metabolite = kegg.compound,
#      metabolic.network = kegg.rpair2,
#      inHouse.compound = inHouse.compound,
#      threads = 3,
#      score.cutoff = 0,
#      path = ".",
#      output.path = file.path(".", "right annotation"),
#      rt.filter = FALSE)
#
#
#
# metA(data = "data2.csv",
#      ms2 = ms2,
#      adduct.table = adduct.table,
#      max.isotope = 4,
#      polarity = "negative",
#      metabolite = kegg.compound,
#      metabolic.network = kegg.rpair2,
#      inHouse.compound = inHouse.compound,
#      threads = 3,
#      score.cutoff = 0,
#      path = ".",
#      output.path = file.path(".", "wrong annotation"),
#      rt.filter = FALSE)

