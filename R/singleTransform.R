# title singleTransform
# description Transform metabolie/gene level information to module/pathway
# level information.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param annotation.result The name of metabolite/gene sample. Column 1 must
# be "name" of metabolite/gene.
# param data.type metabolomics or gene.
# param sample.info The name of sample information. Column 1 is samle.name
# and column 2 is group.
# param group Which group you want to use.
# param trans.to module or pathway.
# param scale Scale sample or not.
# param scale.method The method of scale.
# param species The species.
# param which.peak Select which peak to represent the peak group.
# param polarity The polarity.
# param column The column.
# param method The method to process the module/pathway.
# param path The directory where data are in.
# param output.path The directory you want to output results.
# param metabolic.network kegg.rpair2
# return No.
# export
# singleTransform(annotation.result = "annotation.result.csv",
#                 data.type = "metabolomics",
#                 sample.info = "sample.info.csv",
#                 group,
#                 trans.to = "module",
#                 scale = TRUE,
#                 scale.method = "auto",
#                 species = "dme",
#                 which.peak = "adduct",
#                 polarity = "positive",
#                 column = "hilic",
#                 method = "mean",
#                 path = "/home/jasper/work/metABM/fly pos/DSN identification/W03 vs W30 module test",
#                 output.path = "/home/jasper/work/metABM/fly pos/DSN identification/W03 vs W30 transformation test")

#----------------------------------------------------------------------------
#transform metabolite information to quantitative module information



# annotation.result.pos = sample.pos
# annotation.result.neg = sample.neg
# anno = anno
# polarity = polarity
# group.data = group.data
# data.type = "metabolomics"
# sample.info = sample.info
# group = group
# trans.to = "pathway"
# scale = TRUE
# scale.method = "pareto"
# species = species
# which.peak = "intensity"
# column = column
# method = "mean"
# output.path = file.path(output.path, "Quantitative_information")
# metabolic.network = kegg.rpair2


# singleTransform(annotation.result.pos = sample.pos,
#                 annotation.result.neg = sample.neg,
#                 anno = anno,
#                 polarity = polarity,
#                 group.data = group.data,
#                 data.type = "metabolomics",
#                 sample.info = sample.info,
#                 group = group,
#                 trans.to = "pathway",
#                 scale = TRUE,
#                 scale.method = "pareto",
#                 species = species,
#                 which.peak = "intensity",
#                 column = column,
#                 method = "mean",
#                 output.path = file.path(output.path, "Quantitative_information"),
#                 metabolic.network = kegg.rpair2)



setGeneric(name = "singleTransform",
           def = function(annotation.result.pos = annotation.result.pos,
                          annotation.result.neg = annotation.result.neg,
                          anno,
                          polarity,
                          group.data,
                          data.type = c("metabolomics", "gene"),
                          sample.info = sample.info,
                          group,
                          trans.to = c("module", "pathway"),
                          scale = TRUE,
                          scale.method = c("pareto", "auto"),
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf"),
                          which.peak = c("intensity", "adduct"),
                          column = c("hilic", "rp"),
                          method = c("mean", "sum", "median"),
                          output.path = ".",
                          metabolic.network = kegg.rpair2){


             dir.create(output.path)
             data.type <- match.arg(data.type)
             trans.to <- match.arg(trans.to)
             species <- match.arg(species)
             method <- match.arg(method)
             which.peak <- match.arg(which.peak)
             column <- match.arg(column)

             if(missing(group)){
               warning("You don't set group, so use the default group in sample.info.\n")
               group <- unique(sample.info$group)
             }


             if(any(!is.element(group, unique(sample.info[,2])))){
               idx <- which(!is.element(group, unique(sample.info[,2])))
               stop(paste(group[idx], collapse = " and "),
                    ifelse(length(idx) > 1, " are ", " is "),
                    "not in your sample.info.")
             }

             sample.info <- sample.info[sample.info[,2] %in% group,]

             ##get sample data from annotation.result
             if(polarity == "positive" | polarity == "both"){
               sample.pos <- annotation.result.pos[,match(sample.info[,1], colnames(annotation.result.pos))]
               sample.pos <- as.data.frame(sample.pos)
             }

             if(polarity == "negative" | polarity == "both"){
               sample.neg <- annotation.result.neg[,match(sample.info[,1], colnames(annotation.result.neg))]
               sample.neg <- as.data.frame(sample.neg)
             }


             ##combine sample.pos and sample.neg
             switch(polarity,
                    "positive" = {sample <- sample.pos},
                    "negative" = {sample <- sample.neg},
                    "both" = {sample <- rbind(sample.pos, sample.neg)})



             if(trans.to == 'gene'){
               zero.per <- apply(sample, 1, function(x) {sum(x == 0)/ncol(sample)})
               if(sum(zero.per > 0.5) > 0) {
                 sample <- sample[which(zero.per <= 0.5),,drop = FALSE]
                 cat('There are', sum(zero.per > 0.5),
                     ifelse(data.type == "metabolomics", "peaks", 'genes'),
                     "have more than 50% zero values and have been removed from the sample.")
               }
             }

             ##scale sample
             if(scale){
               switch(scale.method,
                      "auto" = {sample <- t(apply(sample, 1, function(x) {x/sd(x)}))},
                      "pareto" = {sample <- t(apply(sample, 1, function(x) {x/sqrt(sd(x))}))})
             }


             ##begin transform for metabolomics data
             if(data.type == "metabolomics"){
               # peak.name <- showTags2(tags2.after.kegg.matching, slot = "name")
               cat("\n")
               cat("Transform metabolite information to quantitative",
                   ifelse(trans.to == "module", "module", "pathway"),
                   "information.\n")

               ##for each module
               all.id <- unique(unname(unlist(group.data)))

               peak.data <- lapply(all.id, function(x){
                 temp.anno <- anno[which(x == anno$to),]
                 temp.grade <- as.numeric(stringr::str_extract(temp.anno$Confidence, pattern = "[0-9]+"))
                 temp.anno[which(temp.grade == min(temp.grade)),, drop = FALSE]
               })

               names(peak.data) <- all.id

               ##remove the node which have no mapped peak, those node may be
               ## hidden nodes is network
               peak.data <- peak.data[which(unlist(lapply(peak.data, nrow)) != 0)]

               ##for each peak group, select the peak which has the most intense intensity
               if(which.peak == "intensity"){
                 peak.data1 <- lapply(peak.data, function(x){
                   ##only one peak
                   if(nrow(x) == 1) {return(x)}
                   temp.group <- unique(x$group)
                   temp.idx <- lapply(temp.group, function(y) {which(x$group == y)})
                   temp.x <- lapply(temp.idx, function(idx) x[idx, , drop = FALSE])

                   temp.x <- lapply(temp.x, function(z){
                     if(nrow(z) == 1) return(z)
                     # z <- z[z$isotope == "[M]",]
                     temp.peak <- z$name
                     temp.idx <- which.max(apply(sample[match(temp.peak, rownames(sample)),,drop = FALSE], 1, median))
                     z <- z[temp.idx,,drop = FALSE]
                     z
                   })
                   x <- do.call(rbind, temp.x)
                   x
                 })
               }


               if(which.peak == "adduct"){

                 adduct.frequency <- unname(unlist(lapply(peak.data, function(x){
                   x$adduct
                 })))

                 adduct.rank <- names(sort(table(adduct.frequency), decreasing = TRUE))

                 # par(mar = c(5,12,4,2))
                 # temp2 <- barplot(sort(table(adduct.frequency), decreasing = FALSE), horiz = TRUE,
                 #         border = NA, col = "tan1", names.arg = NA, xlab = "Number",
                 #         cex.lab = 1.8, cex.axis = 1.5)
                 # axis(side = 2, at = sort(temp2[,1], TRUE), labels = adduct.rank, las = 1, cex.axis = 1.5)
                 # text(x = sort(table(adduct.frequency)+6, decreasing = FALSE),
                 #      y = temp2[,1], labels = sort(table(adduct.frequency), cex = 1.5))
                 peak.data1 <- lapply(peak.data, function(x){
                   ##only one peak
                   if(nrow(x) == 1) return(x)
                   temp.group <- unique(x$group)
                   temp.idx <- lapply(temp.group, function(y) {which(x$group == y)})
                   temp.x <- lapply(temp.idx, function(idx) x[idx, , drop = FALSE])

                   temp.x <- lapply(temp.x, function(z){
                     if(nrow(z) == 1) return(z)
                     if(any(z$isotope == "[M]")) {z <- z[z$isotope == "[M]",]}
                     temp.adduct <- z$adduct
                     temp.idx <- which.min(match(temp.adduct, adduct.rank))
                     z <- z[temp.idx,,drop = FALSE]
                     z
                   })
                   x <- do.call(rbind, temp.x)
                   x
                 })

               }

               rm(list = "peak.data")
               gc()

               ##remove id which are not in metabolomics.network
               peak.data1 <- peak.data1[which(names(peak.data1) %in% igraph::V(metabolic.network)$name)]
               ###get the subnetwork
               subnetwork <- igraph::subgraph(graph = metabolic.network,
                                              v = names(peak.data1))

               group.number <- unlist(lapply(peak.data1, function(x) {
                 length(unique(x$group))
               }))

               index <- which(group.number > 1)
               old.index <- c(index, 1)
               peak.data2 <- peak.data1
               # round <- 1
               while(length(index) > 0 & (length(old.index) - length(index)) > 0){

                 for(i in index){
                   temp.peak <- peak.data2[[i]]
                   temp.peak.name <- temp.peak$name
                   temp.data <- sample[match(temp.peak.name, rownames(sample)), , drop = FALSE]
                   temp.data <- apply(temp.data, 1, list)
                   temp.name <- names(peak.data1)[i]
                   neighbor <- getNeighbor(metabolite.id = temp.name,
                                           step = 1,
                                           graph = subnetwork)
                   if(is.null(neighbor)) {peak.data2[[i]] <- temp.peak; next()}
                   neighbor <- neighbor$Node.ID
                   temp.idx <- match(names(which(group.number[match(neighbor,
                                                                    names(group.number))] == 1)), names(group.number))
                   if(length(temp.idx) == 0) {peak.data2[[i]] <- temp.peak; next()}

                   neighbor.peak.name <- unlist(lapply(peak.data2[temp.idx], function(x) x$name))
                   neighbor.data <- sample[match(neighbor.peak.name, rownames(sample)), , drop = FALSE]
                   neighbor.data <- apply(neighbor.data, 1, list)

                   ##calculate the correlation
                   temp.correlation <- lapply(temp.data, function(x){
                     x <- as.numeric(x[[1]])
                     abs(unlist(lapply(neighbor.data, function(y) {
                       y <- as.numeric(y[[1]])
                       cor(x, y)
                     })))
                   })

                   temp.correlation.median <- lapply(temp.correlation, median)
                   temp.idx <- which.max(temp.correlation.median)
                   temp.peak <- temp.peak[temp.idx,,drop = FALSE]
                   peak.data2[[i]] <- temp.peak
                 }
                 group.number <- unlist(lapply(peak.data2, function(x) {
                   length(unique(x$group))
                 }))
                 old.index <- index
                 # index <- which(unlist(lapply(peak.data2, nrow)) > 1)
                 index <- which(group.number > 1)
               }

               ###if there are still node have more than two group, those nodes
               ##are detected nodes which are connected with other detected nodes
               ##by hidden nodes, so the peak group which has the biggest intensity
               ## is selected.
               if(length(index) > 0){
                 peak.data2[index] <- lapply(peak.data2[index], function(x){
                   temp.name <- x$name
                   # temp.sample <- sample[match(temp.name, rownames(sample)),,drop = FALSE]
                   temp.score <- as.numeric(anno$score[match(temp.name, anno$name)])
                   temp.idx <- which.max(temp.score)
                   x <- x[temp.idx, , drop = FALSE]
                   x
                 })
               }

               peak.data3 <- do.call(rbind, peak.data2)

               rm(list = c("peak.data1", "subnetwork", "peak.data2"))
               gc()

               ##output node quantative information
               peak.name <- peak.data3$name
               temp.idx <- match(peak.name, rownames(sample))
               node.quantitative.data <- cbind(peak.data3, sample[temp.idx,])

               node.quantitative.data <- node.quantitative.data[,-c(2:3,5:6,8:23)]
               node.quantitative.data <- data.frame(node.quantitative.data[,3],
                                                    node.quantitative.data[,1],
                                                    node.quantitative.data[,2],
                                                    node.quantitative.data[,-c(1:3)],
                                                    stringsAsFactors = FALSE)
               colnames(node.quantitative.data)[1:3] <- c("ID", "peak.name", "annotation.type")

               data("kegg.compound", envir = environment())
               compound.name <- kegg.compound$Name[match(node.quantitative.data$ID, kegg.compound$ID)]
               compound.name <- unlist(lapply(strsplit(x = compound.name, split = ';'), function(x) x[1]))
               rm(list = c("kegg.compound"))

               node.quantitative.data <- data.frame(compound.name,
                                                    node.quantitative.data,
                                                    stringsAsFactors = FALSE)
               node.quantitative.data <- node.quantitative.data[!duplicated(node.quantitative.data$compound.name),,drop=FALSE]

               if(trans.to == "module"){
                 readr::write_csv(as.data.frame(node.quantitative.data),
                                  file.path(output.path, "DNA.node.quantitative.result.csv"))

               }else{
                 readr::write_csv(as.data.frame(node.quantitative.data),
                                  file.path(output.path, "pathway.node.quantitative.result.csv"))
               }

               rm(list = c("node.quantitative.data"))
               gc()

               ##transform metabolite information to group information
               group.sample <- lapply(group.data, function(x){
                 temp.idx <- match(x, peak.data3$to)
                 temp.idx <- temp.idx[!is.na(temp.idx)]
                 if(length(temp.idx) < 3) {return(NULL)}
                 temp.peak.name <- peak.data3$name[temp.idx]
                 temp.sample <- sample[match(temp.peak.name, rownames(sample)),,drop = FALSE]
               })

               names(group.sample) <- names(group.data)
               rm(list = "group.data")
               gc()

               remove.idx <- which(unlist(lapply(group.sample, is.null)))
               if(length(remove.idx > 0)){
                 group.sample <- group.sample[-remove.idx]
               }

               ##get module quantitative inforamtion
               group.quantitative.data <- lapply(group.sample, function(x){
                 x <- apply(x, 2, method)
                 x
               })

               group.quantitative.data <- do.call(rbind, group.quantitative.data)
               group.name <- names(group.sample)
               group.quantitative.data <- data.frame(group.name,
                                                     group.quantitative.data,
                                                     stringsAsFactors = FALSE)

               colnames(group.quantitative.data)[1] <- "name"

               readr::write_csv(as.data.frame(group.quantitative.data),
                                file.path(output.path,
                                          paste("DNA",ifelse(trans.to == "module", "module", "pathway"),
                                                "quantitative.result.csv", sep = ".")))

               rm(list = c("group.sample", "group.quantitative.data"))
               gc()

             }else{
               cat("\n")
               cat("Transform metabolite information to quantitative pathway information.\n")
               ##gene transform
               switch(species,
                      "dme" = {data("dme.gene.kegg.pathway", envir = environment())
                        pathway.data <- dme.gene.kegg.pathway
                        rm(dme.gene.kegg.pathway)
                        gc()},
                      "hsa" = {data("hsa.gene.kegg.pathway", envir = environment())
                        pathway.data <- hsa.gene.kegg.pathway
                        rm(hsa.gene.kegg.pathway)
                        gc()},
                      "mmu" = {data("mmu.gene.kegg.pathway", envir = environment())
                        pathway.data <- mmu.gene.kegg.pathway
                        rm(mmu.gene.kegg.pathway)
                        gc()},
                      "rat" = {data("rat.gene.kegg.pathway", envir = environment())
                        pathway.data <- rat.gene.kegg.pathway
                        rm(rat.gene.kegg.pathway)
                        gc()},
                      "bta" = {data("bta.gene.kegg.pathway", envir = environment())
                        pathway.data <- bta.gene.kegg.pathway
                        rm(bta.gene.kegg.pathway)
                        gc()},
                      "gga" = {data("gga.gene.kegg.pathway", envir = environment())
                        pathway.data <- gga.gene.kegg.pathway
                        rm(gga.gene.kegg.pathway)
                        gc()},
                      "dre" = {data("dre.gene.kegg.pathway", envir = environment())
                        pathway.data <- dre.gene.kegg.pathway
                        rm(dre.gene.kegg.pathway)
                        gc()},
                      "cel" = {data("cel.gene.kegg.pathway", envir = environment())
                        pathway.data <- cel.gene.kegg.pathway
                        rm(cel.gene.kegg.pathway)
                        gc()},
                      "sce" = {data("sce.gene.kegg.pathway", envir = environment())
                        pathway.data <- sce.gene.kegg.pathway
                        rm(sce.gene.kegg.pathway)
                        gc()},
                      "ath" = {data("ath.gene.kegg.pathway", envir = environment())
                        pathway.data <- ath.gene.kegg.pathway
                        rm(ath.gene.kegg.pathway)
                        gc()},
                      "smm" = {data("smm.gene.kegg.pathway", envir = environment())
                        pathway.data <- smm.gene.kegg.pathway
                        rm(smm.gene.kegg.pathway)
                        gc()},
                      "pfa" = {data("pfa.gene.kegg.pathway", envir = environment())
                        pathway.data <- pfa.gene.kegg.pathway
                        rm(pfa.gene.kegg.pathway)},
                      "tbr" = {data("tbr.gene.kegg.pathway", envir = environment())
                        pathway.data <- tbr.gene.kegg.pathway
                        rm(tbr.gene.kegg.pathway)
                        gc()},
                      "eco" = {data("eco.gene.kegg.pathway", envir = environment())
                        pathway.data <- eco.gene.kegg.pathway
                        rm(eco.gene.kegg.pathway)},
                      "ppu" = {data("ppu.gene.kegg.pathway", envir = environment())
                        pathway.data <- ppu.gene.kegg.pathway
                        rm(ppu.gene.kegg.pathway)
                        gc()},
                      "syf" = {data("syf.gene.kegg.pathway", envir = environment())
                        pathway.data <- syf.gene.gene.kegg.pathway
                        rm(syf.gene.gene.kegg.pathway)
                        gc()})



               sample <- as.data.frame(sample)
               pathway.sample <- lapply(pathway.data, function(x){
                 temp.idx <- match(x, rownames(sample))
                 temp.idx <- temp.idx[!is.na(temp.idx)]
                 if(length(temp.idx) == 0) return(NULL)
                 sample[temp.idx,,drop = FALSE]
               })

               remove.idx <- which(unlist(lapply(pathway.sample, is.null)))
               if(length(remove.idx) > 0){
                 pathway.sample <- pathway.sample[-remove.idx]
               }


               pathway.quantitative.data <- lapply(pathway.sample, function(x){
                 x <- apply(x, 2, method)
                 x
               })

               pathway.quantitative.data <- do.call(rbind, pathway.quantitative.data)
               pathway.quantitative.data <- cbind(rownames(pathway.quantitative.data), pathway.quantitative.data)
               colnames(pathway.quantitative.data)[1] <- "name"


               save(pathway.quantitative.data,
                    file = file.path(output.path, "pathway.quantitative.data"), compress = "xz")
               readr::write_csv(as.data.frame(pathway.quantitative.data),
                                file.path(output.path, "pathway.quantitative.data.csv"))

             }

             cat("\n")
             cat("singleTransform is done.\n")
           })
