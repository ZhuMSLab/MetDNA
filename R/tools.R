#---------------------------------------------------------------------------
#' @title uniTest
#' @description Univariate analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param sample The sample abundance. Column is sample and row is peak.
#' @param sample.info The sample.information
#' @param uni.test Test method.
#' @param alternative see ?t.test.
#' @param mu see ?t.test.
#' @param paired see ?t.test.
#' @param var.equal see ?t.test.
#' @param conf.level see ?t.test.
#' @param exact see ?t.test.
#' @param correct see ?t.test.
#' @param conf.int see ?t.test.
#' @param threads How many threads do you want to use? Default is the number of
#' your PC threads - 3.
#' @param ... See ?t.test.
#' @export


setGeneric(name = "uniTest",
           def = function(sample,
                          sample.info,
                          uni.test = c("t","wilcox","anova"),
                          alternative = c("two.sided", "less", "greater"),
                          mu = 0, paired = FALSE, var.equal = FALSE,
                          conf.level = 0.95, exact = NULL, correct = TRUE,
                          conf.int = FALSE,
                          threads = 3,...
           ){
             uni.test <- match.arg(uni.test)
             ###statistic analysis
             sample.info <- sample.info[order(as.character(sample.info[,1])),]
             sample <- sample[,match(sample.info[,1], colnames(sample)), drop = FALSE]
             # sample <- sample[,order(colnames(sample))]
             group <- sample.info$group
             unique.group <- unique(group)
             group1 <- lapply(unique.group, function(x) which(group == x))
             names(group1) <- unique.group
             # cat("Total",nrow(sample),"peaks\n")
             # cat("Univariate test:")

             pbapply::pboptions(type = "timer", style = 1)
             if(uni.test == "t"){
               p.value <- pbapply::pbapply(sample, 1, function(x) {
                 p <- try({
                   t.test(as.numeric(x[group1[[1]]]), as.numeric(x[group1[[2]]]),
                          alternative = alternative, mu = mu, paired = paired,
                          var.equal = var.equal, conf.level = conf.level)$p.value
                   },
                   silent = TRUE)
                 if(class(p) == "try-error"){
                   p <- 1
                 }
                 p
               })



               # for(i in 1:nrow(sample)){
               #   cat(i, " ")
               #   x <- as.numeric(sample[i,])
               #   t.test(as.numeric(x[group1[[1]]]), as.numeric(x[group1[[2]]]),
               #          alternative = alternative, mu = mu, paired = paired,
               #          var.equal = var.equal, conf.level = conf.level)$p.value
               # }

               if(correct) p.value <- p.adjust(p = p.value, method = "fdr")

             }

             if(uni.test == "wilcox"){
               p.value <- pbapply::pbapply(sample, 1, function(x) {
                 p <- try({
                   wilcox.test(as.numeric(x[group1[[1]]]),
                               as.numeric(x[group1[[2]]]), alternative = alternative,
                               mu = mu, paired = paired, exact = exact, correct = correct,
                               conf.int = conf.int, conf.level = conf.level)$p.value
                 },
                 silent = TRUE)
                 if(class(p) == "try-error"){
                   p <- 1
                 }
                 p
               })
             }

             if(uni.test == "anova"){
               p.value <- pbapply::pbapply(sample, 1, function(x){
                 temp.sample <- as.numeric(x)
                 temp.aov <- aov(formula = temp.sample~factor(group))
                 summary(temp.aov)[[1]][1,"Pr(>F)"]
               })
               if(correct) p.value <- p.adjust(p = p.value, method = "fdr")
             }

             p.value <- p.value
           })



#---------------------------------------------------------------------------
#' @title foldChange
#' @description Calculate fold change.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param sample The sample abundance. Column is sample and row is peak.
#' @param sample.info The sample.information
#' @param by.what Use median or mean to calculate fold change.
#' @param group Which group you want to use.
#' @export


setGeneric(name = "foldChange",
           def = function(sample,
                          sample.info,
                          by.what = c("median","mean"),
                          group = c("W03", "W30")
           ){
             sample <- as.data.frame(sample)
             sample.info <- as.data.frame(sample.info)
             by.what <- match.arg(by.what)
             temp.fun <- function(x, by.what = by.what){
               if(by.what == "median"){
                 median(x)
               }else{
                 mean(x)
               }
             }

             ###statistic analysis
             sample.info <- sample.info[order(as.character(sample.info[,1])),]
             sample <- sample[,order(colnames(sample)), drop = FALSE]
             group1 <- sample.info$group
             idx <- lapply(group, function(x) {which(group1 == x)})
             names(idx) <- group

             # cat("Total",nrow(sample),"peaks\n")
             # cat("Calculate fold change:")

             pbapply::pboptions(type = "timer", style = 1)
             fc <- pbapply::pbapply(sample, 1, function(x){
               x <- as.numeric(x)
               temp.fun(x[idx[[2]]], by.what = by.what)/temp.fun(x[idx[[1]]],
                                                                 by.what = by.what)
             })
             fc[is.na(fc)] <- 1
             fc[is.nan(fc)] <- 1
             fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
             return(fc)
           })



# title SXTMTmatch
# description Match two data according to mz and RT.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param data1 First data for matching, first column must be mz
# and seconod column must be rt.
# param data2 Second data for matching, first column must be mz
# and seconod column must be rt.
# param mz.tol mz tol for ms1 and ms2 data matching.
# param rt.tol RT tol for ms1 and ms2 data matching.
# return Return a result which give the matching result of data1 and database.
# export

setGeneric(name = "SXTMTmatch",
           def = function(data1,
                          data2,
                          mz.tol,
                          #rt.tol is relative
                          rt.tol = 30,
                          rt.error.type = c("relative", "abs")){
             rt.error.type <- match.arg(rt.error.type)
             #
             if (nrow(data1) == 0 | nrow(data2) == 0) {
               result <- NULL
               return(result)
             }
             # mz1 <- as.numeric(data1[, 1])
             # rt1 <- as.numeric(data1[, 2])
             info1 <- data1[,c(1,2)]
             info1 <- apply(info1, 1, list)

             mz2 <- as.numeric(data2[, 1])
             rt2 <- as.numeric(data2[, 2])

             result <- pbapply::pblapply(info1, function(x) {
               temp.mz1 <- x[[1]][[1]]
               temp.rt1 <- x[[1]][[2]]
               mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
               if(rt.error.type == "relative"){
                 rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
               }else{
                 rt.error <- abs(temp.rt1 - rt2)
               }

               j <- which(mz.error <= mz.tol & rt.error <= rt.tol)
               if(length(j) == 0){
                 matrix(NA, ncol = 7)
               }else{
                 cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
               }
             })

             if(length(result) == 1){
               result <- cbind(1,result[[1]])
             }else{
               result <- mapply(function(x,y){list(cbind(x,y))},
                                x <- 1:length(info1),
                                y = result)
               result <- do.call(rbind, result)
             }

             result <- matrix(result[which(!apply(result,1,function(x) any(is.na(x)))),], ncol = 8)
             if(nrow(result) == 0) return(NULL)
             colnames(result) <-
               c("Index1",
                 "Index2",
                 "mz1",
                 "mz2",
                 "mz error",
                 "rt1",
                 "rt2",
                 "rt error")
             result <- result
           })





setGeneric(name = "samplePlot",
           def = function(sample,
                          sample.info,
                          file.name,
                          group,
                          output.path = ".",
                          p.value,
                          fc,
                          beeswarm = TRUE,
                          col = c("lightseagreen", "salmon", "orchid4"),
                          pch = 19,
                          xlab = "Group",
                          ylab= "Dysregulated score"){
             sample.info <- as.data.frame(sample.info)
sample.info <- sample.info[sample.info$group %in% group,]
sample <- sample[,match(sample.info[,1], colnames(sample)), drop = FALSE]


group.idx <- lapply(group, function(x){
  which(sample.info$group == x)
})

group.sample <- lapply(group.idx, function(x){sample[,x, drop = FALSE]})
variable.name <- rownames(sample)
if(missing(file.name)) file.name <- variable.name


for(i in 1:nrow(sample)){
  pdf(file = file.path(output.path, paste(file.name[i], ".pdf", sep = "")),
      width = 7, height = 7)
  par(mar = c(5,5,4,2))
  temp.sample <- lapply(group.sample, function(x){as.numeric(x[i,])})

  ylim1 <- min(unlist(temp.sample))
  ylim2 <- max(unlist(temp.sample))
  if(ylim2 >=0){ylim2 <- 1.2*ylim2}else{ylim2 <- 0.8*ylim2}
    temp1 <- boxplot(temp.sample, col = col,
          xlab = xlab, ylab = ylab, names = group, cex.lab = 1.8,
          cex.axis = 1.5,
          main = variable.name[i], cex.main = 1.8,
          ylim = c(ylim1, ylim2), add = FALSE)

    #add something
segments(x0 = 1, y0 = ylim2*0.95, x1 = 2, y1  = ylim2*0.95, lwd = 1.5)

mapply(function(x,y){
  segments(x0 = y-0.05, y0 = max(x)*1.05, x1 = y+0.05, y1  = max(x)*1.05, lwd = 1.5)
},
x = temp.sample,
y = 1:length(temp.sample))


mapply(function(x,y){
  segments(x0 = y, y0 = max(x)*1.05, x1 = y, y1  = ylim2*0.95, lwd = 1.5)
},
x = temp.sample,
y = 1:length(temp.sample))

p <- p.value[i]
if(p < 0.001) label <- "***"
if(0.001 <= p & p < 0.01) label <- "**"
if(0.01 <= p & p < 0.05) label <- "*"
if(0.05 <= p) label <- "NS"

text(x = 1.5, y = ylim2*0.97, cex = 1.5,
     labels = paste(label, unname(round(fc[i],2)), sep = "/FC"))

# text(x = 1.5, y = ylim2*0.97, cex = 1.5,
#      labels = label)


beeswarm::beeswarm(temp.sample,
                   labels = F,add = TRUE,
                   pch = pch)
dev.off()

}


           })




##volcanoPlot
setGeneric(name = "volcanoPlot",
           def = function(p.value,
                          fc,
                          correct = TRUE,
                          p.cutoff = 0.05,
                          fc.cutoff = 1,
                          col = c("grey", "lightseagreen", "salmon"),
                          pch = 19,
                          cex = 1,
                          cex.lab = 1.8,
                          cex.axis = 1.5,
                          xlab = "log2Fold change",
                          ylab = "log10P-Value(adjusted)"){
             if(correct){
               ylab = "log10P-Value(adjusted)"
             }else{
               ylab = "log10P-Value"
             }

             p.value1 <- -log(as.numeric(p.value),10)
             fc1 <- log(as.numeric(fc),2)

             col1 <- rep(NA, length(p.value))
             col1[which(p.value > p.cutoff)] <- col[1]
             col1[which(p.value <= p.cutoff & fc >= fc.cutoff)] <- col[3]
             col1[which(p.value <= p.cutoff & fc <= 1/fc.cutoff)] <- col[2]
             plot(fc1, p.value1, cex = cex, cex.lab = cex.lab,
                  cex.axis = cex.axis, pch = pch,
                  xlab = xlab,
                  ylab = ylab,
                  col = col1)
             abline(h = -log(p.cutoff, 10), lty = 3, lwd = 2,
                    col = "salmon")
             if(fc.cutoff == 1){
               abline(v = log(fc.cutoff, 2), lty = 3, lwd = 2,
                      col = "salmon")
             }else{
               abline(v = log(fc.cutoff, 2), lty = 3, lwd = 2,
                      col = "salmon")
               abline(v = log(1/fc.cutoff, 2), lty = 3, lwd = 2,
                      col = "salmon")
             }


           })





##heatMap


setGeneric(name = "heatMap",
           def = function(sample,
                          sample.info,
                          group,
                          color = c("lightseagreen", "salmon", "orchid4"),
                          scale.method = c("auto", "pareto"),
                          show_rownames = TRUE,
                          show_colnames = FALSE,
                          border_color = NA,
                          fontsize = 10,
                          fontsize_row = 10,
                          fontsize_col = 10,
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          clustering_method = "complete",
                          display_numbers = FALSE
                          ){
             scale.method <- match.arg(scale.method)
             sample.info <- as.data.frame(sample.info)
             sample.info <- sample.info[sample.info$group %in% group,]
             sample <- sample[,match(sample.info$sample.name, colnames(sample)), drop = FALSE]

             group.idx <- lapply(group, function(x){
               which(sample.info$group == x)
             })


             # colnames(sample) <- stringr::str_extract(colnames(score), "W[0-9]{2}")
             if(scale.method == "auto"){
               sample1 <- t(apply(sample, 1, function(x) (x-mean(x))/sd(x)))
             }else{
               sample1 <- t(apply(sample, 1, function(x) (x-mean(x))/sqrt(sd(x))))
             }

             sample1.range <- abs(range(sample1))
             dif <- sample1.range[1] - sample1.range[2]
             if (dif < 0) {
               sample1[sample1 > sample1.range[1]] <- sample1.range[1]
             }
             if (dif > 0) {
               sample1[sample1 < -1 * sample1.range[2]] <- -1 * sample1.range[2]
             }

             annotation_col <- data.frame(Group = factor(c(sample.info$group)))


             rownames(annotation_col) <- sample.info$sample.name


             # Specify colors
             ann_col <- NULL
             for (i in seq_along(group)) {
               ann_col[i] <- color[i]
             }

             ann_colors = list(group = ann_col)
             names(ann_colors[[1]]) <- group


             pheatmap::pheatmap(sample1, colorRampPalette(c("navy", "white", "firebrick"))(50),
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors,
                      show_rownames = show_rownames,
                      show_colnames = show_colnames,
                      border_color = border_color,
                      fontsize = fontsize,
                      fontsize_row = fontsize_row,
                      fontsize_col = fontsize_col,
                      cluster_rows = cluster_rows,
                      cluster_cols = cluster_cols,
                      clustering_distance_rows = clustering_distance_rows,
                      clustering_distance_cols = clustering_distance_cols,
                      clustering_method = clustering_method,
                      display_numbers = display_numbers
                      )
           })




##heatMap
# mse.data <- readr::read_csv(dir()[1])
# mse.data <- as.data.frame(mse.data)
# rownames(mse.data) <- mse.data$X1
# mse.data <- mse.data[,-1]
# pathwayPlot(mse.data = mse.data)
setGeneric(name = "pathwayPlot",
           def = function(mse.data,
                          color = c("lightseagreen", "salmon", "orchid4"),
                          size.range = c(1,5),
                          q.cutoff = 0.05,
                          pch = 19,
                          text = TRUE
           ){

             if(nrow(mse.data) == 0){
               layout(1)
               par(mar = c(5, 5, 4, 2))
               plot(0, 0, col = "white", xaxt = "n", yaxt = "n",
                    xlab = "", ylab = "", cex.main = 1.5)
               text(x = 0, y = 0, labels = 'No pathways are enriched.',
                    cex = 1.5)
             }else{

               mse.data <- mse.data[mse.data$Overlap >0,]
               fdr <- as.numeric(mse.data$p.value)
               overlap <- mse.data$Overlap/mse.data$Pathway.length
               pathway.len <- mse.data$Pathway.length


               fdr1 <- -log(as.numeric(fdr),10)

               col1 <- rep(NA, length(fdr1))
               col1[which(fdr > q.cutoff)] <- color[1]
               col1[which(fdr <= q.cutoff)] <- color[2]

               lm.model <- lm(size.range~range(pathway.len))
               size1 <- coefficients(lm.model)[2]*pathway.len + coefficients(lm.model)[1]

               plot(overlap*100, fdr1, pch = pch,
                    xlab = "Overlap (%)",
                    # ylab = "-log10P-value(adjusted)",
                    ylab = "-log10P-value",
                    col = col1,
                    cex.lab = 1.8, cex.axis = 1.5,
                    cex = size1)

               abline(h = -log(q.cutoff, 10), lty = 3, lwd = 2,
                      col = "salmon")

               temp.idx <- which(fdr < q.cutoff)

               pathway.name <- rownames(mse.data)
               pathway.name <- unlist(lapply(strsplit(pathway.name, split = ";"), function(x) x[1]))

               if(text & length(temp.idx) > 0){
                 # text(x = 100*overlap[temp.idx], y = fdr1[temp.idx], labels = pathway.name[temp.idx])
                 maptools::pointLabel(x = 100*overlap[temp.idx], y = fdr1[temp.idx], pathway.name[temp.idx],
                                      cex = 1)
               }
             }


           })



# setGeneric(name = "pathwayPlot",
#            def = function(mse.data,
#                           color = c("lightseagreen", "salmon", "orchid4"),
#                           size.range = c(1,5),
#                           q.cutoff = 0.05,
#                           pch = 19,
#                           text = TRUE
#            ){
#
#              if(nrow(mse.data) == 0){
#                layout(1)
#                par(mar = c(5, 5, 4, 2))
#                plot(0, 0, col = "white", xaxt = "n", yaxt = "n",
#                     xlab = "", ylab = "", cex.main = 1.5)
#                text(x = 0, y = 0, labels = 'No pathways are enriched.',
#                     cex = 1.5)
#              }else{
#
#                mse.data <- mse.data[mse.data$Overlap >0,]
#                fdr <- as.numeric(mse.data$q.value)
#                overlap <- mse.data$Overlap/mse.data$Pathway.length
#                pathway.len <- mse.data$Pathway.length
#
#
#                fdr1 <- -log(as.numeric(fdr),10)
#
#                col1 <- rep(NA, length(fdr1))
#                col1[which(fdr > q.cutoff)] <- color[1]
#                col1[which(fdr <= q.cutoff)] <- color[2]
#
#                lm.model <- lm(size.range~range(pathway.len))
#                size1 <- coefficients(lm.model)[2]*pathway.len + coefficients(lm.model)[1]
#
#                plot(overlap*100, fdr1, pch = pch,
#                     xlab = "Overlap (%)",
#                     ylab = "-log10P-value(adjusted)",
#                     col = col1,
#                     cex.lab = 1.8, cex.axis = 1.5,
#                     cex = size1)
#
#                abline(h = -log(q.cutoff, 10), lty = 3, lwd = 2,
#                       col = "salmon")
#
#                temp.idx <- which(fdr < q.cutoff)
#
#                pathway.name <- rownames(mse.data)
#                pathway.name <- unlist(lapply(strsplit(pathway.name, split = ";"), function(x) x[1]))
#
#                if(text & length(temp.idx) > 0){
#                  # text(x = 100*overlap[temp.idx], y = fdr1[temp.idx], labels = pathway.name[temp.idx])
#                  maptools::pointLabel(x = 100*overlap[temp.idx], y = fdr1[temp.idx], pathway.name[temp.idx],
#                                       cex = 1)
#                }
#              }
#
#
#            })



setGeneric(name = "pathwayPlot2",
           def = function(pathway.result,
                          color = c("lightseagreen", "salmon", "orchid4"),
                          size.range = c(1,5),
                          p.cutoff = 0.05,
                          pch = 19,
                          text = TRUE
           ){
             pathway.result <- pathway.result[pathway.result$Overlap >0,]
             p1 <- as.numeric(pathway.result$adjusted.p.value.by.pvalue)
             p2 <- as.numeric(pathway.result$adjusted.p.value.by.overlap)

             # overlap <- mse.data$Overlap/mse.data$Pathway.length
             pathway.len <- pathway.result$Pathway.length

             p1 <- -log(as.numeric(p1),10)
             p2 <- -log(as.numeric(p2),10)

             p1[is.infinite(p1)] <- max(p1[!is.infinite(p1)])
             p2[is.infinite(p2)] <- max(p2[!is.infinite(p2)])

             col1 <- rep(NA, length(p1))
             col1[which(p1 > -log(p.cutoff, 10) & p2 > -log(p.cutoff, 10))] <- color[2]
             col1[is.na(col1)] <- color[1]

             lm.model <- lm(size.range~range(pathway.len))
             size1 <- coefficients(lm.model)[2]*pathway.len + coefficients(lm.model)[1]

             plot(p1, p2, pch = pch,
                  xlab = "-log10P-value(adjusted by P-value)",
                  ylab = "-log10P-value(adjusted by overlap)",
                  col = col1,
                  cex.lab = 1.8, cex.axis = 1.5,
                  cex = size1)

             abline(h = -log(p.cutoff, 10), lty = 3, lwd = 2,
                    col = "salmon")

             abline(v = -log(p.cutoff, 10), lty = 3, lwd = 2,
                    col = "salmon")

             temp.idx <- which(p1 > -log(p.cutoff, 10) & p2 > -log(p.cutoff, 10))

             pathway.name <- pathway.result$pathway.name

             if(text & length(temp.idx) > 0){
               # text(x = 100*overlap[temp.idx], y = fdr1[temp.idx], labels = pathway.name[temp.idx])
               maptools::pointLabel(x = p1[temp.idx], y = p2[temp.idx], pathway.name[temp.idx],
                                    cex = 1)
             }

           })




setGeneric(name = "moduleScatterPlot",
           def = function(module.result,
                          color = c("lightseagreen", "salmon", "orchid4"),
                          size.range = c(1,5),
                          p.cutoff = 0.05,
                          pch = 19,
                          text = TRUE
           ){
             # mse.data <- mse.data[mse.data$Overlap >0,]
             p <- as.numeric(module.result$p.value)
             impact <- as.numeric(module.result$Module.impact)
             module.size <- as.numeric(module.result$Module.size)


             log.p <- -log(as.numeric(p),10)
             log.p[is.infinite(log.p)] <- max(log.p[!is.infinite(log.p)])


             col1 <- rep(NA, length(p))
             col1[which(p > p.cutoff)] <- color[1]
             col1[which(p <= p.cutoff)] <- color[2]

             lm.model <- lm(size.range~range(module.size))
             size1 <- coefficients(lm.model)[2]*module.size + coefficients(lm.model)[1]

             plot(impact, log.p, pch = pch,
                  xlab = "Module impact",
                  ylab = "-log10P-value",
                  col = col1,
                  cex.lab = 1.8, cex.axis = 1.5,
                  cex = size1)

             abline(h = -log(p.cutoff, 10), lty = 3, lwd = 2,
                    col = "salmon")

             temp.idx <- which(p < p.cutoff)

             module.name <- module.result$Module.name
             # module.name <- unlist(lapply(strsplit(module.name, split = ";"), function(x) x[1]))

             if(text & length(temp.idx) > 0){
               # text(x = 100*overlap[temp.idx], y = fdr1[temp.idx], labels = pathway.name[temp.idx])
               maptools::pointLabel(x = impact[temp.idx], y = log.p[temp.idx], module.name[temp.idx],
                                    cex = 1)
             }

           })





setGeneric(name = "GetPathwayInfo",
           def = function(mse.data, node.quantitative.data,
                          sample.info,
                          group,
                          species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                      "dre", "cel", "sce", "ath", "smm", "pfa",
                                      "tbr", "eco", "ppu", "syf")){
             species <- match.arg(species)

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

  mse.data <- mse.data[mse.data$Overlap>0,]
  pathway.name.id <- rownames(mse.data)
  pathway.name <- unlist(lapply(strsplit(pathway.name.id, split = ";"), function(x) x[[1]]))
  pathway.id <- unlist(lapply(strsplit(pathway.name.id, split = ";"), function(x) x[[2]]))

  data("kegg.compound", envir = environment())

  for(i in nrow(mse.data)){
  temp.idx <- match(pathway.name.id[i], names(M))
  temp.idx2 <- which(node.quantitative.data$to %in% M[[temp.idx]])

  }


})


setGeneric(name = "dataPlot",
           def = function(data = "data.csv",
                          path = "."){
data <- readr::read_csv(data, col_types = readr::cols(), progress = FALSE)

mz <- as.numeric(data$mz)
rt <- as.numeric(data$rt)

plot(rt, mz, xlab = "Retention time (RT, second)",
     ylab = "Mass to charge ratio (m/z)",
     cex.lab = 1.8, cex.axis = 1.5, pch = 19,
     cex = ifelse(length(rt > 10000), 0.3, 0.5),
     col = "grey")
legend("topleft", legend = paste("Peaks:",length(mz)),
       bty = "n", cex = 1.3)
})




#' @title zipData
#' @description Compress analysis result to zip.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#'@param polarity Acquisition mode.
#'@export


setGeneric(name = "zipData",
           def = function(polarity = c("positive", "negative", "both")){
           polarity <- match.arg(polarity)
           now.path <- getwd()

           if(polarity != "both"){
             dir.create(ifelse(polarity == "positive","POS", 'NEG'))
             now.path1 <- strsplit(x = now.path, split = "/")[[1]]
             new.path <- paste(now.path1[-length(now.path1)], collapse = "/")
               data.name <- c('Analysis_report', "Dysregulated_network_analysis_result",
                                      "MRN_annotation_result",
                                      "MS2_match_result",
                              "MetDNA.parameters.csv")

               file.copy(from = file.path(new.path, ifelse(polarity == "positive","POS", 'NEG'), data.name),
                         to = ifelse(polarity == "positive","POS", 'NEG'),
                         recursive = TRUE, overwrite = TRUE)

               zip::zip(zipfile = ifelse(polarity == "positive","POS.zip", "NEG.zip"),
                        files = ifelse(polarity == "positive", "POS", 'NEG'), recurse = TRUE)
               file.copy(from = ifelse(polarity == "positive","POS.zip", "NEG.zip"), to = "..")
               unlink(x = ifelse(polarity == "positive","POS.zip", "NEG.zip"), recursive = TRUE)
               unlink(x = ifelse(polarity == "positive","POS", "NEG"), recursive = TRUE)
           }else{
             now.path1 <- strsplit(x = now.path, split = "/")[[1]]
             if(tail(now.path1, 1) == "POS" | tail(now.path1, 1) == "NEG"){
               new.path <- paste(now.path1[-length(now.path1)], collapse = "/")
               dir.create("POS and NEG")
               file.copy(file.path(new.path, "POS and NEG",
                                   setdiff(dir(file.path(new.path, "POS and NEG")), "POS and NEG")),
                         to = "POS and NEG", recursive = TRUE, overwrite = TRUE)

             unlink(file.path('POS and NEG',
                              "Dysregulated_network_analysis_result", "intermediate_data"), recursive = TRUE)


             zip::zip(zipfile = "POS and NEG.zip",
                      files = "POS and NEG", recurse = TRUE)

             file.copy(from = "POS and NEG.zip", to = "..")
             unlink(x = "POS and NEG.zip", recursive = TRUE)
             unlink(x = "POS and NEG", recursive = TRUE)

             }
           }

           })



#' @title sendMail
#' @description Sent results to user using email.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param email.address Email address.
#'@param polarity Acquisition mode.
#'@param from The address of sender.
#'@param user.name User name.
#'@param pw Pass word.
#'@param attach Attach or not.
#'@export


setGeneric(name = "sendMail",
           def = function(email.address = "zutaoshen@gmail.com",
                          polarity = c('positive', 'negative', 'both'),
                          from = "kagawabale@163.com",
                          user.name = "kagawabale",
                          pw = "19900521@2016",
                          attach = TRUE){
             now.path <- getwd()
             new.path <-
               paste(strsplit(now.path, split = "/")[[1]][-length(strsplit(now.path, split = "/")[[1]])], collapse = '/')
polarity <- match.arg(polarity)
# from <- "kagawabale@163.com"
to <- email.address

subject <- paste("MetDNA analysis is done", ifelse(polarity != "both", ifelse(polarity == "positive","(Positive mode)", "(Negative mode)"), "(Positive and Negative mode)"))
body <- "<html> <b>Dear user:</b>
         <br/>
         <br/>Thank you for using MetDNA.
         <br/>The analysis of your data is done and the results can be found in attachment.
         <br/>The attachment contains a analysis report (HTML) and a zip file for analysis result.
         <br/>Please don't be hesitate to send email to us (shenxt@sioc.ac.cn) if you have any questions.
         <br/>
         <br/> <b>Your sincerely,</b>
         <br/> <b>XiaotaoShen and Zhengjiaing Zhu<b>
         <br/><b>Website:</b> <a href = \"http://www.zhulab.cn/\"> www.zhulab.cn </a>
         <br/>
         <br/>
         <img src = \"http://a3.qpic.cn/psb?/V12nMOGs2VfuNZ/L0OOyQhz2rQKFWLs0NXqS2UyWhXP8g21XcJXoZDAE80!/b/dD0BAAAAAAAA&bo=sQd2AgAAAAADAOc!&rf=viewer_4\" heigh = '600', width = '600'>
         </html>"
smtp <- list(host.name = "smtp.163.com", port = 465,
             user.name = user.name,
             passwd = pw, ssl=TRUE)

if(polarity != "both"){
  attach.file1 <- file.path(new.path,
                            ifelse(polarity == "positive", "POS", "NEG"),
                            "Analysis_report",
                            "MetDNA.analysis.report.html")

  attach.file2 <- file.path(new.path, ifelse(polarity=="positive", "POS.zip", "NEG.zip"))

  attach.file <- c(attach.file1, attach.file2)
}else{
  attach.file1 <- file.path(new.path,
                            "POS and NEG",
                            "Analysis_report",
                            "MetDNA.analysis.report.html")

  attach.file2 <- file.path(new.path, "POS and NEG.zip")
  # attach.file3 <- file.path(new.path, "POS.zip")
  # attach.file4 <- file.path(new.path, "NEG.zip")

  attach.file <- c(attach.file1, attach.file2)
}

if(attach){
  mailR::send.mail(from = from, to = to, subject = subject, body = body,
                   smtp = smtp, authenticate = TRUE, send = TRUE, html = TRUE,
                   attach.files = attach.file
  )
}else{
  mailR::send.mail(from = from, to = to, subject = subject, body = body,
                   smtp = smtp, authenticate = TRUE, send = TRUE, html = TRUE)
}


           })






#-----------------------------------------------------------------------------
# title getCentrality
# description Get the centrality of all nodes in nwtwork.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param graph The graph.
# param type degree or betweenness.
# return  The centrality of all nodes.


setGeneric(name = "getCentrality",
           def = function(graph,
                          type = c("degree", "betweenness")){
             type = match.arg(type)
             node <- igraph::V(graph)$name
             if(type == "degree"){
               node.degree <- igraph::degree(graph = graph, v = node)
               return(node.degree)
             }

             node.betweenness <- lapply(node, function(x){
               i <- match(x, node)
               if(i == length(node)) return(NULL)
               from <- x
               to <- node[(i+1):length(node)]
               temp  <- igraph::shortest_paths(graph = graph,
                                               from = from, to = to)[[1]]
               temp <- lapply(temp, function(x) {
                 if(length(x) <= 2) return(NULL)
                 names(x)[-c(1,length(x))]
               })
               unlist(temp)
             })

             node.betweenness <- unlist(node.betweenness)
             node.betweenness <- table(node.betweenness)
             node0 <- setdiff(node, names(node.betweenness))
             if(length(node0) == 0) return(node.betweenness)
             add <- rep(0, length(node0))
             names(add) <- node0
             node.betweenness <- c(node.betweenness, add)
             node.betweenness <- node.betweenness[match(node, names(node.betweenness))]
             return(node.betweenness)
           })






#-----------------------------------------------------------------------------
#' @title getNeighbor
#' @description Get the neighbors of one metabolite in graph.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param metabolite.id The metabolite ID.
#' @param step The reaction step.
#' @param graph The kegg.rpair2.
#' @return  The neighbor information.
#' @export

setGeneric(name = "getNeighbor",
           def = function(metabolite.id,
                          step = 1,
                          graph){
             # load("kegg.rpair2")
             # if(!igraph::is.igraph(kegg.rpair1)) {
             #   kegg.rpair2 <- igraph::graph_from_graphnel(kegg.rpair1)
             # }
             if(!metabolite.id %in% igraph::V(graph)$name) return(NULL)

             temp.step <- 1
             temp.id <- metabolite.id
             cum.temp.id <- NULL
             while(temp.step <= step){
               result <-
                 unique(unlist(lapply(temp.id, function(x) names(igraph::neighbors(graph, x)))))
               cum.temp.id <- c(cum.temp.id, temp.id)
               temp.step <- temp.step + 1
               temp.id <- result
             }

             result <- setdiff(result, cum.temp.id)

             result <- sort(result)
             if(length(result) == 0) return(NULL)
             result <-
               data.frame(metabolite.id, result, step,
                          stringsAsFactors = FALSE)
             colnames(result) <- c('Metabolite.ID', 'Node.ID', "Step")
             result <- result
           })



# setGeneric(name = "getGraph",
#            def = function(from, to, edgemode = c("undirected"),
#                           format = c("graphNEL", "igraph"),
#                           weight = NULL){
#              edgemode <- match.arg(edgemode)
#              format <- match.arg(format)
#              if (is.null(weight)) weight <- 1
#              if (format == "graphNEL") {
#                stopifnot(length(from) == length(to))
#                nodes <- unique(c(from, to))
#                ge <- new("graphNEL", nodes = nodes, edgemode = edgemode)
#                g <- graph::addEdge(from, to, ge, weight)
#              }
#              else {
#                if (edgemode == "undirected") {
#                  edgemode <- FALSE
#                }
#                else {
#                  edgemode <- TRUE
#                }
#                g <- igraph::graph.edgelist(cbind(from, to), directed = edgemode)
#              }
#              return(g)
# })




# #####Acc like acc in graph
# setGeneric(name = "findAcc",
#            def = function(object, index, max.step) {
#              nN <- graph::numNodes(object)#node number
#              nNames<- graph::nodes(object)#node name
#              nIndex <- length(index)#index number
#              whN <- match(index, nNames)#index position in node
#              if( any(is.na(whN)) )
#                stop("unmatched node provided")
#
#              rval <- vector("list", length=nIndex)
#              names(rval) <- nNames[whN]
#              for( i in seq_len(nIndex)) {
#                marked <- rep(0, nN)#node marker
#                distv <- rep(0, nN)#node distance
#                names(distv) <- nNames
#                distx <- 1#
#                names(marked) <- nNames
#                nmkd <- 0#
#                marked[index[i]] <- 1#mark index node as 1
#                done <- FALSE
#                while( !done ) {
#                  minds <- nNames[marked==1]#find which node neighbor
#                  for( node in minds) {
#                    avec <- graph::adj(object, node)[[1]]
#                    avec <- avec[marked[avec]==0] #don't mark any already marked
#                    marked[avec] <- 1
#                    distv[avec] <- distx
#                  }
#                  marked[minds] <- 2
#                  distx <- distx+1
#                  newmk <- sum(marked==1)
#                  if( newmk == 0 | distx - 2 == max.step)
#                    done <- TRUE
#                }
#                marked[index[i]] <- 0 ##not the node itself
#                rval[[i]] <- distv[marked==2]
#              }
#              return(rval)
#            })
#
#
#
# ##
# setGeneric(name = "findPath",
#            def = function(object, node1, node2, max.step = 10){
#            step = 1
#            temp.idx <- NULL
#           from = node1
#           result <- list()
#   while(length(temp.idx) == 0 & step <= max.step){
#   temp.result <- findAcc(object, from, max.step = 1)
#   result[[step]] <- temp.result
#   name <- unname(unlist(lapply(temp.result, function(x) {names(x)})))
#   temp.idx <- which(name == node2)
#   step <- step + 1
#   from <- name
#   }
#
#   path.node <- list()
#   i <- step - 1
#   path.node[[i+1]] <- node2
#   before.name <- node2
#   while(!node1 %in% before.name){
#     temp.result <- result[[i]]
#     temp.result <- unlist(temp.result)
#     idx <- lapply(before.name, function(x) grep(x, names(temp.result)))
#     idx <- unlist(idx)
#     before.name <- names(temp.result)[idx]
#     before.name <-
#       unlist(lapply(before.name, function(x) strsplit(x, split = "\\.")[[1]][1]))
#     path.node[[i]] <- before.name
#     i <- i - 1
#   }
#
#   path.node <- unique(unlist(path.node))
#   subgraph <- graph::subGraph(snodes = path.node, object)
#   return(subgraph)
# })


# setGeneric(name = "graph2matrix",
#            def = function(object){
# node.name <- graph::nodes(object)
# edgeL <- graph::edgeL(object)
# idx <- lapply(edgeL, function(x) {x[[1]]})
# pair <- lapply(idx, function(x) {node.name[x]})
# pair <- mapply(function(x,y) {cbind(x, y)}, x = node.name, y = pair)
# pair <- do.call(rbind, pair)
# weight <- unname(unlist(graph::edgeWeights(object)))
# result <- data.frame(pair, weight, stringsAsFactors = FALSE)
# colnames(result) <- c("Node1.ID", "Node2.ID", "weight")
# #remove duplicated pair
# result1 <- apply(result, 1, list)
# result2 <- lapply(result1, function(x) {x[[1]]})
#
# dup <- rep(FALSE, length(result2))
# for(i in seq_along(result2)[-1]){
# before <- result2[1:(i-1)]
# dup[i] <-
#   any(unlist(lapply(before, function(x) {setequal(x, result2[[i]])})))
# }
#
# result3 <- result[!dup,]
# return(result3)
# })




#' @title memoryList
#' @description List the size of objects.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param decreasing Logical. List the object from big to small or not.
#' @param n List how many objects.
#' @return  The size of objects.
#' @export

setGeneric(name = "memoryList",
           def = function(decreasing = TRUE,
                          n = 10){
             object.name <- ls(envir = globalenv())
             cat("Total memory used:")
             temp <- pryr::mem_used()
             print(temp)
             cat("\n")
             if(length(object.name) != 0){
               object.size1 <-
                 sapply(object.name, function(x) {pryr::object_size(get(x))})
               if(decreasing){
                 object.size1 <- object.size1[order(object.size1, decreasing = TRUE)]
               }else{
                 object.size1 <- object.size1[order(object.size1, decreasing = FALSE)]
               }

               for(i in 1:n){
                 temp <- pryr::object_size(get(names(object.size1)[i]))
                 cat(names(object.size1)[i],":")
                 print(temp)
                 cat("\n")
               }
             }
           })






#formula functions
#############------------------------------------------------------------------
#' @title sumFormula
#' @description Combine metabolite and adduct as a new sum formula.
#' If there are no enough element to remove, return NA.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param adduct The adduct of metabolite.
#' @return  A sum formula.
#' @export

setGeneric(name = "sumFormula",
           def = function(formula = "C9H11NO2",
                          adduct = "M-H2O+H"){

             if(is.na(formula)) return(NA)
             if(is.na(adduct)) return(formula)
             if(adduct == "M+" | adduct == "M-"){
               return(formula)
             }

             formula1 <- splitFormula(formula)
             adduct1 <- strsplit(x = adduct, split = "\\-|\\+")[[1]][-1]
             polymer <- as.numeric(gsub(pattern = "M", replacement = "",
                                        strsplit(x = adduct, split = "\\-|\\+")[[1]][1]))
             if (is.na(polymer)) polymer <- 1

             plusorminus <- strsplit(x = adduct, split = "")[[1]]
             plusorminus <- grep("\\+|\\-", plusorminus, value = TRUE)

             formula1$number <- formula1$number * polymer

             adduct1 <- mapply(function(x, y){
               temp <- splitFormula(x)
               temp$number <- temp$number * ifelse(y == "+", 1, -1)
               list(temp)
             },
             x = adduct1,
             y = plusorminus)

             adduct1 <- do.call(rbind, adduct1)

             formula <- rbind(formula1, adduct1)
             rownames(formula) <- NULL

             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               if(any(formula$number < 0)) {
                 return(NA)
               }else{
                 formula$number[formula$number==1] <- "W"
                 formula <- paste(paste(formula$element.name, formula$number, sep = ""), collapse = "")
                 formula <- strsplit(formula, split = "")[[1]]
                 formula[formula == "W"] <- ""
                 formula <- paste(formula, collapse = "")
                 return(formula)
               }
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })

               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })

               formula <- do.call(rbind, formula)
               formula <- formula[formula[,2] != 0,]
               colnames(formula) <- c("element.name", "number")
               if(any(formula$number < 0)) {return(NA)}else{
                 formula$number[formula$number==1] <- "W"
                 formula <- paste(paste(formula$element.name, formula$number, sep = ""), collapse = "")
                 formula <- strsplit(formula, split = "")[[1]]
                 formula[formula == "W"] <- ""
                 formula <- paste(formula, collapse = "")
                 return(formula)
               }
             }
           })

# setGeneric(name = "sumFormula",
#            def = function(formula = "C9H11NO2",
#                           adduct = "M+H"){
#              if(is.na(formula)) return(NA)
#              if(is.na(adduct)) return(formula)
#              if(adduct == "M+" | adduct == "M-"){
#                return(formula)
#              }
#
#              suppressMessages(expr = data(thermo, package = "CHNOSZ"))
#              adduct.species <- c("2ACN", "2Cl", "2H", "2H2O", "2K", "2Na",
#                                  "3ACN", "3H", "3H2O", "3Na", "ACN", "Br",
#                                  "CF3COOH", "CH3OH", "Cl", "H", "H2O", "Hac",
#                                  "HCOOH", "IsoProp", "K", "Na", "NaCOOH", "NH4",
#                                  "DMSO", "CH3CN", "3K", "CH3COO", "F", "NH3",
#                                  "HCOO")
#              adduct.number <- c(2, 2, 2, 2, 2, 2,
#                                 3, 3, 3, 3, 1, 1,
#                                 1, 1, 1, 1, 1, 1,
#                                 1, 1, 1, 1, 1, 1,
#                                 1, 1, 3, 1, 1, 1,
#                                 1)
#
#              adduct.element <- c("C4H6N", "Cl2" ,"H2", "H4O2", "K2", "Na2",
#                                  "C6H9N3", "H3", "H6O3", "Na3", "C2H3N", "Br",
#                                  "C2F3O2H", "CH4O", "Cl", "H", "H2O", "C2H4O2",
#                                  "CO2H2", "C3H8O", "K", "Na", "NaCO2H", "NH4",
#                                  "C2H6OS", "C2H3N", "K3", "C2H3O2", "F", "NH3",
#                                  "CHO2")
#              ###formula has the charge
#              if(adduct == "M+"){return(formula)}
#
#              adduct1 <- strsplit(x = adduct, split = ("\\+|\\-"))[[1]]
#
#              plusorminus <- strsplit(x = adduct, split = "")[[1]]
#              idx1 <- which(plusorminus == "+")
#              idx2 <- which(plusorminus == "-")
#              plusorminus1 <- NULL
#              plusorminus1[idx1] <- "+"
#              plusorminus1[idx2] <- "-"
#              plusorminus1 <- plusorminus1[!is.na(plusorminus1)]
#
#              ## ploymer
#              polymer <- as.numeric(gsub(pattern = "M", replacement = "", x =  adduct1[1]))
#              if (is.na(polymer)) polymer <- 1
#
#              ## adduct information
#              adduct2 <- adduct1[-1]
#
#              ## add polymer
#              formula1 <- CHNOSZ::makeup(unname(formula))
#              formula1 <- CHNOSZ::as.chemical.formula(formula1 * polymer)
#
#              for (i in 1:length(adduct2)) {
#                element <- adduct2[i]
#                if(!is.element(el = element, set = adduct.species)){
#                  stop(element, " is not in the adduct.species")
#                }
#                element <- adduct.element[which(element == adduct.species)]
#
#                ## add some thing
#                if (plusorminus1[i] == "+") {
#                  temp.formula1 <- splitFormula(formula1)
#                  temp.element <- splitFormula(element)
#                  if(nrow(temp.formula1) == 1 & nrow(temp.element) == 1 &
#                     temp.formula1[1,1] == temp.element[1,1]){
#                    formula1 <- CHNOSZ::makeup(formula1) + CHNOSZ::makeup(element)
#                  }else{
#                    temp <- CHNOSZ::i2A(c(formula1, element))
#                    formula1 <- apply(temp, 2, sum)
#                  }
#                  # formula1 <- temp[1,] + temp[2,]
#                  formula1 <- CHNOSZ::as.chemical.formula(formula1)
#                }else{
#                  temp.formula1 <- splitFormula(formula1)
#                  temp.element <- splitFormula(element)
#
#                  if(nrow(temp.formula1) == 1 & nrow(temp.element) == 1 &
#                     temp.formula1[1,1] == temp.element[1,1]){
#                    formula1 <- CHNOSZ::makeup(formula1) - CHNOSZ::makeup(element)
#                  }else{
#                    temp <- CHNOSZ::i2A(c(formula1, element))
#                    formula1 <- temp[1,] - temp[2,]
#                  }
#                  if(any(formula1 < 0)) return(NA)
#                  formula1 <- CHNOSZ::as.chemical.formula(formula1)
#                }
#              }
#              return(formula1)
#
#            })





#-----------------------------------------------------------------------------
#' @title splitFormula
#' @description Split a formula into element and number.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @return  A splited formula.
#' @export

setGeneric(name = "splitFormula",
           def = function(formula = "C9H11NO2"){
             temp.formula <- strsplit(formula, split = "")[[1]]

             number <- NULL
             for(i in 1:length(temp.formula)){
               if(length(grep("[0-9]{1}", temp.formula[i])) == 0){break}
               number[i] <- temp.formula[i]
             }

             if(!is.null(number)) {
               number <- as.numeric(paste(number, collapse = ""))
             }else{
               number <- 1
             }
             ##first select the Na, Cl and so on element
             idx1 <- gregexpr("[A-Z][a-z][0-9]*", formula)[[1]]
             len1 <- attributes(idx1)$match.length
             ##no double element
             if(idx1[1] == -1) {
               double.formula <- matrix(NA, ncol = 2)
               formula1 <- formula
             }else{
               double.letter.element <- NULL
               double.number <- NULL
               remove.idx <- NULL
               for (i in 1:length(idx1)) {
                 double.letter.element[i] <- substr(formula, idx1[i], idx1[i] + len1[i] - 1)
                 if(nchar(double.letter.element[i]) == 2){
                   double.number[i] <- 1
                 }else{
                   double.number[i] <- as.numeric(substr(double.letter.element[i], 3, nchar(double.letter.element[i])))
                 }
                 double.letter.element[i] <- substr(double.letter.element[i], 1, 2)
                 remove.idx <- c(remove.idx, idx1[i] : (idx1[i] + len1[i] - 1))
               }

               double.formula <- data.frame(double.letter.element,
                                            double.number, stringsAsFactors = FALSE)
               formula1 <- strsplit(formula, split = "")[[1]]
               formula1 <- formula1[-remove.idx]
               formula1 <- paste(formula1, collapse = "")
             }

             ## no one element
             if(formula1 == ""){
               one.formula <- matrix(NA, ncol = 2)
             }else{
               idx2 <- gregexpr("[A-Z][0-9]*", formula1)[[1]]
               len2 <- attributes(idx2)$match.length
               one.letter.element <- NULL
               one.number <- NULL
               for (i in 1:length(idx2)) {
                 one.letter.element[i] <- substr(formula1, idx2[i], idx2[i] + len2[i] - 1)
                 if(nchar(one.letter.element[i]) == 1){
                   one.number[i] <- 1
                 }else{
                   one.number[i] <- as.numeric(substr(one.letter.element[i], 2, nchar(one.letter.element[i])))
                 }
                 one.letter.element[i] <- substr(one.letter.element[i], 1, 1)
               }
               one.formula <- data.frame(one.letter.element, one.number,
                                         stringsAsFactors = FALSE)
             }

             colnames(double.formula) <- colnames(one.formula) <- c("element.name","number")
             formula <- rbind(double.formula, one.formula)
             formula <- formula[!apply(formula, 1, function(x) any(is.na(x))),]

             formula <- formula[order(formula$element.name),]
             formula$number <- formula$number * number
             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               return(formula)
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })

               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })

               formula <- do.call(rbind, formula)
               colnames(formula) <- c("element.name", "number")
               return(formula)
             }
           })

#-----------------------------------------------------------------------------
#' @title pasteElement
#' @description Paste formula and element.
#' Combine metabolite and adduct as a new sum formula.
#' If there are no enough element to remove, return NA.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param element The element.
#' @param mode Add or remove a module
#' @export
#' @return  A formula.


pasteElement <- function(formula = "C9H11NO2",
                         element = "H",
                         mode = c("plus", "minus")){

  mode <- match.arg(mode)
  formula <- splitFormula(formula = formula)
  element <- splitFormula(formula = element)


  ## mode = plus
  if(mode == "plus"){
    for (i in 1:nrow(element)){
      temp.name <- as.character(element[i,1])
      temp.number <- as.numeric(element[i,2])
      temp.idx <- match(temp.name, formula[,1])
      if(is.na(temp.idx)) {
        formula <- rbind(formula, element[i,])
      }else{
        formula[temp.idx, 2] <- formula[temp.idx, 2] + temp.number
      }
    }
  }else{
    for (i in 1:nrow(element)){
      temp.name <- as.character(element[i,1])
      temp.number <- as.numeric(element[i,2])
      temp.idx <- match(temp.name, formula[,1])
      if(is.na(temp.idx)) {
        # warning("Formula has no element in adduct!\n")
        return(NA)
      }else{
        formula[temp.idx,2] <- formula[temp.idx,2] - temp.number
        if(formula[temp.idx,2] < 0) {
          # warning("Formula has no enough element in adduct!\n")
          return(NA)
        }
      }
    }
  }

  ###return formula
  formula <- as.data.frame(formula)
  formula <- formula[formula[,2] != 0, , drop = FALSE]
  formula <- c(t(formula))
  formula <- gsub(pattern = " ", replacement = "", x = formula)
  formula <- formula[formula != "1"]
  formula <- paste(formula, collapse = "")
  return(formula)
}




#-----------------------------------------------------------------------------
#' @title checkElement
#' @description Check a formula can add one adduct or not.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param formula The formula of metabolite.
#' @param adduct The adduct.
#' @return  valid, return TRUE; invalid.
#' @export

setGeneric(name = "checkElement",
           def = function(formula = "C9H11NO2",
                          adduct = "M-H2O+H"){
             formula1 <- splitFormula(formula)
             adduct1 <- strsplit(x = adduct, split = "\\-|\\+")[[1]][-1]
             plusorminus <- strsplit(x = adduct, split = "")[[1]]
             plusorminus <- grep("\\+|\\-", plusorminus, value = TRUE)
             if(all(plusorminus == "+")) return(TRUE)

             adduct1 <- mapply(function(x, y){
               temp <- splitFormula(x)
               temp$number <- temp$number * ifelse(y == "+", 1, -1)
               list(temp)
             },
             x = adduct1,
             y = plusorminus)

             adduct1 <- do.call(rbind, adduct1)

             formula <- rbind(formula1, adduct1)
             rownames(formula) <- NULL

             unique.element <- unique(formula$element.name)
             if(length(unique.element) == nrow(formula)){
               if(any(formula$number < 0)) {return(FALSE)}else{return(TRUE)}
             }else{
               formula <- lapply(unique.element, function(x){
                 formula[formula$element.name == x,,drop = FALSE]
               })

               formula <- lapply(formula, function(x){
                 data.frame(unique(x$element.name), sum(x$number),
                            stringsAsFactors = FALSE)
               })

               formula <- do.call(rbind, formula)
               colnames(formula) <- c("element.name", "number")
               if(any(formula$number < 0)) {return(FALSE)}else{return(TRUE)}
             }
           })



# setGeneric(name = "checkElement",
#            def = function(formula = "C9H11NO2",
#                           adduct = "M-H2O+H"){
#              #
#              if(is.na(formula)) return(NA)
#              if(is.na(adduct)) return(formula)
#              if(adduct == "M+" | adduct == "M-"){
#                return(formula)
#              }
#              suppressMessages(expr = data(thermo, package = "CHNOSZ"))
#              adduct.species <- c("2ACN", "2Cl", "2H", "2H2O", "2K", "2Na",
#                                  "3ACN", "3H", "3H2O", "3Na", "ACN", "Br",
#                                  "CF3COOH", "CH3OH", "Cl", "H", "H2O", "Hac",
#                                  "HCOOH", "IsoProp", "K", "Na", "NaCOOH", "NH4",
#                                  "DMSO", "CH3CN", "3K", "CH3COO", "F", "NH3",
#                                  "HCOO")
#              adduct.number <- c(2, 2, 2, 2, 2, 2,
#                                 3, 3, 3, 3, 1, 1,
#                                 1, 1, 1, 1, 1, 1,
#                                 1, 1, 1, 1, 1, 1,
#                                 1, 1, 3, 1, 1, 1,
#                                 1)
#
#              adduct.element <- c("C4H6N", "Cl2" ,"H2", "H4O2", "K2", "Na2",
#                                  "C6H9N3", "H3", "H6O3", "Na3", "C2H3N", "Br",
#                                  "C2F3O2H", "CH4O", "Cl", "H", "H2O", "C2H4O2",
#                                  "CO2H2", "C3H8O", "K", "Na", "NaCO2H", "NH4",
#                                  "C2H6OS", "C2H3N", "K3", "C2H3O2", "F", "NH3",
#                                  "CHO2")
#
#              adduct1 <- strsplit(x = adduct, split = ("\\+|\\-"))[[1]]
#              plusorminus <- strsplit(x = adduct, split = "")[[1]]
#              idx1 <- which(plusorminus == "+")
#              idx2 <- which(plusorminus == "-")
#              plusorminus1 <- NULL
#              plusorminus1[idx1] <- "+"
#              plusorminus1[idx2] <- "-"
#              plusorminus1 <- plusorminus1[!is.na(plusorminus1)]
#
#              ## ploymer
#              polymer <- as.numeric(gsub(pattern = "M", replacement = "", x =  adduct1[1]))
#              if (is.na(polymer)) polymer <- 1
#
#              ## adduct information
#              adduct2 <- adduct1[-1]
#              #
#              #   ## add polymer
#              #   formula1 <- splitFormula(formula)
#              #   formula1[,2] <- formula1[,2] * polymer
#              #   formula1 <- c(t(formula1))
#              #   formula1 <- gsub(pattern = " ", replacement = "", x = formula1)
#              #   formula1 <- formula1[formula1 != "1"]
#              #   formula1 <- paste(formula1, collapse = "")
#
#
#              ## add polymer
#              formula1 <- CHNOSZ::makeup(unname(formula))
#              formula1 <- CHNOSZ::as.chemical.formula(formula1 * polymer)
#
#              for (i in 1:length(adduct2)) {
#                element <- adduct2[i]
#                element <- adduct.element[which(element == adduct.species)]
#                ## add some thing
#                if (plusorminus1[i] == "+") {
#                  temp.formula1 <- splitFormula(formula1)
#                  temp.element <- splitFormula(element)
#                  if(nrow(temp.formula1) == 1 & nrow(temp.element) == 1 &
#                     temp.formula1[1,1] == temp.element[1,1]){
#                    formula1 <- CHNOSZ::makeup(formula1) + CHNOSZ::makeup(element)
#                  }else{
#                    temp <- CHNOSZ::i2A(c(formula1, element))
#                    formula1 <- apply(temp, 2, sum)
#                  }
#                  # formula1 <- temp[1,] + temp[2,]
#                  formula1 <- CHNOSZ::as.chemical.formula(formula1)
#                }else{
#                  temp.formula1 <- splitFormula(formula1)
#                  temp.element <- splitFormula(element)
#
#                  if(nrow(temp.formula1) == 1 & nrow(temp.element) == 1 &
#                     temp.formula1[1,1] == temp.element[1,1]){
#                    formula1 <- CHNOSZ::makeup(formula1) - CHNOSZ::makeup(element)
#                  }else{
#                    temp <- CHNOSZ::i2A(c(formula1, element))
#                    formula1 <- temp[1,] - temp[2,]
#                  }
#                  if(any(formula1 < 0)) return(FALSE)
#                  formula1 <- CHNOSZ::as.chemical.formula(formula1)
#                }
#              }
#              return(TRUE)
#            })






# title changeTags: Chages tags (matrix or data.frame) to S4 class peakInfo
# description changeTages is used to change tags information (matrix or
# data frame) to peakInfo (S4 class).
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param tags The tags information.
# param ms2 The ms2.
# param adduct.table The adduct table.
# param weight.mz mz weight.
# param weight.rt RT weight
# param weight.dp Dop product weight
# param mz.tol The mz tolerance.
# param rt.tol The RT tolerance.
# param dp.tol The dp tolerance.
# param ... other parameters.
# return Tags2 data.
# export

setGeneric(name = "changeTags",
           def = function(tags = tags,
                          ms2 = ms2,
                          adduct.table = adduct.table,
                          weight.mz = 0.25,
                          weight.rt = 0.25,
                          weight.dp = 0.5,
                          mz.tol = 25,
                          rt.tol = 30,
                          dp.tol = 0.5,
                          ...) {
             #score formula
             ##calculate score for peak, using mz error, rt error
             ##and int error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)

             temp.x <- c(dp.tol, 1)
             temp.y <- c(0, 1)
             lm.reg.dp <- lm(temp.y~temp.x)


             #load ms2 information
             if(!is.null(ms2)){
               ms2 <- ms2
               ms2.name <-
                 unname(unlist(lapply(ms2, function(x) {x[[1]]["NAME",]})))
             }
             #change space and nan to NA
             tags[which(tags == "", arr.ind = TRUE)] <- NA
             # tags[which(is.nan(tags), arr.ind = TRUE)] <- NA
             #change tags to a list, a element is a peak.
             tags1 <- apply(tags, 1, list)
             #change tags1 to peakInfo S4 class

             tags2 <- lapply(tags1, function(x){
               temp <- new(Class = "peakInfo",
                           name = stringr::str_trim(as.character(x[[1]]["Peak.name"])),
                           mz = as.numeric(x[[1]]["mz"]),
                           rt = as.numeric(x[[1]]["rt"]),
                           ##add a empty ms2
                           ms2 = data.frame(),
                           annotation = list()
               )
               if(!is.na(x[[1]][["KEGG.ID"]])){
                 kegg.id <- x[[1]]["KEGG.ID"]
                 kegg.id <- strsplit(kegg.id, split = ";")[[1]]
                 kegg.id[kegg.id == "NA"] <- NA
                 if(any(!is.na(kegg.id))){
                   adduct <- strsplit(x[[1]]["adduct"], split = ";")[[1]]
                   ms2.sim <-
                     as.numeric(strsplit(x[[1]]["ms2.sim"], split = ";")[[1]])
                   rt.error <-
                     as.numeric(strsplit(x[[1]]["rt.error"], split = ";")[[1]])

                   ##for the identification whose rt error is NA, make it as
                   ##rt.tol
                   rt.error[is.na(rt.error)] <- rt.tol

                   mz.error <-
                     as.numeric(strsplit(x[[1]]["mz.error"], split = ";")[[1]])
                   formula <- strsplit(x[[1]]["Formula"], split = ";")[[1]]
                   temp.idx <- which(!is.na(kegg.id))

                   kegg.id <- kegg.id[temp.idx]
                   adduct <- adduct[temp.idx]
                   ms2.sim <- ms2.sim[temp.idx]
                   rt.error <- rt.error[temp.idx]
                   mz.error <- mz.error[temp.idx]
                   formula <- formula[temp.idx]


                   for(i in 1:length(kegg.id)){
                     score1 <- coefficients(lm.reg.mz)[1] +
                       coefficients(lm.reg.mz)[2] * mz.error[i]
                     score2 <- coefficients(lm.reg.rt)[1] +
                       coefficients(lm.reg.rt)[2] * rt.error[i]
                     score3 <- coefficients(lm.reg.dp)[1] +
                       coefficients(lm.reg.dp)[2] * ms2.sim[i]
                     score <-
                       score1*weight.mz + score2*weight.rt + score3 * weight.dp

                     temp@annotation[[i]] <-
                       list(type = "seed",
                            From = NA,
                            From.peak = NA,
                            to = stringr::str_trim(kegg.id[i]),
                            step = NA,
                            level = 1,
                            as.seed = FALSE,
                            as.seed.round = NA,
                            isotope = "[M]",
                            adduct = stringr::str_trim(adduct[i]),
                            charge =
                              as.numeric(adduct.table$charge[match(adduct[i], adduct.table$name)]),
                            Formula = stringr::str_trim(formula[i]),
                            mz.error = as.numeric(mz.error[i]),
                            rt.error = as.numeric(rt.error[i]),
                            int.error = NA,
                            ms2.sim = as.numeric(ms2.sim[i]),
                            score = unname(score)
                       )
                   }


                 }
               }
               temp
             })

             ##if has ms2, add ms2 information
             if(!is.null(ms2)){
               tags2 <- lapply(tags2, function(x){
                 x@ms2 <- data.frame(ms2[[match(x@name, ms2.name)]]$spec,
                                     stringsAsFactors = FALSE)
                 x
               })
             }
             # tags2 <-
             #   lapply(tags2, function(x) {filterPeak(object = x, score.thr = 0)})

             ##if score of annotation is less than 0, then change it to 0.001
             tags2 <- lapply(tags2, function(x){
               if(length(x@annotation) == 0) return(x)
                annotation <- x@annotation
                annotation <- lapply(annotation, function(y){
                if(y$score < 0) y$score <- 0.001
                y
                })
                x@annotation <- annotation
                x
             })

             tags2 <- tags2
           }
)


#####---------------------------------------------------------------
# title isotope2peakInfo
# description Add isotope result into tags2 data.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param isotopes.result The isotope result.
# param mz.tol The mz tolerance.
# param rt.tol The RT tolerance.
# param int.tol The intensity ratio tolerance.
# param weight.mz mz weight.
# param weight.rt RT weight.
# param weight.int Intesnity ratio weight.
# param peak.info The tags2 data.
# param peak.idx The index of seeds.
# param anno.idx The index of annotation
# param annotation.idx The index of peaks which are annotated.
# return Tags2 data.
# export

setGeneric(name = "isotope2peakInfo",
           def = function(isotopes.result,
                          mz.tol = 25,
                          rt.tol = 3,
                          int.tol = 500,
                          weight.mz = 0.45,
                          weight.rt = 0.45,
                          weight.int = 0.1,
                          peak.info,
                          peak.idx,
                          anno.idx,
                          annotation.idx = c(1:length(peak.info))){
             #
             if (is.null(isotopes.result)) return(peak.info)
             index <- annotation.idx[isotopes.result[,"peakIndex"]]
             seed <- peak.info[[peak.idx]]
             seed.name <- seed@name
             #score formula
             ##calculate score for peak, using mz error, rt error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)

             temp.x <- c(0, int.tol)
             temp.y <- c(1, 0)
             lm.reg.int <- lm(temp.y~temp.x)

             for(i in 1:length(index)){
               if(length(peak.info[[index[i]]]@annotation) == 0){
                 k = 1
                 peak.info[[index[i]]]@annotation <- list(
                   list(type = NA,
                        From = NA,
                        From.peak = NA,
                        to = NA,
                        step = NA,
                        level = NA,
                        as.seed = FALSE,
                        as.seed.round = NA,
                        isotope = NA,
                        adduct = NA,
                        charge = NA,
                        Formula = NA,
                        mz.error = NA,
                        rt.error = NA,
                        int.error = NA,
                        ms2.sim = NA,
                        score = NA)
                 )
               }else{
                 k = length(peak.info[[index[i]]]@annotation) + 1
                 peak.info[[index[i]]]@annotation <-
                   c(peak.info[[index[i]]]@annotation,
                     list(peak.info[[index[i]]]@annotation[[1]]))
               }
               peak.info[[index[i]]]@annotation[[k]]$type = "isotopeAnnotation"
               peak.info[[index[i]]]@annotation[[k]]$From =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$From.peak =
                 unname(seed.name)
               peak.info[[index[i]]]@annotation[[k]]$to =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$step = NA
               peak.info[[index[i]]]@annotation[[k]]$level =
                 unname(seed@annotation[[anno.idx]]$level + 1)
               peak.info[[index[i]]]@annotation[[k]]$as.seed = FALSE
               peak.info[[index[i]]]@annotation[[k]]$as.seed.round = NA
               peak.info[[index[i]]]@annotation[[k]]$isotope =
                 unname(isotopes.result$isotopes[i])
               peak.info[[index[i]]]@annotation[[k]]$adduct =
                 unname(seed@annotation[[anno.idx]]$adduct)
               peak.info[[index[i]]]@annotation[[k]]$charge =
                 unname(seed@annotation[[anno.idx]]$charge)
               peak.info[[index[i]]]@annotation[[k]]$Formula =
                 unname(seed@annotation[[anno.idx]]$Formula)
               peak.info[[index[i]]]@annotation[[k]]$mz.error =
                 unname(isotopes.result$mzError.ppm[i])
               peak.info[[index[i]]]@annotation[[k]]$rt.error =
                 unname(isotopes.result$rtError.s[i])
               peak.info[[index[i]]]@annotation[[k]]$int.error =
                 unname(isotopes.result$IntensityRatioError..[i])
               peak.info[[index[i]]]@annotation[[k]]$ms2.sim = NA

               score1 <- coefficients(lm.reg.mz)[1] +
                 coefficients(lm.reg.mz)[2] * isotopes.result$mzError.ppm[i]
               score2 <- coefficients(lm.reg.rt)[1] +
                 coefficients(lm.reg.rt)[2] * isotopes.result$rtError.s[i]
               score3 <- coefficients(lm.reg.int)[1] +
                 coefficients(lm.reg.int)[2] * isotopes.result$IntensityRatioError..[i]
               score <- score1*weight.mz + score2*weight.rt + score3*weight.int

               peak.info[[index[i]]]@annotation[[k]]$score = unname(score)
             }
             peak.info <- peak.info
           })





#####---------------------------------------------------------------
# title adduct2peakInfo
# description Add adduct result into tags2 data.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param adduct.result The adduct result.
# param mz.tol The mz tolerance.
# param rt.tol The RT tolerance.
# param weight.mz mz weight.
# param weight.rt RT weight.
# param peak.info The tags2 data.
# param peak.idx The index of seeds.
# param anno.idx The index of annotation
# param annotation.idx The index of peaks which are annotated.
# return Tags2 data.
# export
setGeneric(name = "adduct2peakInfo",
           def = function(adduct.result,
                          mz.tol = 25,
                          rt.tol = 3,
                          weight.mz = 0.8,
                          weight.rt = 0.2,
                          peak.info,
                          peak.idx,
                          anno.idx,
                          annotation.idx = c(1:length(peak.info))){
             if(is.null(adduct.result)) return(peak.info)
             index <- annotation.idx[adduct.result[,"peakIndex"]]
             seed <- peak.info[[peak.idx]]

             #------------------------------------------------------------------
             ##seed information
             seed.name <- seed@name
             seed.from.name <- seed@annotation[[anno.idx]]$From.peak
             seed.from <- seed@annotation[[anno.idx]]$From

             ##find the peaks which annotate the seed. For example, if seed
             #annotatte
             #if seed.from.name or seed.form is NA, this means this seed is the
             #raw seed for annotation
             if(!is.na(seed.from.name) | !is.na(seed.from)){
               #index1 is the index peak which have annotation
               index1 <-
                 index[which(showTags2(peak.info[index],
                                       slot = "annotation.len") > 0)]
               if(length(index1) == 0){
                 index <- index
               }else{
                 #remove.i is the index of annotated peaks which should not be
                 #annotated by this seed
                 remove.i <- NULL
                 for(i in 1:length(index1)){
                   name <- peak.info[[index1[i]]]@name
                   # id <- peak.info[[index1[i]]]@annotation[[1]]$to
                   id <-
                     unlist(lapply(peak.info[[210]]@annotation, function(x) x$to))
                   as.seed <-
                     unlist(lapply(peak.info[[210]]@annotation, function(x) x$as.seed))
                   id <- id[as.seed]
                   if(length(id) == 0) id <- "no"
                   if(seed.from.name == name & any(seed.from == id)){
                     remove.i[i] <- index1[i]
                   }
                 }
                 remove.i <- remove.i[!is.na(remove.i)]
                 index <- setdiff(index, remove.i)
               }
             }
             if(length(index) == 0) return(peak.info)
             #------------------------------------------------------------------

             #score formula
             ##calculate score for peak, using mz error, rt error
             ##and int error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)


             for(i in 1:length(index)){
               if(length(peak.info[[index[i]]]@annotation) == 0){
                 k = 1
                 peak.info[[index[i]]]@annotation <- list(
                   list(type = NA,
                        From = NA,
                        From.peak = NA,
                        to = NA,
                        step = NA,
                        level = NA,
                        as.seed = FALSE,
                        as.seed.round = NA,
                        isotope = NA,
                        adduct = NA,
                        charge = NA,
                        Formula = NA,
                        mz.error = NA,
                        rt.error = NA,
                        int.error = NA,
                        ms2.sim = NA,
                        score = NA)
                 )
               }else{
                 k = length(peak.info[[index[i]]]@annotation) + 1
                 peak.info[[index[i]]]@annotation <-
                   c(peak.info[[index[i]]]@annotation,
                     list(peak.info[[index[i]]]@annotation[[1]]))
               }
               peak.info[[index[i]]]@annotation[[k]]$type = "adductAnnotation"
               peak.info[[index[i]]]@annotation[[k]]$From =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$From.peak =
                 unname(seed.name)
               peak.info[[index[i]]]@annotation[[k]]$to =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$step = NA
               peak.info[[index[i]]]@annotation[[k]]$level =
                 unname(seed@annotation[[anno.idx]]$level + 1)
               peak.info[[index[i]]]@annotation[[k]]$as.seed = FALSE
               peak.info[[index[i]]]@annotation[[k]]$as.seed.round = NA
               peak.info[[index[i]]]@annotation[[k]]$isotope =
                 '[M]'
               peak.info[[index[i]]]@annotation[[k]]$adduct =
                 unname(adduct.result$adducts[i])
               peak.info[[index[i]]]@annotation[[k]]$charge =
                 unname(seed@annotation[[anno.idx]]$charge)
               peak.info[[index[i]]]@annotation[[k]]$Formula =
                 unname(seed@annotation[[anno.idx]]$Formula)
               peak.info[[index[i]]]@annotation[[k]]$mz.error =
                 unname(adduct.result$mzError.ppm[i])
               peak.info[[index[i]]]@annotation[[k]]$rt.error =
                 unname(adduct.result$rtError.s[i])
               peak.info[[index[i]]]@annotation[[k]]$int.error = NA
               peak.info[[index[i]]]@annotation[[k]]$ms2.sim = NA

               score1 <- coefficients(lm.reg.mz)[1] +
                 coefficients(lm.reg.mz)[2] * adduct.result$mzError.ppm[i]
               score2 <- coefficients(lm.reg.rt)[1] +
                 coefficients(lm.reg.rt)[2] * adduct.result$rtError.s[i]

               score <- score1*weight.mz + score2*weight.rt

               peak.info[[index[i]]]@annotation[[k]]$score = unname(score)
             }
             peak.info <- peak.info
           })

#####---------------------------------------------------------------
# title metabolite2peakInfo
# description Add annotation information form metAnnotation into
# peakInfo.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param metabolite.result metabolite.result from metAnnotation
# param mz.tol m/z tolerrance
# param rt.tol RT tolerrance
# param dp.tol Dot product tolerrance
# param weight.mz m/z weight
# param weight.rt RT weight
# param weight.dp Dot product weight
# param peak.info Peak info data
# param peak.idx The index of seed peaks.
# param anno.idx The index of seed annotations.
# param metabolite KEGG compound database
# param annotation.idx The index of peaks which are used for annotation
# return The tags2 data.
# export

setGeneric(name = "metabolite2peakInfo",
           def = function(metabolite.result,
                          mz.tol = 25,
                          rt.tol = 30,
                          dp.tol = 0.5,
                          weight.mz = 0.25,
                          weight.rt = 0.25,
                          weight.dp = 0.5,
                          peak.info,
                          peak.idx = 333,
                          anno.idx = 1,
                          metabolite = metabolite,
                          annotation.idx = c(1:length(peak.info))){
             #
             ##check parameters
             if(!is.list(peak.info)) stop("peak.info provided is invalid.")

             if(is.null(metabolite.result)) return(peak.info)

             #index is the index of peaks which are annotated by seed
             index <- annotation.idx[metabolite.result[,"peakIndex"]]
             #seed is the seed which annotate peaks
             seed <- peak.info[[peak.idx]]

             #------------------------------------------------------------------------------
             ##seed information
             seed.name <- seed@name
             seed.from.name <- seed@annotation[[anno.idx]]$From.peak
             seed.from <- seed@annotation[[anno.idx]]$From

             ##find the peaks which annotate the seed. For example, if seed
             #annotatte
             #if seed.from.name or seed.form is NA, this means this seed is the
             #raw seed for annotation
             if(!is.na(seed.from.name) & !is.na(seed.from)){
               #index1 is the index peak which have annotation
               index1 <-
                 index[which(showTags2(peak.info[index], slot = "annotation.len") > 0)]
               if(length(index1) == 0){
                 index <- index
               }else{
                 #remove.i is the index of annotated peaks which should not be
                 #annotated by this seed
                 remove.i <- NULL
                 for(i in 1:length(index1)){
                   name <- peak.info[[index1[i]]]@name
                   # id <- peak.info[[index1[i]]]@annotation[[1]]$to
                   id <-
                     unlist(lapply(peak.info[[index1[i]]]@annotation, function(x) x$to))
                   as.seed <-
                     unlist(lapply(peak.info[[index1[i]]]@annotation, function(x) x$as.seed))
                   id <- id[as.seed]
                   if(length(id) == 0) id <- "no"
                   if(seed.from.name == name & any(seed.from == id)){
                     remove.i[i] <- index1[i]
                   }
                 }
                 remove.i <- remove.i[!is.na(remove.i)]
                 remove.i <- unique(remove.i)
                 if(length(remove.i) != 0) {
                   metabolite.result <- metabolite.result[-which(remove.i == index),]
                   if(nrow(metabolite.result) != 0) {
                     index <- metabolite.result$peakIndex}else{
                       index <- NULL
                     }
                 }
               }
             }
             if(length(index) == 0) return(peak.info)
             #------------------------------------------------------------------
             #
             #score formula
             ##calculate score for peak, using mz error, rt error
             ##and int error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)


             temp.x <- c(dp.tol, 1)
             temp.y <- c(0, 1)
             lm.reg.dp <- lm(temp.y~temp.x)

             for(i in 1:length(index)){
               if(length(peak.info[[index[i]]]@annotation) == 0){
                 k = 1
                 peak.info[[index[i]]]@annotation <- list(
                   list(type = NA,
                        From = NA,
                        From.peak = NA,
                        to = NA,
                        step = NA,
                        level = NA,
                        as.seed = FALSE,
                        as.seed.round = NA,
                        isotope = NA,
                        adduct = NA,
                        charge = NA,
                        Formula = NA,
                        mz.error = NA,
                        rt.error = NA,
                        int.error = NA,
                        ms2.sim = NA,
                        score = NA))
               }else{
                 k = length(peak.info[[index[i]]]@annotation) + 1
                 peak.info[[index[i]]]@annotation <-
                   c(peak.info[[index[i]]]@annotation,
                     list(peak.info[[index[i]]]@annotation[[1]]))
               }

               ##remove the annotation which is from the peak that it annotate
               ##for example, if peak A is metabolte 1, and the it annotate peak B as
               ##metabolite 2, and then Peak B is used as seed and annotate peak A, so to
               #peak A, all the annotation from Peak B (metabolite 1) should be removed

               peak.info[[index[i]]]@annotation[[k]]$type = "metAnnotation"
               peak.info[[index[i]]]@annotation[[k]]$From =
                 unname(seed@annotation[[anno.idx]]$to)
               peak.info[[index[i]]]@annotation[[k]]$From.peak =
                 unname(seed.name)
               peak.info[[index[i]]]@annotation[[k]]$to =
                 unname(metabolite.result$peakID[i])
               peak.info[[index[i]]]@annotation[[k]]$step =
                 unname(metabolite.result$step[i])
               peak.info[[index[i]]]@annotation[[k]]$level =
                 unname(seed@annotation[[anno.idx]]$level + 1)
               peak.info[[index[i]]]@annotation[[k]]$as.seed = FALSE
               peak.info[[index[i]]]@annotation[[k]]$as.seed.round = NA
               peak.info[[index[i]]]@annotation[[k]]$isotope = '[M]'
               peak.info[[index[i]]]@annotation[[k]]$adduct =
                 unname(metabolite.result$adduct[i])
               peak.info[[index[i]]]@annotation[[k]]$charge =
                 unname(metabolite.result$charge[i])
               peak.info[[index[i]]]@annotation[[k]]$Formula =
                 unname(metabolite$Formula[match(metabolite.result$peakID[i], metabolite$ID)])
               peak.info[[index[i]]]@annotation[[k]]$mz.error =
                 unname(metabolite.result$mzError.ppm[i])
               peak.info[[index[i]]]@annotation[[k]]$rt.error =
                 unname(metabolite.result$rtError[i])
               peak.info[[index[i]]]@annotation[[k]]$int.error = NA
               peak.info[[index[i]]]@annotation[[k]]$ms2.sim =
                 unname(metabolite.result$dotProduct[i])



               score1 <- coefficients(lm.reg.mz)[1] +
                 coefficients(lm.reg.mz)[2] * metabolite.result$mzError.ppm[i]
               score2 <- coefficients(lm.reg.rt)[1] +
                 coefficients(lm.reg.rt)[2] * metabolite.result$rtError[i]

               score3 <- coefficients(lm.reg.dp)[1] +
                 coefficients(lm.reg.dp)[2] * metabolite.result$dotProduct[i]

               score <- score1*weight.mz + score2*weight.rt + score3 * weight.dp
               peak.info[[index[i]]]@annotation[[k]]$score = unname(score)
             }
             peak.info <- peak.info
           })




#####---------------------------------------------------------------
# title kegg2peakInfo
# description Add annotation information form kegg matching result (mz and rt)
# author Xiaotao Shen
#  \email{shenxt@@sioc.ac.cn}
# param kegg.result KEGG.result from kegg matching
# param mz.tol m/z tolerance
# param rt.tol RT tolerance
# param dp.tol Dot product tolerance
# param weight.mz m/z weight
# param weight.rt RT weight
# param weight.dp Dot product weight
# param peak.info Peak info data
# param kegg.compound KEGG compound database
# param annotation.idx The index of peaks which are used for annotation
# return The tags2 data.
# export



setGeneric(name = "kegg2peakInfo",
           def = function(kegg.result,
                          mz.tol = 25,
                          rt.tol = 30,
                          dp.tol = 0.5,
                          weight.mz = 0.5,
                          weight.rt = 0.5,
                          weight.dp = 0,
                          peak.info,
                          kegg.compound,
                          annotation.idx = c(1:length(peak.info))){
             #
             ##check parameters
             if(!is.list(peak.info)) stop("peak.info provided is invalid.")

             if(is.null(kegg.result)) return(peak.info)

             #index is the index of peaks which are annotated by seed
             index <- annotation.idx[kegg.result[,"peakIndex"]]
             #score formula
             ##calculate score for peak, using mz error, rt error
             ##and int error
             temp.x <- c(0, mz.tol)
             temp.y <- c(1, 0)
             lm.reg.mz <- lm(temp.y~temp.x)

             temp.x <- c(0, rt.tol)
             temp.y <- c(1, 0)
             lm.reg.rt <- lm(temp.y~temp.x)

             temp.x <- c(dp.tol, 1)
             temp.y <- c(0, 1)
             lm.reg.dp <- lm(temp.y~temp.x)

             for(i in 1:length(index)){
               if(length(peak.info[[index[i]]]@annotation) == 0){
                 k = 1
                 peak.info[[index[i]]]@annotation <- list(
                   list(type = NA,
                        From = NA,
                        From.peak = NA,
                        to = NA,
                        step = NA,
                        level = NA,
                        as.seed = FALSE,
                        as.seed.round = NA,
                        isotope = NA,
                        adduct = NA,
                        charge = NA,
                        Formula = NA,
                        mz.error = NA,
                        rt.error = NA,
                        int.error = NA,
                        ms2.sim = NA,
                        score = NA))
               }else{
                 k = length(peak.info[[index[i]]]@annotation) + 1
                 peak.info[[index[i]]]@annotation <-
                   c(peak.info[[index[i]]]@annotation,
                     list(peak.info[[index[i]]]@annotation[[1]]))
               }

               peak.info[[index[i]]]@annotation[[k]]$type = "keggMatching"
               peak.info[[index[i]]]@annotation[[k]]$From = NA
               peak.info[[index[i]]]@annotation[[k]]$From.peak = NA
               peak.info[[index[i]]]@annotation[[k]]$to =
                 unname(kegg.result$peakID[i])
               peak.info[[index[i]]]@annotation[[k]]$step = NA
               peak.info[[index[i]]]@annotation[[k]]$level = NA
               peak.info[[index[i]]]@annotation[[k]]$as.seed = FALSE
               peak.info[[index[i]]]@annotation[[k]]$as.seed.round = NA
               peak.info[[index[i]]]@annotation[[k]]$isotope = '[M]'
               peak.info[[index[i]]]@annotation[[k]]$adduct =
                 unname(kegg.result$adduct[i])
               peak.info[[index[i]]]@annotation[[k]]$charge =
                 unname(kegg.result$charge[i])
               peak.info[[index[i]]]@annotation[[k]]$Formula =
                 unname(kegg.compound$Formula[match(kegg.result$peakID[i], kegg.compound$ID)])
               peak.info[[index[i]]]@annotation[[k]]$mz.error =
                 unname(kegg.result$mzError.ppm[i])
               peak.info[[index[i]]]@annotation[[k]]$rt.error =
                 unname(kegg.result$rtError[i])
               peak.info[[index[i]]]@annotation[[k]]$int.error = NA
               peak.info[[index[i]]]@annotation[[k]]$ms2.sim = NA



               score1 <- coefficients(lm.reg.mz)[1] +
                 coefficients(lm.reg.mz)[2] * kegg.result$mzError.ppm[i]
               score2 <- coefficients(lm.reg.rt)[1] +
                 coefficients(lm.reg.rt)[2] * kegg.result$rtError[i]

               score3 <- 0

               score <- score1*weight.mz + score2*weight.rt + score3 * weight.dp
               peak.info[[index[i]]]@annotation[[k]]$score = unname(score)
             }
             peak.info <- peak.info
           })




###0---------------------------------------------------------
# removeAnnotation <- function(peak.info,
#                              wrong.id,
#                              wrong.name,
#                              from.round = 0){
#   #
#   #
#   #
#   for(i in 1:length(wrong.id)){
#     index <- findAnnotation(peak.info = peak.info,
#                             name = wrong.name[i],
#                             id = wrong.id[i],
#                             from.round = from.round)
#     if(!is.null(index) & length(index) != 0){
#       for (j in 1:length(index)){
#         temp.idx <- index[[j]]
#         for(k in 1:length(temp.idx)){
#           peak.info[[as.numeric(names(temp.idx[k]))]]@annotation <-
#             peak.info[[as.numeric(names(temp.idx[k]))]]@annotation[-temp.idx[[k]]]
#         }
#
#       }
#
#     }else{
#       next
#     }
#   }
#   return(peak.info)
# }


#
# findAnnotation: Find the peak indexs and annotation indexs of
# round x seed i (level is x + 1) and all the annotation peaks from the seed.
# description findAnnotation
# param peak.info tags2 class which contains peakInfo class.
# param name The name of the seed.
# param id The id of the seed.
# param from.round The round of the peak used as seed.
# author Xiaotao Shen
# export
#
# setGeneric(
#   name = "findAnnotation",
#   def = function(peak.info,
#                  name = "M102T561",
#                  id = "C00576",
#                  from.round = 1) {
#     #
#     from <- from.round
#     ##find the biggest round of peak.info
#     level <-  lapply(peak.info, function(x) {
#       annotation <- x@annotation
#       if (length(annotation) != 0) {
#         level <- unlist(lapply(annotation, function(x) {
#           x$level
#         }))
#         level
#       } else{
#         NA
#       }
#     })
#     level <- unlist(level)
#     level <- level[!is.na(level)]
#     max.level <- max(level)
#     max.round = max.level - 1
#     #
#     #     if (from.round > max.round) {
#     #       stop("No round ", from.round, " annotation")
#     #     }
#
#
#     return.list <- list()
#     ###find the seed's peak index and annotation index
#     peak.name <- showTags2(peak.info, slot = "name")
#     peak.index <- which(peak.name == name)
#     annotation <- peak.info[[peak.index]]@annotation
#     annotation.index <-
#       which(unlist(lapply(annotation, function(x) {x$to})) == id)[1]
#
#     annotation.index <- list(annotation.index)
#     names(annotation.index) <- peak.index
#     return.list[[1]] <- annotation.index
#     names(return.list)[1] <- paste("round", from.round, sep = "")
#
#     if(from.round == max.round){
#       # return(return.list)
#       return(NULL)
#     }
#
#     from.id <- id
#     from.name <- name
#     for (round.number in (from.round+1):max.round) {
#       # cat(round.number);cat(" ")
#       index <- which(showTags2(tags2 = peak.info, slot = "annotation.len")  > 0)
#       result <- lapply(peak.info[index], function(x) {
#         annotation <- x@annotation
#         anno.from <- unlist(lapply(annotation, function(x) {x$From}))
#         anno.from.peak <- unlist(lapply(annotation, function(x) {x$From.peak}))
#         anno.round <- unlist(lapply(annotation, function(x) {x$level})) - 1
#         anno.to <- unlist(lapply(annotation, function(x) {x$to}))
#         temp.idx <-
#           which(anno.from %in% from.id & anno.from.peak %in% from.name & anno.round == round.number)
#         if(length(temp.idx) == 0){
#           NA
#         }else{
#           annotation.index <- temp.idx
#           annotation.id <- anno.to[temp.idx]
#           idx_id <- c(list(annotation.index), list(annotation.id))
#           return(idx_id)
#         }
#
#       })
#       names(result) <- index
#       result <- result[!is.na(result)]
#       if(length(result) == 0){
#         break
#       }else{
#         idx <- lapply(result, function(x) {x[[1]]})
#         names(idx) <- names(result)
#         idx <- list(idx)
#         return.list <- c(return.list, idx)
#         names(return.list)[length(return.list)] <-
#           paste("round", round.number, sep = "")
#
#         from.id <- unlist(lapply(result, function(x) {x[[2]]}))
#         from.name <- showTags2(tags2 = peak.info, slot = "name")[as.numeric(names(from.id))]
#       }
#     }
#     return.list <- return.list[-1]
#     return(return.list)
#   }
# )


###----------------------------------------------------------------------------
#####---------------------------------------------------------------
#' @title showTags2
#' @description Show information of tags2 data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param tags2 The tags2 data.
#' @param slot The inforamtion.
#' @return The information of tags2 data.
#' @export

setGeneric(name = "showTags2",
           def = function(tags2,
                          slot = c("name", "mz", "rt",
                                   "annotation.len", "annotation.id",
                                   "annotation.from", "annotation.from.peak",
                                   "annotation.score", "annotation.level",
                                   "annotation.step", "seed.number",
                                   "annotation.type","annotation.formula",
                                   "as.seed")){
             slot <- match.arg(slot)
             if(slot == "name"){
               name <- unlist(lapply(tags2, function(x) {x@name}))
               return(name)
             }

             if(slot == "mz"){
               mz <- unlist(lapply(tags2, function(x) {x@mz}))
               return(mz)
             }

             if(slot == "rt"){
               rt <- unlist(lapply(tags2, function(x) {x@rt}))
               return(rt)
             }


             annotation.len <- unlist(lapply(tags2, function(x) {length(x@annotation)}))

             if(slot == "annotation.len"){
               return(annotation.len)
             }

             index <- which(annotation.len > 0)
             if(length(index) == 0) {cat('No peaks have annotation.\n')
             }

             if(slot == "annotation.id"){
               id <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$to})
               })
               names(id) <- index
               return(id)
             }

             if(slot == "annotation.type"){
               type <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$type})
               })
               names(type) <- index
               return(type)
             }


             if(slot == "annotation.from"){
               from <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$From})
               })
               names(from) <- index
               return(from)
             }

             if(slot == "annotation.formula"){
               formula <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$Formula})
               })
               names(formula) <- index
               return(formula)
             }


             if(slot == "annotation.from.peak"){
               from.peak <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$From.peak})
               })
               names(from.peak) <- index
               return(from.peak)
             }

             if(slot == "annotation.score"){
               score <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$score})
               })
               names(score) <- index
               return(score)
             }

             if(slot == "annotation.level"){
               level <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$level})
               })
               names(level) <- index
               return(level)
             }

             if(slot == "annotation.step"){
               step <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 lapply(annotation, function(x) {x$step})
               })
               names(step) <- index
               return(step)
             }

             if(slot == "seed.number"){
               seed.number <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 sum(unlist(lapply(annotation, function(x) {x$as.seed})))
               })
               seed.number <- unlist(seed.number)
               names(seed.number) <- index
               return(seed.number)
             }


             if(slot == "as.seed"){
               as.seed <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 unlist(lapply(annotation, function(x) {x$as.seed}))
               })


               seed.round <- lapply(tags2[index], function(x) {
                 annotation <- x@annotation
                 unlist(lapply(annotation, function(x) {x$as.seed.round}))
               })

               names(as.seed) <- index
               names(seed.round) <- index
               result <- list(as.seed, seed.round)
               names(result) <- c("seed", "round")
               return(result)
             }


           })



#####---------------------------------------------------------------
# title trans2Matrix
# description Transform tags2 data to matrix.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param tags2 The tags2 data.
# param base Sort according to annotation or peak.
# return The matrix of tags2 data.
# export

setGeneric(name = "trans2Matrix",
           def = function(tags2,
                          base = c("annotation", "peak")) {

             base = match.arg(base)
             result <- lapply(tags2, function(x) {
               name <- x@name
               mz <- x@mz
               rt <- x@rt
               annotation <- x@annotation
               if(length(annotation) == 0){
                 NULL
               }else{
                 # annotation <- lapply(annotation, function(x) {x <- x[names(x) != "addcut"]})
                 annotation <- do.call(rbind, lapply(annotation, unlist))
                 if(any(colnames(annotation) == "int.ratio.error") & any(colnames(annotation) == "int.error")){
                   annotation[,"int.error"] <- annotation[,"int.ratio.error"]
                   annotation <- annotation[,c(1:17)]
                   colnames(annotation)[15] <- "int.error"
                 }
                 if(any(colnames(annotation) == "int.ratio.error") & all(colnames(annotation) != "int.error")){
                   colnames(annotation)[15] <- "int.error"
                 }

                 data.frame(name, mz, rt, annotation,
                            stringsAsFactors = FALSE)
               }

             })
             result <- do.call(what = "rbind", args = result)
             if(base == "annotation"){
               result <- result[order(result$to),]
             }else{
               result <- result[order(result$name),]
             }
             return(result)
           }
)


#####---------------------------------------------------------------
#' @title tags2Match
#' @description Find annotation in tags2 data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param query The ID of annotation.
#' @param tags2 The tags2 data.
#' @param is.seed Is seed or not.
#' @return The index of annotation.
#' @export

setGeneric(name = "tags2Match",
           def = function(query = c("M102T561", "C00576"),
                          tags2,
                          is.seed = FALSE){
             peak.name <- query[1]
             peak.id <- query[2]
             idx1 <- match(peak.name, showTags2(tags2, slot = "name"))
             id <- unlist(showTags2(tags2[idx1], slot = "annotation.id")[[1]])
             if(is.seed){
               idx2 <-
                 which(showTags2(tags2[idx1], slot = "as.seed")[[1]][[1]] & peak.id == id)

             }else{
               idx2 <- grep(peak.id, id)
             }
             names(idx1) <- "peak.index"
             names(idx2) <- rep("annotation.index", length(idx2))
             return(c(idx1, idx2))
           })



#####---------------------------------------------------------------
#' @title compareAnnotation
#' @description compareAnnotation.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param tags2_1 The tags2 data.
#' @param tags2_2 The tags2 data.
#' @return The comparesion result.
#' @export


setGeneric(name = "compareAnnotation",
           def = function(tags2_1, tags2_2){
             peak.name1 <- showTags2(tags2_1, slot = "name")
             peak.name2 <- showTags2(tags2_2, slot = "name")
             if(any(peak.name1 != peak.name2)) {stop("No same peaks")}
             annotation1 <- lapply(tags2_1, function(x) {x@annotation})
             annotation2 <- lapply(tags2_2, function(x) {x@annotation})

             annotation1 <- lapply(annotation1, function(x) {if(length(x) > 0) {
               temp1 <- lapply(x, function(y) {unlist(y[c(1,4,9,10,12,17)])})
               temp2 <- do.call(rbind, temp1)
               temp3 <- apply(temp2, 2, function(x) {paste(x, collapse = ";")})
               temp3
             }else{
               temp3 <- rep(NA, 6)
               names(temp3) <- c("type", "to", "isotope","adduct","Formula", "score")
               temp3
             }})


             annotation2 <- lapply(annotation2, function(x) {if(length(x) > 0) {
               temp1 <- lapply(x, function(y) {unlist(y[c(1,4,9,10,12,17)])})
               temp2 <- do.call(rbind, temp1)
               temp3 <- apply(temp2, 2, function(x) {paste(x, collapse = ";")})
               temp3
             }else{
               temp3 <- rep(NA, 6)
               names(temp3) <- c("type", "to", "isotope","adduct","Formula", "score")
               temp3
             }})

             annotation1 <- do.call(what = rbind, args = annotation1)
             annotation2 <- do.call(what = rbind, args = annotation2)

             annotation1 <- cbind(peak.name1, annotation1)
             annotation2 <- cbind(peak.name2, annotation2)

             colnames(annotation1)[1] <- colnames(annotation2)[2] <- "peak.name"

             result <- data.frame(annotation1[,1], annotation1[,2], annotation2[,2],
                                  annotation1[,3], annotation2[,3],
                                  annotation1[,4], annotation2[,4],
                                  annotation1[,5], annotation2[,5],
                                  annotation1[,6], annotation2[,6],
                                  annotation1[,7], annotation2[,7],
                                  stringsAsFactors = FALSE)
             colnames(result) <- c("peak.name", "Type1", "Type2", "ID1", "ID2",
                                   "Isotope1","Isotope2",
                                   "Adduct1", "Adduct2","Formula1", "Formula2",
                                   "Score1", 'Score2')

             id1 <- result$ID1
             id2 <- result$ID2
             formula1 <- result$Formula1
             formula2 <- result$Formula2

             Note <- mapply(function(x1, y1, x2, y2) {
               if(any(is.na(c(x1,y1)))){
                 note <- "No annotation"
                 note
               }else{
                 x1 <- strsplit(x1, split = ";")[[1]]
                 y1 <- strsplit(y1, split = ";")[[1]]
                 x2 <- strsplit(x2, split = ";")[[1]]
                 y2 <- strsplit(y2, split = ";")[[1]]
                 if(length(intersect(x1, y1)) == 0 & length(intersect(x2,y2)) == 0) {
                   note <- "Wrong annotation"
                 }

                 if(length(intersect(x1, y1)) == 0 & length(intersect(x2,y2)) != 0) {
                   note <- "Isomer annotation"
                 }

                 if(length(intersect(x1, y1)) != 0) {
                   note <- "Right annotation"
                 }
                 note

               }
             },

             x1 = id1,
             y1 = id2,
             x2 = formula1,
             y2 = formula2)

             result <- data.frame(result, Note, stringsAsFactors = FALSE)
             cat(names(table(Note)))
             cat("\n")
             cat(table(Note))
             return(result)

           })



#####---------------------------------------------------------------
#' @title tags2Result
#' @description Transform tags2 data to matrix like ms2Annotation.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param tags2 The tags2 data.
#' @param score.cutoff The score cutoff.
#' @param candidate.num The number of candidates.
#' @return The matrix of tags2 data.
#' @export


setGeneric(name = "tags2Result",
           def = function(tags2,
                          score.cutoff = 0,
                          candidate.num = 3000){
             cat("Transform tags2 to matrix:\n")
             pbapply::pboptions(type = "timer", style = 1)
             result <- pbapply::pblapply(tags2, function(x) {
               name <- x@name
               mz <- x@mz
               rt <- x@rt
               annotation <- x@annotation
               if(length(annotation) == 0){
                 # c(name, mz, rt, rep(NA, 17))
                 result <- data.frame(name, mz, rt,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,##17 NAs
                            stringsAsFactors = FALSE)
                 colnames(result) <- c("name", "mz", "rt", "type", "From", "From.peak", "to", "step",
                 "level", "as.seed", "as.seed.round", "isotope", "adduct", "charge", "Formula",
                 "mz.error", "rt.error", "int.error", "ms2.sim", "score")
                 result
               }else{
                 type <- unlist(lapply(annotation, function(x) x$type))
                 From <- unlist(lapply(annotation, function(x) x$From))
                 From.peak <- unlist(lapply(annotation, function(x) x$From.peak))
                 to <- unlist(lapply(annotation, function(x) x$to))
                 step <- unlist(lapply(annotation, function(x) x$step))
                 level <- unlist(lapply(annotation, function(x) x$level))
                 as.seed <- unlist(lapply(annotation, function(x) x$as.seed))
                 as.seed.round <- unlist(lapply(annotation, function(x) x$as.seed.round))
                 isotope <- unlist(lapply(annotation, function(x) x$isotope))
                 adduct <- unlist(lapply(annotation, function(x) x$adduct))
                 charge <- unlist(lapply(annotation, function(x) x$charge))
                 Formula <- unlist(lapply(annotation, function(x) x$Formula))
                 mz.error <- unlist(lapply(annotation, function(x) x$mz.error))
                 rt.error <- unlist(lapply(annotation, function(x) x$rt.error))
                 int.error <- unlist(lapply(annotation, function(x) x$int.error))
                 if(is.null(int.error)) {int.error <- unlist(lapply(annotation, function(x) x$int.ratio.error))}
                 ms2.sim <- unlist(lapply(annotation, function(x) x$ms2.sim))
                 score <- unlist(lapply(annotation, function(x) x$score))
                 result <- data.frame(name, mz, rt, type, From, From.peak, to, step,
                                      level, as.seed,
                                      as.seed.round, isotope, adduct, charge, Formula,
                                      mz.error, rt.error, int.error, ms2.sim,
                                      score, stringsAsFactors = FALSE)
                 result <- result[!duplicated(result),]
                 result <- result[order(result$score, decreasing = TRUE),]##order result according to score
                 result <- result[result$score >= score.cutoff,]
                 if(nrow(result) > candidate.num) result <- result[1:candidate.num,]
                 if(nrow(result) == 0){return(c(name, mz, rt, rep(NA, 17)))}
                 type <- paste(result$type, collapse = ";")
                 From <- paste(result$From, collapse = ";")
                 From.peak <- paste(result$From.peak, collapse = ";")
                 to <- paste(result$to, collapse = ";")
                 step <- paste(result$step, collapse = ";")
                 level <- paste(result$level, collapse = ";")
                 as.seed <- paste(result$ as.seed, collapse = ";")
                 as.seed.round <- paste(result$as.seed.round, collapse = ";")
                 as.seed.round <- paste(result$as.seed.round, collapse = ";")
                 isotope <- paste(result$isotope, collapse = ";")
                 adduct <- paste(result$adduct, collapse = ";")
                 charge <- paste(result$charge, collapse = ";")
                 Formula <- paste(result$Formula, collapse = ";")
                 mz.error <- paste(result$mz.error, collapse = ";")
                 rt.error <- paste(result$rt.error, collapse = ";")
                 int.error <- paste(result$int.error, collapse = ";")
                 ms2.sim <- paste(result$ms2.sim, collapse = ";")
                 score <- paste(result$score, collapse = ";")
                 result <- data.frame(name, mz, rt, type, From, From.peak, to, step,
                                      level, as.seed,
                                      as.seed.round, isotope, adduct, charge, Formula,
                                      mz.error, rt.error, int.error, ms2.sim,
                                      score, stringsAsFactors = FALSE)
                 rm(list = c("type", "From", "From.peak", "to", "step", "level",
                             "as.seed", "as.seed.round", "isotope", "adduct",
                             "charge", "Formula", "mz.error", "rt.error",
                             "int.error", "ms2.sim", "score"))
                 return(result)
               }
             })

             #########################fix bugs
             # for(i in 1:length(tags2)){
             #  cat(i);cat(" ")
             #   x <- tags2[[i]]
             #   name <- x@name
             #   mz <- x@mz
             #   rt <- x@rt
             #   annotation <- x@annotation
             #   if(length(annotation) == 0){
             #     c(name, mz, rt, rep(NA, 17))
             #   }else{
             #     type <- unlist(lapply(annotation, function(x) x$type))
             #     From <- unlist(lapply(annotation, function(x) x$From))
             #     From.peak <- unlist(lapply(annotation, function(x) x$From.peak))
             #     to <- unlist(lapply(annotation, function(x) x$to))
             #     step <- unlist(lapply(annotation, function(x) x$step))
             #     level <- unlist(lapply(annotation, function(x) x$level))
             #     as.seed <- unlist(lapply(annotation, function(x) x$as.seed))
             #     as.seed.round <- unlist(lapply(annotation, function(x) x$as.seed.round))
             #     isotope <- unlist(lapply(annotation, function(x) x$isotope))
             #     adduct <- unlist(lapply(annotation, function(x) x$adduct))
             #     charge <- unlist(lapply(annotation, function(x) x$charge))
             #     Formula <- unlist(lapply(annotation, function(x) x$Formula))
             #     mz.error <- unlist(lapply(annotation, function(x) x$mz.error))
             #     rt.error <- unlist(lapply(annotation, function(x) x$rt.error))
             #     int.error <- unlist(lapply(annotation, function(x) x$int.error))
             #     if(is.null(int.error)) {int.error <- unlist(lapply(annotation, function(x) x$int.ratio.error))}
             #     ms2.sim <- unlist(lapply(annotation, function(x) x$ms2.sim))
             #     score <- unlist(lapply(annotation, function(x) x$score))
             #     result <- data.frame(name, mz, rt, type, From, From.peak, to, step,
             #                          level, as.seed,
             #                          as.seed.round, isotope, adduct, charge, Formula,
             #                          mz.error, rt.error, int.error, ms2.sim,
             #                          score, stringsAsFactors = FALSE)
             #     result <- result[!duplicated(result),]
             #     result <- result[result$score > score.cutoff,]
             #     if(nrow(result) == 0){return(c(peak.name, peak.mz, peak.rt, rep(NA, 17)))}
             #     type <- paste(result$type, collapse = ";")
             #     From <- paste(result$From, collapse = ";")
             #     From.peak <- paste(result$From.peak, collapse = ";")
             #     to <- paste(result$to, collapse = ";")
             #     step <- paste(result$step, collapse = ";")
             #     level <- paste(result$level, collapse = ";")
             #     as.seed <- paste(result$ as.seed, collapse = ";")
             #     as.seed.round <- paste(result$as.seed.round, collapse = ";")
             #     as.seed.round <- paste(result$as.seed.round, collapse = ";")
             #     isotope <- paste(result$isotope, collapse = ";")
             #     adduct <- paste(result$adduct, collapse = ";")
             #     charge <- paste(result$charge, collapse = ";")
             #     Formula <- paste(result$Formula, collapse = ";")
             #     mz.error <- paste(result$mz.error, collapse = ";")
             #     rt.error <- paste(result$rt.error, collapse = ";")
             #     int.error <- paste(result$int.error, collapse = ";")
             #     ms2.sim <- paste(result$ms2.sim, collapse = ";")
             #     score <- paste(result$score, collapse = ";")
             #     result <- data.frame(name, mz, rt, type, From, From.peak, to, step,
             #                          level, as.seed,
             #                          as.seed.round, isotope, adduct, charge, Formula,
             #                          mz.error, rt.error, int.error, ms2.sim,
             #                          score, stringsAsFactors = FALSE)
             #     rm("result")
             #     rm(list = c("type", "From", "From.peak", "to", "step", "level",
             #                 "as.seed", "as.seed.round", "isotope", "adduct",
             #                 "charge", "Formula", "mz.error", "rt.error",
             #                 "int.error", "ms2.sim", "score"))
             #     # return(result)
             #   }
             # }

             ##-----------------------------------------------------------

             rm(list = "tags2")
             result <- do.call(rbind, result)
             result <- as.data.frame(result)
             result <- result
           })


#####---------------------------------------------------------------
# title tags2Result1
# description Transform tags2 data to matrix like ms2Annotation.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param tags2 The tags2 data.
# param score.cutoff The score cutoff.
# param threads The number of threads.
# return The matrix of tags2 data.
# export

setGeneric(name = "tags2Result1",
           def = function(tags2,
                          score.cutoff = 0,
                          threads = parallel::detectCores()-3){
             cat("tags2 to matrix:\n")
             # pbapply::pboptions(type = "timer", style = 1)
             temp.fun <- function(x, score.cutoff){
               name <- x@name
               mz <- x@mz
               rt <- x@rt
               annotation <- x@annotation
               if(length(annotation) == 0){
                 c(name, mz, rt, rep(NA, 17))
               }else{
                 type <- unlist(lapply(annotation, function(x) x$type))
                 From <- unlist(lapply(annotation, function(x) x$From))
                 From.peak <- unlist(lapply(annotation, function(x) x$From.peak))
                 to <- unlist(lapply(annotation, function(x) x$to))
                 step <- unlist(lapply(annotation, function(x) x$step))
                 level <- unlist(lapply(annotation, function(x) x$level))
                 as.seed <- unlist(lapply(annotation, function(x) x$as.seed))
                 as.seed.round <- unlist(lapply(annotation, function(x) x$as.seed.round))
                 isotope <- unlist(lapply(annotation, function(x) x$isotope))
                 adduct <- unlist(lapply(annotation, function(x) x$adduct))
                 charge <- unlist(lapply(annotation, function(x) x$charge))
                 Formula <- unlist(lapply(annotation, function(x) x$Formula))
                 mz.error <- unlist(lapply(annotation, function(x) x$mz.error))
                 rt.error <- unlist(lapply(annotation, function(x) x$rt.error))
                 int.error <- unlist(lapply(annotation, function(x) x$int.error))
                 ms2.sim <- unlist(lapply(annotation, function(x) x$ms2.sim))
                 score <- unlist(lapply(annotation, function(x) x$score))
                 result <- data.frame(name, mz, rt, type, From, From.peak, to, step,
                                      level, as.seed,
                                      as.seed.round, isotope, adduct, charge, Formula,
                                      mz.error, rt.error, int.error, ms2.sim,
                                      score, stringsAsFactors = FALSE)
                 result <- result[!duplicated(result),]
                 result <- result[result$score >= score.cutoff,]
                 if(nrow(result) == 0){return(c(peak.name, peak.mz, peak.rt, rep(NA, 17)))}
                 type <- paste(result$type, collapse = ";")
                 From <- paste(result$From, collapse = ";")
                 From.peak <- paste(result$From.peak, collapse = ";")
                 to <- paste(result$to, collapse = ";")
                 step <- paste(result$step, collapse = ";")
                 level <- paste(result$level, collapse = ";")
                 as.seed <- paste(result$ as.seed, collapse = ";")
                 as.seed.round <- paste(result$as.seed.round, collapse = ";")
                 as.seed.round <- paste(result$as.seed.round, collapse = ";")
                 isotope <- paste(result$isotope, collapse = ";")
                 adduct <- paste(result$adduct, collapse = ";")
                 charge <- paste(result$charge, collapse = ";")
                 Formula <- paste(result$Formula, collapse = ";")
                 mz.error <- paste(result$mz.error, collapse = ";")
                 rt.error <- paste(result$rt.error, collapse = ";")
                 int.error <- paste(result$int.error, collapse = ";")
                 ms2.sim <- paste(result$ms2.sim, collapse = ";")
                 score <- paste(result$score, collapse = ";")
                 result <- data.frame(name, mz, rt, type, From, From.peak, to, step,
                                      level, as.seed,
                                      as.seed.round, isotope, adduct, charge, Formula,
                                      mz.error, rt.error, int.error, ms2.sim,
                                      score, stringsAsFactors = FALSE)
                 return(result)
               }
             }

             result <-
               BiocParallel::bplapply(X = tags2,
                                      temp.fun,
                                      BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                        progressbar = TRUE),
                                      score.cutoff = score.cutoff)

             result <- do.call(rbind, result)
             result

           })



#####---------------------------------------------------------------
# title result2Tags
# description Filter tags2 data according to result.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param result The result.
# param tags2 The tags2 data.
# return The tags2 data.
# export

setGeneric(name = "result2Tags",
           def = function(result,
                          tags2 = tags2){

             peak.name <- showTags2(tags2, slot = "name")
             unique.name <- unique(result$name)
             ##remove the annotation of peaks which are not in result
             idx <- which(!peak.name %in% unique.name)
             tags2[idx] <- lapply(tags2[idx], function(x) {
               x@annotation <- list()
               x
             })

             for(i in 1:length(unique.name)){
               temp.name <- unique.name[i]
               temp.idx1 <- which(temp.name == peak.name)
               temp.idx2 <- which(result$name == temp.name)
               temp.result <- result[temp.idx2,]
               remain.annotation <- apply(temp.result, 1, list)
               remain.annotation <- lapply(remain.annotation, function(x) {
                 paste(x[[1]][c("type", "From", "From.peak", "to", "step", "level", "as.seed",
                                "as.seed.round","isotope", "adduct", "charge", "Formula",
                                "mz.error", "rt.error","int.error","ms2.sim", "score")],
                       collapse = ";")
               })
               remain.annotation <- unlist(remain.annotation)

               annotation <- tags2[[temp.idx1]]@annotation
               annotation <- lapply(annotation, function(x) {
                 x <- unlist(x)
                 paste(x[c("type", "From", "From.peak", "to", "step", "level", "as.seed",
                           "as.seed.round","isotope", "adduct", "charge", "Formula",
                           "mz.error", "rt.error","int.error","ms2.sim", "score")],
                       collapse = ";")
               })
               annotation <- unlist(annotation)

               remain.idx <- which(annotation %in% remain.annotation)
               tags2[[temp.idx1]]@annotation <- tags2[[temp.idx1]]@annotation[remain.idx]
             }
             tags2
           })





setGeneric(name = "removeIsotopeFromResult",
           def = function(result){
             ##temp function
             temp.fun <- function(x,y){
               unlist(mapply(function(x,y){
                 if(all(is.na(y))) return(x)
                 x <- strsplit(x, split = ";")[[1]]
                 x <- x[-y]
                 if(length(x) == 0) return(NA)
                 paste(x, collapse = ";")
               },
               x = x,
               y = y))
             }

             #-----------------------------

             result <- as.data.frame(result)
             isotope <- result$isotope

             idx <- lapply(isotope, function(x){
               if(is.na(x)) return(NA)
               x <- strsplit(x, split = ";")[[1]]
               temp.idx <- which(x != "[M]")
               if(length(temp.idx) == 0) return(NA)
               temp.idx
             })

             result$type <- temp.fun(result$type, y = idx)
             result$From <- temp.fun(result$From, y = idx)
             result$From.peak <- temp.fun(result$From.peak, y = idx)
             result$to <- temp.fun(result$to, y = idx)
             result$step <- temp.fun(result$step, y = idx)
             result$level <- temp.fun(result$level, y = idx)
             result$as.seed <- temp.fun(result$as.seed, y = idx)
             result$as.seed.round <- temp.fun(result$as.seed.round, y = idx)
             result$isotope <- temp.fun(result$isotope, y = idx)
             result$adduct <- temp.fun(result$adduct, y = idx)
             result$charge <- temp.fun(result$charge, y = idx)
             result$Formula <- temp.fun(result$Formula, y = idx)
             result$mz.error <- temp.fun(result$mz.error, y = idx)
             result$rt.error <- temp.fun(result$rt.error, y = idx)
             result$int.error <- temp.fun(result$int.error, y = idx)
             result$ms2.sim <- temp.fun(result$ms2.sim, y = idx)
             result$score <- temp.fun(result$score, y = idx)

             result <- result
           })



#---------------------------------------------------------------------------
#' @title removeIsotope
#' @description remove isotope from result.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param result The MRN annotation table form MetDNA.
#' @export

setGeneric(name = "removeIsotope",
           def = function(result){
             ##temp function
             temp.fun <- function(x,y){
               unlist(mapply(function(x,y){
                 if(all(is.na(y))) return(x)
                 x <- strsplit(x, split = ";")[[1]]
                 x <- x[-y]
                 if(length(x) == 0) return(NA)
                 paste(x, collapse = ";")
               },
               x = x,
               y = y))
             }

             #-----------------------------

             result <- as.data.frame(result)
             isotope <- result$isotope

             idx <- lapply(isotope, function(x){
               if(is.na(x)) return(NA)
               x <- strsplit(x, split = ";")[[1]]
               temp.idx <- which(x != "[M]")
               if(length(temp.idx) == 0) return(NA)
               temp.idx
             })

             result$Annotation.type <- temp.fun(result$Annotation.type, y = idx)
             result$annotated.from.ID <- temp.fun(result$annotated.from.ID, y = idx)
             result$annotated.from.peak <- temp.fun(result$annotated.from.peak, y = idx)
             result$ID <- temp.fun(result$ID, y = idx)
             result$compound.name <- temp.fun(result$compound.name, y = idx)
             result$isotope <- temp.fun(result$isotope, y = idx)
             result$adduct <- temp.fun(result$adduct, y = idx)
             result$Formula <- temp.fun(result$Formula, y = idx)
             result$score <- temp.fun(result$score, y = idx)
             result$peak.group <- temp.fun(result$peak.group, y = idx)
             result$confidence <- temp.fun(result$confidence, y = idx)
             result <- result
           })



#----------------------------------------------------------------------------
# load("tags2.after.redundancy.remove")
# load("ms2")
# tags2 <- tags2.after.redundancy.remove
# type = "jpg"
# seed.col = "salmon"
# neighbor.col = "lightseagreen"
# no.matched.col = "grey"
# width = 1200
# height = 400
# is.include.precursor = TRUE
# is.tune.ms2.exp = TRUE
# is.tune.ms2.lib = FALSE
# path = "."
# ms2 = ms2
# data("kegg.compound")
# kegg.compound = kegg.compound
#
# dir.create("Seed_Neighbor")
# plotSeedNeighbor(tags2 = tags2.after.redundancy.remove,
#                  type = "jpg",
#                  path = "Seed_Neighbor",
#                  ms2 = ms2,
#                  kegg.compound = kegg.compound)

setGeneric(name = "plotSeedNeighbor",
           def = function(tags2,
                          type = c("pdf", "jpg"),
                          seed.col = "salmon",
                          neighbor.col = "lightseagreen",
                          no.matched.col = "grey",
                          width = ifelse(type == "pdf", 20, 1200),
                          height = ifelse(type == "pdf", 7, 400),
                          is.include.precursor = TRUE,
                          is.tune.ms2.exp = TRUE,
                          is.tune.ms2.lib = FALSE,
                          path = ".",
                          ms2 = ms2,
                          kegg.compound = kegg.compound,
                          ...){
             type <- match.arg(type)

             ms2.name <- unname(unlist(lapply(ms2, function(x){
               x[[1]][1,1]
             })))

             index <- which(showTags2(tags2 = tags2, slot = "annotation.len") > 0)
             if(length(index) == 0) return("No annotaion")
             # tags2 <- tags2[index]
             annotation.type <- showTags2(tags2 = tags2, slot = "annotation.type")

             index <- unname(which(unlist(lapply(annotation.type, function(x){
               any(unlist(x) == "metAnnotation")
             }))))

             annotation.type <- annotation.type[index]

             peak.index <- as.numeric(names(annotation.type))
             annotation.index <- lapply(annotation.type, function(x){
               which(x == "metAnnotation")
             })


#
#              x <- unlist(mapply(function(x, y){
#              unlist(lapply(tags2[[x]]@annotation[y], function(z) z$From.peak))
#              },
#              x = peak.index,
#              y = annotation.index))
#              x <- unique(x)



             # x <- NULL
             for(i in 1:length(peak.index)){
               # cat(match(i, peak.index)); cat(" ")
               cat(i); cat(" ")
              temp.neighbor <- tags2[[peak.index[i]]]
              temp.annotation.index <- annotation.index[[i]]
              neighbor.peak.name <- temp.neighbor@name
              neighbor.mz <- temp.neighbor@mz

              for(j in temp.annotation.index){
                # cat(j); cat(" ")
              temp.annotation <- temp.neighbor@annotation[[j]]
              neighbor.id <- temp.annotation$to
              neighbor.compound.name <-
                strsplit(kegg.compound$Name[match(neighbor.id, kegg.compound$ID)], split = ";")[[1]][1]
              neighbor.dp <- temp.annotation$ms2.sim

              seed.peak.name <- temp.annotation$From.peak
              seed.mz <- tags2[[match(seed.peak.name, showTags2(tags2, slot = "name"))]]@mz
              seed.id <- temp.annotation$From
              seed.compound.name <- strsplit(kegg.compound$Name[match(seed.id, kegg.compound$ID)], split = ";")[[1]][1]

              seed.info <- c(seed.peak.name, seed.id, seed.compound.name, seed.mz)
              neighbor.info <- c(neighbor.peak.name, neighbor.id, neighbor.compound.name, neighbor.mz, neighbor.dp)

              seed.ms2 <- ms2[[match(seed.peak.name, ms2.name)]][[2]]
              neighbor.ms2 <- ms2[[match(neighbor.peak.name, ms2.name)]][[2]]

              if(is.null(seed.ms2) | is.null(neighbor.ms2)){
                # x <- c(x, seed.peak.name)
                next()
              }else{

                if(type == "pdf"){
                  pdf(file.path(path, paste("Seed",seed.peak.name, seed.id,
                                            "Neighbor", neighbor.peak.name, neighbor.id, "pdf", sep = '.')),
                      width = width,
                      height = height)
                }else{
                  jpeg(file.path(path, paste("Seed",seed.peak.name, seed.id,
                                             "Neighbor", neighbor.peak.name, neighbor.id, "jpeg", sep = '.')),
                       width = width,
                       height = height)
                }

                par(mar = c(5,5,4,2))
                matchPlot(seed.info = seed.info, neighbor.info = neighbor.info,
                          seed.ms2 = seed.ms2, neighbor.ms2 = neighbor.ms2,
                          seed.col = seed.col, neighbor.col = neighbor.col,
                          no.matched.col = no.matched.col)
                dev.off()
              }
              }
             }

})







matchPlot <- function(seed.info,
                      neighbor.info,
                      seed.ms2,
                      neighbor.ms2,
                      seed.col = seed.col,
                      neighbor.col = neighbor.col,
                      no.matched.col = no.matched.col
){
  seed.ms2[,2] <- seed.ms2[,2]/max(seed.ms2[,2])
  neighbor.ms2[,2] <- neighbor.ms2[,2]/max(neighbor.ms2[,2])

  x.max <- max(as.numeric(seed.info[4]), as.numeric(neighbor.info[4]))
  x.min <- min(seed.ms2[,1], neighbor.ms2[,1])

  note <- annotateMS2(spec = seed.ms2, matched.spec = neighbor.ms2)
  col1 <- rep(no.matched.col, length(note))
  col1[note=="matched"] <- seed.col
  seed.ms2 <- data.frame(seed.ms2, note, stringsAsFactors = FALSE)

  note <- annotateMS2(spec = neighbor.ms2, matched.spec = seed.ms2)
  col2 <- rep(no.matched.col, length(note))
  col2[note=="matched"] <- neighbor.col
  neighbor.ms2 <- data.frame(neighbor.ms2, note,
                             stringsAsFactors = FALSE)

  plot(0, xlim = c(x.min, x.max), ylim = c(-1,1), col = "white",
       xlab = "m/z", ylab = "Relative intensity",
       cex.lab = 1.8, cex.axis = 1.5)
  abline(h = 0)
  abline(v = min(as.numeric(seed.info[4]), as.numeric(neighbor.info[4])) + 1,
         lwd = 1.5, col = "orchid4")
  legend("topright",
         legend = paste("No matched range\n", ">=mz",
                        min(as.numeric(seed.info[4]), as.numeric(neighbor.info[4]))),
         bty = "n", text.col = 'orchid4', cex = 1.2)

  points(seed.ms2[,1], seed.ms2[,2], type = "h", col = col1)
  points(neighbor.ms2[,1], -neighbor.ms2[,2], type = "h", col = col2)

  points(seed.ms2[seed.ms2[,3] == "matched",1],
         seed.ms2[seed.ms2[,3] == "matched",2], type = "p", col = seed.col, pch = 19)
  points(neighbor.ms2[neighbor.ms2[,3] == "matched",1],
         -neighbor.ms2[neighbor.ms2[,3] == "matched",2], type = "p", col = neighbor.col, pch = 19)


  maptools::pointLabel(x = x.min, y = 1,
                       labels = paste("Seed",
                                      "\nPeak name:", seed.info[1],
                                      "\nKEGG ID:", seed.info[2],
                                      "\nCompound name:", seed.info[3]),
                       cex = 1.2)

  maptools::pointLabel(x = x.min, y = -1,
                       labels = paste("Neighbor",
                                      "\nPeak name:", neighbor.info[1],
                                      "\nKEGG ID:", neighbor.info[2],
                                      "\nCompound name:", neighbor.info[3],
                                      "\nDot product:", round(as.numeric(neighbor.info[5]), 2)),
                       cex = 1.2)

}



annotateMS2 <- function(spec, matched.spec, ppm.ms2match = 30){
  note <- rep(NA, nrow(spec))
  for(i in 1:length(note)){
    note[i] <-
      ifelse(any(abs(as.numeric(spec[i,1]) - as.numeric(matched.spec[,1]))*10^6/ifelse(as.numeric(spec[i,1])>=400, as.numeric(spec[i,1]), 400) < ppm.ms2match), "matched", "no")
  }
  note
}




#-----------------------------------------------------------------------------
##Note, Now the score cutoff is not open for users
setGeneric(name = "getAnnotationResult",
           def = function(tags2,
                          p.value,
                          correct = TRUE,
                          foldchange,
                          tags.result2,
                          kegg.compound = kegg.compound,
                          candidate.num = 3000,
                          score.cutoff = 0.4){

             temp <- tags2Result(tags2 = tags2,
                                 score.cutoff = score.cutoff,
                                 candidate.num = candidate.num)

             colnames(temp)[7] <- "ID"

             confidence <- lapply(temp$name, function(x){
               temp.idx <- which(tags.result2$name == x)
               if(length(temp.idx) == 0) return(NA)
               paste(tags.result2$Confidence[temp.idx], collapse = ";")
             })

             peak.group <- lapply(temp$name, function(x){
               temp.idx <- which(tags.result2$name == x)
               if(length(temp.idx) == 0) return(NA)
               paste(tags.result2$group[temp.idx], collapse = ";")
             })

             confidence <- unlist(confidence)
             peak.group <- unlist(peak.group)

             id <- temp$ID
             compound.name <- unlist(lapply(id, function(x){
               if(is.na(x)){
                 return(NA)
               }else{
                 x <- strsplit(x = x, split = ";")[[1]]
                 temp.name <- kegg.compound$Name[match(x, kegg.compound$ID)]
                 temp.name <- unlist(lapply(strsplit(temp.name, split = ";"), function(x) x[1]))
                 paste(temp.name, collapse = ";")
               }
             }))

             temp <- data.frame(temp, compound.name, peak.group,
                                confidence, stringsAsFactors = FALSE)

             if(!missing(p.value) & !missing(foldchange)){
               temp <- data.frame(temp, p.value, foldchange,
                                  stringsAsFactors = FALSE)
               if(correct){colnames(temp)[24] <- "p.value.adjusted"}
             }

             temp <- temp[,-c(8, 9, 10, 11, 14, 16, 17, 18, 19)]
             colnames(temp)[c(5,6)] <- c("annotated.from.ID", "annotated.from.peak")
             colnames(temp)[4] <- "Annotation.type"
             # temp$type <- lapply(temp.type, function(x){
             #   if(is.na(x)) return(NA)
             #   x <- strsplit(x, split = ";")[[1]]
             # })

             temp <- data.frame(temp[,c(1:7, 12)], temp[,-c(1:7, 12)],
                                stringsAsFactors = FALSE)
             temp <- temp

})




##--------------------------------------------------------------------------
setGeneric(name = "getPathway", def = function(species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                                             "dre", "cel", "sce", "ath", "smm", "pfa",
                                                                   "tbr", "eco", "ppu", "syf"),
                                               type = c("gene", "metabolite")){

 type <- match.arg(type)
 species <- match.arg(species)

 if(type == "metabolite"){
   switch(species,
          "hsa" = {data("hsa.kegg.pathway", envir = environment())
            pathway = hsa.kegg.pathway},
          "dme" = {data("dme.kegg.pathway", envir = environment())
            pathway = dme.kegg.pathway},
          "mmu" = {data("mmu.kegg.pathway", envir = environment())
            pathway = mmu.kegg.pathway},
          "rat" = {data("rat.kegg.pathway", envir = environment())
            pathway = rat.kegg.pathway},
          "bta" = {data("bta.kegg.pathway", envir = environment())
            pathway = bta.kegg.pathway},
          "gga" = {data("gga.kegg.pathway", envir = environment())
            pathway = gga.kegg.pathway},
          "dre" = {data("dre.kegg.pathway", envir = environment())
            pathway = dre.kegg.pathway},
          "cel" = {data("cel.kegg.pathway", envir = environment())
            pathway = cel.kegg.pathway},
          "sce" = {data("sce.kegg.pathway", envir = environment())
            pathway = sce.kegg.pathway},
          "ath" = {data("ath.kegg.pathway", envir = environment())
            pathway = ath.kegg.pathway},
          "smm" = {data("smm.kegg.pathway", envir = environment())
            pathway = smm.kegg.pathway},
          "pfa" = {data("pfa.kegg.pathway", envir = environment())
            pathway = pfa.kegg.pathway},
          "tbr" = {data("tbr.kegg.pathway", envir = environment())
            pathway = tbr.kegg.pathway},
          "eco" = {data("eco.kegg.pathway", envir = environment())
            pathway = eco.kegg.pathway},
          "ppu" = {data("ppu.kegg.pathway", envir = environment())
            pathway = ppu.kegg.pathway},
          "syf" = {data("syf.kegg.pathway", envir = environment())
            pathway = syf.kegg.pathway}
   )
 }


 if(type == "gene"){
   switch(species,
          "hsa" = {data("hsa.gene.kegg.pathway", envir = environment())
            pathway = hsa.gene.kegg.pathway},
          "dme" = {data("dme.gene.kegg.pathway", envir = environment())
            pathway = dme.gene.kegg.pathway},
          "mmu" = {data("mmu.gene.kegg.pathway", envir = environment())
            pathway = mmu.gene.kegg.pathway},
          "rat" = {data("rat.gene.kegg.pathway", envir = environment())
            pathway = rat.gene.kegg.pathway},
          "bta" = {data("bta.gene.kegg.pathway", envir = environment())
            pathway = bta.gene.kegg.pathway},
          "gga" = {data("gga.gene.kegg.pathway", envir = environment())
            pathway = gga.gene.kegg.pathway},
          "dre" = {data("dre.gene.kegg.pathway", envir = environment())
            pathway = dre.gene.kegg.pathway},
          "cel" = {data("cel.gene.kegg.pathway", envir = environment())
            pathway = cel.gene.kegg.pathway},
          "sce" = {data("sce.gene.kegg.pathway", envir = environment())
            pathway = sce.gene.kegg.pathway},
          "ath" = {data("ath.gene.kegg.pathway", envir = environment())
            pathway = ath.gene.kegg.pathway},
          "smm" = {data("smm.gene.kegg.pathway", envir = environment())
            pathway = smm.gene.kegg.pathway},
          "pfa" = {data("pfa.gene.kegg.pathway", envir = environment())
            pathway = pfa.gene.kegg.pathway},
          "tbr" = {data("tbr.gene.kegg.pathway", envir = environment())
            pathway = tbr.gene.kegg.pathway},
          "eco" = {data("eco.gene.kegg.pathway", envir = environment())
            pathway = eco.gene.kegg.pathway},
          "ppu" = {data("ppu.gene.kegg.pathway", envir = environment())
            pathway = ppu.gene.kegg.pathway},
          "syf" = {data("syf.gene.kegg.pathway", envir = environment())
            pathway = syf.gene.kegg.pathway}
   )
 }

 pathway <- pathway



})






##---------------------------------------------------------------------------
setGeneric(name = "filterAnnotationByEnrichment",
           def = function(module,
                          index,
                          anno1,
                          anno2,
                          result2,
                          tags2,
                          match.result2,
                          mz.tol = 25,
                          rt.tol = 30){

             annotation1 <- vector(mode = "list", length = length(index))
             annotation2 <- vector(mode = "list", length = length(index))


             for(i in index){
               # cat(i); cat(" ")
               temp.module <- module[[i]]
               anno1.result <- lapply(anno1, function(x){
                 if(is.na(x)){return(NA)}
                 x <- strsplit(x, split = ";")[[1]]
                 idx <- which(x %in% temp.module)
                 if(length(idx) == 0) return(NA)
                 return(paste(x[idx], collapse = ";"))
               })
               anno1.result <- unlist(anno1.result)

               anno2.result <- lapply(anno2, function(x){
                 if(is.na(x)){return(NA)}
                 x <- strsplit(x, split = ";")[[1]]
                 idx <- which(x %in% temp.module)
                 if(length(idx) == 0) return(NA)
                 return(paste(x[idx], collapse = ";"))
               })

               anno2.result <- unlist(anno2.result)

               annotation1[[i]] <- anno1.result
               annotation2[[i]] <- anno2.result
             }


             annotation1 <- do.call(cbind, annotation1)
             annotation2 <- do.call(cbind, annotation2)


             peak.name <- result2$name

             annotation1 <- data.frame(peak.name,
                                       unlist(anno1), annotation1,
                                       stringsAsFactors = FALSE)
             annotation2 <- data.frame(peak.name,
                                       unlist(anno2), annotation2,
                                       stringsAsFactors = FALSE)

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


             #new anno should be from prefer.annotation1 and prefer.annotation2
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
                                        paste('Module', 1:(ncol(annotation1)-2), sep = ""))
             colnames(annotation2) <- c("Peak.name","KEGG.annotation",
                                        paste('Module', 1:(ncol(annotation2)-2), sep = ""))

             annotation1 <- data.frame(annotation1[,1:2], perfect.annotation1,
                                       annotation1[,-c(1:2)],
                                       stringsAsFactors = FALSE)
             annotation2 <- data.frame(annotation2[,1:2], perfect.annotation2,
                                       annotation2[,-c(1:2)],
                                       stringsAsFactors = FALSE)
             colnames(annotation1)[3] <- colnames(annotation2)[3] <- "Annotation"


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
                 mzError.ppm <- abs(temp.match.result$Accurate.mass - peakMz)*10^6/ifelse(peakMz>=400,peakMz,400)
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


               data("kegg.compound", envir = environment())
               for(i in 1:length(kegg.result)){
                 tags2 <- kegg2peakInfo(kegg.result = kegg.result[[i]],
                                        mz.tol = mz.tol,
                                        rt.tol = rt.tol,
                                        weight.mz = 0.5, weight.rt = 0.5,weight.dp = 0,
                                        peak.info = tags2,
                                        kegg.compound = kegg.compound)
               }
             }else{
               tags2 <- NULL
             }


             return.result <- list(tags2, anno, anno1, anno2, annotation1, annotation2)
             names(return.result) <- c("tags2", "anno", "anno1", "annot2", "annotation1", "annotation2")
             return.result <- return.result
           })





#------------------------------------------------------------------------------
setGeneric(name = "removePeak",
           def = function(sample,
                          sample.info,
                          by = c("mv", "zero",),
                          mz.cutoff = 0.5,
                          any = FALSE){
             sample.name <- colnames(sample)
             group <- NULL
             group[grep("\\.01", colnames(sample))] <- "Cancer"
             group[grep("\\.11", colnames(sample))] <- "Benign"
             sample.info <- data.frame(sample.name, group,
                                       stringsAsFactors = FALSE)

           })


#
# sample.info <- readr::read_table2(dir()[1])
# sample.info <- as.data.frame(sample.info)




#------------------------------------------------------------------------------
setGeneric(name = "getPostfix", def = function(x){
  unlist(lapply(strsplit(x = x, split = "\\."),  function(y) {
    if(length(y) == 1) return(NA)
    y[[length(y)]]}
    ))
})




#------------------------------------------------------------------------------
errorDisplay <- function(expr,
                         error.info = "error"){
  process.result <- try(expr,
                        silent = TRUE)
  if(class(process.result) == "try-error"){
    cat(error.info, "\n")
    process.result <- process.result
  }else{
    process.result <- "Right"
  }
}



#------------------------------------------------------------------------------
setGeneric(name = "changeFile", def = function(path){
file.rename(from = file.path(path, "volcano.plot.pdf"), to = file.path(path, "1Volcano.plot.pdf"))
  file.rename(from = file.path(path, "1Volcano.plot.pdf"), to = file.path(path, "Volcano.plot.pdf"))
  file.rename(from = file.path(path, "pathway.heatmap.pdf"), to = file.path(path, "Pathway.heatmap.pdf"))
  file.rename(from = file.path(path, "pathway.heatmap.pdf"), to = file.path(path, "1pathway.heatmap.pdf"))
  file.rename(from = file.path(path, "1pathway.heatmap.pdf"), to = file.path(path, "Pathway.heatmap.pdf"))
  file.rename(from = file.path(path, "Quantitative_information", "pathway.node.quantitative.result.csv"),
              to = file.path(path, "Quantitative_information", "Quantitative.pathway.metabolite.result.csv"))
  file.copy(from = file.path(path, "Quantitative_information", "Quantitative.pathway.metabolite.result.csv"),
            to = file.path(path, "Quantitative.pathway.metabolite.result.csv"))

  file.rename(from = file.path(path, "Quantitative_information", "pathway.quantitative.result.csv"),
              to = file.path(path, "Quantitative_information", "Quantitative.pathway.result.csv"))
  file.copy(from = file.path(path, "Quantitative_information", "Quantitative.pathway.result.csv"),
            to = file.path(path, "Quantitative.pathway.result.csv"))

  unlink(x = file.path(path, "Quantitative_information"), recursive = TRUE)
})





#####removeTags.result
setGeneric(name = "removeTagsResult",
           function(tags.result,
                    candidate.num = 5){

             tags.result <- lapply(unique(tags.result$name),
                                   function(x){
                              temp.idx <- which(tags.result$name == x)
                              temp.data <- tags.result[temp.idx,,drop = FALSE]
                              if(nrow(temp.data) > candidate.num) temp.data <- temp.data[1:candidate.num,]
                              temp.data
             })

            tags.result <- do.call(rbind, tags.result)
            return(tags.result)
})

