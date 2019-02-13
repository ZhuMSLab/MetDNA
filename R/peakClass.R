#' @title peakInfo
#' @description The Peakinfo class is the S4 class for Peak annotation
#' information.
#' @slot name Peak name.
#' @slot mz Peak m/z.
#' @slot rt Peak rentention time (RT).
#' @slot ms2 MS/MS spectrum of peak.
#' @slot annotation Annotation information of peak. It is a list contains
#' different annotation result.
#' @name peakInfo-Class
#' @exportClass peakInfo
#' @author Xiaotao Shen
#' @export
setClass("peakInfo",
         representation(name = "character",
                        mz = "numeric",
                        rt = "numeric",
                        ms2 = "data.frame",
                        annotation = "list")
         # prototype(name = "M102T561",
         #           mz = 102.0897,
         #           rt = 560.593,
         #           ms2 = data.frame(
         #             mz = c(41.03722,42.03253,43.01675,43.04068,
         #                    43.05172,44.04867,56.04833,57.88836,
         #                    58.06352,58.40606,59.07152,61.02475,
         #                    70.05883,72.04043,74.09432,87.06190,
         #                    102.08545,219.80564,219.90548,220.13019,
         #                    220.52996,222.31249,223.11005,223.85028),
         #             intensity = c(73,122,261,49,49,73,24,49,5113,24,
         #                           381,24,24,24,24,24,24,24,24,
         #                           24,24,24,24,24),stringsAsFactors = FALSE),
         #           annotation = list(
         #             list(type = "seed",
         #                  From = NA,
         #                  From.peak = NA,
         #                  to = "C00576",
         #                  step = NA,
         #                  level = 1,
         #                  as.seed = FALSE,
         #                  as.seed.round = NA,
         #                  isotope = "[M]",
         #                  addcut = "M+H",
         #                  charge = 1,
         #                  Formula = "C5H11NO",
         #                  mz.error = 0,
         #                  rt.error = 0,
         #                  int.error = NA,
         #                  ms2.sim = 1,
         #                  score = 1
         #                  )
         #           )
         #           )
)


#-----------------------------------------------------------------------------
setMethod(f = "show",
          signature   = "peakInfo",
          definition = function(object) {
            cat("Peak information\n")
            cat("--------------\n")
            cat("Peak name:", object@name)
            cat("\n")
            cat("mz:", object@mz)
            cat("\n")
            cat("rt:", object@rt)
            cat("\n")
            cat("MS/MS spectrum:", ifelse(nrow(object@ms2) > 0, "Yes\n", "No\n"))
            # cat(ifelse(nrow(object@ms2) > 0, "Yes\n", "No\n"))
            # format(object@ms2)
            cat('\n')
            cat("--------------\n")
            cat("Annotation information\n")
            cat("--------------\n")
            if(length(object@annotation) == 0){
              cat("No annotation.")
            }else{
              for(i in 1:length(object@annotation)){
                cat("\n")
                cat("Annotation:",i)
                cat("\n")
                cat("Type:", object@annotation[[i]]$type)
                cat("\n")
                cat("From:", object@annotation[[i]]$From)
                cat("\n")
                cat("From which peak:", object@annotation[[i]]$From.peak)
                cat("\n")
                cat("Step:", object@annotation[[i]]$step)
                cat("\n")
                cat("Annotation result (KEGG ID):", object@annotation[[i]]$to)
                cat("\n")
                cat("Annotation level:", object@annotation[[i]]$level)
                cat("\n")
                cat("Is seed or not:", object@annotation[[i]]$as.seed)
                cat("\n")
                cat("In which round as seed:", object@annotation[[i]]$as.seed.round)
                cat("\n")
                cat("Isotpe:", object@annotation[[i]]$isotope)
                cat("\n")
                cat("Adduct:", object@annotation[[i]]$adduct)
                cat("\n")
                cat("Charge:", object@annotation[[i]]$charge)
                cat("\n")
                cat("Formula:", object@annotation[[i]]$Formula)
                cat("\n")
                cat("m/z error (ppm):", object@annotation[[i]]$mz.error)
                cat("\n")
                cat("RT error:", object@annotation[[i]]$rt.error)
                cat("\n")
                cat("Intensity ratio error:", object@annotation[[i]]$int.error)
                cat("\n")
                cat("MS/MS similarity (dot product):", object@annotation[[i]]$ms2.sim)
                cat("\n")
                cat("Score:", object@annotation[[i]]$score)
                cat("\n")
              }
            }

          }
)




setGeneric(
  name = "getInfo",
  def  = function(object,...) {
    standardGeneric("getInfo")
  }
)


setMethod(f = "getInfo",
          signature  = "peakInfo",
          definition = function(object) {
            name <- object@name
            mz <- object@mz
            rt <- object@rt
            annotation <- object@annotation
            if(length(annotation) == 1){
              annotation <- unlist(annotation)
            }else{
              idx.max.score <- which.max(unlist(lapply(annotation, function(x) {x$score})))
              annotation <- unlist(annotation[[idx.max.score]])
              names(annotation)[14] <- "Score"
            }
            result <- c(name, mz, rt, annotation)
            names(result)[1:3] <- c("name", "mz", "rt")
            return(result)
          }
)







setGeneric(
  name = "filterPeak",
  def  = function(object,...) {
    standardGeneric("filterPeak")
  }
)

# title filterPeak
# description Filter peakInfo data.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param object peakInfo data.
# param score.thr Score cutoff.
# return peakInfo data.


setMethod(f = "filterPeak",
          signature  = "peakInfo",
          definition = function(object,
                                score.thr = 0) {
            annotation <- object@annotation
            if(length(annotation) != 0){
              score <- unlist(lapply(annotation, function(x) {x$score}))
              annotation <- annotation[order(score, decreasing = TRUE)]
              score <- unlist(lapply(annotation, function(x) {x$score}))
              if(length(which(score >= score.thr)) > 0){
                annotation <- annotation[which(score >= score.thr)]
                object@annotation <- annotation
              }else{
                object@annotation <- list()
              }

            }
            return(object)
          }
)



#' @title equalPeakinfo
#' @description Compare two peakInfo data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object1 peakInfo.
#' @param object2 peakInfo.
#' @return TRUE or FALSE.
#' @export

setGeneric(name = "equalPeakinfo",
           def = function(object1,
                          object2) {
             name1 <- object1@name
             mz1 <- object1@mz
             rt1 <- object1@rt
             annotation1 <- object1@annotation

             name2 <- object2@name
             mz2 <- object2@mz
             rt2 <- object2@rt
             annotation2 <- object2@annotation

             if(length(annotation1) != length(annotation2)){
               return(FALSE)
             }

             if(name1 == name2 & mz1 == mz2 & rt1 == rt2){
               if(length(annotation1) == 0){
                 return(TRUE)
               }else{
                 anno1 <- unlist(annotation1)
                 anno2 <- unlist(annotation2)
                 anno1 <- anno1[!is.na(anno1)]
                 anno2 <- anno2[!is.na(anno2)]
                 if(all(anno1 == anno2)){
                   return(TRUE)
                 }else{
                   return(FALSE)
                 }

               }
             }else{
               return(FALSE)
             }

           })


#' @title equalPeakinfo1
#' @description Compare two tags2 data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param x The tags2 data.
#' @param y The tags2 data.
#' @return TRUE or FALSE.
#' @export


setGeneric(name = "equalPeakinfo1",
           def = function(x,y){
             result <-
               mapply(function(x,y){equalPeakinfo(x, y)}, x, y)
             result
           })







# dataInfo: Data information
#
# description The data inforamtion  class is the S4 class for metabolomics
# information.
# slot name Peak name.
# slot mz Peak m/z.
# slot rt Peak rentention time (RT).
# slot ms2 MS/MS spectrum of peak.
# slot annotation Annotation information of peak. It is a list contains
# different annotation result.
# name dataInfo-Class
# exportClass dataInfo
# author Xiaotao Shen
# export

# setClass("dataInfo",
#          representation(sample.name = "character",
#                         peak.name = "character",
#                         sample.group = "list",
#                         sample.batch = "list",
#                         sample.data = "data.frame",
#                         sample.pheno = "data.frame",
#                         peak.tags = "data.frame")
# )
#
#
#
# setGeneric(name = "readData",
#            def = function(sample = "sample.csv",
#                           tags = "tags.csv",
#                           info = "info.csv",
#                           pheno = NULL){
#              sample <- readr::read_csv(sample, progress = FALSE, col_types = readr::cols())
#              tags <- readr::read_csv(tags, progress = FALSE, col_types = readr::cols())
#              info <- readr::read_csv(info, progress = FALSE, col_types = readr::cols())
#
#              if(is.null(pheno)){
#                sample.pheno <- data.frame()
#              }else{
#                sample.pheno <- readr::read_csv(pheno, progress = FALSE, col_types = readr::cols())
#              }
#
#              sample.name1 <- colnames(sample)
#              sample.name2 <- info$sample.name
#              if(any(sort(sample.name1) != sort(sample.name2))) {
#                stop("Sample names from sample are different with sample names
#                     which are provided in info.")
#              }
#
#              sample <- sample[,order(sample.name2)]
#              sample.name <- sample.name2
#              peak.name <- tags$Peak.name
#
#              sample.group <- info$group
#              unique.sample.group <- unique(sample.group)
#              sample.group <-
#                lapply(unique.sample.group, function(x) {sample.name[which(info$group == x)]})
#              names(sample.group) <- unique.sample.group
#
#              sample.batch <- info$batch
#              unique.sample.batch <- unique(sample.batch)
#              sample.batch <-
#                lapply(unique.sample.batch, function(x) {sample.name[which(info$batch == x)]})
#              names(sample.batch) <- unique.sample.batch
#
#              sample.data <- as.data.frame(sample)
#
#              peak.tags <- as.data.frame(tags)
#
#              dataInfo <- new(Class = "dataInfo",
#                              sample.name = sample.name,
#                              peak.name = peak.name,
#                              sample.group = sample.group,
#                              sample.batch = sample.batch,
#                              sample.data = sample.data,
#                              sample.pheno = sample.pheno,
#                              peak.tags = peak.tags)
#              })
#
#
#
# setMethod(f = "show",
#           signature   = "dataInfo",
#           definition = function(object) {
#             cat("The data has",nrow(object@sample.data), "peaks and", ncol(object@sample.data), "samples\n")
#             cat("(1) Sample information\n")
#             cat("--------------\n")
#             cat("Sample has", length(object@sample.group), "groups:\n")
#             cat(names(object@sample.group))
#             cat("\n")
#             cat("Sample has", length(object@sample.batch), "batches.\n")
#
#             cat("Sample names are",
#                 object@sample.name[ifelse(length(object@sample.name) > 5,
#                                           list(c(1:5)),
#                                           list(1:length(object@sample.name)))[[1]]])
#             if(length(object@sample.name) > 5){
#               cat("...")
#             }
#             cat("\n")
#             if(nrow(object@sample.pheno) == 0){
#               cat("There are not sample pheno information.\n")
#             }else{
#               cat("The sample pheno contains:\n")
#               cat(colnames(sample.pheno))
#             }
#             cat("\n")
#             cat("(2) Peak information\n")
#             cat("--------------\n")
#             cat("The peak information contains:\n")
#             cat(colnames(peak.tags))
#           }
# )
#
# setGeneric(
#   name = "tTest",
#   def  = function(object,...) {
#     standardGeneric("tTest")
#   }
# )
#
# setMethod(f = "tTest",
#           signature = "dataInfo",
#           definition = function(object,
#                                 group = c("W03", "W30"),
#                                 adjust = TRUE,
#                                 method = c("holm", "hochberg",
#                                            "hommel", "bonferroni", "BH", "BY",
#                                            "fdr", "none"),
#                                 alternative = c("two.sided", "less", "greater"),
#                                 mu = 0, paired = FALSE, var.equal = FALSE,
#                                 conf.level = 0.95, ...){
#             method = match.arg(method)
#             alternative = match.arg(alternative)
#             sample <- object@sample.data
#             group <- object@sample.group[group]
#             group1 <- sample[,group[[1]]]
#             group2 <- sample[,group[[2]]]
#
#             group1 <- apply(group1, 1, list)
#             group1 <- lapply(group1, function(x) x[[1]])
#
#             group2 <- apply(group2, 1, list)
#             group2 <- lapply(group2, function(x) x[[1]])
#
#             p.value <- mapply(function(x, y) {t.test(x, y)$p.value}, x = group1, y = group2)
#             if(adjust){
#               p.value <- p.adjust(p = p.value, method = method)
#             }
#             peak.name <- object@peak.tags$Peak.name
#             names(p.value) <- peak.name
#             return(p.value)
#           })




# setGeneric(
#   name = "foldChange",
#   def  = function(object,...) {
#     standardGeneric("foldChange")
#   }
# )
#
# setMethod(f = "foldChange",
#           signature = "dataInfo",
#           definition = function(object,
#                                 group = c("W03", "W30"),
#                                 method = c("median", "mean")){
#             method = match.arg(method)
#             sample <- object@sample.data
#             group <- object@sample.group[group]
#             group1 <- sample[,group[[1]]]
#             group2 <- sample[,group[[2]]]
#
#             group1 <- apply(group1, 1, method)
#             group2 <- apply(group2, 1, method)
#
#             fc <- group2/group1
#             fc[is.nan(fc)] <- 1
#             fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
#             fc[fc==0] <- min(fc[fc!=0])
#
#             peak.name <- object@peak.tags$Peak.name
#             names(fc) <- peak.name
#             return(fc)
#           })
