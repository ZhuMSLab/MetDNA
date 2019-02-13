
#' @title moduleInfo
#' @description The moduleInfo class is the S4 class for module.
#' @slot name Peak name.
#' @slot mz Peak m/z.
#' @slot ms2 MS/MS spectrum of peak.
#' @slot annotation Annotation information of peak. It is a list contains
#' different annotation result.
#' @name peakInfo-Class
#' @exportClass peakInfo
#' @author Xiaotao Shen
#' @export
setClass("moduleInfo",
         representation(Module.name = "character",
                        Module.size = "numeric",
                        Module.impact = "numeric",
                        p.value = "numeric",
                        Detected.number = "numeric",
                        Hidden.number = "numeric",
                        All.metabolite.id = "character",
                        Detected.ID = "character",
                        Hidden.ID = "character",
                        All.metabolite.name = "character",
                        Detected.metabolite.name = "character",
                        Hidden.metabolite.name = "character",
                        msea = "data.frame")
         # prototype(Module.name = "module1",
         # Module.size = 46,
         # Module.impact = 0.6470588,
         # p.value = 0.00158,
         # Detected.number = 27,
         # Hidden.number = 19,
         # All.metabolite.id = c("C00031","C00185","C00092","C00095","C00103",
         # "C00198","C00207","C00208","C00221","C00243","C00252","C00561",
         # "C00267","C00345", "C00530", "C01172","C00689","C00844","C01113",
         # "C01594","C02013","C05394","C06217","C01083","C16143","C02659",
         # "C03107","C02655", "C01236","C02995","C04534","C02048","C01093",
         # "C06219","C02964","C05403","C08325","C05396","C01870","C06218",
         # "C06187","C13636", "C15548","C16827","C01518","C06468"),
         # Detected.ID = c("C00031","C00185","C00095","C00198","C00208"
         # "C00221","C00243","C00252","C00267","C01594","C02013","C05394"
         # "C06217","C01083","C16143","C03107","C02655","C02048","C01093"
         # "C06219","C02964","C05403","C08325","C01870","C15548","C01518",
         # "C06468"),
         # Hidden.ID = c("C00092","C00103","C00207","C00561","C00345","C00530"
         # "C01172","C00689","C00844","C01113","C02659","C01236","C02995",
         # "C04534","C05396","C06218","C06187","C13636","C16827"),
         # All.metabolite.name = "character",
         # Detected.metabolite.name = "character",
         # Hidden.metabolite.name = "character",
         #            msea= list()
         #           )
)


setMethod(f = "show",
          signature   = "moduleInfo",
          definition = function(object) {
            cat("Module name:", object@Module.name, "\n", sep = " ")
            cat("Module size:", object@Module.size, "\n", sep = " ")
            cat("Module impact:", object@Module.impact, "\n", sep = " ")
            cat("Module p value:", object@p.value, "\n", sep = " ")
            cat("Detected metabolite number:", object@Detected.number, "\n", sep = " ")
            cat("Hidden metabolite number:", object@Hidden.number, "\n", sep = " ")
            cat("All metabolite:",
                object@All.metabolite.id[1:5][!is.na(object@All.metabolite.id[1:5])],
                "...\n", sep = " ")
            cat("Detected metabolite:",
                object@Detected.ID[1:5][!is.na(object@Detected.ID[1:5])],
                "...\n", sep = " ")
            cat("Hidden metabolite:",
                object@Hidden.ID[1:5][!is.na(object@Hidden.ID[1:5])],
                "...\n", sep = " ")
            if(nrow(object@msea) == 0){
              cat("No pathway are enriched.")
            }else{
              cat("---------------------------------------------------------\n")
              cat("Metabolite sets enrichment result\n")
              print(object@msea)
            }

          }
)




#-----------------------------------------------------------------------------
#' @title moduleplot
#' @description Barplot of pathway enrichment result.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object moduleInfo object.
#' @param cex.label The font size of label.
#' @param mar margin size.
#' @param n The number of pathway.
#' @param main Main of plot.
#' @param show p.value or voerlap.
#' @param ... see ?barplot
#' @export

setGeneric(
  name = "moduleplot",
  def  = function(object,
                  cex.label = 1,
                  mar,
                  n,
                  main,
                  show = c("p", "overlap"),
                  ...) {
    msea <- object@msea
    module.name <- object@Module.name
    if(nrow(msea) == 0) {
      layout(1)
      par(mar = c(5, 5, 4, 2))
      plot(0, 0, col = "white", xaxt = "n", yaxt = "n",
           xlab = "", ylab = "", main = module.name, cex.main = 1.5)
      text(x = 0, y = 0, labels = 'No pathways are enriched.',
           cex = 1.5)
    }else{
      show <- match.arg(show)
      module.name <- object@Module.name
      overlap <- msea$Overlap/msea$Pathway.length
      msea$Overlap <- overlap


      switch(show,
             "p" = {msea <- msea[order(msea$q.value),]},
             "overlap" = {msea <- msea[order(msea$Overlap, decreasing = TRUE),]})

      if(nrow(msea) == 0) {
        layout(1)
        par(mar = c(5, 5, 4, 2))
        plot(0, 0, col = "white", xaxt = "n", yaxt = "n",
             xlab = "", ylab = "", main = module.name, cex.main = 1.5)
        text(x = 0, y = 0, labels = 'No pathways are enriched.',
             cex = 1.5)
      }else{
        if(missing(n)) {
          msea <- msea
        }else{
          # msea <- msea[order(msea$p.value),]
          msea <- msea[1:n,]
          msea <- msea[!is.na(msea[,1]),]
        }

        if(show == "p"){
          msea <- msea[order(msea$q.value, decreasing = TRUE),]
          pathway.name <- rownames(msea)
          pathway.name <- unname(unlist(sapply(pathway.name, function(x){
            strsplit(x, split = ";")[[1]][1]
          })))
          q.value <- msea$q.value
          q.value <- -log(q.value, 10)
          # overlap <- msea$Overlap
          max.len <- max(nchar(pathway.name))
          left.margin <- max.len * cex.label / 2.4

          if(missing(mar)) {
            par(mar = c(5,left.margin,4,2))
          }else{
            par(mar = mar)
          }

          color <- NULL
          color[which(q.value >= 1.3)] <- "salmon"
          color[which(q.value < 1.3)] <- "black"
          temp <- barplot(q.value, horiz = TRUE, border = NA,
                          main = ifelse(missing(main), module.name, main),
                          col = color, cex.lab = 1.5, cex.axis = 1,
                          xlab = "-log10P-value(adjusted)", cex.main = 1.5, ...)
          idx1 <- which(q.value >= 1.3)
          idx2 <- which(q.value < 1.3)
          axis(side = 2, at = temp[idx1], labels = pathway.name[idx1],
               las = 2, cex.axis = cex.label,
               col.axis = "firebrick", tick = FALSE)
          axis(side = 2, at = temp[idx2], labels = pathway.name[idx2],
               las = 2, cex.axis = cex.label,
               col.axis = "black", tick = FALSE)

        }else{
          msea <- msea[order(msea$Overlap, decreasing = FALSE),]
          pathway.name <- rownames(msea)
          pathway.name <- unname(unlist(sapply(pathway.name, function(x){
            strsplit(x, split = ";")[[1]][1]
          })))
          overlap <- msea$Overlap

          max.len <- max(nchar(pathway.name))
          left.margin <- max.len * cex.label / 2.4

          if(missing(mar)) {
            par(mar = c(5,left.margin,4,2))
          }else{
            par(mar = mar)
          }

          color <- NULL
          color[which(overlap >= 0.1)] <- "salmon"
          color[which(overlap < 0.1)] <- "black"
          temp <- barplot(overlap*100, horiz = TRUE, border = NA,
                          main = ifelse(missing(main), module.name, main),
                          col = color, cex.lab = 1.5, cex.axis = 1,
                          xlab = "Overlap (%)", cex.main = 1.5, ...)
          idx1 <- which(overlap >= 0.1)
          idx2 <- which(overlap < 0.1)
          axis(side = 2, at = temp[idx1], labels = pathway.name[idx1],
               las = 2, cex.axis = cex.label,
               col.axis = "firebrick", tick = FALSE)
          axis(side = 2, at = temp[idx2], labels = pathway.name[idx2],
               las = 2, cex.axis = cex.label,
               col.axis = "black", tick = FALSE)
        }

      }
      # abline(v = 1.3, lty = 2, col = "tomato", lwd = 1.5)
    }
    }


)



