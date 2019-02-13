
#' @title readAnnotation
#' @description Read the data from ms2Annotation and seperate it to
#' tags and sample
#' @param data The file name (csv) from ms2Annotation
#' @param path The work directory
#' @param rt.filter Filter annotation according or not.
#' @param rt.tol The tolerance of RT (\%).
#' @param inHouse.compound in house compound database.
#' @author Xiaotao Shen
#' @export

setGeneric(name = "readAnnotation",
           def = function(data = "data.csv",
                          path = ".",
                          rt.filter = FALSE,
                          rt.tol = 30,
                          inHouse.compound = inhouse.compound){

             options(warn = -1)
             data <- readr::read_csv(file.path(path, data),
                                     col_types = readr::cols(),
                                     progress = TRUE)
             data <- as.data.frame(data)

             tags.idx <- match(c("name", "mz", "rt", "nhits.reverse",
                                 "hits.reverse", "nhits.forward",
                                 "hits.forward"), colnames(data))
             tags <- data[,tags.idx]
             sample <- data[,-tags.idx]

             rm(list = "data")
             gc()

             hit.reverse <- tags$hits.reverse
             idx1 <- which(!is.na(hit.reverse))

             hit.reverse1 <- hit.reverse[idx1]
             peak.rt <- tags$rt[idx1]
             peak.mz <- tags$mz[idx1]

             rm(list = c("hit.reverse"))
             gc()

             hit.reverse2 <- lapply(hit.reverse1,
                                    function(x) {strsplit(x, split = ";")[[1]]})

             rm(list = c("hit.reverse1"))
             gc()

             labid2 <- lapply(hit.reverse2, function(x) {
               y <- stringr::str_extract(string = x,
                                         pattern = "labid\\{[a-zA-Z0-9]+\\}")
               y <- gsub(pattern = "labid\\{", replacement = "", x = y)
               y <- gsub(pattern = "\\}", replacement = "", x = y)
               y
             })

             adduct2 <- stringr::str_extract_all(string = hit.reverse2,
                                                 pattern = "adduct\\{[\\(0-9a-zA-Z\\+\\-]+")
             adduct2 <-
               lapply(adduct2, function(x) {
                 gsub(pattern = "adduct\\{", replacement = "", x = x)})
             adduct2 <-
               lapply(adduct2, function(x) {
                 gsub(pattern = "\\(", replacement = "", x = x)})

             # ##remove M- adduct
             # mapply(function(x, y){
             # temp.index <- which(y != "M-")
             # if(length(temp.index) == 0) return(NA)
             # temp.index
             # },
             # x = labid2,
             # y = adduct2)

             ##remove wrong annotation, for examle n01170 (C18H22N2) can't
             ## -H2O

             right.idx <- mapply(function(x, y) {
               temp.idx <- match(x, inHouse.compound$Lab.ID)
               temp.formula <- inHouse.compound$Formula[temp.idx]
               temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
                                       x = temp.formula,
                                       y = y)
               which(!is.na(temp.formula1))
             },
             x = labid2,
             y = adduct2)

             ##fix bugs
             # for(i in 1:length(labid2)){
             #   cat(i); cat(" ")
             #   x <- labid2[[i]]
             #   y <- adduct2[[i]]
             #   temp.idx <- match(x, inHouse.compound$Lab.ID)
             #   temp.formula <- inHouse.compound$Formula[temp.idx]
             #   temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
             #                           x = temp.formula,
             #                           y = y)
             #
             # }


             remove.idx <- which(unlist(lapply(right.idx, length)) == 0)
             if(length(remove.idx) != 0){
             right.idx <- right.idx[-remove.idx]
             hit.reverse2 <- hit.reverse2[-remove.idx]
             labid2 <- labid2[-remove.idx]
             adduct2 <- adduct2[-remove.idx]
             idx1 <- idx1[-remove.idx]
             }

             peak.rt <- tags$rt[idx1]
             peak.mz <- tags$mz[idx1]
#------------------------------------------------------------------------------
             ##fix bugs
             # for(i in 1:length(labid2)){
             #   cat(i); cat(" ")
             #   x <- labid2[[i]]
             #   y <- adduct2[[i]]
             #   temp.idx <- match(x, inHouse.compound$Lab.ID)
             #   temp.formula <- inHouse.compound$Formula[temp.idx]
             #   temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
             #                           x = temp.formula,
             #                           y = y)
             # }
#------------------------------------------------------------------------------


             hit.reverse2 <- mapply(function(x, y) {
               x[y]
             },
             x = hit.reverse2,
             y = right.idx)

             labid2 <- mapply(function(x, y) {
               x[y]
             },
             x = labid2,
             y = right.idx)

             adduct2 <- mapply(function(x, y) {
               x[y]
             },
             x = adduct2,
             y = right.idx)

             rm(list = c("adduct2"))
             gc()

             if(!rt.filter){
               standard.rt <- lapply(labid2, function(x) {
                 rep(NA, length(x))
               })

               rt.error2 <- mapply(function(x,y) {
                 error <- abs(y - x)*100/y
               }, x = peak.rt, y = standard.rt)


               index <- mapply(function(x,y) {
                 idx <- abs(y - x)*100/y
                 idx[is.na(idx)] <- 0
                 idx <- sapply(idx, function(x) {ifelse(x > rt.tol, FALSE, TRUE)})
                 idx
               }, x = peak.rt, y = standard.rt)
             }else{
               standard.rt <- lapply(labid2, function(x) {
                 # temp.idx <- match(x, rt.all$`in house ID`)
                 temp.idx <- match(x, inHouse.compound$Lab.ID)
                 # temp.rt <- rt.all$`RT （min）`[temp.idx]
                 temp.rt <- inHouse.compound$RT[temp.idx]
                 temp.rt
               })

               rt.error2 <- mapply(function(x,y) {
                 error <- abs(y - x)*100/y
               }, x = peak.rt, y = standard.rt)

               index <- mapply(function(x,y) {
                 idx <- abs(y - x)*100/y
                 idx[is.na(idx)] <- rt.tol+1
                 idx <- sapply(idx, function(x) {ifelse(x > rt.tol, FALSE, TRUE)})
                 idx
               }, x = peak.rt, y = standard.rt)
             }

             rm(list = c("peak.rt"))
             gc()

             hit.reverse3 <-
               mapply(function(x,y){x[y]}, x = hit.reverse2, y = index)

             rm(list = c("hit.reverse2"))
             gc()

             hit.reverse4 <-
               lapply(hit.reverse3, function(x) {paste(x, collapse = ";")})
             hit.reverse4 <- unlist(hit.reverse4)

             rm(list = c("hit.reverse3"))
             gc()

             rt.error3 <- mapply(function(x,y){x[y]}, x = rt.error2, y = index)
             rm(list = c("rt.error2"))
             gc()
             rt.error4 <-
               lapply(rt.error3, function(x) {paste(x, collapse = ";")})
             rt.error4 <- unlist(rt.error4)
             rm(list = c("rt.error3"))
             gc()

             ms2.sim3 <-
               stringr::str_extract_all(string = hit.reverse4,
                                        pattern = "score\\{0\\.[0-9]+\\}|score\\{1\\}")
             ms2.sim3 <-
               lapply(ms2.sim3, function(x) {
                 gsub(pattern = "score\\{", replacement = "", x = x)})
             ms2.sim3 <-
               lapply(ms2.sim3, function(x) {
                 gsub(pattern = "\\}", replacement = "", x = x)})
             ms2.sim4 <-
               lapply(ms2.sim3, function(x) {paste(x, collapse = ";")})
             ms2.sim4 <- unlist(ms2.sim4)

             rm(list = c("ms2.sim3"))
             gc()

             adduct3 <- stringr::str_extract_all(string = hit.reverse4,
                                                 pattern = "adduct\\{[\\(0-9a-zA-Z\\+\\-]+")
             adduct3 <-
               lapply(adduct3, function(x) {
                 gsub(pattern = "adduct\\{", replacement = "", x = x)})
             adduct3 <-
               lapply(adduct3, function(x) {
                 gsub(pattern = "\\(", replacement = "", x = x)})
             adduct4 <- lapply(adduct3, function(x) {paste(x, collapse = ";")})
             adduct4 <- unlist(adduct4)

             name3 <-
               stringr::str_extract_all(string = hit.reverse4,
                                        pattern = "name\\{[^\\{]+\\}")

             rm(list = c("hit.reverse4"))
             gc()

             name3 <- lapply(name3, function(x) {
               gsub(pattern = "name\\{", replacement = "", x = x)})
             name3 <- lapply(name3, function(x) {
               gsub(pattern = "\\}", replacement = "", x = x)})
             name4 <- lapply(name3, function(x) {paste(x, collapse = ";")})
             name4 <- unlist(name4)
             rm(list = "name3")
             gc()

             labid3 <- mapply(function(x,y){x[y]}, x = labid2, y = index)
             rm(list = "labid2")
             gc()
             labid4 <- lapply(labid3, function(x) {paste(x, collapse = ";")})
             labid4 <- unlist(labid4)

             ##kegg.id
             KEGG.ID3 <- lapply(labid3, function(x){
               inHouse.compound$KEGG.ID[match(x, inHouse.compound$Lab.ID)]
             })

             ##if is more than 1 ID, remain first ID
             KEGG.ID3 <- lapply(KEGG.ID3, function(x){
               id.len <- nchar(x)
               id.len[is.na(id.len)] <- 6
               temp.idx <- which(id.len > 6)
               if(length(temp.idx) > 0){
                 x[temp.idx] <-
                   unlist(lapply(strsplit(x[temp.idx], split = "/"), function(x){
                     x[[1]]}
                   ))}
               x
             })

             KEGG.ID4 <- lapply(KEGG.ID3, function(x) {
               paste(x, collapse = ";")
             })
             KEGG.ID4 <- unlist(KEGG.ID4)
             KEGG.ID4[grep('NA', KEGG.ID4)] <- NA
             rm(list = c("KEGG.ID3"))
             gc()
             ###calculate mz error
             accurate.mass3 <- mapply(function(x, y) {
               temp.idx <- match(x, inHouse.compound$Lab.ID)
               temp.formula <- inHouse.compound$Formula[temp.idx]
               temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
                                       x = temp.formula,
                                       y = y)
               a.mass <- sapply(temp.formula1, function(x) {
                 Rdisop::getMass(Rdisop::getMolecule(formula = x))})
               a.mass
             },
             x = labid3,
             y = adduct3)

             rm(list = c("adduct3"))
             gc()
             # for(i in 1:length(labid3)){
             #   cat(i); cat(" ")
             # x <- labid3[[i]]
             # y <- adduct3[[i]]
             # temp.formula <- inHouse.compound$Formula[match(x, inHouse.compound$lab.ID)]
             # temp.formula1 <- mapply(function(x, y){sumFormula(x,y)},
             #                         x = temp.formula,
             #                         y = y)
             # a.mass <- sapply(temp.formula1, function(x) {
             #   Rdisop::getMass(Rdisop::getMolecule(formula = x))})
             # }


             mz.error3 <- mapply(function(x, y) {
               if(is.list(x)) x <- NA
               # error <- abs(x - y)*10^6/y
               error <- abs(x - y)*10^6/ifelse(y>=400,y,400)
               unname(error)
             },
             x = accurate.mass3,
             y = peak.mz)

             rm(list = c("peak.mz"))
             gc()
             ##should be fixed
             mz.error3 <- lapply(mz.error3, function(x){
               x[!is.na(x)] <- 10
               x
             })

             mz.error4 <-
               lapply(mz.error3, function(x) {paste(x, collapse = ";")})
             mz.error4 <- unlist(mz.error4)
             mz.error4[mz.error4=="NA"] <- ""


             #formula
             Formula3 <- lapply(labid3, function(x){
               inHouse.compound$Formula[match(x, inHouse.compound$Lab.ID)]
             })

             rm(list = c("labid3"))
             gc()

             Formula4 <- lapply(Formula3,  function(x) {
               paste(x, collapse = ";")
             })
             Formula4 <- unlist(Formula4)

             tags <- as.data.frame(tags)
             tags <- tags[,c("name","mz","rt")]
             colnames(tags) <- c("Peak.name","mz","rt")
             ms2.sim <- rep(NA, nrow(tags))
             adduct <- rep(NA, nrow(tags))
             rt.error <- rep(NA, nrow(tags))
             mz.error <- rep(NA, nrow(tags))
             Metabolite.name <- rep(NA, nrow(tags))
             KEGG.ID <- rep(NA, nrow(tags))
             Formula <- rep(NA, nrow(tags))
             isotope <- rep(NA, nrow(tags))
             labid <- rep(NA, nrow(tags))

             ms2.sim[idx1] <- ms2.sim4
             adduct[idx1] <- adduct4
             rt.error[idx1] <- rt.error4
             mz.error[idx1] <- mz.error4
             Metabolite.name[idx1] <- name4
             KEGG.ID[idx1] <- KEGG.ID4
             Formula[idx1] <- Formula4
             labid[idx1] <- labid4
             isotope[idx1] <- "[M]"

             rm(list = c("labid4", "adduct4", "ms2.sim4", "rt.error4",
                         "mz.error4", "name4", "KEGG.ID4", "Formula4"))
             gc()

             tags <- data.frame(tags, adduct, isotope, ms2.sim, rt.error, mz.error,
                                Formula, Metabolite.name, KEGG.ID, labid,
                                stringsAsFactors = FALSE)
             rm(list = c("adduct", "isotope", "ms2.sim", "rt.error", "mz.error",
                         "Formula", "Metabolite.name", "KEGG.ID", "labid"))
             gc()

             tags[which(tags == "", arr.ind = TRUE)] <- NA

             ##change tags to tags2 style
             result <- list(tags, sample)
             names(result) <- c("tags", "sample")
             result <- result
           })






################################################################################

setGeneric(name = "removeByRT",
           def = function(data,
                          rt.data,
                          rt.tol = 10,#second
                          direction = c("less", "bigger")
                          ){
             direction <- match.arg(direction)
             tags.idx <- match(c("name", "mz", "rt", "nhits.reverse",
                                 "hits.reverse", "nhits.forward",
                                 "hits.forward"), colnames(data))
             tags <- data[,tags.idx]
             sample <- data[,-tags.idx]

             #####reverse
             hit.reverse <- tags$hits.reverse
             idx.reverse <- which(!is.na(hit.reverse))

             new.hit.reverse <- hit.reverse[idx.reverse]
             reverse.rt <- tags$rt[idx.reverse]

             new.hit.reverse <- lapply(new.hit.reverse,
                                    function(x) {strsplit(x, split = ";")[[1]]})

             reverse.labid <- lapply(new.hit.reverse, function(x) {
               y <- stringr::str_extract(string = x,
                                         pattern = "labid\\{[a-zA-Z0-9]+\\}")
               y <- gsub(pattern = "labid\\{", replacement = "", x = y)
               y <- gsub(pattern = "\\}", replacement = "", x = y)
               y
             })

             standard.reverse.rt <- lapply(reverse.labid, function(x){
               rt.data$RT[match(x, rt.data$labID)]
             })


             reverse.rt.error <- mapply(function(x, y){
               abs(x - y)
             },
             x = reverse.rt,
             y = standard.reverse.rt)


             if(direction == "less"){
               temp.index <- lapply(reverse.rt.error, function(x){
                 x[is.na(x)] <- rt.tol + 100
                 which(x <= rt.tol)
               })
             }else{
               temp.index <- lapply(reverse.rt.error, function(x){
                 x[is.na(x)] <- rt.tol - 100
                 which(x > rt.tol)
               })
             }



             new.hit.reverse <- mapply(function(x, y){
             if(length(y) == 0) return(NA)
             x[y]
             },
             x = new.hit.reverse,
             y = temp.index)

             new.hit.reverse <- unlist(lapply(new.hit.reverse, function(x){
             paste(x, collapse = ";")
             }))

             new.hit.reverse[new.hit.reverse == "NA"] <- NA
             hit.reverse[idx.reverse] <- new.hit.reverse



             ##############################################################
             #####forward
             hit.forward <- tags$hits.forward
             idx.forward <- which(!is.na(hit.forward))

             new.hit.forward <- hit.forward[idx.forward]
             forward.rt <- tags$rt[idx.forward]

             new.hit.forward <- lapply(new.hit.forward,
                                       function(x) {strsplit(x, split = ";")[[1]]})

             forward.labid <- lapply(new.hit.forward, function(x) {
               y <- stringr::str_extract(string = x,
                                         pattern = "labid\\{[a-zA-Z0-9]+\\}")
               y <- gsub(pattern = "labid\\{", replacement = "", x = y)
               y <- gsub(pattern = "\\}", replacement = "", x = y)
               y
             })

             standard.forward.rt <- lapply(forward.labid, function(x){
               rt.data$RT[match(x, rt.data$labID)]
             })


             forward.rt.error <- mapply(function(x, y){
               abs(x - y)
             },
             x = forward.rt,
             y = standard.forward.rt)


             if(direction == "less"){
               temp.index <- lapply(forward.rt.error, function(x){
                 x[is.na(x)] <- rt.tol + 100
                 which(x <= rt.tol)
               })
             }else{
               temp.index <- lapply(forward.rt.error, function(x){
                 x[is.na(x)] <- rt.tol - 100
                 which(x > rt.tol)
               })
             }


             new.hit.forward <- mapply(function(x, y){
               if(length(y) == 0) return(NA)
               x[y]
             },
             x = new.hit.forward,
             y = temp.index)

             new.hit.forward <- unlist(lapply(new.hit.forward, function(x){
               paste(x, collapse = ";")
             }))

             new.hit.forward[new.hit.forward == "NA"] <- NA
             hit.forward[idx.forward] <- new.hit.forward


             tags$hits.reverse <- hit.reverse
             tags$hits.forward <- hit.forward

             data <- data.frame(tags, sample, stringsAsFactors = FALSE)
             return(data)
})





################################################################################

setGeneric(name = "removeByDP",
           def = function(data,
                          dp.cutoff = 0.8,#second
                          direction = c("less", "bigger")
           ){
             direction <- match.arg(direction)
             tags.idx <- match(c("name", "mz", "rt", "nhits.reverse",
                                 "hits.reverse", "nhits.forward",
                                 "hits.forward"), colnames(data))
             tags <- data[,tags.idx]
             sample <- data[,-tags.idx]

             #####reverse
             hit.reverse <- tags$hits.reverse

             reverse.dp <- lapply(hit.reverse, function(x){
               if(is.na(x)) return(NA)
               temp.score <- stringr::str_extract_all(string = x, pattern = "score\\{[0-9\\.]{0,30}\\}")
               unlist(lapply(temp.score, function(x){
                 as.numeric(gsub(pattern = "score\\{|\\}", replacement = "", x = x))
               }))
             })



             new.hit.reverse <- lapply(hit.reverse,
                                       function(x) {strsplit(x, split = ";")[[1]]})



             reverse.index <- lapply(reverse.dp, function(x){
               if(is.na(x)) return(NA)
               if(direction == "less"){
               temp.idx <- which(x <= dp.cutoff)
               }else{
                 temp.idx <- which(x >= dp.cutoff)
               }
               if(length(temp.idx)==0) return(NA)
               temp.idx
             })

             new.hit.reverse <- mapply(function(x, y){
               if(is.na(y)) return(NA)
               x[y]
             },
             x = new.hit.reverse,
             y = reverse.index)

             new.hit.reverse <- unlist(lapply(new.hit.reverse, function(x){
               paste(x, collapse = ";")
             }))

             new.hit.reverse[new.hit.reverse == "NA"] <- NA



             ##############################################################
             #####forward
             hit.forward <- tags$hits.forward

             forward.dp <- lapply(hit.forward, function(x){
               if(is.na(x)) return(NA)
               temp.score <- stringr::str_extract_all(string = x, pattern = "score\\{[0-9\\.]{0,30}\\}")
               unlist(lapply(temp.score, function(x){
                 as.numeric(gsub(pattern = "score\\{|\\}", replacement = "", x = x))
               }))
             })



             new.hit.forward <- lapply(hit.forward,
                                       function(x) {strsplit(x, split = ";")[[1]]})



             forward.index <- lapply(forward.dp, function(x){
               if(is.na(x)) return(NA)
               if(direction == "less"){
                 temp.idx <- which(x <= dp.cutoff)
               }else{
                 temp.idx <- which(x >= dp.cutoff)
               }
               if(length(temp.idx)==0) return(NA)
               temp.idx
             })

             new.hit.forward <- mapply(function(x, y){
               if(is.na(y)) return(NA)
               x[y]
             },
             x = new.hit.forward,
             y = forward.index)

             new.hit.forward <- unlist(lapply(new.hit.forward, function(x){
               paste(x, collapse = ";")
             }))

             new.hit.forward[new.hit.forward == "NA"] <- NA


             tags$hits.reverse <- new.hit.reverse
             tags$hits.forward <- new.hit.forward

             data <- data.frame(tags, sample, stringsAsFactors = FALSE)
             return(data)
           })


