# title isotopeAnnotation
# description Find the isotopes of a known metabolite according to rt, mz
#  and intensity.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param id The KEGG ID of metabolite.
# param formula The formula of metabolite.
# param adduct The adduct of metabolite.
# param charge The charge of metabolite.
# param mz The mz of metabolite.
# param rt The rt of metabolite.
# param int The mean/median intensity of metabolite accross the samples.
# param peak.mz The mz of peaks.
# param peak.rt The rt of peaks.
# param peak.int The mean/median intensity of peaks accross the samples.
# param cor The correlations between metabolite and peaks.
# param rt.tol The rt tolerance (second).
# param mz.tol The mz tolerance.
# param cor.tol The cor tolerance.
# param int.tol The intensity ratio tolerance.
# param max.isotope The max number of isotopes.
# return  Isotope annotation information of peaks.



##test data
# setwd('/home/jasper/work/MetDNA/isotopeAnnotation')
# id = NA
# formula = "C21H27N7O14P2"
# adduct = "M+H"
# charge = 1
# mz = 664.1139
# rt = 812.143
# int = 205109.4
# data <- readr::read_csv("annotation.result.csv")
# data <- as.data.frame(data)
# sample <- data[,grep("Sample", colnames(data))]
# tags <- data[,-grep("Sample", colnames(data))]
# peak.mz <- tags$mz
# peak.rt <- tags$rt
# peak.int <- apply(sample, 1, median)
# cor <- rep(1, nrow(tags))
# rt.tol = 5
# or.tol = 0
# int.tol = 500
# max.isotope = 4
#
# system.time(a <- isotopeAnnotation(peak.mz = peak.mz, peak.rt = peak.rt, peak.int = peak.int, cor = cor))
# system.time(b <- isotopeAnnotation1(peak.mz = peak.mz, peak.rt = peak.rt, peak.int = peak.int, cor = cor))





#backup of old version
setGeneric(name = "isotopeAnnotation",
           def = function(# metabolite information
             id = NA,
             formula = "C21H27N7O14P2",
             adduct = "M+H",
             charge = 1,
             mz = 664.1139,
             rt = 812.143,
             int = 205109.4,
             ## peak information
             peak.mz,
             peak.rt,
             peak.int,
             cor,
             ## other parameters
             rt.tol = 5,
             mz.tol = 25,
             cor.tol = 0,
             int.tol = 500,
             max.isotope = 4){

             ##Nicotinamide adenine dinucleotide (NAD) as the example
             formula1 <- sumFormula(formula = formula,adduct = adduct)
             ###should be fix latter
             if(is.na(formula1)) formula1 <- formula
             molecule <- Rdisop::getMolecule(formula = formula1,
                                             z = charge,
                                             maxisotopes = max.isotope + 1)
             isotopes <- t(Rdisop::getIsotope(molecule = molecule))
             rownames(isotopes) <-
               c("[M]",paste("[M","+",c(1:(nrow(isotopes)-1)),"]", sep = ""))
             isotopes <- data.frame(isotopes, rownames(isotopes),
                                    stringsAsFactors = FALSE)
             colnames(isotopes) <- c("mz", "intensity", "isotope")
             accurate.mz <- mz
             peak.int <- peak.int/int
             ##if int is 0
             peak.int[is.na(peak.int)] <- 1
             peak.int[is.nan(peak.int)] <- 1
             peak.int[is.infinite(peak.int)] <- 100000

             ###rt and cor filtering
             rt.error <- abs(rt - peak.rt)
             index1 <- which(rt.error <= rt.tol & cor >= cor.tol)
             if(length(index1) == 0) return(NULL)
             #all the peaks are filtered using rt and names as 1
             peak.mz1 <- peak.mz[index1]
             peak.rt1 <- peak.rt[index1]
             peak.int1 <- peak.int[index1]
             peak.int1 <- peak.int[index1]
             cor1 <- cor[index1]

             #iso.idx is the index of peak for isotopes
             iso.idx <- NULL
             mz.error <- NULL
             int.error <- NULL
             iso <- NULL
             correlation <- NULL

             for(i in 2:nrow(isotopes)){
               ##mz and intensity infromation
               temp.mz <- as.numeric(isotopes[i,"mz"])
               temp.int <-
                 as.numeric(isotopes[i,"intensity"])/as.numeric(isotopes[1,"intensity"])

               ## calculate error
               # peak.mz.error <- abs(temp.mz - peak.mz1)*10^6/temp.mz
               peak.mz.error <- abs(temp.mz - peak.mz1)*10^6/ifelse(temp.mz>=400,temp.mz,400)
               peak.int.error <-
                 abs(temp.int - peak.int1)*100/temp.int

               ##has peak matched the mz tolerance
               idx <- which(peak.mz.error <= mz.tol)
               ###idx=0, no peaks matched
               if(length(idx) == 0) {
                 iso.idx[i] <- NA
                 mz.error[i] <- NA
                 int.error[i] <- NA
                 correlation[i] <- NA
                 iso[i] <- NA
               }
               ## one peak matched, it is
               if(length(idx) == 1) {
                 iso.idx[i] <- idx
                 mz.error[i] <- peak.mz.error[idx]
                 int.error[i] <- peak.int.error[idx]
                 correlation[i] <- cor1[idx]
                 iso[i] <- isotopes[i,3]
               }
               ## more than one matched, see the intensity ratio error
               if(length(idx) > 1) {
                 idx <- idx[which.min(peak.int.error[idx])]
                 iso.idx[i] <- idx
                 mz.error[i] <- peak.mz.error[idx]
                 int.error[i] <- peak.int.error[idx]
                 correlation[i] <- cor1[idx]
                 iso[i] <- isotopes[i,3]
               }
             }
             #
             ##--------------------------------------------------------------------------
             if(all(is.na(iso.idx))) {
               return(NULL)
             }else{
               iso.idx <- iso.idx[!is.na(iso.idx)]
               index2 <- index1[iso.idx]
               mz.error <- mz.error[!is.na(mz.error)]
               int.error <- int.error[!is.na(int.error)]
               correlation <- correlation[!is.na(correlation)]
               rt.error <- rt.error[index2]
               iso <- iso[!is.na(iso)]
               iso.info <- data.frame(index2, rt.error, mz.error,
                                      int.error, correlation,
                                      iso, stringsAsFactors = FALSE)
             }
             colnames(iso.info) <- c("peakIndex", "rtError.s",
                                     "mzError.ppm", "IntensityRatioError.%",
                                     "correlation", "isotopes")
             ##remove peaks which have large intensity ratio error. For M+1 isotope, the
             ##tolerance is set as 30%, and for other isotopes, the tolerance is set as
             ##100%.
             remove.idx1 <-
               which(iso.info$isotopes == "[M+1]" & iso.info$"IntensityRatioError.%" > 500)

             remove.idx2 <-
               which(iso.info$isotopes != "[M+1]" & iso.info$"IntensityRatioError.%" > 500)

             remove.idx <- c(remove.idx1, remove.idx2)
             if(length(remove.idx) != 0){
               iso.info <- iso.info[-c(remove.idx1, remove.idx2),]
             }
             if(nrow(iso.info) == 0) {return(NULL)}
             iso.info <- data.frame(id, iso.info, stringsAsFactors = FALSE)
             rm(list = c("peak.mz", "peak.rt", "peak.int", "cor"))
             gc()
             iso.info <- iso.info
           })










#-------------------------------------------------------------------------------
#' @title adductAnnotation
#' @description adduct annotation
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#'@param id The kegg ID of metabolite.
#'@param formula The formula of metabolite.
#'@param adduct The adduct of metabolite.
#'@param polarity The mode of data.
#'@param mz The mz of metabolite.
#'@param rt The rt of metabolite.
#'@param adduct.table The adduct table.
#'@param peak.mz The mz of peaks.
#'@param peak.rt The rt of peaks.
#'@param cor The correlations between metabolite and peaks.
#'@param mz.tol The mz tol.
#'@param rt.tol The rt tol (second).
#'@param cor.tol The correlation tol.
#'@return  Isotope annotation information of peaks.
#'@export

#test data
# setwd('/home/jasper/work/MetDNA/adductAnnotation')
# id = NA
# formula = "C6H14N4O2"
# adduct = "M+H"
# polarity = "positive"
# mz = 175.118
# rt = 961.6225
# data <- readr::read_csv("annotation.result.csv")
# data <- as.data.frame(data)
# sample <- data[,grep("Sample", colnames(data))]
# tags <- data[,-grep("Sample", colnames(data))]
# peak.mz <- tags$mz
# peak.rt <- tags$rt
# cor <- rep(1, nrow(tags))
# mz.tol = 25
# rt.tol = 3
# cor.tol = 0
# load("adduct.table.hilic.rda")
# adduct.table <- adduct.table.hilic
#
# system.time(a <- adductAnnotation(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor))
# system.time(b <- adductAnnotation1(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor))


# system.time(lapply(1:100, function(x) adductAnnotation(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor)))
# system.time(lapply(1:100, function(x) adductAnnotation1(adduct.table = adduct.table, peak.mz = peak.mz, peak.rt = peak.rt, cor = cor)))

setGeneric(name = "adductAnnotation",
           def = function(
             # metabolite information
             id = NA,
             formula = "C6H14N4O2",
             adduct = "M+H",
             polarity = c("positive", "negative"),
             mz = 175.118,
             rt = 961.6225,
             adduct.table,
             ## peak information
             peak.mz,
             peak.rt,
             cor,
             ## other parameters
             mz.tol = 25,
             rt.tol = 3,
             cor.tol = 0.5){

             polarity = match.arg(polarity)

             if(polarity == "positive") {
               adduct.table <- adduct.table[adduct.table[,"mode"] == "positive",]
             }else{
               adduct.table <- adduct.table[adduct.table[,"mode"] == "negative",]
             }

             ###remove metabolite adduct
             idx <- which(adduct.table[,1] == adduct)
             if(length(idx) == 1){
               adduct.table <- adduct.table[-idx,]
             }

             if(adduct != "M+") {
               adduct.table <- adduct.table[which(adduct.table[,1] != "M+"),]
             }
             adduct.name <- as.character(adduct.table[,1])
             adduct.charge <- as.numeric(adduct.table[,3])

             ##remove invalid
             remain.index <- which(sapply(adduct.name, function(x) {
               checkElement(formula, x)}))

             # for(i in 1:length(adduct.name)){
             #   checkElement(formula, adduct.name[i])
             # }

             if(length(remain.index) == 0) return(NULL)

             adduct.name <- adduct.name[remain.index]
             adduct.charge <- adduct.charge[remain.index]

             ##formula1 is the new adduct formula
             formula1 <- sapply(adduct.name, function(x) {sumFormula(formula,x)})

             # ##if the metabolite don't have enough element to remove, formula is NA.
             charge <- adduct.charge[!is.na(formula1)]
             formula1 <- formula1[!is.na(formula1)]

             ##get the accurate mass for different adduct format
             accurate.mass <- sapply(formula1, function(x) {Rdisop::getMolecule(x)$exactmass})
             accurate.mz <- accurate.mass/charge
             ##adduct.info is the adduct table of adduct
             adduct.info <- data.frame(accurate.mz, rt, charge,
                                       names(accurate.mz),
                                       formula, formula1,
                                       stringsAsFactors = FALSE)
             colnames(adduct.info) <-
               c("mz", "rt","charge", "adduct" , "Formula", "Adduct.Formula")

             ###rt mz, and cor filtering
             rt.error <- lapply(adduct.info$rt, function(x){
               abs(x - peak.rt)
             })

             mz.error <- lapply(adduct.info$mz, function(x){
               # abs(x - peak.mz)*10^6/x
               abs(x - peak.mz)*10^6/ifelse(x>=400,x,400)
             })

             index1 <- mapply(function(x, y){
               temp.idx <- which(x <= rt.tol & y <= mz.tol)
             },
             x = rt.error,
             y = mz.error)


             rm(list = c("peak.mz", "peak.rt", "cor"))
             gc()
             if(all(unlist(lapply(index1, length)) == 0)) return(NULL)

             add.info <- mapply(function(temp.idx, temp.mz.error, temp.rt.error){
               if(length(temp.idx) == 0) return(list(rep(NA, 4)))
               if(length(temp.idx) == 1) {
                 list(c(temp.idx, temp.rt.error[temp.idx], temp.mz.error[temp.idx],
                        1))
               }else{
                 temp.idx <- temp.idx[which.min(temp.rt.error[temp.idx])][1]
                 list(c(temp.idx, temp.rt.error[temp.idx], temp.mz.error[temp.idx],
                        1))
               }

             },
             temp.idx = index1,
             temp.mz.error = mz.error,
             temp.rt.error = rt.error)

             rm(list = c("mz.error", "rt.error"))
             gc()

             add.info <- do.call(rbind, add.info)
             add.info <- data.frame(add.info, adduct.info$adduct,
                                    stringsAsFactors = FALSE)




             colnames(add.info) <-
               c("peakIndex", "rtError.s", "mzError.ppm", "correlation", "adducts")

             add.info <- add.info[which(apply(add.info, 1,
                                              function(x) all(!is.na(x)))),,drop = FALSE]
             add.info <- data.frame(id, add.info, stringsAsFactors = FALSE)

             add.info <- add.info

           })








########################------------------------------------------------------
#'@title metAnnotation
#'@description Annotate peak table from one metabolite.
#'@author Xiaotao Shen
#'\email{shenxt@@sioc.ac.cn}
#'@param metabolite.name The metabolite name.
#'@param metabolite.id The metabolite ID.
#'@param formula The formula of metabolite.
#'@param adduct The adduct of metabolite.
#'@param polarity The polarity.
#'@param mz The mz of metabolite.
#'@param rt The RT of metabolite.
#'@param peak.mz The mz of all peaks.
#'@param peak.rt The RT of all peaks.
#'@param cor The correlation of metabolite between all peaks.
#'@param ms2 The ms2 data of peak table.
#'@param mz.tol The mz tolerance.
#'@param rt.tol The RT tolerance for metabolite annotation. (\%)
#'@param cor.tol The cor tolerance.
#'@param dp.tol The tolerance of dot product.
#'@param step The reaction step.
#'@param metabolite The kegg compound database.
#'@param metabolic.network kegg.rpair2
#'@param adduct.table Adduct table.
#'@return  The metabolite annotation information.
#'@export

##test data
# setwd('/home/jasper/work/MetDNA/metAnnotation')
# metabolite.name = "M175T962"
# metabolite.id = "C00062"
# formula = "C6H14N4O2"
# adduct = "M+H"
# polarity = "positive"
# mz = 175.118
# rt = 961.6225
# ## peak information
# load("ms2")
# ## other parameters
# mz.tol = 25
# #r tol is relative(%)
# rt.tol = 30
# cor.tol = 0
# dp.tol = 0.5
# step = 1
# load("kegg.compound.rda")
# load("kegg.rpair2.rda")
# metabolite <- kegg.compound
# metabolic.network <- kegg.rpair2
#
# data <- readr::read_csv("annotation.result.csv")
# data <- as.data.frame(data)
# sample <- data[,grep("Sample", colnames(data))]
# tags <- data[,-grep("Sample", colnames(data))]
# peak.mz <- tags$mz
# peak.rt <- tags$rt
# names(peak.rt) <- names(peak.mz) <- tags$name
# cor <- rep(1, nrow(tags))
# load("adduct.table.hilic.rda")
# adduct.table <- adduct.table.hilic
#
#
#
#
# system.time(a <- metAnnotation(polarity = "positive",
#                                 peak.mz = peak.mz, peak.rt = peak.rt,
#                                 cor = cor, ms2 = ms2, metabolite = kegg.compound,
#                                 metabolic.network = kegg.rpair2,
#                                 adduct.table = adduct.table))
#
# a <- a[order(a$peakName),]
#
# system.time(b <- metAnnotation1(polarity = "positive",
#                                peak.mz = peak.mz, peak.rt = peak.rt,
#                                cor = cor, ms2 = ms2, metabolite = kegg.compound,
#                                metabolic.network = kegg.rpair2,
#                                adduct.table = adduct.table))

setGeneric(name = "metAnnotation",
           def = function(metabolite.name = "M175T962",
                          metabolite.id = "C00062",
                          formula = "C6H14N4O2",
                          adduct = "M+H",
                          polarity = c("positive", "negative"),
                          mz = 175.118,
                          rt = 961.6225,
                          ## peak information
                          peak.mz,
                          peak.rt,
                          cor,
                          ms2,
                          ## other parameters
                          mz.tol = 25,
                          #r tol is relative(%)
                          rt.tol = 30,
                          cor.tol = 0,
                          dp.tol = 0.5,
                          step = 1,
                          metabolite,
                          metabolic.network,
                          adduct.table){
             polarity <- match.arg(polarity)
             adduct.table <- adduct.table[adduct.table[,"mode"]==polarity,]

             met.result <- getNeighbor(metabolite.id = metabolite.id,
                                       step = step,
                                       graph = metabolic.network)
             if(is.null(met.result)) return(NULL)
             ##should be fixed
             if(nrow(met.result) > 100){met.result <- met.result[1:100,]}

             ##remove the metabolite which has no standard formula
             remove.idx <-
               which(is.na(metabolite$Exact.mass[match(met.result$Node.ID, metabolite$ID)]))
             if(length(remove.idx) != 0) met.result <- met.result[-remove.idx,,drop = FALSE]
             if(nrow(met.result) == 0) return(NULL)

             ### add RT information to met.result
             temp.idx <- match(met.result$Node.ID, metabolite$ID)

             met.result <-
               data.frame(met.result,
                          metabolite[temp.idx,c("RT", "attribute")],
                          stringsAsFactors = FALSE)

             # met.result <-
             #   data.frame(met.result,
             #              metabolite[temp.idx,c("Amide23.RT", "Amide23.RT.attribute")],
             #              stringsAsFactors = FALSE)
             # colnames(met.result)[4:5] <- c("RT", "attribute")
             # met.result$RT <- met.result$RT*60


             ##only the peaks with MS2 are remainded
             peak.name <- names(peak.rt)
             raw.peak.name <- peak.name
             ms2.name <- unname(unlist(lapply(ms2, function(x) x[[1]][1,1])))
             temp.index <- match(ms2.name, peak.name)
             if(length(temp.index) == 0) return(NULL)

             peak.mz <- peak.mz[temp.index]
             peak.rt <- peak.rt[temp.index]
             peak.name <- peak.name[temp.index]
             cor <- cor[temp.index]

             # rm(list = c("peak.mz", "peak.rt", "cor", "peak.name"))

             ###using RT to remove some peaks
             temp <- sapply(met.result$RT, function(x) {abs(peak.rt - x)*100/(x)})
             temp.index <- which(apply(temp, 1, function(x) {any(x <= rt.tol)}))
             temp.index2 <- which(apply(temp, 2, function(x) {any(x <= rt.tol)}))

             if(length(temp.index) == 0 | length(temp.index2) == 0) return(NULL)

             ###peak.mz and rt and cor are peak after removing using RT
             peak.mz <- peak.mz[temp.index]
             peak.rt <- peak.rt[temp.index]
             cor <- cor[temp.index]
             peak.name <- peak.name[temp.index]

             met.result <- met.result[temp.index2,,drop = FALSE]

             # rm(list = c("peak.mz1", "peak.rt1", "cor1", "peak.name1", "met.result"))

             ## add adduct to neighbor
             met.result <- apply(met.result, 1, list)

             new.met.result <- lapply(met.result, function(x){
               x <- x[[1]]
               node.step <- x["Step"]
               met.id <- x["Metabolite.ID"]
               node.id <- x["Node.ID"]
               node.rt <- x["RT"]
               node.rt.attribute <- x["attribute"]

               met.mz <-
                 as.numeric(metabolite[,"Exact.mass"][match(node.id, metabolite[,"ID"])])
               node.formula <-
                 as.character(metabolite[,"Formula"][match(node.id, metabolite[,"ID"])])
               node.adduct <- as.character(adduct.table[,1])

               temp.idx <- which(sapply(node.adduct, function(x) {checkElement(node.formula, x)}))
               node.adduct <- as.character(adduct.table[temp.idx,1])
               node.polymer <- as.numeric(adduct.table[temp.idx,2])
               node.charge <- as.numeric(adduct.table[temp.idx,"charge"])
               node.mz <- (met.mz*node.polymer + adduct.table[temp.idx,"massdiff"])/node.charge

               temp <-
                 data.frame(node.step, met.id, node.id, node.adduct,
                            node.charge, node.mz,  node.rt,
                            node.rt.attribute, stringsAsFactors = FALSE)
               temp
             })


             new.met.result <- do.call(rbind, new.met.result)

             # rm(met.result2)

             ##use mz and rt to remove peaks and nodes
             mz.error <- lapply(new.met.result$node.mz, function(x){
               # abs(as.numeric(x) - peak.mz)*10^6/x
               abs(as.numeric(x) - peak.mz)*10^6/ifelse(x>=400,x,400)
             })

             rt.error <- lapply(new.met.result$node.rt, function(x){
               abs(as.numeric(x) - peak.rt)*100/as.numeric(x)
             })

             index <- mapply(function(x, y){
               which(x <= mz.tol & y <= rt.tol)
             },
             x = mz.error,
             y = rt.error)

             if(all(unlist(lapply(index, length)) == 0)) return(NULL)

             new.met.result <- new.met.result[which(unlist(lapply(index, length)) != 0),]
             index <- index[which(unlist(lapply(index, length)) != 0)]
             temp.idx <- unique(unlist(index))

             peak.mz <- peak.mz[temp.idx]
             peak.rt <- peak.rt[temp.idx]
             cor <- cor[temp.idx]
             peak.name <- peak.name[temp.idx]


             ###the last calculation
             mz.error <- lapply(new.met.result$node.mz, function(x){
               # abs(as.numeric(x) - peak.mz)*10^6/x
               abs(as.numeric(x) - peak.mz)*10^6/ifelse(x>=400,x,400)
             })

             rt.error <- lapply(new.met.result$node.rt, function(x){
               abs(as.numeric(x) - peak.rt)*100/as.numeric(x)
             })


             metabolite.ms2 <- ms2[[match(metabolite.name, ms2.name)]]

             peak.ms2 <- ms2[match(names(peak.rt), ms2.name)]


             stru.sim <-
               mapply(function(x,y){
                 if(y >= mz){
                   as.numeric(IdentifyFeature(x$spec, metabolite.ms2$spec,
                                              direction = "reverse"))
                 }else{
                   as.numeric(IdentifyFeature(metabolite.ms2$spec, x$spec,
                                              direction = "reverse"))
                 }

               },
               peak.ms2,
               peak.mz
               )

             stru.sim[is.nan(stru.sim)] <- 0
             stru.sim[is.na(stru.sim)] <- 0


             dp <- vector(mode = "list", length = length(mz.error))
             dp <- lapply(dp, function(x){
               x <- stru.sim
               x
             })


             index <- mapply(function(x, y, z){
               which(x <= mz.tol & y <= rt.tol & z >= dp.tol)
             },
             x = mz.error,
             y = rt.error,
             z = dp)

             if(all(unlist(lapply(index, length)) == 0)) return(NULL)
             # new.met.result <- new.met.result[which(unlist(lapply(index, length)) != 0),]
             new.met.result <- apply(new.met.result, 1, list)
             # index <- index[which(unlist(lapply(index, length)) != 0)]
             # names(index) <- new.met.result$node.id

             # peak.name <- names(peak.mz)
             metabolite.info <- mapply(function(temp.idx, temp.mz.error, temp.rt.error,
                                                temp.met.result, temp.dp){
               if(length(temp.idx) == 0) {
                 return(list(rep(NA, 13)))
               }else{
                 list(data.frame(temp.idx, temp.mz.error[temp.idx], temp.rt.error[temp.idx],
                                 1, temp.dp[temp.idx],
                                 matrix(rep(temp.met.result[[1]], length(temp.idx)), nrow = length(temp.idx), byrow = TRUE),
                                 stringsAsFactors = FALSE))
               }
             },
             temp.idx = index,
             temp.mz.error = mz.error,
             temp.rt.error = rt.error,
             temp.met.result = new.met.result,
             temp.dp = dp)


             metabolite.info <- do.call(rbind, metabolite.info)
             metabolite.info <- metabolite.info[which(apply(metabolite.info, 1,
                                                            function(x) all(!is.na(x)))),,drop = FALSE]

             colnames(metabolite.info) <-
               c("peakIndex","mzError.ppm", "rtError","correlation", "dotProduct",
                 "step","node1ID","peakID","adduct","charge",
                 "theoreticalMZ","RT","attribute")

             peakIndex <- metabolite.info$peakIndex
             peakName <- names(peak.mz)[peakIndex]
             peakMz <- peak.mz[peakIndex]
             peakRT <- peak.rt[peakIndex]
             peakIndex <- match(peakName, raw.peak.name)

             metabolite.info$peakIndex <- peakIndex

             metabolite.info <- data.frame(peakName, peakMz, peakRT, metabolite.info,
                                           stringsAsFactors = FALSE)

             metabolite.info <- metabolite.info[,c(1,4,2,3,5:16)]

             metabolite.info <- metabolite.info

           })




