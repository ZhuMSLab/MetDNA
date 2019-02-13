#' @title metABM
#' @description Annotate peaks based on metabolic reaction network (MRN).
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param annotation.result The name of annotation result from ms2Annotation.
#' @param ms2.data The name of ms2 data from ms2Annotation.
#' @param prefer.adduct Which adduct you want to use for RT prediction.
#' @param use.default.md Use default molecular descriptors for RT prediction or not.
#' @param column hilic or rp.
#' @param polarity positive or negative.
#' @param threads The number of threads.
#' @param path The work directory.
#' @param output.path The directory for output results.
#' @param max.isotope The number of isotope peaks
#' @param mz.tol The m/z tolerance. Default is 25 ppm.
#' @param rt.tol1 The RT tolerance for isotope and adduct annotation. (second)
#' Default is 3.
#' @param rt.tol2 The RT tolerance for metabolite annotation. (\%) Default
#' is 30\%.
#' @param cor.tol The correlation tolerance.
#' @param int.tol The intensity ratio tolerance of isotope annotation.
#' @param dp.tol The tolerance of dot product.
#' @param max.step The max number of reaction step.
#' @param score.cutoff Score cutoff of annotations.
#' @param remain Remain some seeds as validation or not.
#' @param remain.per The percentage of remained seeds. Default is 50\%.
#' @param seed.neighbor.match.plot Output MS/MS match plot of seed and neighbor or not.
#' @param candidate.num How many candidates for peaks are outputted. Default is 3.
#' @return Annotation result.
#' @export
#'

# annotation.result = "annotation.result.csv"
# ms2.data = "ms2"
# prefer.adduct = "M+H"
# column = "hilic"
# polarity = "positive"
# threads = 3
# path = "."
# output.path = file.path(path, "metABM")
# max.isotope = 4
# mz.tol = 25
# rt.tol1 = 3
# rt.tol2 = 30
# cor.tol = 0
# int.tol = 500
# dp.tol = 0.5
# max.step = 3
# score.cutoff = 0
# remain = FALSE
# remain.per = 0.5


setGeneric(name = "metABM",
           def = function(annotation.result = "ms2.match.annotation.result.csv",
                          ms2.data = "ms2",
                          prefer.adduct = c("all", "M+H", "M+Na", "M-H"),
                          use.default.md = TRUE,
                          column = c("hilic", "rp"),
                          polarity = c("positive", "negative"),
                          threads = 3,
                          path = file.path(".", "MS2_match_result"),
                          output.path = file.path(".", "MRN_annotation_result"),
                          max.isotope = 4,
                          mz.tol = 25,
                          rt.tol1 = 3,
                          rt.tol2 = 30,
                          cor.tol = 0,
                          int.tol = 500,
                          dp.tol = 0.5,
                          max.step = 3,
                          score.cutoff = 0,
                          remain = FALSE,
                          remain.per = 0.5,
                          seed.neighbor.match.plot = FALSE,
                          candidate.num = 3
                          ){

             path <- path[1]
             output.path <- output.path[1]
             dir.create(output.path)
             dir.create(file.path(output.path, "intermediate_data"))
             ##parameters
             column <- match.arg(column)
             polarity <- match.arg(polarity)
             prefer.adduct <- match.arg(prefer.adduct)

# suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
##check for data
file <- c(dir(path), dir(file.path(path, "intermediate_data")))
need.file <- c(annotation.result, ms2.data)

file.check <- which(!need.file %in% file)

if(length(file.check) > 0) {
  stop(paste(need.file[file.check], collapse = " & "),
       " are not in the directory.")
}


#---------------------------------------------------------------------------
##read annotation.result from ms2Annotation
data("inHouse.compound", envir = environment())
data("kegg.compound", envir = environment())
cat("Read annotation result from ms2Annotation.\n")

# temp <- readr::read_csv(file.path(path, annotation.result), progress = FALSE, col_types = readr::cols())
# file.copy(from = file.path(path, annotation.result), to = file.path(output.path, "intermediate_data",annotation.result))
# file.copy(from = file.path(path, "intermediate_data","sample.info.csv"), to = file.path(output.path, "intermediate_data"))


#read data
data <- readAnnotation(data = annotation.result, rt.filter = FALSE,
                     inHouse.compound = inHouse.compound,
                     path = path)
# sample <- data[[2]]

##RT prediction
#---------------------------------------------------------------------------
data("inHouse.compound.md", envir = environment())
data("kegg.compound.md", envir = environment())
cat("\n")
cat("Construct RT prediction model.\n")

rt.result <- rtPrediction(data = data,
                          prefer.adduct = prefer.adduct,
                          threads = threads,
                          inHouse.compound.md = inHouse.compound.md,
                          kegg.compound.md = kegg.compound.md,
                          use.default.md = use.default.md,
                          column = column)

save(rt.result, file = file.path(output.path, "intermediate_data", "rt.result"), compress = "xz")

rm(list = c("data"))
gc()


##add prediction RT to inhouse compound and KEGG compound
inhouse.rt <- rt.result[[1]]
kegg.rt <- rt.result[[2]]

idx <- match(inHouse.compound$Lab.ID, row.names(inhouse.rt))


inHouse.compound <- data.frame(inHouse.compound, inhouse.rt[idx,],
                               stringsAsFactors = FALSE)

##Add the predicted RT to inHouse compound and KEGG compound
inHouse.compound$RT[is.na(inHouse.compound$RT)] <-
  median(inHouse.compound$RT[!is.na(inHouse.compound$RT)])

idx <- match(kegg.compound$ID, row.names(kegg.rt))

kegg.compound <- data.frame(kegg.compound, kegg.rt[idx,],
                            stringsAsFactors = FALSE)

kegg.compound$RT[is.na(kegg.compound$RT)] <-
  median(kegg.compound$RT[!is.na(kegg.compound$RT)])

rm(list = c("inHouse.compound.md", "kegg.compound.md", "rt.result"))
gc()

#==============================================================================
##Metabolite annotation
load(file.path(path, "intermediate_data", ms2.data))
ms2.data <- ms2
# save(ms2.data, file = file.path(output.path, "intermediate_data","ms2"), compress = "xz")
data("kegg.rpair2", envir = environment())

if(column == "hilic") {
  data("adduct.table.hilic", envir = environment())
  adduct.table <- adduct.table.hilic
  rm(list = c("adduct.table.hilic"))
  gc()
}

if(column == "rp") {
  data("adduct.table.rp", envir = environment())
  adduct.table <- adduct.table.rp
  rm(list = c("adduct.table.rp"))
  gc()
  }

adduct.table <- adduct.table[adduct.table$mode == polarity,]


metA(data = annotation.result,
     ms2 = ms2.data,
     adduct.table = adduct.table,
     max.isotope = max.isotope,
     polarity = polarity,
     mz.tol = mz.tol,
     rt.tol1 = rt.tol1,
     rt.tol2 = rt.tol2,
     cor.tol = cor.tol,
     int.tol = int.tol,
     dp.tol = dp.tol,
     max.step = max.step,
     metabolite = kegg.compound,
     metabolic.network = kegg.rpair2,
     inHouse.compound = inHouse.compound,
     threads = threads,
     score.cutoff = score.cutoff,
     remain = remain,
     remain.per = remain.per,
     path = path,
     output.path = output.path,
     seed.neighbor.match.plot = seed.neighbor.match.plot,
     candidate.num = candidate.num)

rm(list = c("inHouse.compound", "kegg.rpair2"))
gc()

cat("\n")
cat("metABM is done.\n")
})









# title metA
# description Annotate peak table from seeds.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param data The data (csv) from ms2Annotation.
# param ms2 The ms2 data of peak table.
# param adduct.table Adduct table.
# param max.isotope The number of isotope peaks.
# param polarity The polarity.
# param mz.tol The mz tolerance.
# param rt.tol1 The RT tolerance for isotope and adduct annotation. (second)
# param rt.tol2 The RT tolerance for metabolite annotation. (\%)
# param cor.tol The cor tolerance.
# param int.tol The intensity ratio tolerance.
# param dp.tol The tolerance of dot product.
# param max.step The max number of reaction step.
# param remain Remain some seeds as validation or not.
# param remain.per The percentage of remained seeds.
# param metabolite The kegg compound database.
# param metabolic.network kegg.rpair2.
# param inHouse.compound The inhouse compound database with predicted RT.
# param threads The number of threads.
# param score.cutoff Score cutoff.
# param path The directory.
# param output.path The directory of outputing results.
# return  The tags data with annotation result.




# data = "data2.csv"
# ms2 = ms2
# adduct.table = adduct.table
# max.isotope = 4
# polarity = "positive"
# mz.tol = 25
# rt.tol1 = 3
# rt.tol2 = 30
# cor.tol = 0
# int.tol = 500
# dp.tol = 0
# max.step = 3
# remain = FALSE
# remain.per = 0.5
# metabolite = kegg.compound
# metabolic.network = kegg.rpair2
# inHouse.compound = inHouse.compound
# threads = 3
# score.cutoff = 0
# path = "."
# output.path = file.path(".", "wrong annotation")
# rt.filter = FALSE

setGeneric(name = "metA", def = function(data = "data.csv",
                                         ms2,
                                         adduct.table,
                                         max.isotope = 4,
                                         polarity = c("positive", "negative"),
                                         mz.tol = 25,
                                         rt.tol1 = 3,
                                         rt.tol2 = 30,
                                         cor.tol = 0,
                                         int.tol = 500,
                                         dp.tol = 0.5,
                                         max.step = 3,
                                         remain = FALSE,
                                         remain.per = 0.5,
                                         metabolite,
                                         metabolic.network = kegg.rpair2,
                                         inHouse.compound,
                                         threads = 3,
                                         score.cutoff = 0,
                                         path = ".",
                                         output.path = ".",
                                         rt.filter = TRUE,
                                         seed.neighbor.match.plot = FALSE,
                                         candidate.num = 3){


  dir.create(output.path)
  dir.create(file.path(output.path, "intermediate_data"))
  polarity <- match.arg(polarity)
  adduct.table <- adduct.table[adduct.table$mode == polarity,]
  #thermo data is necessary for CHNOSZ
  # suppressMessages(expr = data(thermo, package = "CHNOSZ", envir = environment()))
  ##data is imported from XCMS processed MS2 data, which is used to get the
  ##identification result and sample abundance information

  # cat("\n")
  # cat("Read annotation result from ms2Annotation.\n")
  data <- readAnnotation(data = data, path = path, rt.tol = rt.tol2,
                         rt.filter = rt.filter,
                         inHouse.compound = inHouse.compound)
  sample <- data[[2]]
  peak.int <- apply(sample, 1, median)
  rm(list = "sample")
  gc()

  tags <- data[[1]]
  tags[which(tags=="NA", arr.ind = TRUE)] <- NA
  rm(list = "data")
  gc()

  ms2.name <- unname(unlist(lapply(ms2, function(x) {x[[1]]["NAME",]})))
  ##tags is transformed in to tags2 type data
  cat("\n")
  cat("Transform tags to tags2.\n")


  ####use the tags2 data if it is in the directory
  if(any(dir(path) == "tags2")) {
    load(file.path(path, "tags2"))
    cat("=============================IMPORTANT NOTE=========================\n")
    cat("===                      Use the tags2 from local               ====\n")
    cat("====================================================================\n")
  }else{
    tags2 <- changeTags(tags = tags,
                        ms2 = ms2,
                        weight.mz = 0.25,
                        weight.rt = 0.25,
                        weight.dp = 0.5,
                        mz.tol = 25,
                        rt.tol = 30,
                        dp.tol = 0.5,
                        adduct.table = adduct.table)
  }

  rm(list = c("tags"))
  gc()
  #save raw tags2 data in local
  # raw.tags2 <- tags2
  # save(raw.tags2, file = file.path(output.path, "raw.tags2"))
  # rm(list = "raw.tags2")
  ##remain some seeds as validation dataset
  if(remain){
    cat("\n")
    cat(remain.per*100,"% of seeds are remained as validation data.\n")
    idx <-
      which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)
    save(idx, file = file.path(output.path,"idx"), compress = "xz")
    if(any(dir(output.path)=="remain.idx")){
      load(file.path(output.path, "remain.idx"))
    }else{
      remain.idx <- lapply(1:30, function(x) sort(sample(idx, round(remain.per*length(idx)))))
      remain.idx <- remain.idx[[sample(1:30, 1)]]
      # remain.idx <- idx[sample(1:length(idx), round(remain.per@*length(idx)))]
      save(remain.idx, file = file.path(output.path,"remain.idx"), compress = "xz")
    }

    ##remove annotation of reamin idx in tags2
    tags2[remain.idx] <- lapply(tags2[remain.idx], function(x){
      x@annotation <- list()
      x
    })
  }

  #get the initinal seeds for annotation
  idx <-
    which(unlist(lapply(tags2, function(x) {length(x@annotation)})) > 0)

  ##add seed information to tags2

  tags2 <- seedLabel(tags2 = tags2, round = 1,
                     seed.idx = idx, score.cutoff = score.cutoff)


  # tags2 <- seedLabel(tags2 = tags2, round = 1,
  #                    seed.idx = idx, score.cutoff = 0)

  seed.annotation.number <- sum(unlist(showTags2(tags2, slot = "as.seed")[[1]]))
  seed.peak.number <- length(idx)
  ##first round is 1

  if(any(dir(file.path(output.path, "intermediate_data")) == "tags2.after.annotation")){
    load(file.path(output.path,"intermediate_data/tags2.after.annotation"))
    new.tags <- tags2.after.annotation
    rm(list = c("tags2.after.annotation"))
  }else{
    round <- 1
    while(length(idx) > 0){
      cat("\n")
      cat("Round", round, "annotation\n")
      cat("Seed peak number:", seed.peak.number)
      cat("\n")
      cat("Seed annotation number:", seed.annotation.number)
      cat("\n")
      cat("Memory used: ")
      temp <- pryr::mem_used()
      print(temp)
      cat("-------------------------------------------------------------------\n")
      ##Ronund 1 annotation
      ##parameter preparation for metIden


      new.tags <- metIden(
        peak.int = peak.int,
        tags2 = tags2,
        seed.idx = idx,
        max.isotope = max.isotope,
        polarity = polarity,
        rt.tol1 = rt.tol1,
        rt.tol2 = rt.tol2,
        mz.tol = mz.tol,
        cor.tol = cor.tol,
        int.tol = int.tol,
        dp.tol = dp.tol,
        metabolite = metabolite,
        metabolic.network = metabolic.network,
        max.step = max.step,
        ms2 = ms2,
        round = round,
        remove.index = NULL,
        iso.annotation = TRUE,
        add.annotation = TRUE,
        met.annotation = TRUE,
        threads = threads,
        adduct.table = adduct.table
      )

      ##remove the metAnnotation which is from the peak itself
      annotation.type <- showTags2(new.tags, slot = "annotation.type")
      annotation.idx <- lapply(annotation.type, function(x){
        which(unlist(x) == "metAnnotation")
      })

      annotation.idx <- annotation.idx[which(unlist(lapply(annotation.idx, length))!=0)]

      peak.index <- as.numeric(names(annotation.idx))

      # count <- NULL
      for(i in 1:length(peak.index)){
        # cat("\n");cat("i:");cat(i);cat("\n")
        temp.tags2 <- new.tags[[peak.index[i]]]
        temp.annotation.idx <- annotation.idx[[i]]
        for(j in temp.annotation.idx){
          # cat(j);cat(" ")
          temp.annotation <- temp.tags2@annotation[[j]]
          if(temp.annotation$From.peak ==  temp.tags2@name & temp.annotation$type == "metAnnotation"){
            temp.tags2@annotation[[j]] <- list()
            # count <- c(count, list(c(i,j)))
          }else{
            next()
          }
        }
        temp.tags2@annotation <- temp.tags2@annotation[unlist(lapply(temp.tags2@annotation, length)) !=0]
        if(length(temp.tags2@annotation) == 0){temp.tags2@annotation <- list()}
        new.tags[[peak.index[i]]] <- temp.tags2
      }

      rm(list = c("annotation.type", "annotation.idx", "peak.index"))


      ##-------------------------------------------------------------------
      # temp.result <- trans2Matrix(tags2 = new.tags, base = "annotation")
      # temp.result <- temp.result[temp.result$isotope == "[M]",]
      # temp.result <- temp.result[temp.result$type != "seed",]
      # temp.mz <- as.numeric(temp.result$mz)
      # temp.formula <- temp.result$Formula
      # temp.adduct <- temp.result$adduct
      # ###
      # temp.formula1 <- mapply(function(x,y){sumFormula(x,y)},x = temp.formula, y = temp.adduct)
      # temp.mass <- sapply(temp.formula1, function(x) Rdisop::getMass(Rdisop::getMolecule(x)))
      # temp.mz.error <- abs(temp.mz - temp.mass)*10^6/temp.mass
      # temp.mz.error[is.na(temp.mz.error)] <- 100
      # sum(temp.mz.error > 25)
      # idx1 <- which(temp.mz.error > 100)
      # table(temp.result$type[idx1])
      # temp.result[idx1,-c(10,11,14,17,18,19,20)]
      # temp.mz.error[idx1]
      #------------------------------------------------------------------
      ##find seed for next round annotation
      #
      #round should be added 1
      round <- round + 1
      ##find new seed index for next round annotation
      annotation.idx <- which(showTags2(new.tags,slot = "annotation.len") > 0)
      ##temp.idx is index of seeds which have at least one annotation that is not
      ##seed before
      temp.idx <-
        which(unlist(lapply(showTags2(new.tags[annotation.idx], slot = "as.seed")[[1]], function(x) {any(!x)})))

      idx <- annotation.idx[temp.idx]

      ##seed criteria: 1, not isotope; 2, not seed before; 3, score bigger than
      ##cutoff
      idx.len <- unlist(lapply(new.tags[idx], function(x) {
        annotation <- x@annotation
        annotation.isotope <- unlist(lapply(annotation, function(x) {x$isotope}))
        annotation.id <- unlist(lapply(annotation, function(x) {x$to}))
        annotation.as.seed <- unlist(lapply(annotation, function(x) {x$as.seed}))
        annotation.score <- unlist(lapply(annotation, function(x) {x$score}))
        temp.idx <-
          which(annotation.isotope == "[M]" &
                  !duplicated(annotation.id) &
                  !annotation.as.seed &
                  annotation.score >= score.cutoff)
        length(temp.idx)
      }))

      idx <- idx[which(idx.len > 0)]

      if(length(idx) == 0) {break()}
      ##add seed information to tags2
      new.tags <- seedLabel(tags2 = new.tags, round = round,
                            seed.idx = idx, score.cutoff = score.cutoff)

      tags2 <- new.tags
      rm(list = "new.tags")
      gc()

      #need be fixed
      seed.annotation.number <-
        unlist(showTags2(tags2, slot = "as.seed")[[2]])
      seed.annotation.number <- seed.annotation.number[!is.na(seed.annotation.number)]
      seed.annotation.number <- sum(seed.annotation.number == round)
      seed.peak.number <- length(idx)
    }
    ##-----------------------------------------------------------------------
    ###-----------------------------------------------------------------------

    ###order annotation for each peak according to score
    new.tags <-
      lapply(new.tags,
             function(x) {filterPeak(object = x, score.thr = 0)})

    ##save new.tags2 after annotation
    tags2.after.annotation <- new.tags

    save(tags2.after.annotation,
         file = file.path(output.path, "intermediate_data","tags2.after.annotation"), compress = "xz")
    rm(tags2.after.annotation)
    gc()
  }

  ##begin to give confidence to annotation of peaks
  result <- trans2Matrix(tags2 = new.tags, base = "peak")
  # result3 <- result
  # write.table(result3, "Annotation.of.peaks.xlsx", row.names = FALSE, sep = "\t")

  #####give the confidence level of each ID
  cat('\n')
  cat("Assign confidence.\n")
  unique.id <- unique(result$to)

  pbapply::pboptions(type="timer", style = 1)

  # id.result <- pbapply::pblapply(unique.id, function(x) {
    id.result <- lapply(unique.id, function(x) {
    temp.id <- x
    temp.idx <- which(result$to == temp.id)
    temp.result <- result[temp.idx,]
    temp.rt <- temp.result$rt
    ##group peaks according to RT
    rt.class <- groupRT(rt = temp.rt, rt.tol = rt.tol1)
    temp.result <- lapply(rt.class, function(x) {temp.result[x,]})

    temp.result <-  lapply(temp.result, function(x) {
      x1 <- data.frame()
      unique.peak.name <- unique(x$name)

      for(j in seq_along(unique.peak.name)){
        temp.idx <- which(x$name == unique.peak.name[j])
        if(length(temp.idx) == 1){
          Note <- 1
          x1 <- rbind(x1, data.frame(x[temp.idx,], Note, stringsAsFactors = FALSE))
          # colnames(x1)[ncol(x1)] <- "Note"
        }else{
          temp.x <- x[temp.idx,]
          if(any(temp.x$type == "seed")){
            temp.x <- temp.x[which(temp.x$type == "seed")[1],]
          }else{
            temp.x <- temp.x[order(temp.x$score, decreasing = TRUE),]
            temp.x <- temp.x[1,]
          }
          Note <- length(temp.idx)
          x1 <- rbind(x1, data.frame(temp.x, Note, stringsAsFactors = FALSE))
          # colnames(x1)[ncol(x1)] <- "Note"
        }
      }

      colnames(x1)[ncol(x1)] <- "Note"
      x1
    })
    temp.confidence <- unlist(lapply(temp.result, function(x) confidence(x, polarity = polarity)))
    #
    # ##if the metabolite has more than groups, and one group confidence is grade4,
    # ##then the group with the grade 4 is removed.
    # if(length(unique(temp.confidence)) > 1 & any(temp.confidence == "grade4")){
    #   temp.result <- temp.result[-which(temp.confidence == "grade4")]
    #   temp.confidence <- temp.confidence[-which(temp.confidence == "grade4")]
    # }

    temp.result <- mapply(function(x,
                                   Confidence,
                                   group) {
      result <-
        data.frame(x,
                   Confidence,
                   paste(x$to, group, sep = "_"),
                   stringsAsFactors = FALSE)
      colnames(result)[ncol(result)] <- "group"
      list(result)
    },
    x = temp.result,
    Confidence = temp.confidence,
    group = c(1:length(temp.result)))

    temp.result <- do.call(rbind, temp.result)
    temp.result
  })

  result1 <- do.call(rbind, id.result)
  # save(result1, file = file.path(path, "result1"))
  ##remove redundancy
  cat("\n")
  cat("Remove redundancy.\n")
  tags.result <- removeRedundancy(result = result1, path = output.path, polarity = polarity)
  load(file.path(output.path, "redun"))
  redun1 <- redun
  save(redun1, file = file.path(output.path, "intermediate_data","redun1"), compress = "xz")
  unlink(x = file.path(output.path, "redun"), recursive = TRUE)
  save(tags.result, file = file.path(output.path, "intermediate_data","tags.result"), compress = "xz")

  ##filter new.tags according to annotation.result
  new.tags <- result2Tags(result = tags.result, tags2 = new.tags)
  tags2.after.redundancy.remove <- new.tags
  save(tags2.after.redundancy.remove,
       file = file.path(output.path, "intermediate_data","tags2.after.redundancy.remove"), compress = "xz")

  ##change tags2.after.redundancy.remove to csv like ms2Annotation

  cat("\n")
  data("kegg.compound", envir = environment())

  temp <- getAnnotationResult(tags2 = tags2.after.redundancy.remove,
                              tags.result2 = tags.result,
                              kegg.compound = kegg.compound,
                              candidate.num = candidate.num)

  cat("\n")
  readr::write_csv(temp, file.path(output.path, "MRN.annotation.result.csv"))

  ##output MS/MS matching spectra
  if(seed.neighbor.match.plot){
    dir.create(file.path(output.path, "Seed_Neighbor_MS2_match_spectra"))

    cat('Plot spectra match results (Seed vs Neighbor).\n')
    # data("kegg.compound", envir = environment())
    # load(file.path(output.path, "intermediate_data", "ms2"))
    plotSeedNeighbor(tags2 = tags2.after.redundancy.remove,
                     type = "pdf",
                     path = file.path(output.path, "Seed_Neighbor_MS2_match_spectra"),
                     ms2 = ms2,
                     kegg.compound = kegg.compound)
    cat("\n")
  }


  rm(list = c("tags2.after.redundancy.remove", "new.tags", "temp"))
  gc()
})


#-----------------------------------------------------------------------------
# title removeRedundancy
# description Remove annotation and peak redundancy.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param result The result from metA.
# param path The directory.
# return  A result after removing redundancy.
# export


setGeneric(name = "removeRedundancy",
           def = function(result,
                          polarity = c("positive", "negative"),
                          path = "."){
             polarity <- match.arg(polarity)
             #redundancy
             redun <- list()
             redun[[1]] <- calculateRedundancy(result = result)
             count <- 1
             delta.redun <- 1
             cat("Round: ")
             while(delta.redun > 0.01 & count <= 5){
               cat(count); cat(" ")
               unique.id <- unique(result$to)
               result <- lapply(unique.id, function(temp.id){
                 temp.idx <- which(result$to == temp.id)
                 temp.result <- result[temp.idx,]
                 temp.group <- unique(temp.result$group)
                 temp.idx <- lapply(temp.group, function(x) which(temp.result$group == x))
                 temp.result <- lapply(temp.idx, function(x) temp.result[x,])
                 ##confidence
                 temp.confidence <- unlist(lapply(temp.result, function(x) confidence(x, polarity = polarity)))

                 ##if the metabolite has more than two groups, and one group confidence is grade4,
                 ##then the group with the grade 4 is removed.
                 if(length(unique(temp.confidence)) > 1 & any(temp.confidence == "grade4")){
                   temp.result <- temp.result[-which(temp.confidence == "grade4")]
                   temp.confidence <- temp.confidence[-which(temp.confidence == "grade4")]
                 }

                 temp.result <- mapply(function(x, Confidence) {
                   # result <- data.frame(x, Confidence, stringsAsFactors = FALSE)
                   result <- x
                   result$Confidence <- Confidence
                   result
                   list(result)
                 },
                 x = temp.result,
                 Confidence = temp.confidence)

                 temp.result <- do.call(rbind, temp.result)
                 temp.result
               })

               result <- do.call(rbind, result)
               ##remove the one peak vs many annotation
               ##if one peak match more than two annotations, ony the annotation which has
               ## the biggest grade is remained
               unique.peak <- unique(result$name)

               remain.idx <- lapply(unique.peak, function(temp.peak){
                 temp.idx <- which(temp.peak == result$name)
                 if(length(temp.idx) == 1) return(temp.idx)
                 temp.result <- result[temp.idx,]
                 temp.confidence <- temp.result$Confidence
                 if(length(unique(temp.confidence)) > 1){
                   temp.grade <- as.numeric(substr(temp.confidence, 6,6))
                   remain.idx <- temp.idx[which(temp.grade == temp.grade[which.min(temp.grade)])]
                   return(remain.idx)
                 }
                 return(temp.idx)
               })

               remain.idx <- unlist(remain.idx)

               result <- result[remain.idx,]

               result <- result
               ###give confidence to result again
               group <- unique(result$group)

               result <- lapply(group, function(temp.group){
                 temp.result <- result[result$group == temp.group,]
                 new.confidence <- confidence(temp.result, polarity = polarity)
                 temp.result$Confidence <- new.confidence
                 temp.result
               })

               result <- do.call(rbind, result)
               redun <- c(redun, list(calculateRedundancy(result = result)))
               delta.redun <- abs(mean(redun[[length(redun)-1]]) - mean(redun[[length(redun)]]))
               count <- count + 1
             }
             cat("\n")
             save(redun, file = file.path(path, "redun"), compress = "xz")
             result
           })



#-----------------------------------------------------------------------------
# title calculateRedundancy
# description Calculate redundancy of result from metA.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param result The result from metA.
# return  peak and metabolite redundancy.
# export

setGeneric(name = "calculateRedundancy",
           def = function(result){
             unique.name <- unique(result$name)
             unique.group <- unique(result$group)
             unique.id <- unique(result$to)
             peak.redun <- nrow(result)/length(unique.name)
             id.redun <- length(unique.group)/length(unique.id)
             redundancy <- c(peak.redun, id.redun)
             names(redundancy) <- c("Peak.redundancy", "Metabolite.redundancy")
             redundancy
           })



#-----------------------------------------------------------------------------
# title groupRT
# description Group peaks according to RT.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param rt The RTs of peaks.
# param rt.tol The tolerance of RT.
# return  The index of peaks in each group.
# export

setGeneric(name = "groupRT",
           def = function(rt,
                          rt.tol = 3) {
             rt <- round(rt, 2)
             breaks <- seq(min(rt), max(rt) + rt.tol, by = rt.tol)
             w <- hist(rt, breaks = breaks, plot = FALSE)
             counts <- w$counts
             idx <- which(counts != 0)
             interval.range <- data.frame(breaks[idx],breaks[idx+1], stringsAsFactors = FALSE)
             interval.range <-
               data.frame(interval.range, c(1:nrow(interval.range)), stringsAsFactors = FALSE)
             colnames(interval.range) <- c("From", "To", "Module")

             rt.class <-
               sapply(rt, function(x) {
                 a <- x - round(interval.range$From, 2)
                 b <- round(interval.range$To, 2) - x
                 idx <- which(a >= 0 & b >= 0)
                 if(length(idx) > 0) {idx <- idx[1]}
                 idx
               })

             rt.class <-
               lapply(sort(unique(rt.class)), function(x) {
                 unname(which(rt.class == x))
               })
             return(rt.class)
           }
)


#-----------------------------------------------------------------------------
# title confidence
# description Assign confidence to metabolite ID group.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param x The metabolite ID group.
# return  The confidence of metabolite ID group.



setGeneric(name = "confidence",
           def = function(x,
                          polarity = c("positive", "negative")){
             polarity <- match.arg(polarity)

             if(polarity == "positive"){
               prefer.adduct <- c("M+H", "M+Na", "M+NH4")
             }else{
               prefer.adduct <- c("M-H", "M+CH3COO", "M+Cl")
             }

             adduct <- x$adduct
             isotope <- x$isotope
             type <- x$type
             if(any(type == "seed")) return("grade1")
             temp.idx1 <- which(isotope %in% c("[M+1]","[M+2]","[M+3]","[M+4]"))
             if(length(temp.idx1) > 0) return("grade2")
             if(length(temp.idx1) == 0){
               temp.idx2 <- which(adduct %in% prefer.adduct)
               if(length(temp.idx2) > 0) {return("grade3")}
               if(length(temp.idx2) == 0) {return("grade4")}
             }
           })


#-----------------------------------------------------------------------------
# title seedLabel
# description Label seed in tags2 data.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param tags2 The tags2 data.
# param round The annotation round.
# param seed.idx The index of seeds.
# param score.cutoff The score cutoff.
# return  The tags2 data with labeled seeds.

setGeneric(name = "seedLabel",
           def = function(tags2,
                          round,
                          seed.idx,
                          score.cutoff = 0){
             tags2[seed.idx] <- lapply(tags2[seed.idx], function(x) {
               annotation <- x@annotation
               annotation.isotope <-
                 unlist(lapply(annotation, function(x) {x$isotope}))
               annotation.id <-
                 unlist(lapply(annotation, function(x) {x$to}))
               annotation.as.seed <-
                 unlist(lapply(annotation, function(x) {x$as.seed}))
               annotation.score <-
                 unlist(lapply(annotation, function(x) {x$score}))
               temp.idx <-
                 which(annotation.isotope == "[M]" &
                         !duplicated(annotation.id) &
                         !annotation.as.seed &
                         annotation.score >= score.cutoff)
               if(length(temp.idx) != 0){
                 annotation[temp.idx] <-
                   lapply(annotation[temp.idx], function(x) {x$as.seed <- TRUE;x})
                 annotation[temp.idx] <-
                   lapply(annotation[temp.idx], function(x) {x$as.seed.round <- round;x})
                 x@annotation <- annotation
               }else{
                 x@annotation <- annotation
               }
               x
             })
             return(tags2)
           })






# load("new.tags")
# load("raw.tags2")
# load("remain.idx")
#
# new.tags[[remain.idx[1]]]
# raw.tags2[[remain.idx[1]]]
#
#
# result1 <- tags2Result(tags2 = new.tags[remain.idx], score.cutoff = 0)
# result2 <- tags2Result(tags2 = raw.tags2[remain.idx], score.cutoff = 0)
#
#
# to1 <- result1$to
# to2 <- result2$to
#
# to1 <- sapply(to1, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# to1 <- unname(to1)
#
#
# to2 <- sapply(to2, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# to2 <- unname(to2)
#
#
#
# formula1 <- result1$Formula
# formula2 <- result2$Formula
#
# formula1 <- sapply(formula1, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# formula1 <- unname(formula1)
#
#
# formula2 <- sapply(formula2, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# formula2 <- unname(formula2)
#
#
#
# note <- mapply(function(x,y,z,w){
#   if(is.na(x)) return("no")
#   x <- strsplit(x, split = ";")[[1]]
#   y <- strsplit(y, split = ";")[[1]]
#   z <- strsplit(z, split = ";")[[1]]
#   w <- strsplit(w, split = ";")[[1]]
#   temp <- intersect(x, y)
#   if(length(temp) > 0) return("right")
#   if(length(temp) == 0) {
#     temp <- intersect(z, w)
#     if(length(temp) > 0) return("isomer")
#     if(length(temp) == 0) return("wrong")
#   }
# },
# x = to1,
# y = to2,
# z = formula1,
# w = formula2)
#
# note <- unname(note)
#
# annotation.result <- data.frame(to1, to2, note, stringsAsFactors = FALSE)
# write.csv(annotation.result, "annotation.result.csv")
#
# par(mar = c(5,5,4,2))
# pie(c(91*3, 106*3-91*3), col = c("lightseagreen", "salmon"), border = 'white',
#     labels = c("Annotation", "No annotation"), clockwise = T, cex = 1.5, radius = 1.0)
# par(new = T)
# pie(1, col = "white", border = NA, labels = "", radius = 0.6)
#
#
#
#
# pie(c(76+73+77, 14, 33), col = c("orchid", "lightseagreen", "firebrick1"), border = 'white',
#     labels = c("Right", "Isomer", "Wrong"), clockwise = T, cex = 1.5, radius = 1.0)
# par(new = T)
# pie(1, col = "white", border = NA, labels = "", radius = 0.6)
#
#
#




# idx1 <- which(result2_old$type == "seed")
# name1 <- result2_old$name[idx1]
# a <- result2_old[idx1,]
# b <- result2[idx1,]
#
#
#
# index1 <- which(showTags2(raw.tags2, "annotation.len") > 0)
# result1 <- tags2Result(tags2 = raw.tags2[index1], score.cutoff = 0)
# result2 <- tags2Result(tags2 = tags2[index1], score.cutoff = 0)
#
#
# to1 <- result1$to
# to2 <- result2$to
#
# to1 <- sapply(to1, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# to1 <- unname(to1)
#
#
# to2 <- sapply(to2, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# to2 <- unname(to2)
#
#
#
# formula1 <- result1$Formula
# formula2 <- result2$Formula
#
# formula1 <- sapply(formula1, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# formula1 <- unname(formula1)
#
#
# formula2 <- sapply(formula2, function(x) {
#   if(is.na(x)) return(x)
#   temp <- strsplit(x, split = ";")[[1]]
#   temp <- unique(temp)
#   temp <- paste(temp, collapse = ";")
#   temp
# })
#
# formula2 <- unname(formula2)
#
#
#
# note <- mapply(function(x,y,z,w){
#   if(is.na(x)) return("no")
#   x <- strsplit(x, split = ";")[[1]]
#   y <- strsplit(y, split = ";")[[1]]
#   z <- strsplit(z, split = ";")[[1]]
#   w <- strsplit(w, split = ";")[[1]]
#   temp <- intersect(x, y)
#   if(length(temp) > 0) return("right")
#   if(length(temp) == 0) {
#     temp <- intersect(z, w)
#     if(length(temp) > 0) return("isomer")
#     if(length(temp) == 0) return("wrong")
#   }
# },
# x = to1,
# y = to2,
# z = formula1,
# w = formula2)
#
# note <- unname(note)
#
#
# annotation.result <- data.frame(to1, formula1, to2, formula2, note,stringsAsFactors = FALSE)
# write.csv(annotation.result, "annotation.result.csv")
#
#
# temp <- mapply(function(x, y){
#   x <- strsplit(x, split = ";")[[1]]
#   y <- strsplit(y, split = ";")[[1]]
#   len1 <- unname(length(x))
#   len2 <- unname(length(y))
#   len3 <- unname(length(intersect(x, y)))
#   return(list(c(len1, len2, len3)))
# },
# x = to1,
# y = to2)
#
#
#
# temp <- do.call(rbind, temp)
#
# par(mar = c(5,5,4,2))
# plot(temp[,1], temp[,2], xlab = "Annotation number (Measured RT)",
#      ylab = "Annotation number (Predicted RT)",
#      cex.lab = 1.8, cex.axis = 1.5, pch = 19)
# par(xpd = FALSE)
# abline(0,1,lty =2, lwd = 1.5, col = "tomato")
#
# plot(temp[,4], temp[,1], xlab = "Peak index", ylab = "Annotation number",
#      cex.lab = 1.8, cex.axis = 1.5,
#      pch = 19)
#
# points(temp[,4], temp[,2], pch = 15, col = 'lightseagreen')
#
#
# x <- barplot(table(temp[,2]), border = NA, cex = 1.5, cex.lab = 1.8,cex.axis = 1.5,
#              xlab = "Annotation number", ylab = "Peak number",
#              col = c("lightseagreen", "salmon", "tan1","black" ,"orchid4", rep("black", 2)))
# par(xpd = TRUE)
# text(x = x[,1], y = table(temp[,2])+5, labels = table(temp[,2]), cex = 1.3)
#
#
# pie(table(temp[,2]), clockwise = TRUE, border = 'white',
#     col = c("lightseagreen", "salmon", "tan1","black","orchid4", rep("black", 2)))
# par(new = T)
# pie(1, col = "white", border = "white", labels = "", radius = 0.5)
#
#
# temp <- as.data.frame(temp)
# temp <- data.frame(temp, c(1:nrow(temp)), stringsAsFactors = FALSE)
# colnames(temp) <- c("v1", "v2", "v3", "v4")
# library(ggplot2)
# ggplot(temp, aes(x = v4, y =v2))+geom_point()+
#   coord_polar()
#
#
# temp1 <- temp[temp[,1] == 1,]
# table(temp1[,2])
# temp2 <- temp[temp[,1] == 2,]
# temp3 <- temp[temp[,1] == 3,]
# temp4 <- temp[temp[,1] == 4,]
# temp5 <- temp[temp[,1] == 5,]
# temp6 <- temp[temp[,1] == 6,]
# temp7 <- temp[temp[,1] == 7,]
#
# barplot(c(166,4), border = NA, col = c("lightseagreen", "salmon"))
#
# pie(table(temp1[,2]), clockwise = TRUE, border = 'white',
#     col = c("lightseagreen", "salmon", "tan1","black","orchid4", rep("black", 2)),
#     cex = 1.8)
#
# pie(table(temp7[,2]), clockwise = TRUE, border = 'white',
#     col = c("lightseagreen", "salmon", "tan1","black","orchid4", rep("black", 2)),
#     cex = 1.8)
# par(new = T)
# pie(1, col = "white", border = "white", labels = "", radius = 0.5)
#
# table(temp1[,2])
# table(temp2[,2])
# table(temp3[,2])
# table(temp3[,2])
# pie(c(231, 9), clockwise = TRUE, border = 'white',
#     col = c("lightseagreen", "salmon", "tan1","black","orchid4", rep("black", 2)),
#     cex = 1.8, labels = c("Completely same", "Contains"), radius = 1)
# par(new = T)
# pie(1, col = "white", border = "white", labels = "", radius = 0.5)
#
# barplot()


# title metIden
# description Annotate peak table from seeds.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param peak.int The median/mean intensity of peaks.
# param tags2 The tags2 data of all peaks.
# param seed.idx The index of seeds.
# param max.isotope The number of isotope peaks.
# param polarity The polarity.
# param rt.tol1 The RT tolerance for isotope and adduct annotation. (second)
# param rt.tol2 The RT tolerance for metabolite annotation (\%).
# param mz.tol The mz tolerance.
# param cor.tol The cor tolerance.
# param int.tol The intensity ratio tolerance.
# param dp.tol The tolerance of dot product.
# param metabolite The kegg compound database.
# param metabolic.network kegg.rpair2
# param max.step The max number of reaction step.
# param ms2 The ms2 data of peak table.
# param round The round of annotation.
# param remove.index The index of peaks which you dont want to annotate.
# param iso.annotation Isotope annotation or not.
# param add.annotation Adduct annotation or not.
# param met.annotation Metabolite annotation or not.
# param threads The number of threads.
# param adduct.table Adduct table.
# return  The tags data with annotation result.

# setwd("/home/jasper/work/MetDNA/metIden")
# load("peak.int")
# load("tags2")
# load("idx")
# load("kegg.compound.rda")
# metabolite <- kegg.compound
# load("kegg.rpair2.rda")
# metabolic.network <- kegg.rpair2
# load('adduct.table')
# load("ms2")
#
# max.isotope = 4
# polarity = "positive"
# #rt.tol1 is absolute
# rt.tol1 = 3
# #rt.tol2 is relative
# rt.tol2 = 30
# mz.tol = 25
# cor.tol = 0
# #for isotope annotation
# int.tol = 500
# dp.tol = 0.5
# max.step = 3
# round <- 1
# # cor.matrix,
# remove.index = NULL
# iso.annotation = TRUE
# add.annotation = TRUE
# met.annotation = TRUE
# threads = 3
# seed.idx <- idx

setGeneric(name = "metIden",
           def = function(peak.int,
                          tags2,
                          seed.idx,
                          max.isotope = 4,
                          polarity = c("positive", "negative"),
                          #rt.tol1 is absolute
                          rt.tol1 = 3,
                          #rt.tol2 is relative
                          rt.tol2 = 30,
                          mz.tol = 25,
                          cor.tol = 0,
                          #for isotope annotation
                          int.tol = 500,
                          dp.tol = 0.5,
                          metabolite,
                          metabolic.network,
                          max.step = 3,
                          ms2,
                          round,
                          # cor.matrix,
                          remove.index = NULL,
                          iso.annotation = TRUE,
                          add.annotation = TRUE,
                          met.annotation = TRUE,
                          threads = 3,
                          adduct.table){
             #

             polarity <- match.arg(polarity)
             adduct.table <- adduct.table[adduct.table$mode == polarity,]
             # load the database
             ##ms2 is the ms2 data base, ms2.name is the name of ms2 spectra
             # ms2 <- ms2
             ms2.name <- unname(unlist(lapply(ms2, function(x) {x[[1]]["NAME",]})))

             ##tags2 is the S4 class peakInfo for tags
             ##peak.mz and peak.rt is the mz and rt of all peaks
             peak.name <- unname(unlist(lapply(tags2, function(x) x@name)))
             peak.mz <- unname(unlist(lapply(tags2, function(x) x@mz)))
             peak.rt <- unname(unlist(lapply(tags2, function(x) x@rt)))
             names(peak.mz) <- names(peak.rt) <- peak.name

             ##isotope annotation
             if(iso.annotation){

               # save(seed.idx, file = "seed.idx")
               # save(tags2, file = "tags2")
               # save(annotation.idx, file = "annotation.idx")
               cat("\n")
               cat("Isotope annotation.\n")
               annotation.idx <- 1:length(tags2)
               if(!is.null(remove.index)) {
                 annotation.idx <- setdiff(annotation.idx, remove.index)
               }
               pbapply::pboptions(type="timer", style = 1)
               # isotope.result <- pbapply::pblapply(seed.idx, function(i){
                 isotope.result <- lapply(seed.idx, function(i){
                 seed <- tags2[[i]]
                 seed.as.seed.round <- showTags2(list(seed), slot = "as.seed")[[2]][[1]]
                 seed.i <- which(seed.as.seed.round == round)
                 if(length(seed.i)==0){
                   return(NULL)
                 }else{
                   temp.result <- lapply(seed.i,function(j){
                     query.name <- seed@name
                     query.id <- seed@annotation[[j]]$to
                     query.charge <- seed@annotation[[j]]$charge
                     query.level <- seed@annotation[[j]]$level
                     query.formula <- stringr::str_trim(seed@annotation[[j]]$Formula)
                     query.adduct <- seed@annotation[[j]]$adduct
                     query.mz <- seed@mz
                     query.rt <- seed@rt
                     query.int <- peak.int[i]
                     peak.mz.temp <- peak.mz[annotation.idx]
                     peak.rt.temp <- peak.rt[annotation.idx]
                     peak.int.temp <- peak.int[annotation.idx]

                     # cor.temp <- cor.matrix[annotation.idx, anno.i]
                     cor.temp <- rep(1, length(annotation.idx))
                     isotopes.result <- isotopeAnnotation(id = query.id,
                                                          formula = query.formula,
                                                          adduct = query.adduct,
                                                          charge = query.charge,
                                                          mz = query.mz,
                                                          rt = query.rt,
                                                          int = query.int,
                                                          peak.mz = peak.mz.temp,
                                                          peak.rt = peak.rt.temp,
                                                          peak.int = peak.int.temp,
                                                          cor = cor.temp,
                                                          rt.tol = rt.tol1,
                                                          mz.tol = mz.tol,
                                                          cor.tol = cor.tol,
                                                          int.tol = int.tol,
                                                          max.isotope = max.isotope)
                     isotopes.result
                   })
                   temp.result
                 }
               })




               ###fix bugs
               # load("seed.idx")
               # load("tags2")
               # load("peak.int")
               # load("peak.mz")
               # load("peak.rt")
               # load("annotation.idx")
               #
               # for (i in seed.idx){
               #   cat(i);cat(" ")
               #   seed <- tags2[[i]]
               #   seed.as.seed.round <- showTags2(list(seed), slot = "as.seed")[[2]][[1]]
               #   seed.i <- which(seed.as.seed.round == round)
               #   if(length(seed.i)==0){
               #     return(NULL)
               #   }else{
               #     temp.result <- lapply(seed.i,function(j){
               #       query.name <- seed@name
               #       query.id <- seed@annotation[[j]]$to
               #       query.charge <- seed@annotation[[j]]$charge
               #       query.level <- seed@annotation[[j]]$level
               #       query.formula <- stringr::str_trim(seed@annotation[[j]]$Formula)
               #       query.adduct <- seed@annotation[[j]]$adduct
               #       query.mz <- seed@mz
               #       query.rt <- seed@rt
               #       query.int <- peak.int[i]
               #       peak.mz.temp <- peak.mz[annotation.idx]
               #       peak.rt.temp <- peak.rt[annotation.idx]
               #       peak.int.temp <- peak.int[annotation.idx]
               #
               #       # cor.temp <- cor.matrix[annotation.idx, anno.i]
               #       cor.temp <- rep(1, length(annotation.idx))
               #       isotopes.result <- isotopeAnnotation(id = query.id,
               #                                            formula = query.formula,
               #                                            adduct = query.adduct,
               #                                            charge = query.charge,
               #                                            mz = query.mz,
               #                                            rt = query.rt,
               #                                            int = query.int,
               #                                            peak.mz = peak.mz.temp,
               #                                            peak.rt = peak.rt.temp,
               #                                            peak.int = peak.int.temp,
               #                                            cor = cor.temp,
               #                                            rt.tol = rt.tol1,
               #                                            mz.tol = mz.tol,
               #                                            cor.tol = cor.tol,
               #                                            int.tol = int.tol,
               #                                            max.isotope = max.isotope)
               #       isotopes.result
               #     })
               #     temp.result
               #   }
               # }




               ##add isotope.result to tags2
               anno.idx <- showTags2(tags2[seed.idx], slot = "as.seed")[[2]]
               anno.idx <- lapply(anno.idx, function(x) which(x == round))
               for(j in 1:length(seed.idx)){
                 # cat(j); cat(" ")
                 temp.result <- isotope.result[[j]]
                 for(k in 1:length(temp.result)){
                   tags2 <- isotope2peakInfo(isotopes.result = temp.result[[k]],
                                             mz.tol = mz.tol,
                                             rt.tol = rt.tol1,
                                             int.tol = int.tol,
                                             weight.mz = 0.45,
                                             weight.rt = 0.45,
                                             weight.int = 0.1,
                                             peak.info = tags2,
                                             peak.idx = seed.idx[j],
                                             anno.idx = anno.idx[[j]][k],
                                             annotation.idx = annotation.idx)
                 }

               }

               rm(list = "isotope.result")
               gc()
             }
             #------------------------------------------------------------------
             ##adduct annotation
             if(add.annotation){
               cat("\n")
               cat("Adduct annotation.\n")

               temp.fun <- function(index, tags2, peak.mz, peak.rt,
                                    polarity, rt.tol, mz.tol,
                                    cor.tol, adduct.table,
                                    annotation.idx,
                                    round){
                 annotation.idx <- 1:length(tags2)
                 library(Rcpp, warn.conflicts = FALSE, quietly = TRUE)
                 # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
                 temp.tags2 <- tags2[index]
                 rm(tags2)
                 gc()
                 pbapply::pboptions(type="timer", style = 1)
                 # result <- pbapply::pblapply(temp.tags2, function(x){
                   result <- lapply(temp.tags2, function(x){
                   seed <- x
                   seed.as.seed.round <-
                     showTags2(list(seed), slot = "as.seed")[[2]][[1]]
                   seed.i <- which(seed.as.seed.round == round)
                   peak.mz.temp <- peak.mz[annotation.idx]
                   peak.rt.temp <- peak.rt[annotation.idx]
                   cor.temp <- rep(1, length(annotation.idx))
                   rm(list=c("peak.mz", "peak.rt"))
                   gc()
                   if(length(seed.i)==0){
                     return(NULL)
                   }else{
                     temp.result <- lapply(seed.i,function(j){
                       query.name <- seed@name
                       query.id <- seed@annotation[[j]]$to
                       query.charge <- seed@annotation[[j]]$charge
                       query.level <- seed@annotation[[j]]$level
                       query.formula <- stringr::str_trim(seed@annotation[[j]]$Formula)
                       query.adduct <- seed@annotation[[j]]$adduct
                       query.mz <- seed@mz
                       query.rt <- seed@rt

                       adduct.result <-
                         adductAnnotation(id = query.id,
                                          formula = query.formula,
                                          adduct = query.adduct,
                                          polarity = polarity,
                                          mz = query.mz,
                                          rt = query.rt,
                                          peak.mz = peak.mz.temp,
                                          peak.rt = peak.rt.temp,
                                          cor = cor.temp,
                                          rt.tol = rt.tol1,
                                          mz.tol = mz.tol,
                                          adduct.table = adduct.table,
                                          cor.tol = cor.tol)
                       adduct.result
                     })
                     temp.result
                   }
                 })
                 result
               }


               cl <- snow::makeCluster(threads, type = "SOCK")

               snow::clusterExport(cl, list = list("showTags2",
                                                   "adductAnnotation",
                                                   "sumFormula",
                                                   "checkElement",
                                                   "splitFormula","pasteElement"))
               nc <- length(cl)
               options(warn = -1)
               if(length(seed.idx) < threads){
                 threads <-  1
               }
               ichunks <- split(seed.idx, 1:threads)
               # library(Rcpp)
               adduct.result <-
                 snow::clusterApply(cl = cl, ichunks,
                                    fun = temp.fun,
                                    tags2 = tags2,
                                    peak.mz = peak.mz,
                                    peak.rt = peak.rt,
                                    polarity = polarity,
                                    rt.tol = rt.tol1,
                                    mz.tol = mz.tol,
                                    cor.tol = cor.tol,
                                    round = round,
                                    adduct.table = adduct.table,
                                    annotation.idx = annotation.idx)
               snow::stopCluster(cl)

               adduct.result1 <- adduct.result[[1]]
               if(threads > 1){
                 for(i in 2:length(adduct.result)){
                   adduct.result1 <- c(adduct.result1, adduct.result[[i]])
                 }
               }

               adduct.result1 <- adduct.result1[order(unlist(ichunks))]
               rm(list = "adduct.result")
               gc()

               ##add adduct.result to tags2
               for(j in 1:length(seed.idx)){
                 temp.result <- adduct.result1[[j]]
                 for(k in 1:length(temp.result)){
                   tags2 <- adduct2peakInfo(adduct.result = temp.result[[k]],
                                            mz.tol = mz.tol,
                                            rt.tol = rt.tol1,
                                            weight.mz = 0.8,
                                            weight.rt = 0.2,
                                            peak.info = tags2,
                                            peak.idx = seed.idx[j],
                                            anno.idx = anno.idx[[j]][k],
                                            annotation.idx = annotation.idx)

                 }
               }
               rm(list = "adduct.result1")
               gc()
             }


             # temp.result <- trans2Matrix(tags2 = tags2, base = "annotation")
             # temp.result <- temp.result[temp.result$isotope == "[M]",]
             # temp.mz <- as.numeric(temp.result$mz)
             # temp.formula <- temp.result$Formula
             # temp.formula1 <- mapply(function(x,y){sumFormula(x, y)},
             #                         x = temp.formula, y = temp.adduct)
             #
             # # temp.mass <- sapply(temp.formula1, function(x) Rdisop::getMass(getMolecule(x)))
             # temp.mass <- NULL
             # for(i in 1:length(temp.formula1)){
             #   cat(i);cat(" ")
             #   temp.mass[i] <- Rdisop::getMass(Rdisop::getMolecule(temp.formula1[i]))
             # }
             #
             #
             # temp.adduct <- temp.result$adduct
             # temp.id <- temp.result$to
             #
             # temp.mass <-
             # temp.shift <- as.numeric(adduct.table$massdiff[match(temp.adduct, adduct.table$name)])
             # temp.em <- temp.mass + temp.shift
             # temp.mz.error <- abs(temp.mz - temp.em)*10^6/temp.em
             # temp.mz.error[is.na(temp.mz.error)] <- 100
             # sum(temp.mz.error > 25)
             # temp.idx <- which(temp.mz.error > 25)
             # table(temp.result$type[temp.idx])
             # temp.result[temp.idx,]



             #------------------------------------------------------------------
             ##metabolite annotation
             if(met.annotation){
               cat("\n")
               cat("Metabolite annotation.\n")

               temp.fun <- function(index,
                                    tags2,
                                    peak.mz,
                                    peak.rt,
                                    polarity,
                                    adduct.table,
                                    annotation.idx,
                                    round,
                                    mz.tol,
                                    rt.tol,
                                    cor.tol,
                                    dp.tol,
                                    max.step,
                                    met.database,
                                    metabolic.network,
                                    ms2){
                 annotation.idx <- 1:length(tags2)

                 library(Rcpp, quietly = TRUE,
                         logical.return = FALSE, warn.conflicts = FALSE)
                 # suppressMessages(data(thermo, package = "CHNOSZ", envir = environment()))
                 temp.tags2 <- tags2[index]
                 rm(tags2)
                 gc()
                 result <- lapply(temp.tags2, function(x){
                   seed <- x
                   seed.as.seed.round <-
                     showTags2(list(seed), slot = "as.seed")[[2]][[1]]
                   seed.i <- which(seed.as.seed.round == round)
                   peak.mz.temp <- peak.mz[annotation.idx]
                   peak.rt.temp <- peak.rt[annotation.idx]
                   cor.temp <- rep(1, length(annotation.idx))
                   rm(list=c("peak.mz", "peak.rt"))
                   gc()
                   if(length(seed.i) == 0){
                     return(NULL)
                   }else{
                     temp.result <- lapply(seed.i,function(j){
                       query.name <- seed@name
                       query.id <- seed@annotation[[j]]$to
                       query.charge <- seed@annotation[[j]]$charge
                       query.level <- seed@annotation[[j]]$level
                       query.formula <- stringr::str_trim(seed@annotation[[j]]$Formula)
                       query.adduct <- seed@annotation[[j]]$adduct
                       query.mz <- seed@mz
                       query.rt <- seed@rt
                       metabolite.result <- NULL
                       reaction.step <- 1
                       while(is.null(metabolite.result) & reaction.step <= max.step){
                         metabolite.result <-
                           metAnnotation(metabolite.name = query.name,
                                         metabolite.id = query.id,
                                         formula = query.id,
                                         adduct = query.adduct,
                                         polarity = polarity,
                                         mz = query.mz,
                                         rt = query.rt,
                                         peak.mz = peak.mz.temp,
                                         peak.rt = peak.rt.temp,
                                         ms2 = ms2,
                                         cor = cor.temp,
                                         mz.tol = mz.tol,
                                         rt.tol = rt.tol2,
                                         cor.tol = cor.tol,
                                         step = reaction.step,
                                         metabolite = met.database,
                                         metabolic.network = metabolic.network,
                                         adduct.table = adduct.table,
                                         dp.tol = dp.tol)

                         reaction.step <- reaction.step + 1
                       }
                       metabolite.result
                     })
                     temp.result
                   }
                 })
                 result
               }


               cl <- snow::makeCluster(threads, type = "SOCK")
               # cl <- snow::makeCluster(2, type = "SOCK")
               snow::clusterExport(cl, list = list("showTags2", "metAnnotation",
                                                   "checkElement", "sumFormula",
                                                   "splitFormula","pasteElement",
                                                   "getNeighbor",
                                                   "IdentifyFeature",
                                                   "GetMatchResult",
                                                   "GetSpec2Match",
                                                   "MatchFromTemp",
                                                   "MatchSpec",
                                                   "GetDiffMZppm",
                                                   "GetWeightedInt",
                                                   "GetDotProduct"))

               nc <- length(cl)
               options(warn = -1)
               if(length(seed.idx) < threads){
                 threads <-  1
               }
               ichunks <- split(seed.idx, 1:threads)
               # library(Rcpp)
               system.time(metabolite.result <-
                             snow::clusterApply(cl = cl, ichunks,
                                                fun = temp.fun,
                                                tags2 = tags2,
                                                peak.mz = peak.mz,
                                                peak.rt = peak.rt,
                                                rt.tol = rt.tol2,
                                                mz.tol = mz.tol,
                                                cor.tol = cor.tol,
                                                dp.tol = dp.tol,
                                                max.step = max.step,
                                                polarity = polarity,
                                                adduct.table = adduct.table,
                                                annotation.idx = annotation.idx,
                                                round = round,
                                                met.database = metabolite,
                                                metabolic.network = metabolic.network,
                                                ms2 = ms2))



               snow::stopCluster(cl)

               metabolite.result1 <- metabolite.result[[1]]
               if(threads > 1){
                 for(i in 2:length(metabolite.result)){
                   metabolite.result1 <- c(metabolite.result1, metabolite.result[[i]])
                 }
               }

               metabolite.result1 <- metabolite.result1[order(unlist(ichunks))]

               # a <- lapply(metabolite.result1, function(x) {do.call(rbind, x)})
               #
               # b <- do.call(rbind, a)
               #
               # b.mz <- b$peakMz
               # b.formula <- kegg.compound$Formula[match(b$peakID, kegg.compound$ID)]
               # b.mass <- as.numeric(kegg.compound$Exact.mass[match(b$peakID, kegg.compound$ID)])
               # b.shift <- as.numeric(adduct.table$massdiff[match(b$adduct, adduct.table$name)])
               # b.em <- b.mass + b.shift
               # b.mz.error <- abs(b.mz - b.em)*10^6/b.em
               # sum(b.mz.error > 25)
               ##add metabolite.result to tags2
               for(j in 1:length(seed.idx)){
                 temp.result <- metabolite.result1[[j]]
                 for(k in 1:length(temp.result)){
                   tags2 <- metabolite2peakInfo(metabolite.result = temp.result[[k]],
                                                mz.tol = mz.tol,
                                                rt.tol = rt.tol2,
                                                dp.tol = dp.tol,
                                                weight.mz = 0.25,
                                                weight.rt = 0.25,
                                                weight.dp = 0.5,
                                                peak.info = tags2,
                                                peak.idx = seed.idx[j],
                                                anno.idx = anno.idx[[j]][k],
                                                metabolite = metabolite,
                                                annotation.idx = annotation.idx)
                 }
               }
               rm(list = "metabolite.result1", "metabolite.result")
               gc()
             }
             new.tags <- tags2
             rm(list = "tags2")
             gc()
             new.tags <- new.tags
           }
)













