 # title processTags2
 # description Process the tags2 data after KEGG annotations are added,
 #  including isotope annotation, confidence assignment and redundancy removel.
 # author Xiaotao Shen
 #  \email{shenxt@@sioc.ac.cn}
 # param new.tags2 The tags2 data after KEGG annotations are added.
 # param sample Sample information.
 # param mz.tol The mz tolerance of isotope annotation.
 # param rt.tol The RT tolerance of isotope annotation.
 # param cor.tol The correlation tolerance of isotope annotation.
 # param int.tol The intensity ratio tolerance of isotope annotation.
 # param max.isotope The max number of isotope.
 # param path The directory.
 # return  The tags data with annotation result.
#
# new.tags2 = new.tags.neg
# sample = sample.neg
# mz.tol = mz.tol
# rt.tol = rt.tol1
# cor.tol = cor.tol
# int.tol = int.tol
# max.isotope = max.isotope
# polarity = "negative"
# path = output.path


setGeneric(name = "processTags2",
           def = function(new.tags2,
                          sample,
                          mz.tol = 25,
                          rt.tol = 3,
                          cor.tol = 0,
                          int.tol = 500,
                          max.isotope = 4,
                          polarity = c("positive", "negative"),
                          path = "."){
             polarity <- match.arg(polarity)
             peak.int <- apply(sample, 1, mean)
             rm(sample)
             gc()
  index <- which(showTags2(tags2 = new.tags2, slot = "annotation.len") > 0)
  type <-
    unname(unlist(lapply(showTags2(new.tags2[index], slot = 'annotation.type'), function(x) x[[1]])))
  index <- index[which(type == "keggMatching")]
  if(length(index) == 0) {
    tags.result2 <- NA
    result <- list(new.tags2, list(tags.result2))
    names(result) <- c("new.tags", "tags.result")
    return(result)
    }
  ##find isotopes for kegg matching annotation
  cat("\n")
  cat("Isotopes annotation.\n")
  cat("Seed number:",length(index))
  cat("\n")
  cat("Memory used:")
  temp <- pryr::mem_used()
  print(temp)
  cat("--------------------------------------------------------------------\n")

  # cat("Peak:\n")

  peak.mz <- unname(unlist(lapply(new.tags2, function(x) x@mz)))
  peak.rt <- unname(unlist(lapply(new.tags2, function(x) x@rt)))
  annotation.idx <- 1:length(new.tags2)

  pbapply::pboptions(type="timer", style = 1)

  isotope.result <- pbapply::pblapply(index, function(x) {
    seed <- new.tags2[[x]]
    seed.i <- 1:length(seed@annotation)
      temp.result <- lapply(seed.i,function(j){
        query.name <- seed@name
        query.id <- seed@annotation[[j]]$to
        query.charge <- seed@annotation[[j]]$charge
        query.level <- seed@annotation[[j]]$level
        query.formula <- stringr::str_trim(seed@annotation[[j]]$Formula)
        query.adduct <- seed@annotation[[j]]$adduct
        query.mz <- seed@mz
        query.rt <- seed@rt
        query.int <- peak.int[x]
        peak.mz.temp <- peak.mz[annotation.idx]
        peak.rt.temp <- peak.rt[annotation.idx]
        peak.int.temp <- peak.int[annotation.idx]

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
                                             rt.tol = rt.tol,
                                             mz.tol = mz.tol,
                                             cor.tol = cor.tol,
                                             int.tol = int.tol,
                                             max.isotope = max.isotope)
        isotopes.result
      })
      temp.result
  })

  rm(list = c('peak.mz', "peak.rt", "peak.int", "peak.mz.temp",
              "peak.rt.temp", "peak.int.temp"))
  gc()

  anno.idx <- lapply(index, function(x) {
    seed <- new.tags2[[x]]
  1:length(seed@annotation)
  })

  ##add isotope.result to tags2
  for(j in 1:length(index)){
    # cat(j); cat("\n")
    temp.result <- isotope.result[[j]]
    for(k in 1:length(temp.result)){
      # cat(k); cat(" ")
      new.tags2 <- isotope2peakInfo(isotopes.result = temp.result[[k]],
                                mz.tol = mz.tol,
                                rt.tol = rt.tol,
                                int.tol = int.tol,
                                weight.mz = 0.45,
                                weight.rt = 0.45,
                                weight.int = 0.1,
                                peak.info = new.tags2,
                                peak.idx = index[j],
                                anno.idx = anno.idx[[j]][[k]],
                                annotation.idx = annotation.idx)
    }
  }

  rm(list = "isotope.result", "annotation.idx")
  gc()
  # save(new.tags2, file = file.path(path,"new.tags2"))
  ##
  ##begin to give confidence to annotation of peaks
  result <- trans2Matrix(tags2 = new.tags2, base = "peak")
  # result3 <- result
  # write.table(result3, "Annotation.of.peaks.xlsx", row.names = FALSE, sep = "\t")

  #####give the confidence level of each ID
  cat('\n')
  cat("Assign confidence.\n")
  unique.id <- unique(result$to)

  pbapply::pboptions(type="timer", style = 1)

  id.result <- pbapply::pblapply(unique.id, function(x) {
    temp.id <- x
    temp.idx <- which(result$to == temp.id)
    temp.result <- result[temp.idx,]
    temp.rt <- temp.result$rt
    ##group peaks according to RT
    rt.class <- groupRT(rt = temp.rt, rt.tol = rt.tol)
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
  ##remove redundancy
  # cat("\n")
  # cat("Remove redundancy.\n")
  # tags.result2 <- removeRedundancy(result = result1, path = path, polarity = polarity)
  # load(file.path(path, "redun"))
  # redun2 <- redun
  # save(redun2, file = file.path(path, "redun2"), compress = "xz")
  # unlink(x = file.path(path, "redun"), recursive = TRUE)

  ##filter new.tags2 according to annotation.result
  new.tags2 <- result2Tags(result = result1, tags2 = new.tags2)
  result <- list(new.tags2, result1)
  names(result) <- c("new.tags", "tags.result")
  result
})



