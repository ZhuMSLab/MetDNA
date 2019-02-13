#' @title ms2Annotation
#' @description Annotate peak using MS/MS spectra in in-house database.
#' @author Xiaotao Shen, Yandong Yin
#' \email{shenxt@@sioc.ac.cn}
#' @param ms1.file The name of ms1 peak table. Column 1 is "name", Column 2 is
#' "mz" and column is "rt".
#' @param instrument The instrument you used to acquire data. "AgilentQTOF",
#' "SciexTripleTOF", "OtherQTOF" and "ThermoOrbitrap" are supported.
#' @param ms2.file The vector of names of ms2 files. MS2 file must be mzXML.
#' @param sample.info The name of sample.info.
#' @param mz.tol mz tol for ms1 and ms2 data matching.
#' @param rt.tol RT tol for ms1 and ms2 data matching.
#' @param dp.cutoff The cutoff of dot product. Default is 0.8.
#' @param polarity The polarity of mode.
#' @param path Directory.
#' @param output.path The directory for outputing results.
#' @param ms2.type "mgf", "msp" or "mzXML".
#' @param column The column.
#' @param ce The collision energy.
#' @param ms2.match.plot Output MS2 match plot or not.
#' @return Return the annotation result.
#' @export


# ms2Annotation(ms2.type = "msp", polarity = "negative")


setGeneric(name = "ms2Annotation",
           def = function(ms1.file = "data.csv",
                          instrument = c("AgilentQTOF", "SciexTripleTOF",
                                         "OtherQTOF", "ThermoOrbitrap"),
                          ms2.file,
                          sample.info = "sample.info.csv",
                          mz.tol = 25,
                          rt.tol = 10,
                          dp.cutoff = 0.8,
                          polarity = c("positive", "negative"),
                          path = ".",
                          output.path = file.path(path, "MS2_match_result"),
                          ms2.type = c("mgf", "mzXML", "msp"),
                          column = c("hilic", "rp"),
                          ce = c("30", "10", "20", "35,15", "40", "50"),
                          ms2.match.plot = TRUE
           ){


             polarity <- match.arg(polarity)
             column <- match.arg(column)
             ms2.type <- match.arg(ms2.type)
             instrument <- match.arg(instrument)
             ce <- match.arg(ce)

             output.path <- output.path[1]
             dir.create(output.path)
             dir.create(file.path(output.path, "intermediate_data"))
             ##check ms1.file and ms2.file
             file <- dir(path)
             if(all(file != ms1.file)){stop("There is not ms1.file in you directory.")}

             if(length(grep(ms2.type, file)) == 0) {
               stop("There are not ms2.file in you directory.")
             }

             if(missing(ms2.file)){
               ms2.file <- grep(ms2.type, file, value = TRUE)
             }
             ##
             # temp.idx <- gregexpr("/", output.path)[[1]][length(gregexpr("/", output.path)[[1]])]
             # if(length(temp.idx) == 0){
             # output.name <- output.path
             # output.path <- file.path(path, output.name)
             # }
             #
             # round <- 1
             # while(file.exists(output.path)){
             # output.path <- paste(output.path, round, sep = "_")
             # round <- round + 1
             # }

             # any(dir(path) == output.path)


             # sample.info <- read.csv(file.path(path, "sample.info.csv"),
             #                         stringsAsFactors = FALSE)
             file.copy(from = file.path(path, sample.info),
                       to = file.path(output.path, "intermediate_data"))

             if(polarity == "positive") pol <- "pos"
             if(polarity == "negative") pol <- "neg"


             #load data
             if(instrument == "AgilentQTOF"){
               data("zhuMetlib", envir = environment())
               library = zhuMetlib
               ce <- as.character(as.numeric(ce) + 10)
             }

             if(instrument == "SciexTripleTOF"){
               data("zhuMetlib", envir = environment())
               library = zhuMetlib
             }

             if(instrument == "OtherQTOF"){
               data("zhuMetlib", envir = environment())
               library = zhuMetlib
               ce <- as.character(as.numeric(ce) + 10)
             }

             if(instrument == "ThermoOrbitrap"){
               data("orbitrapMetlib", envir = environment())
               library = orbitrapMetlib
             }

             #------------------------------------------------------------------

             if(polarity == "positive" & column == "hilic"){
               data("hilic.pos", envir = environment())
               adduct <- hilic.pos
             }

             if(polarity == "positive" & column == "rp"){
               data("rp.pos", envir = environment())
               adduct <- rp.pos
             }

             if(polarity == "negative" & column == "hilic"){
               data("hilic.neg", envir = environment())
               adduct <- hilic.neg
             }

             if(polarity == "negative" & column == "rp"){
               data("rp.neg", envir = environment())
               adduct <- rp.neg
             }


             #read MS1 data and get MS2 data
             result <- combineMS1MS2(ms1.file = ms1.file, ms2.file = ms2.file,
                                     ms2.type = ms2.type,
                                     mz.tol = mz.tol, rt.tol = rt.tol,
                                     path = path,
                                     output.path = file.path(output.path, "intermediate_data"))

             ms1 <- as.data.frame(result[[1]])
             ms2 <- result[[2]]

             ##construct spec for MetID
             ms2.name <- unname(unlist(lapply(ms2, function(x) x[[1]][1,])))
             ms1.name <- ms1$name
             temp.idx <- match(ms2.name , ms1.name)
             spec <- vector(mode = "list", length = nrow(ms1))
             names(spec) <- ms1.name
             spec[temp.idx] <- ms2

             spec <- lapply(spec, function(x){
               if(is.null(x)) return(x)
               x[[2]]
             })

             colnames(ms1)[2] <- "mzmed"
             temp.path <- file.path(output.path, "MS2_match_spectra")
             dir.create(temp.path)
             annotation.result <- MetID(info = ms1,
                                        spec = spec,
                                        dp.cutoff = dp.cutoff,
                                        lib = library,
                                        pol = pol,
                                        ce = ce,
                                        adduct = adduct,
                                        d.out = temp.path,
                                        lc = 'HILIC',
                                        ms2.match.plot = ms2.match.plot)

             colnames(annotation.result)[2] <- "mz"
             # save(annotation.result, file = file.path(path, "annotation.result"))
             readr::write_csv(annotation.result,
                              file.path(output.path, "ms2.match.annotation.result.csv"))
             annotation.result <- annotation.result
           })




#------------------------------------------------------------------------------
# title combineMS1MS2
# description Combine MS1 and MS2 data.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param ms1.file The name of ms1 peak table. Column 1 is "name", Column 2 is
# "mz" and column 3 is "rt".
# param ms2.file The vector of names of ms2 files. MS2 file must be mzXML.
# param mz.tol mz tol for ms1 and ms2 data matching.
# param rt.tol RT tol for ms1 and ms2 data matching.
# param path Directory.
# param ms2.type The type of MS2 file, default is mzXML.
# return Return ms1 and ms2 data.


setGeneric(name = "combineMS1MS2",
           function(ms1.file, ms2.file,
                    mz.tol = 25,
                    rt.tol = 10,
                    path = ".",
                    output.path,
                    ms2.type = c("mgf", "mzXML", "msp")){
             #

             ms2.type <- match.arg(ms2.type)
             ##read ms1 data (peak table)
             cat("Read MS1 data.\n")
             cat("\n")
             ms1 <- readr::read_csv(file.path(path, ms1.file), progress = FALSE, col_types = readr::cols())
             ms1.name <- ms1$name
             ms1.mz <- ms1$mz
             ms1.rt <- ms1$rt
             ms1.info <- data.frame(ms1.mz, ms1.rt, ms1.name,
                                    stringsAsFactors = FALSE)

             ##read ms2 data (msp)
             cat("Read MS2 data.\n")
             if(ms2.type == "mzXML"){
               ms2 <- readMZXML(file = file.path(path, ms2.file))
             }

             if(ms2.type == "mgf"){
               ms2 <- readMGF(file = file.path(path, ms2.file))
             }


             if(ms2.type == "msp"){
               ms2 <- readMSP(file = file.path(path, ms2.file))
             }

             ms2.info <- lapply(ms2, function(x) x[[1]])
             ms2.info <- do.call(rbind, ms2.info)
             ms2.info <- as.data.frame(ms2.info)

             ##if ms2.info is from metAnalyzer
             if(length(grep("[A-Ma-z]", ms2.info[1,2])) > 0){
               if(length(grep("POS|NEG", ms1.info[1,3]))>0){
                 ms1.name <- gsub("[_POS|_NEG]", "", ms1.info[,3])
               }
              temp.idx <- match(ms2.info[,2], ms1.name)
              ms2.info[,2] <- ms1.info$ms1.rt[temp.idx]
              ms2.info[,1] <- ms1.info$ms1.mz[temp.idx]
              ms2.info$rt[is.na(ms2.info$rt)] <- 0
              ms2.info$mz[is.na(ms2.info$mz)] <- 0

              ##
              temp <- apply(ms2.info, 1, list)
              temp <- lapply(temp, unlist)

              ms2 <- mapply(function(x, y){
              x[[1]] <- y
              list(x)
              },
              x = ms2,
              y = temp)
             }

             if(max(ms2.info[,2], na.rm = TRUE) < 60){
               ms2.info[,2] <- ms2.info[,2] * 60
             }

             ##match ms1 and ms2
             if(ms2.type == "msp"){
               cat("\n")
               cat("Match MS1 and MS2 data.(mz tolerance ",
                   1," ppm, RT tolerance ", 1, " second).\n", sep = "")
               mz.tol <- 1
               rt.tol <- 1
             }else{
               cat("\n")
               cat("Match MS1 and MS2 data.(mz tolerance ",
                   mz.tol," ppm, RT tolerance ", rt.tol, " second).\n", sep = "")
             }
             match.result <- SXTMTmatch(data1 = ms2.info, data2 = ms1.info,
                                        mz.tol = mz.tol, rt.tol = rt.tol,
                                        rt.error.type = "abs")


             ms2.name <- rep(NA, nrow(ms2.info))
             ms2.name[match.result[,"Index1"]] <- ms1.info[match.result[,"Index2"],"ms1.name"]
             ms2.name <- as.character(ms2.name)

             ##remove no matched MS2
             remove.idx <- which(is.na(ms2.name))
             if(length(remove.idx) > 0){
               ms2.name <- ms2.name[-remove.idx]
               ms2 <- ms2[-remove.idx]
             }

             ###remove duplicated MS2 spectrum, if one peak has more than 1 ms2 spectrum,
             ###select the biggest of the sum intensity of top 10 fragments.
             unique.name <- unique(ms2.name)
             cat("\n")
             cat("Select the most intense MS2 spectrum for one peak.\n")
             remain.idx <- pbapply::pblapply(unique.name, function(name){
               temp.idx <- which(ms2.name == name)
               if(length(temp.idx) == 1) return(temp.idx)
               temp.ms2 <- ms2[temp.idx]
               temp.ms2 <- lapply(temp.ms2, function(x) x[[2]])
               temp.int <- lapply(temp.ms2, function(x) {
                 x <- as.data.frame(x)
                 x <- x[order(x[,2], decreasing = TRUE),]
                 if(nrow(x) >= 10) {sum(x[1:10,2])}else{sum(x[,2])}
               })
               temp.int <- unlist(temp.int)
               temp.idx <- temp.idx[which.max(temp.int)]
               temp.idx
             })

             remain.idx <- sort(unlist(remain.idx))

             ms2 <- ms2[remain.idx]
             ms2.name <- ms2.name[remain.idx]


             ms2 <- mapply(function(info, name){
               NAME <- name
               PRECURSORMZ <- info[[1]][1]
               PRECURSORRT <- info[[1]][2]

               info[[1]] <- data.frame(NAME, PRECURSORMZ,
                                       PRECURSORRT, stringsAsFactors = FALSE)
               info[[1]] <- t(info[[1]])
               # info[[1]] <- matrix(c(NAME, PRECURSORMZ, PRECURSORRT), ncol = 1)
               rownames(info[[1]]) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT')
               list(info)
             },
             info = ms2,
             name = as.character(ms2.name))

             save(ms2, file = file.path(output.path, "ms2"), compress = "xz")
             result <- list(ms1, ms2)
             names(result) <- c("ms1", "ms2")
             return(result)
           })








#---------------------------------------------------------------------------
setGeneric(name = "ConcensusSpec",
           def = function(ms2.all, ppm.merge = 30, mz.bin.min = 0.004,
                          ring.mz.diff.thr = 0.3, ring.int.rel.thr = 0.2,
                          minfrac.vote = 0.5){
             # merge fragments within the same m/z bin across spectra
             ms2.merged <- MergeFragments(ms2.all, ppm.merge, mz.bin.min)
             # quanlity control: removing ring effect and low intensity spectra
             ms2.merged <- RemoveRingEffect(ms2.merged,
                                            mz.diff.thr = ring.mz.diff.thr,
                                            int.rel.thr = ring.int.rel.thr)
             # get concensus spectra with peak voting
             ms2.final <- VoteSpectra(ms2.merged, ms2.all,
                                      minfrac.vote = minfrac.vote,
                                      ppm = ppm.merge,
                                      mz.bin.min = mz.bin.min)
           })



#---------------------------------------------------------------------------
RemoveRingEffect <- function(spec, mz.diff.thr = 0.3, int.rel.thr = 0.2) {
  spec <- spec[order(spec[, 'mz']), , drop = FALSE]
  nr.ring <- nrow(spec) + 1
  mz <- spec[, 'mz']

  mz.diff <- diff(mz)
  idx.mzdiff <- which(mz.diff <= mz.diff.thr)
  if (length(idx.mzdiff) == 0) {
    return(spec)
  }

  nr.ring.possible <- unique(c(idx.mzdiff, idx.mzdiff + 1))
  while (TRUE) {
    idx.int.max <- which.max(spec[nr.ring.possible, 2])
    nr.int.max <- nr.ring.possible[idx.int.max]
    int.thr <- spec[nr.int.max, 2] * int.rel.thr

    mz.diff <- abs(mz[nr.ring.possible[-idx.int.max]] - mz[nr.int.max])
    int <- spec[nr.ring.possible[-idx.int.max], 2]
    nr.ring <- append(nr.ring, nr.ring.possible[-idx.int.max][which(mz.diff <= mz.diff.thr & int <= int.thr)])
    nr.ring.possible <- nr.ring.possible[!nr.ring.possible %in% c(nr.ring, nr.int.max)]
    if (length(nr.ring.possible) == 0) {
      break
    }
  }

  return(spec[-nr.ring, , drop = FALSE])
}


#---------------------------------------------------------------------------
MergeFragments <- function(ms2.all, ppm = 30, mz.bin.min = 0.004) {
  spec.all <- do.call(rbind, ms2.all)
  spec.all <- spec.all[order(spec.all[, 'mz']), , drop = FALSE]

  idx.left <- seq(nrow(spec.all))

  spec.merged <- {}
  while (length(idx.left) > 0 ) {
    idx <- tail(idx.left, 1)
    mz <- spec.all[idx, 'mz']
    mz.range <- c(-1, 0) * max(prod(mz, ppm, 1e-6), mz.bin.min) + mz
    idx.range <- idx.left[spec.all[idx.left, 'mz'] >= mz.range[1] &
                            spec.all[idx.left, 'mz'] <= mz.range[2]]
    spec.tmp <- sapply(c('mz', 'intensity'), function(x) {
      quantile(spec.all[idx.left[idx.range], x], 0.5)
    })
    spec.merged <- rbind(spec.merged, spec.tmp)
    idx.left <- idx.left[-idx.range]
  }
  colnames(spec.merged) <- c('mz', 'intensity')
  rownames(spec.merged) <- NULL
  spec.merged <- spec.merged[order(spec.merged[, 'mz']), , drop = FALSE]

  return(spec.merged)
}



#---------------------------------------------------------------------------
VoteSpectra <- function(ms2.merged, ms2.all, num.smp,
                        ms2.find.method = c('concensus', 'combined'),
                        minfrac.vote = 0.25,
                        ppm = 30,
                        mz.bin.min = 0.004) {
  ms2.find.method <- match.arg(ms2.find.method)

  num.contained <- sapply(ms2.merged[, 'mz'], function(mz) {
    mz.range <- mz + max(prod(mz, ppm, 1e-6), mz.bin.min) * c(-1, 1)
    is.contained <- sapply(ms2.all, function(spec) {
      any(spec[, 'mz'] >= mz.range[1] & spec[, 'mz'] <= mz.range[2])
    })
    sum(is.contained)
  })
  if (ms2.find.method == 'combined') {
    ms2.merged[num.contained/num.smp >= minfrac.vote, , drop = FALSE]
  } else {
    ms2.merged[num.contained/length(ms2.all) >= minfrac.vote, , drop = FALSE]
  }
}







#---------------------------------------------------------------------------
setGeneric(name = "MetID", def = function(info, spec, pol, ce = '30',
                                          lib,
                                          dp.cutoff = 0.8,
                                          adduct,
                                          d.out,
                                          lc = 'HILIC',
                                          ms2.match.plot = TRUE){


  rerun.ms2 <- FALSE
  is.plotMSMS <- FALSE
  dt.peaktable <- info
  # colnames(dt.peaktable) <- c('mzmed', colnames(info)[-1])
  colnames(dt.peaktable)[2] <- "mzmed"
  dt.peaktable$mzmed <- as.numeric(dt.peaktable$mzmed)
  # dt.peaktable <- cbind('name' = seq(nrow(dt.peaktable)), dt.peaktable)
  # fn.lib <- 'zhuMetlib.RData'
  # fn.lib <- 'zhuMetlib.RData'
  lc <- 'HILIC'
  # fn.adduct.lc <- switch(pol,
  #                        'pos' = paste0(lc, '_POS.csv'),
  #                        'neg' = paste0(lc, '_NEG.csv'))
  # fp.adduct.lc <- file.path('.', fn.adduct.lc)
  # fp.adduct.lc <- adduct.table.hilic
  # lib <- LoadData(fn.lib)
  lib.meta <- lib$meta[[pol]][[ce]]
  lib.name <- apply(lib.meta, 1, function(info) {
    nr <- which(lib$meta$compound[, 'labid'] == as.character(info['labid']))
    lib$meta$compound[nr, 'name']
  })
  lib.meta$name <- lib.name

  lib.spec <- lib$compound[[pol]]

  # adduct <- read.csv(fp.adduct.lc,
  #                    stringsAsFactors = FALSE)
  lib.mz <- t(sapply(lib.meta[, 'mz'], function(mz) {
    mz <- as.numeric(mz)
    apply(adduct, 1, function(info.adduct) {
      x <- gsub('\\(', '',
                gsub('M.*', '', info.adduct['adduct']))
      xm <- ifelse(x == '', 1, as.numeric(x))
      xm * mz + as.numeric(info.adduct['mz'])
    })
  }))
  colnames(lib.mz) <- adduct$adduct
  idx.peak.match <- unname(which(!sapply(spec, is.null)))
  # ShowSeperateLines('Identifying metabolites with library ...')
  cat("\n")
  cat('Identify metabolites with MS/MS library.\n')
  # match library reverse and forward
  for (direction in c('reverse', 'forward')) {
    cat("\n")
    cat('Search ', direction, '.\n', sep = "")
    # fn.skip <- paste0('Search_', direction, '.Rda')
    #
    pbapply::pboptions(style = 1)
    id.peaks.list <- pbapply::pblapply(idx.peak.match, ParIdentifyFeatureMet,
                                       dt.peaktable,
                                       # spec[idx.peak.match],
                                       spec,
                                       lib.meta, lib.spec, lib.mz,
                                       ce,
                                       search.scope = 'ms1.based',
                                       ppm.ms1match = 25,
                                       ppm.ms2match = 35,
                                       cutoff = dp.cutoff,
                                       top = 10, # report top 10 satisfactories
                                       top.below = 5, # if no satisfactory, report top 5
                                       is.tune.ms2.exp = TRUE,
                                       is.tune.ms2.lib = FALSE,
                                       mz.range.ms2 = NULL,
                                       is.include.precursor = ifelse(pol == 'neg', TRUE, FALSE),
                                       weight.mz = 0,
                                       weight.int = 1,
                                       int.ms2.min.relative = 0.01,
                                       is.apply.ms2.min.relative = TRUE,
                                       noise.ms2 = 10,
                                       snthr = 3,
                                       ppm.sanity.check = 100,
                                       is.sanity.check = FALSE,
                                       direction = direction)
    # save(id.peaks.list, file = fn.skip)


    num.hits <- rep(0, nrow(dt.peaktable))
    num.hits[idx.peak.match] <- sapply(id.peaks.list, function(id) {
      if (all(is.na(id))) {
        0
      } else {
        nrow(id)
      }
    })

    info.hits <- rep('', nrow(dt.peaktable))
    info.hits[idx.peak.match] <- sapply(id.peaks.list, function(id){
      if (!all(is.na(id)) & length(id) > 0) {
        paste(paste('score{', round(id$score, 6), '}',
                    'adduct{', id$adduct, '}',
                    'name{', id$name, '}',
                    'labid{', id$labid, '}', sep = ''),
              collapse = ';')
      } else {
        ''
      }
    })
    if (all(num.hits == 0)) {
      return()
    }

    # if (is.plotMSMS) {
    # d.spec <- file.path('MetLibMatch/MSMSfigures', d.out)
    d.spec <- d.out
    dir.create(d.spec, recursive = T)

    info.peak.plot <- data.frame('nr.peaktable' = idx.peak.match[which(num.hits[idx.peak.match] > 0)],
                                 'idx.peaklist' = which(num.hits[idx.peak.match] > 0),
                                stringsAsFactors = FALSE)

    if(ms2.match.plot){
      cat("\n")
      cat('Plot spectra match results.\n')
      PlotIDResults(id.peaks.list, info.peak.plot,
                    spec, dt.peaktable,
                    lib.spec,
                    lib.meta,
                    ce,
                    polarity = pol,
                    direction = direction,
                    d.spec = d.spec,
                    width = 20, height = 7,
                    is.include.precursor = FALSE,
                    is.tune.ms2.exp = TRUE,
                    is.tune.ms2.lib = FALSE,
                    mz.range.ms2 = NULL)
    }
    # }
    assign(paste0('id.peaks.info.', direction),
           cbind('nhits' = num.hits, 'hits' = info.hits))
  }

  id.peaks.info <- cbind(id.peaks.info.reverse, id.peaks.info.forward)
  colnames(id.peaks.info) <- c('nhits.reverse', 'hits.reverse', 'nhits.forward', 'hits.forward')
  # if(!file.exists('MetLibMatch')) dir.create('MetLibMatch')
  # fn <- file.path('MetLibMatch', paste0(gsub('/', '_', d.out), '.csv'))
  # write.csv(cbind(dt.peaktable, id.peaks.info), fn, row.names = FALSE)
  annotation.result <- data.frame(dt.peaktable, id.peaks.info,
                                  stringsAsFactors = FALSE)
  return(annotation.result)
})


#---------------------------------------------------------------------------
setGeneric(name = "ParIdentifyFeatureMet",
           def = function(idx.pk,
                          dt.peaktable,
                          ms.assigned,
                          lib.meta,
                          lib.spec,
                          lib.mz,
                          ce,
                          search.scope,
                          ppm.ms1match = 30,
                          ppm.ms2match = 30,
                          cutoff = 0.6,
                          top = 10, # report top 10 satisfactories
                          top.below = 5, # if no satisfactory, report top 5
                          is.tune.ms2.exp = FALSE,
                          is.include.precursor = is.include.precursor,
                          ...){
             # ...:
             # is.tune.ms2.lib = is.tune.ms2.lib,
             # weight.mz = weight.mz,
             # weight.int = weight.int,
             # int.ms2.min.relative = int.ms2.min.relative,
             # is.apply.ms2.min.relative = is.apply.ms2.min.relative,
             # noise.ms2 = noise.ms2,
             # snthr = 3,
             # ppm.sanity.check = 100,
             # is.sanity.check = FALSE,
             # direction = direction
             # suppressMessages(require(MetMatch))
             # cat(idx.pk, '\n')
             pk.precursor <- dt.peaktable[idx.pk, ]
             pk.mz <- pk.precursor$mzmed
             pk.spec <- ms.assigned[[idx.pk]]

             if (is.tune.ms2.exp) {
               pk.spec <- TuneMS2(pk.spec, pk.mz, is.include.precursor = is.include.precursor,...)
             }
             if (length(pk.spec) == 0) {
               return(NA)
             }
             pk.mz.range <- GetRangePPM(pk.mz, ppm.ms1match)

             switch(search.scope,
                    'ms1.based' = {
                      idx.lib.match <- apply(lib.mz, 2, function(mz.col) {
                        which(mz.col >= pk.mz.range[1] & mz.col <= pk.mz.range[2])
                      })

                      idx.remove <- which(sapply(idx.lib.match, length) == 0)
                      if (length(idx.remove) > 0) {
                        idx.lib.match <- idx.lib.match[-idx.remove]
                      }
                    },
                    'all' = {
                      idx.lib.match <- list(seq(nrow(lib.meta)))
                    })

             if (length(idx.lib.match) == 0) {
               return(NA)
             }
             #
             id.list.adduct <- lapply(seq(idx.lib.match), function(idx) {
               idx.match <- idx.lib.match[[idx]]
               ids <- lib.meta[idx.match, 'labid']
               mz.lib <- lib.meta[idx.match, 'mz']
               id.list <- lapply(seq_along(ids), function(idx) {
                 id <- ids[idx]
                 lib.spec.match <- lib.spec[[id]][[ce]]
                 if (!is.include.precursor) {
                   mz.cutoff <- GetRangePPM(mz.lib[idx], 20)[2]
                   lib.spec.match <- lib.spec.match[lib.spec.match[, 'mz'] <= mz.cutoff, ,
                                                    drop = FALSE]
                 }
                 IdentifyFeature(pk.spec, lib.spec.match,
                                 ppm.ms2match = ppm.ms2match,
                                 is.include.precursor = is.include.precursor,...)
               })
               #

               sc <- do.call(c, id.list)
               sc[is.nan(sc) | is.na(sc)] <- 0
               id.result <- data.frame(sc, lib.meta[idx.match, c(1, 5),
                                                    drop = FALSE], stringsAsFactors = FALSE)
               id.result$adduct <- rep(names(idx.lib.match[idx]), length(ids))
               id.result
             })
             id.result <- do.call(rbind, id.list.adduct)
             switch(as.character(ncol(id.result)),
                    '3' = {
                      colnames(id.result) <- c('score', 'labid', 'name')
                    },
                    '4' = {
                      colnames(id.result) <- c('score', 'labid', 'name', 'adduct')
                    })

             id.result.cutoff <- id.result[id.result[, 'score'] >= cutoff, , drop = FALSE]
             id.result.cutoff <- id.result.cutoff[order(id.result.cutoff[, 'score'], decreasing = TRUE), , drop = FALSE]
             id.result.cutoff[, 'score'] <- round(id.result.cutoff[, 'score'], 4)

             if (length(id.result.cutoff) > 0) {
               id.result.output <- head(id.result.cutoff, top)
             } else {
               id.result.output <- head(id.result, top.below)
             }

             if (nrow(id.result.output) == 0) {
               id.result.output <- NA
             }

             return(id.result.output)
           })



#---------------------------------------------------------------------------

setGeneric(name = "PlotIDResults",
           def = function(id.peaks.list, info.peak.plot, ms.assigned, dt.peaktable,
                          lib.spec, lib.meta, ce, polarity = c("positive", "negative"),
                          direction = c("reverse", "forward"), d.spec = "ms2Result/MSMSfigures",
                          width = 20, height = 7, is.include.precursor = TRUE, is.tune.ms2.exp = TRUE,
                          is.tune.ms2.lib = FALSE,
                          ...){
             polarity <- match.arg(polarity)
             direction <- match.arg(direction)
             col.plot <- c(lib = "salmon", exp = "lightseagreen", filtered = "gray")
             for (i in rev(seq(nrow(info.peak.plot)))) {
               cat(i, " ")
               apply(id.peaks.list[[info.peak.plot$idx.peaklist[i]]],
                     1, function(r.id) {
                       id <- as.character(r.id["labid"])
                       score <- round(as.numeric(r.id["score"]), 3)
                       spec.lib <- spec.lib.all <- lib.spec[[id]][[ce]]
                       if (!is.include.precursor) {
                         mz.precursor <- as.numeric(lib.meta[lib.meta[,
                                                                      "labid"] == id, "mz"])
                         nr.remove <- which(spec.lib[, "mz"] == mz.precursor)
                         if (length(nr.remove) > 0) {
                           spec.lib <- spec.lib[-nr.remove, , drop = FALSE]
                         }
                       }
                       if (is.tune.ms2.lib) {
                         spec.lib <- TuneMS2(spec.lib, mz.precursor,
                                             is.include.precursor = is.include.precursor,
                                             ...)
                       }
                       spec.lib.filtered <- spec.lib.all[which(!spec.lib.all[,
                                                                             "mz"] %in% spec.lib[, "mz"]), , drop = FALSE]
                       nr.peaktable <- info.peak.plot$nr.peaktable[i]
                       spec.exp <- spec.exp.all <- ms.assigned[[nr.peaktable]]
                       if (is.tune.ms2.exp) {
                         spec.exp <- TuneMS2(spec.exp, dt.peaktable[nr.peaktable,
                                                                    "mzmed"], is.include.precursor = is.include.precursor,
                                             ...)
                       }
                       spec.exp.filtered <- spec.exp.all[which(!spec.exp.all[,
                                                                             "mz"] %in% spec.exp[, "mz"]), , drop = FALSE]
                       spec2match <- GetSpec2Match(spec.exp, spec.lib,
                                                   direction = direction)
                       nr.matched <- which(spec2match$exp[, "intensity"] >
                                             0 & spec2match$lib[, "intensity"] > 0)
                       spec.matched <- lapply(spec2match, function(spec) {
                         spec[nr.matched, , drop = FALSE]
                       })
                       d.plot <- file.path(d.spec, paste(dt.peaktable[nr.peaktable,
                                                                      "name"], direction, sep = "_"))
                       dir.create(d.plot)
                       cmpd.replaced <- paste(sapply(strsplit(r.id["name"],
                                                              "")[[1]], function(x) {
                                                                switch(x, `:` = "_", `/` = "-", x)
                                                              }), collapse = "")
                       fn.plot <- switch(as.character(length(r.id)),
                                         `3` = file.path(d.plot, paste(score, ",",
                                                                       cmpd.replaced, ".pdf", sep = "")), `4` = file.path(d.plot,
                                                                                                                          paste(score, ",", cmpd.replaced, ",", r.id["adduct"],
                                                                                                                                ".pdf", sep = "")))
                       range.mz <- range(c(spec.lib.all[, "mz"], spec.exp.all[,
                                                                              "mz"]))
                       range.int <- c(-1, 1)
                       pdf(file = fn.plot, height = 7, width = 20,
                           family = "mono")

                       par(mar = c(5,5,4,2))

                       plot(range.mz, range.int, type = "n", main = r.id["name"],
                            xlab = "m/z", ylab = "Relative intensity",  cex.lab = 1.5,
                            cex.axis = 1.3, cex.main = 1.5)

                       abline(h = 0, col = "black")
                       ref.lib <- max(spec.lib.all[, "intensity"])
                       points(NormalizeSpec(spec.lib, ref.lib, "down"),
                              type = "h", col = col.plot["lib"])
                       if (nrow(spec.lib.filtered) > 0) {
                         points(NormalizeSpec(spec = spec.lib.filtered,
                                              ref = ref.lib, pos = "down"), type = "h",
                                col = col.plot["filtered"])
                       }
                       ref.exp <- max(spec.exp.all[, "intensity"])
                       points(NormalizeSpec(spec.exp, ref.exp, "top"),
                              type = "h", col = col.plot["exp"])
                       if (nrow(spec.exp.filtered) > 0) {
                         points(NormalizeSpec(spec = spec.exp.filtered,
                                              ref = ref.exp, pos = "top"), type = "h",
                                col = col.plot["filtered"])
                       }
                       points(NormalizeSpec(spec.matched$lib, ref.lib,
                                            "down"), type = "p", pch = 20, col = col.plot["lib"])
                       points(NormalizeSpec(spec.matched$exp, ref.exp,
                                            "top"), type = "p", pch = 20, col = col.plot["exp"])
                       if (ce == "spec") {
                         legend("bottomleft", legend = c(paste("Name:",
                                                               r.id["name"]), paste("Polarity:", polarity)),
                                pch = NA, bty = "n")
                       }
                       else {
                         legend("bottomleft", legend = c(paste("LabID:",
                                                               r.id["labid"]), paste("Polarity:", polarity),
                                                         paste("CE:", ce)), pch = NA, bty = "n", cex = 1.3)
                       }
                       legend("topleft", cex = 1.3, legend = c(paste("Score:",
                                                                     score), paste("Matched peaks:", c("data",
                                                                                                       "library"))), col = c("white", col.plot[c("exp",
                                                                                                                                                 "lib")]), pch = list(NA, 20, 20), bty = "n")
                       dev.off()
                     })
             }
           })




#-----------------------------------------------------------------------------
setGeneric(name = "removeNoise",
           def = function(spec, mz.tol = 30){
             spec <- matrix(spec, ncol = 2)
             colnames(spec) <- c("mz", "intensity")
             if(nrow(spec) == 1) return(spec)
             spec <- spec[order(spec[,1]),]
             mz <- as.numeric(spec[,1])

             new.spec <- NULL
             i = 1
             while(i < length(mz)){
               # cat(i); cat(" ")
               temp.mz <- mz[i]
               mz.error <- abs(temp.mz - mz)*10^6/ifelse(temp.mz >= 400, temp.mz, 400)
               temp.idx <- which(mz.error <= mz.tol)
               temp.spec <- spec[temp.idx,]
               if(length(temp.idx) == 1) {
                 new.spec <- rbind(new.spec, temp.spec)
                 i <- max(temp.idx) + 1
                 next()
               }
               temp.mz <- median(temp.spec[,1])
               temp.int <- max(temp.spec[,2])
               new.spec <- rbind(new.spec, c(temp.mz, temp.int))
               # spec <- spec[-temp.idx,]
               i <- max(temp.idx) + 1
             }

             row.names(new.spec) <- NULL
             colnames(new.spec) <- c("mz", "intensity")
             return(new.spec)
})




#-----------------------------------------------------------------------------
# title readMGF
# description Read mgf data.
# author Xiaotao Shen, Yandong Yin
#  \email{shenxt@@sioc.ac.cn}
# param file The vector of names of ms2 files. MS2 file must be mgf.
# return Return ms2 data.
# export

setGeneric(name = "readMGF",
           def = function(file){
             pbapply::pboptions(style = 1)
             # cat("Reading MS2 data (step1)\n")
             # mgf.data.list <- pbapply::pblapply(file, ListMGF)
             ms2 <- pbapply::pblapply(file, function(mgf.data) {
               mgf.data <- ListMGF(mgf.data)
               # nl.spec <- grep('^\\d', mgf.data)
               nl.spec <- lapply(mgf.data, function(x) grep('^\\d', x))
               info.mz <- lapply(mgf.data, function(x) grep('^PEPMASS', x, value = T))
               info.rt <- lapply(mgf.data, function(x) grep('^RTINSECONDS', x, value = T))

               info.mz <- unlist(info.mz)
               #for orbitrap data, the intensity of precursor ion should be removed
               info.mz <- unlist(lapply(strsplit(x = info.mz, split = " "), function(x) x[1]))
               info.mz <- as.numeric(gsub(pattern = "\\w+=", "", info.mz))
               info.rt <- unlist(info.rt)
               info.rt <- as.numeric(gsub(pattern = "\\w+=", "", info.rt))

               spec <- mapply(function(x, y){
                 do.call(rbind, strsplit(x[y], split = " "))
               },
               x = mgf.data,
               y = nl.spec)

               spec <- lapply(spec, function(x){
                 temp <- cbind(as.numeric(x[,1]),as.numeric(x[,2]))
                 temp <- matrix(temp, ncol = 2)
                 if(nrow(temp) > 0) temp <- temp[temp[,2] >= max(temp[,2])*0.01,]
                 temp <- matrix(temp, ncol = 2)
                 colnames(temp) <- c("mz", "intensity")
                 temp
               })

               ms2 <- mapply(function(x,y,z){
                 info <- c(y, z)
                 names(info) <- c("mz", "rt")
                 spectrum <- as.matrix(x)
                 temp <- list(info, spectrum)
                 names(temp) <- c("info", "spec")
                 list(temp)
               },
               x = spec,
               y = info.mz,
               z = info.rt)

               ms2

             })


             spec.info <- ms2[[1]]
             if(length(ms2) > 1){
               for(i in 2:length(ms2)){
                 spec.info <- c(spec.info, ms2[[i]])
               }
             }

             remove.idx <- which(unlist(lapply(spec.info, function(x) nrow(x[[2]]))) == 0)
             if(length(remove.idx) != 0) spec.info <- spec.info[-remove.idx]
             # ##remove noise
             # cat("\n")
             # cat("Remove noise of MS/MS spectra...\n")
             # spec.info <- pbapply::pblapply(spec.info, function(x){
             #   temp.spec <- x[[2]]
             #   temp.spec <- removeNoise(temp.spec)
             #   x[[2]] <- temp.spec
             #   x
             # })

             spec.info <- spec.info
           })


#----------------------------------------------------------------------------
setGeneric(name = "ListMGF",
           def = function(file){
             mgf.data <- readLines(file)
             nl.rec.new <- 1
             idx.rec <- 1
             rec.list <- list()
             for(nl in 1:length(mgf.data))
             {
               if(mgf.data[nl]=="END IONS")
               {
                 rec.list[idx.rec] <- list(Compound = mgf.data[nl.rec.new : nl])
                 nl.rec.new <- nl + 1
                 idx.rec <- idx.rec + 1
               }
             }
             rec.list
})

setGeneric(name = "CheckInRange", def = function(targets, range){
  targets >= range[1] & targets <= range[2]
})





#-----------------------------------------------------------------------------
# title readMZXML
# description Read mzXML data.
# author Xiaotao Shen
#  \email{shenxt@@sioc.ac.cn}
# param file The vector of names of ms2 files. MS2 file must be mzXML.
# return Return ms2 data.
# export

setGeneric(name = "readMZXML", function(file){
  # cat("Open MS2 file.\n")
  mzxml.data <- lapply(file, function(x) {
    mzR::openMSfile(x)
  })
  # cat("Extract MS2 information\n")
  mzxml.info <- lapply(mzxml.data, function(x){
    mzR::header(x)
  })
  # cat("Extract MS2 spectrum\n")
  pbapply::pboptions(type = "timer", style = 1)
  mzxml.peak <- pbapply::pblapply(mzxml.data, function(x){
    mzR::peaks(x)
  })


  ms2 <- mapply(function(info, peak){
    ms2.idx <- which(info$msLevel == 2)
    ms2.info <- info[ms2.idx, c("precursorMZ", "retentionTime")]
    ms2.info <- apply(ms2.info, 1, list)
    ms2.spec <- peak[ms2.idx]

    temp.ms2 <- mapply(function(x, y){
      names(x[[1]]) <- c("mz", "rt")
      colnames(y) <- c("mz", "intensity")
      if(nrow(y) > 0) y <- y[y[,2] >= max(y[,2])*0.01,]
      y <- matrix(y, ncol = 2)
      temp <- list(x[[1]], y)
      names(temp) <- c("info", "spec")
      list(temp)
    }, x = ms2.info,
    y = ms2.spec)

    list(temp.ms2)
  },
  info = mzxml.info,
  peak = mzxml.peak)

  spec.info <- ms2[[1]]
  if(length(ms2) > 1){
    for(i in 2:length(ms2)){
      spec.info <- c(spec.info, ms2[[i]])
    }
  }


  remove.idx <- which(unlist(lapply(spec.info, function(x) nrow(x[[2]]))) == 0)
  if(length(remove.idx) != 0) spec.info <- spec.info[-remove.idx]

  # ##remove noise
  # cat("\n")
  # cat("Remove noise of MS/MS spectra...\n")
  # spec.info <- pbapply::pblapply(spec.info, function(x){
  # temp.spec <- x[[2]]
  # temp.spec <- removeNoise(temp.spec)
  # x[[2]] <- temp.spec
  # x
  # })

  spec.info <- spec.info
})



#---------------------------------------------------------------------------
#title ReadMSP
#aliases MSP file reader
#description  Read a MSP file and return a list of spectra for all feature
# with feature information
#param file path of the msp file
setGeneric('readMSP', function(file) {
  msp.data <- readLines(file)
  # n.tot <- length(msp.data)
  n.null <- which(msp.data == '')

  temp.idx1 <- c(1, n.null[-length(n.null)])
  temp.idx2 <- n.null - 1

  temp.idx <- data.frame(temp.idx1, temp.idx2,
                         stringsAsFactors = FALSE)
  temp.idx <- apply(temp.idx, 1, list)

  temp.idx <- lapply(temp.idx, unlist)

  # n.spec <- which(grepl('^\\d', msp.data))
  # n.info <- seq(n.tot)[-c(n.spec, n.null)]

  pbapply::pboptions(style = 1)
  info.spec <- pbapply::pblapply(temp.idx, function(idx) {

    temp.msp.data <- msp.data[idx[1]:idx[2]]

    temp.msp.data <- temp.msp.data[temp.msp.data != ""]
    info.idx <- grep("[A-Za-z]", temp.msp.data)
    temp.info <- temp.msp.data[info.idx]
    temp.info <- strsplit(temp.info, split = ":")
    temp.info <- do.call(rbind, temp.info)
    temp.info <- data.frame(temp.info,
                            stringsAsFactors = FALSE)
    temp.info[,2] <- stringr::str_trim(temp.info[,2])
    colnames(temp.info) <- rownames(temp.info) <- NULL
    rownames(temp.info) <- temp.info[,1]
    temp.info <- temp.info[,-1,drop = FALSE]

    temp.spec <- temp.msp.data[-info.idx]

    if(length(temp.spec) != 0){
      if(length(grep(" ", temp.spec[1])) == 1){
        temp.spec <- strsplit(temp.spec, split = ' ')
      }

      if(length(grep("\t", temp.spec[1])) == 1){
        temp.spec <- strsplit(x = temp.spec, split = "\t")
      }

      temp.spec <- do.call(rbind, temp.spec)
      temp.spec <- data.frame(temp.spec,
                              stringsAsFactors = FALSE)
      colnames(temp.spec) <- c('mz', 'intensity')
      rownames(temp.spec) <- NULL
      temp.spec$mz <- as.numeric(as.character(temp.spec$mz))
      temp.spec$intensity <- as.numeric(temp.spec$intensity)
      temp.spec <- temp.spec[temp.spec$intensity != 0,]
    }else{
      temp.spec <- NULL
    }

    list('info' = temp.info,
         'spec' = temp.spec)
  })

  mz.idx <- grep("[Mm][Zz]", rownames(info.spec[[1]][[1]]))
  rt.idx <- grep("Time|TIME|time|RT|rt|Rt", rownames(info.spec[[1]][[1]]))

  ##fix bug in msp data from metAnalyzer
  if(length(rt.idx)==0){
    cat("The msp data are from MetAnalyzer software.\n")
    rt.idx <- grep("NAME|Name|name", rownames(info.spec[[1]][[1]]))
    ##rt.idx is the name of peak
    info.spec <- lapply(info.spec, function(x){
      info <- x[[1]]
      mz <- as.numeric(info[mz.idx, 1])
      rt <- as.character(info[rt.idx, 1])
      info <- c(mz, rt)
      names(info) <- c("mz", "rt")
      x[[1]] <- info
      x
    })
  }else{
    info.spec <- lapply(info.spec, function(x){
      info <- x[[1]]
      mz <- as.numeric(info[mz.idx, 1])
      rt <- as.numeric(info[rt.idx, 1])
      info <- c(mz, rt)
      names(info) <- c("mz", "rt")
      x[[1]] <- info
      x
    })
  }


  remove.idx <- which(unlist(lapply(info.spec, function(x) is.null(x[[2]]))))
  if(length(remove.idx) > 0){
    info.spec <- info.spec[-remove.idx]
  }

  info.spec <- info.spec
})



#------------------------------------------------------------------------------
  setGeneric(name = "WriteMSP",
             def = function(info, fn.pre, spec.all){
               fn.save <- paste0(fn.pre, '_spectra.msp')
               #

               sink(fn.save)
               for (idx in seq(nrow(info))) {
                 if (!is.null(spec.all[[idx]])) {
                   if (nrow(spec.all[[idx]]) > 0) {
                     mz <- info[idx, 'Mass']
                     spec <- spec.all[[idx]]
                     cat('IDX: ', idx, '\n', sep = '')
                     cat('PRECURSORMZ: ', mz, '\n', sep = '')
                     cat('Num Peaks: ', nrow(spec), '\n', sep = '')
                     for (nr in seq(nrow(spec))) {
                       cat(paste(spec[nr, ], collapse = ' '), '\n', sep = '')
                     }
                     cat('\n')
                   }
                 }
               }
               sink()
  })


#------------------------------------------------------------------------------
setGeneric(name = "GetPpmRange",
           def = function(ref, tol){
             ref + c(-1, 1) * ref * tol * 1e-6
})



setGeneric(name = "LoadData", def = function(file, keep.name = FALSE, env){
  if (missing(env))
    env <- new.env()
  b <- load(file, envir = env)
  if (keep.name | length(b) > 1) {
    r <- lapply(b, function(b1) env[[b1]])
    names(r) <- b
    r
  }
  else {
    env[[b]]
  }
})

