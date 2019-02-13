###IdentifyFeature
#'@title IdentifyFeature
#'@description IdentifyFeature
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param spec.exp spec.exp.
#'@param spec.lib spec.lib.
#'@param is.include.precursor is.include.precursor.
#'@param is.tune.ms2.lib is.tune.ms2.lib.
#'@param direction direction.
#'@param ... other parameters
#'@export

IdentifyFeature <- function(spec.exp,
                            spec.lib,
                            is.include.precursor = FALSE,
                            is.tune.ms2.lib = FALSE,
                            direction = "forward",
                            ...) {
  # ()
  if (is.tune.ms2.lib) {
    spec.lib <- TuneMS2(spec.lib, is.include.precursor = is.include.precursor, ...)
  }
  if (length(spec.lib) == 0) {
    return(NA)
  }


  return(GetMatchResult(spec.exp, spec.lib, direction = direction, ...))
}

#-----------------------------------------------------------------------------
###
###GetMatchResult
#'@title GetMatchResult
#'@description GetMatchResult
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param spec.exp spec.exp
#'@param spec.lib spec.lib
#'@param weight.int weight.int
#'@param weight.mz weight.mz
#'@param ppm.ms2match ppm.ms2match
#'@param mz.ppm.thr mz.ppm.thr
#'@param ppm.sanity.check ppm.sanity.check
#'@param is.sanity.check is.sanity.check
#'@param direction direction
#'@param ... other parameters
#'@export
#'

GetMatchResult <- function(spec.exp, spec.lib,
                           weight.int = 1,
                           weight.mz = 0,
                           ppm.ms2match = 30,
                           mz.ppm.thr = 400,
                           ppm.sanity.check = 100,
                           is.sanity.check = FALSE,
                           direction = direction,
                           ...) {
  # ()
  if (is.sanity.check) {
    switch(direction,
           'forward' = {
             if (any(c(GetDiffMZppm(spec.exp[, 'mz']), GetDiffMZppm(spec.lib[, 'mz'])) <= ppm.sanity.check)) {
               stop('Difference between m/z is too small!!')
             }
           },
           'reverse' = {
             if (any(GetDiffMZppm(spec.lib) <= ppm.sanity.check)) {
               stop('Difference between m/z is too small!!')
             }
           },
           stop('Error setup for parameter: direction!!!'))
  }
  #
  spec2match <- GetSpec2Match(spec.exp, spec.lib,
                              ppm.ms2match = ppm.ms2match,
                              mz.ppm.thr = mz.ppm.thr,
                              direction = direction)
  int.weighted.pk  <- GetWeightedInt(spec2match$exp, weight.mz, weight.int)
  int.weighted.lib <- GetWeightedInt(spec2match$lib, weight.mz, weight.int)
  match.score <- GetDotProduct(int.weighted.pk, int.weighted.lib)
  attr(match.score, 'spec') <- spec2match
  attr(match.score, 'spec.compared') <- cbind('exp' = int.weighted.pk,
                                              'lib' = int.weighted.lib)
  return(match.score)
}



#-----------------------------------------------------------------------------
###GetSpec2Match
#'@title GetSpec2Match
#'@description GetSpec2Match
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param spec.exp spec.exp
#'@param spec.lib spec.lib
#'@param ppm.ms2match ppm.ms2match
#'@param mz.ppm.thr mz.ppm.thr
#'@param direction direction
#'@export
#'

GetSpec2Match <- function(spec.exp, spec.lib,
                          ppm.ms2match = 30,
                          mz.ppm.thr = 400,
                          direction = c('reverse', 'forward')) {
  #
  direction = match.arg(direction)
  mz.pool   <- sort(c(spec.exp[, 'mz'], spec.lib[, 'mz']))
  spec.temp <- cbind('mz' = mz.pool, 'intensity' = 0)

  spec.exp.temp  <- MatchFromTemp(spec.exp, spec.temp)
  spec.lib.temp <- MatchFromTemp(spec.lib, spec.temp)

  # combine nearby peaks
  pk.spec  <- MatchSpec(spec.exp.temp,  ppm.ms2match = ppm.ms2match, mz.ppm.thr = mz.ppm.thr)
  lib.spec <- MatchSpec(spec.lib.temp, ppm.ms2match = ppm.ms2match, mz.ppm.thr = mz.ppm.thr)

  if (direction == 'reverse') {
    idx.rm <- which(lib.spec[, 'intensity'] == 0)
    if (length(idx.rm > 0)) {
      pk.spec  <- pk.spec[-idx.rm, , drop = FALSE]
      lib.spec <- lib.spec[-idx.rm, , drop = FALSE]
    }
  }

  return(list('exp' = pk.spec, 'lib' = lib.spec))
}



#-----------------------------------------------------------------------------
###MatchFromTemp
#'@title MatchFromTemp
#'@description MatchFromTemp
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param spec spec
#'@param temp temp
#'@export
#'

MatchFromTemp <- function(spec, temp) {
  temp[match(spec[, 'mz'], temp[, 'mz']), 'intensity'] <- spec[, 'intensity']
  temp
}




#-----------------------------------------------------------------------------
###MatchSpec
#'@title MatchSpec
#'@description MatchSpec
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param spec spec
#'@param ppm.ms2match ppm.ms2match
#'@param mz.ppm.thr mz.ppm.thr
#'@export
#'


MatchSpec <- function(spec, ppm.ms2match = 30, mz.ppm.thr = 400) {
  while (TRUE) {
    mz.diff.ppm <- GetDiffMZppm(spec[, 'mz'], mz.ppm.thr = mz.ppm.thr)
    idx <- which(mz.diff.ppm < ppm.ms2match)
    if (length(idx) > 0) {
      i <- tail(idx, 1)
      j <- which.max(spec[c(i, i + 1), 'intensity'])
      spec[i, 'intensity'] <- spec[i + j - 1, 'intensity']
      i2 <- i + 1
      spec[i, 'mz'] <- spec[i2, 'mz']
      spec <- spec[-i - 1, , drop = FALSE]
    } else {
      break
    }
  }
  return(spec)
}



#-----------------------------------------------------------------------------
#'@title GetDiffMZppm
#'@description GetDiffMZppm
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param mz mz
#'@param mz.ppm.thr mz.ppm.thr
#'@export

GetDiffMZppm <- function(mz, mz.ppm.thr = NULL) {
  mz.diff <- diff(mz) / mz[-1] * 1e6
  if (!is.null(mz.ppm.thr)) {
    idx <- which(mz[-1] <= mz.ppm.thr)
    mz.diff[idx] <- mz.diff[idx] * mz[-1] / mz.ppm.thr
  }
  mz.diff
}



#-----------------------------------------------------------------------------
#'@title GetWeightedInt
#'@description GetWeightedInt
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param spec spec
#'@param weight.mz weight.mz
#'@param weight.int weight.int
#'@export
#'


GetWeightedInt <- function(spec, weight.mz = 0, weight.int = 1) {
  return(spec[, 'mz'] ^ weight.mz * spec[, 'intensity'] ^ weight.int)
}


#-----------------------------------------------------------------------------
#'@title GetDotProduct
#'@description GetDotProduct
#'@author Hao Li, Yandong Yin, Kai Weng
#'\email{shenxt@@sioc.ac.cn}
#'@param x x
#'@param y y
#'@export

GetDotProduct <- function(x, y) {
  # return(sum(x * y) ^ 2 / (sum(x ^ 2) * sum(y ^ 2)))
  return(sum(x * y) / sqrt((sum(x ^ 2) * sum(y ^ 2))))
}


#-----------------------------------------------------------------------------

TuneMS2 <- function(spec,
                    mz.precursor,
                    is.include.precursor = TRUE,
                    snthr = 3,
                    noise.ms2 = 3,
                    mz.range.ms2 = NULL,
                    int.ms2.min.abs,
                    int.ms2.min.relative = 0.01,
                    is.apply.ms2.min.relative = TRUE,
                    is.check.sanity = TRUE,
                    int.check.sanity = 50,
                    ppm.precursor.filter = 20, ...) {
  # re-orgnize spec with increasing mz
  spec <- spec[order(spec[, 'mz']), , drop = FALSE]
  # define the lowest non-noise signal
  if (missing(int.ms2.min.abs)) {
    int.ms2.min.abs <- noise.ms2 * snthr
  }
  #
  # provent noise over estimation
  if (is.check.sanity & int.ms2.min.abs > int.check.sanity) {
    stop('Noise estimation is too high!')
  }

  # considering precursor ion
  if (missing(mz.precursor)) {
    mz.precursor <- max(spec[, 'mz'])
  }
  mz.precursor.range <- GetRangePPM(mz.precursor, ppm.precursor.filter)
  idx.mz.precursor.range <- ifelse(is.include.precursor, 2, 1)
  mz.cutoff <- mz.precursor.range[idx.mz.precursor.range]
  spec <- spec[spec[,'mz'] < mz.cutoff, , drop = FALSE]

  if (nrow(spec) == 0) {
    return()
  }

  if (!is.null(mz.range.ms2)) {
    nr.keep <- which(spec[, 'mz'] >= mz.range.ms2[1] &
                       spec[, 'mz'] <= mz.range.ms2[2])
    if (length(nr.keep) > 0) {
      spec <- spec[nr.keep, , drop = FALSE]
    }
    else {
      return()
    }
  }

  if ((max(spec[, 'intensity']) * int.ms2.min.relative) < int.ms2.min.abs) {
    is.warning.lowspec <- TRUE
  }

  # discarding low intensity spec (1% highest int and int.ms2.min.abs)
  int.cutoff <- max(max(spec[, 'intensity']) * int.ms2.min.relative,
                    int.ms2.min.abs)
  spec <- spec[spec[, 'intensity'] >= int.cutoff, , drop = FALSE]
  if (nrow(spec) == 0) {
    return()
  }

  # discarding ring effects
  spec <- RemoveRingEffect(spec)
}


#-----------------------------------------------------------------------------

RemoveRingEffect <- function(spec, mz.diff.thr = 0.3, int.rel.thr = 0.2) {
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


#-----------------------------------------------------------------------------

NormalizeSpec <- function(spec, ref = max(abs(spec[, 2])), pos = 'top') {
  spec[, 2] <- spec[, 2] / ref * ifelse(pos == 'down', -1, 1)
  return(spec)
}



#-----------------------------------------------------------------------------

GetRangePPM <- function(data, ppm) {
  t(sapply(data, function(dt) dt * (1 + c(-1, 1) * ppm * 1e-6)))
}
#-----------------------------------------------------------------------------

# PlotIDResults <- function(id.peaks.list, info.peak.plot,
#                           ms.assigned, dt.peaktable,
#                           lib.spec,
#                           lib.meta,
#                           ce,
#                           polarity = c('positive', 'negative'),
#                           direction = c('reverse', 'forward'),
#                           d.spec = 'ms2Result/MSMSfigures',
#                           width = 20, height = 7,
#                           is.include.precursor = TRUE,
#                           is.tune.ms2.exp = TRUE,
#                           is.tune.ms2.lib = FALSE,
#                           ...) {
#   polarity <- match.arg(polarity)
#   direction <- match.arg(direction)
#   col.plot <- c( 'lib' = 'red', 'exp' = 'blue', 'filtered' = 'gray')
#
#   for (i in rev(seq(nrow(info.peak.plot)))) {
#     cat(i); cat(" ")
#     apply(id.peaks.list[[info.peak.plot$idx.peaklist[i]]], 1, function(r.id) {
#       id <- as.character(r.id['labid'])
#       score <- round(as.numeric(r.id['score']), 3)
#       spec.lib <- spec.lib.all <- lib.spec[[id]][[ce]]
#       if (!is.include.precursor) {
#         mz.precursor <- as.numeric(lib.meta[lib.meta[, 'labid'] == id, 'mz'])
#         nr.remove <- which(spec.lib[, 'mz'] == mz.precursor)
#         spec.lib <- spec.lib[-nr.remove, , drop = FALSE]
#       }
#       if (is.tune.ms2.lib) {
#         spec.lib <- TuneMS2(spec.lib, mz.precursor,
#                             is.include.precursor = is.include.precursor,
#                             ...)
#       }
#
#       spec.lib.filtered <- spec.lib.all[
#         which(!spec.lib.all[, 'mz'] %in% spec.lib[, 'mz']), ,
#         drop = FALSE]
#
#       nr.peaktable <- info.peak.plot$nr.peaktable[i]
#       spec.exp <- spec.exp.all <- ms.assigned[[nr.peaktable]]
#       if (is.tune.ms2.exp) {
#         spec.exp <- TuneMS2(spec.exp, dt.peaktable[nr.peaktable, 'mzmed'],
#                             is.include.precursor = is.include.precursor,
#                             ...)
#       }
#       spec.exp.filtered <- spec.exp.all[
#         which(!spec.exp.all[, 'mz'] %in% spec.exp[, 'mz']), ,
#         drop = FALSE]
#
#       spec2match <- GetSpec2Match(spec.exp, spec.lib, direction = direction)
#
#       nr.matched <- which(spec2match$exp[, 'intensity'] > 0 &
#                             spec2match$lib[, 'intensity'] > 0)
#       spec.matched <- lapply(spec2match, function(spec) {
#         spec[nr.matched, , drop = FALSE]
#       })
#
#       d.plot <- file.path(d.spec, paste(dt.peaktable[nr.peaktable, 'name'], direction, sep = '_'))
#       dir.create(d.plot)
#       cmpd.replaced <- paste(
#         sapply(strsplit(r.id['name'], '')[[1]], function(x) {
#           # switch(x, ':' = '：', '/' = '／', x)
#           switch(x, ':' = '_', '/' = '-', x)
#         }),
#         collapse = '')
#
#       fn.plot <- switch(as.character(length(r.id)),
#                         '3' = file.path(d.plot,
#                                         paste(score, ',',
#                                               cmpd.replaced, '.pdf',
#                                               sep = ''))
#                         ,
#                         '4' = file.path(d.plot,
#                                         paste(score, ',',
#                                               cmpd.replaced, ',',
#                                               r.id['adduct'], '.pdf',
#                                               sep = ''))
#       )
#       range.mz <- range(c(spec.lib.all[, 'mz'], spec.exp.all[, 'mz']))
#       range.int <- c(-1, 1)
#
#       pdf(file = fn.plot, height = 7, width = 20, family = 'mono')
#       plot(range.mz, range.int, type = 'n', main = r.id['name'],
#            xlab = 'm/z', ylab = 'Relative intensity')
#       abline(h = 0, col = 'black')
#
#       ref.lib <- max(spec.lib.all[, 'intensity'])
#       points(NormalizeSpec(spec.lib, ref.lib, 'down'),
#              type = 'h', col = col.plot['lib'])
#       if (nrow(spec.lib.filtered) > 0) {
#         points(NormalizeSpec(spec = spec.lib.filtered,
#                              ref = ref.lib,
#                              pos = 'down'),
#                type = 'h', col = col.plot['filtered'])
#       }
#       ref.exp <- max(spec.exp.all[, 'intensity'])
#
#       points(NormalizeSpec(spec.exp, ref.exp, 'top'),
#              type = 'h', col = col.plot['exp'])
#       if (nrow(spec.exp.filtered) > 0) {
#         points(NormalizeSpec(spec = spec.exp.filtered,
#                              ref = ref.exp,
#                              pos = 'top'),
#                type = 'h', col = col.plot['filtered'])
#       }
#
#       points(NormalizeSpec(spec.matched$lib, ref.lib, 'down'),
#              type = 'p', pch = 20, col = col.plot['lib'])
#       points(NormalizeSpec(spec.matched$exp, ref.exp, 'top'),
#              type = 'p', pch = 20, col = col.plot['exp'])
#       if (ce == 'spec') {
#         legend('bottomleft',
#                legend = c(paste('Name:', r.id['name']),
#                           paste('Polarity:', polarity)),
#                pch = NA, bty = 'n')
#       } else {
#         legend('bottomleft',
#                legend = c(paste('LabID:', r.id['labid']),
#                           paste('Polarity:', polarity),
#                           paste('CE:', ce)),
#                pch = NA, bty = 'n')
#       }
#
#       legend('topleft',
#              legend = c(paste('Score:', score),
#                         paste('Matched peaks:', c('data', 'library'))),
#              col = c('white', col.plot[c('exp', 'lib')]),
#              pch = list(NA, 20, 20),
#              bty = 'n')
#       dev.off()
#     })
#   }
# }
