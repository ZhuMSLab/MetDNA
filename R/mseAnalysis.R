#' @title mseAnalysis
#' @description Metabolite sets enrichment analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param metabolite.id The vector of metabolite IDs.
#' @param species The species. "hsa" is Homo sapiens (human),
#' "dme" is Drosophlia melanogaster (fruit fly),
#' "mmu" is Mus musculus (mouse), "rat" is Rattus norvegicus (rat),
#' "bta" is Bos taurus (cow), "gga" is Gallus domesticus (chicken),
#' "dre" is Danio rerio (zebrafish), "cel" is Caenorharomyces elegans (nematode),
#'  "sce" is Saccharomyces cerevisaiae (yeast), "ath" is Arabidopsis thaliana (thale cress),
#'  "smm" is Schistosoma mansoni, "pfa" is Plasmodum falciparum 3D7,
#'  "tbr" is Trypanosoma brucei, "eco" is Escherichia coli K-12 MG1655,
#' "ppu" is Pseudomonas putida KT2440, "syf" is Synechococcus elongatus.
#' @param test.method Hypergeometric or fisher test.
#' @return The MSE analysis result.
#' @export


mseAnalysis <- function (metabolite.id,
                         species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                     "dre", "cel", "sce", "ath", "smm", "pfa",
                                     "tbr", "eco", "ppu", "syf"),
                         test.method = c("hypergeometric", "fisher")){
  #
  metabolite.id <- unique(metabolite.id)
  species <- match.arg(species)
  test.method <- match.arg(test.method)


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

  metabolite.id <- metabolite.id[which(metabolite.id %in% unique(unlist(M)))]
  if(length(metabolite.id) == 0) return(NULL)

  ALL <- unname(unique(unlist(M)))
  ALL <- as.character(as.matrix(ALL))
  SIG <- as.character(as.matrix(metabolite.id))
  num_all <- length(ALL)
  num_sig <- length(SIG)

  Lall0 <- unlist(lapply(M, function(x) {
    length(intersect(x, ALL))
  }))

  Lall <- Lall0

  Lsig <- unlist(lapply(M, function(x) {
    length(intersect(x, SIG))
  }))

  remove.idx <- which(unname(Lall0, length) == 0)

  if(length(remove.idx) > 0){
    Lall <- Lall0[-remove.idx]
    Lsig <- Lsig[-remove.idx]
  }

  # error handling
  if (length(Lall) < 2){
    # stop function
    return(NULL)
  }

  P<-NaN
  system.time(for (i in 1:length(Lall)){
    # ------------------------------------
    #Generating 2?~2 table
    # -------------------------------------
    a1 <- Lsig[i]# significant and including pathway
    a2 <- Lall[i]-Lsig[i]# not significant and including pathway
    a3 <- length(SIG)-a1# significant and not including pathway
    a4 <- (length(ALL)-length(SIG))-a2# not significant and not including pathway

    if(test.method == "hypergeometric"){
      P[i] <- phyper(q = a1 - 1, m = Lall[i], n = num_all - Lall[i],
                     k = num_sig, lower.tail = FALSE)
    }else{
      tab <- t(matrix(c(a1,a2,a3,a4),2))
      # ----------------------------------
      # Fisher's exact test
      # ----------------------------------
      check <- tryCatch({
        resfish <- fisher.test(tab, alternative="greater")
      }, error = function(e){
        NA
      })

      if(class(check) != "htest"){
        P[i] <- 1
      }else{
        resfish <- fisher.test(tab, alternative="greater")
        P[i] <- resfish$p.value
      }
    }


  })



  # -----------------------
  #q-value
  # -----------------------
  Q <- p.adjust(P, method="BH")
  FDR <- p.adjust(P, method="fdr")
  # --------------------------------------------------------
  #significant metabolites for metabolite set
  # --------------------------------------------------------
  # LES <- NaN
  # for (i in 1:ncol(Lsig)){
  #   les <- SIG[Lsig[,i]==1]
  #   LES[i] <- list(les)
  #
  # }
  # names(LES) <- colnames(Lsig)

  # ----------------------
  #Result
  # ----------------------
  PQ <- cbind(P,Q, FDR)
  rownames(PQ) <- colnames(Lsig)
  colnames(PQ) <- c("p.value","q.value", "FDR")

  ##calculate the impact of pathway
  info <- lapply(M, function(module) {
    overlap.number <- length(intersect(module, metabolite.id))
    pathway.number <- length(module)
    c(pathway.number, overlap.number)
  })

  info <- do.call(rbind, info)
  colnames(info) <- c("Pathway.length", "Overlap")

  info <- data.frame(PQ, info, stringsAsFactors = FALSE)

  info <- info[order(info[,1]),]
  #     RES <- list(PQ, LES)
  #     names(RES) <- c("Result of MSEA(ORA)","significant metabolites")

  # -------------------
  #Return
  # -------------------
  info <- info
}


setlabel <- function (M_ID, M){
  temp <- lapply(M_ID, function(x){
    unlist(lapply(M, function(y){
      sum(is.element(x, y))
    }))
  })
  temp <-do.call(rbind, temp)
  rownames(temp) <- M_ID
  temp <- temp
}





#-------------------------------------------------------------------------------
mseAnalysis1 <- function (metabolite.id,
                         M,
                         test.method = c("hypergeometric", "fisher")){
  #
  metabolite.id <- unique(metabolite.id)
  test.method <- match.arg(test.method)

  metabolite.id <- metabolite.id[which(metabolite.id %in% unique(unlist(M)))]
  if(length(metabolite.id) == 0) return(NULL)

  ALL <- unname(unique(unlist(M)))
  ALL <- as.character(as.matrix(ALL))
  SIG <- as.character(as.matrix(metabolite.id))
  num_all <- length(ALL)
  num_sig <- length(SIG)

  Lall0 <- unlist(lapply(M, function(x) {
    length(intersect(x, ALL))
  }))

  Lall <- Lall0

  Lsig <- unlist(lapply(M, function(x) {
    length(intersect(x, SIG))
  }))

  remove.idx <- which(unname(Lall0, length) == 0)

  if(length(remove.idx) > 0){
    Lall <- Lall0[-remove.idx]
    Lsig <- Lsig[-remove.idx]
  }

  # error handling
  if (length(Lall) < 2){
    # stop function
    return(NULL)
  }

  P<-NaN
  system.time(for (i in 1:length(Lall)){
    # ------------------------------------
    #Generating 2?~2 table
    # -------------------------------------
    a1 <- Lsig[i]# significant and including pathway
    a2 <- Lall[i]-Lsig[i]# not significant and including pathway
    a3 <- length(SIG)-a1# significant and not including pathway
    a4 <- (length(ALL)-length(SIG))-a2# not significant and not including pathway

    if(test.method == "hypergeometric"){
      P[i] <- phyper(q = a1 - 1, m = Lall[i], n = num_all - Lall[i],
                     k = num_sig, lower.tail = FALSE)
    }else{
      tab <- t(matrix(c(a1,a2,a3,a4),2))
      # ----------------------------------
      # Fisher's exact test
      # ----------------------------------
      check <- tryCatch({
        resfish <- fisher.test(tab, alternative="greater")
      }, error = function(e){
        NA
      })

      if(class(check) != "htest"){
        P[i] <- 1
      }else{
        resfish <- fisher.test(tab, alternative="greater")
        P[i] <- resfish$p.value
      }
    }


  })



  # -----------------------
  #q-value
  # -----------------------
  Q <- p.adjust(P, method="BH")
  FDR <- p.adjust(P, method="fdr")
  # --------------------------------------------------------
  #significant metabolites for metabolite set
  # --------------------------------------------------------
  # LES <- NaN
  # for (i in 1:ncol(Lsig)){
  #   les <- SIG[Lsig[,i]==1]
  #   LES[i] <- list(les)
  #
  # }
  # names(LES) <- colnames(Lsig)

  # ----------------------
  #Result
  # ----------------------
  PQ <- cbind(P,Q, FDR)
  rownames(PQ) <- colnames(Lsig)
  colnames(PQ) <- c("p.value","q.value", "FDR")

  ##calculate the impact of pathway
  info <- lapply(M, function(module) {
    overlap.number <- length(intersect(module, metabolite.id))
    pathway.number <- length(module)
    c(pathway.number, overlap.number)
  })

  info <- do.call(rbind, info)
  colnames(info) <- c("Pathway.length", "Overlap")

  info <- data.frame(PQ, info, stringsAsFactors = FALSE)

  info <- info[order(info[,1]),]
  #     RES <- list(PQ, LES)
  #     names(RES) <- c("Result of MSEA(ORA)","significant metabolites")

  # -------------------
  #Return
  # -------------------
  info <- info
}
