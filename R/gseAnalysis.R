#' @title gseAnalysis
#' @description Gene sets enrichment analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param gene.id The vector of gene IDs.
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
#' @return  The MSE analysis result.
#' @export


gseAnalysis <- function (gene.id,
                         species = c("hsa","dme", "mmu", "rat", "bta", "gga",
                                     "dre", "cel", "sce", "ath", "smm", "pfa",
                                     "tbr", "eco", "ppu", "syf"),
                         test.method = c("hypergeometric", "fisher")){
  gene.id <- unique(gene.id)
  species <- match.arg(species)
  test.method <- match.arg(test.method)


  # data(paste(species, "gene.kegg.pathway", sep = "."), envir = environment())

  switch(species,
         "hsa" = {data("hsa.gene.kegg.pathway", envir = environment())
           M = hsa.gene.kegg.pathway},
         "dme" = {data("dme.gene.kegg.pathway", envir = environment())
           M = dme.gene.kegg.pathway},
         "mmu" = {data("mmu.gene.kegg.pathway", envir = environment())
           M = mmu.gene.kegg.pathway},
         "rat" = {data("rat.gene.kegg.pathway", envir = environment())
           M = rat.gene.kegg.pathway},
         "bta" = {data("bta.gene.kegg.pathway", envir = environment())
           M = bta.gene.kegg.pathway},
         "gga" = {data("gga.gene.kegg.pathway", envir = environment())
           M = gga.gene.kegg.pathway},
         "dre" = {data("dre.gene.kegg.pathway", envir = environment())
           M = dre.gene.kegg.pathway},
         "cel" = {data("cel.gene.kegg.pathway", envir = environment())
           M = cel.gene.kegg.pathway},
         "sce" = {data("sce.gene.kegg.pathway", envir = environment())
           M = sce.gene.kegg.pathway},
         "ath" = {data("ath.gene.kegg.pathway", envir = environment())
           M = ath.gene.kegg.pathway},
         "smm" = {data("smm.gene.kegg.pathway", envir = environment())
           M = smm.gene.kegg.pathway},
         "pfa" = {data("pfa.gene.kegg.pathway", envir = environment())
           M = pfa.gene.kegg.pathway},
         "tbr" = {data("tbr.gene.kegg.pathway", envir = environment())
           M = tbr.gene.kegg.pathway},
         "eco" = {data("eco.gene.kegg.pathway", envir = environment())
           M = eco.gene.kegg.pathway},
         "ppu" = {data("ppu.gene.kegg.pathway", envir = environment())
           M = ppu.gene.kegg.pathway},
         "syf" = {data("syf.gene.kegg.pathway", envir = environment())
           M = syf.gene.kegg.pathway}
         )

  # if(species == "hsa") {data("hsa.gene.kegg.pathway", envir = environment()); M = hsa.gene.kegg.pathway}
  # if(species == "dme") {data("dme.gene.kegg.pathway", envir = environment()); M = dme.gene.kegg.pathway}
  # if(species == "mmu") {data("mmu.gene.kegg.pathway", envir = environment()); M = mmu.gene.kegg.pathway}
  # if(species == "rat") {data("rat.gene.kegg.pathway", envir = environment()); M = rat.gene.kegg.pathway}

  ##filter pathway (M) according to overlap
  # result <- lapply(M, function(x){
  #   map.number <- length(intersect(x, gene.id))
  #   pathway.number <- length(x)
  #   map.ratio <- map.number/pathway.number
  #   c(pathway.number, map.number, map.ratio)
  # })
  #
  #
  # result <- as.data.frame(do.call(rbind, result))
  # colnames(result) <- c("pathway.number", "map.number", "map.ratio")
  #
  # remove.idx <- unname(which(result[,1] == 1 | result[,2] == 0 | result[,3] < 0.02))
  # if(length(remove.idx) > 0) {
  #   M <- M[-remove.idx]
  #   result <- result[-remove.idx,,drop = FALSE]
  # }
  #
  # if(length(M) == 0) return(NULL)
  # if(all(result[,2] == 1)) return(NULL)

  ALL <- unname(unique(unlist(M)))
  ALL <- as.character(as.matrix(ALL))

  ##remove the ID which is not in the ALL
  gene.id <- gene.id[which(is.element(gene.id, ALL))]

  if(length(gene.id) == 0) return(NULL)

  SIG <- as.character(as.matrix(gene.id))
  num_all <- length(ALL)
  num_sig <- length(SIG)

  # --------------------------------------------------------------
  #Generating label matrix for detected metabolites
  # --------------------------------------------------------------

  Lall0 <- setlabel(ALL,M)

  # delete metabolite set
  Lall <- Lall0[,colSums(Lall0)!=0, drop = FALSE]

  # error handling
  if (ncol(Lall) < 2){
    # stop function
    return(NULL)
  }

  # --------------------------------------------------------------
  #Generating label matrix for significant metabolites
  # --------------------------------------------------------------

  Lsig <- setlabel(SIG,M)
  Lsig <- Lsig[,colSums(Lall0)!=0, drop = FALSE]

  l <- colSums(Lall0)!=0

  # -------------------------------
  #Calculating  ORA
  # -------------------------------

  P<-NaN;
  for (i in 1:sum(l)){

    # ------------------------------------
    #Generating 2?~2 table
    # -------------------------------------
    a1 <- sum(Lsig[,i])# significant and including pathway
    a2 <- sum(Lall[,i])-sum(Lsig[,i])# not significant and including pathway
    a3 <- length(SIG)-a1# significant and not including pathway
    a4 <- (length(ALL)-length(SIG))-a2# not significant and not including pathway

    tab <- t(matrix(c(a1,a2,a3,a4),2));

    if(test.method == "hypergeometric"){
      P[i] <- phyper(q = a1 - 1, m = sum(Lall[,i]), n = num_all - sum(Lall[,i]),
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

  }

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
    overlap.number <- length(intersect(module, gene.id))
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
