# ######code of pathifier
# ".getmeasuredgenesinpathway"<- function(syms,allgenes){
#   l<-length(syms)
#   pathways<-vector('list',l)
#   for (i in 1:l) {
#     n<-length(syms[[i]])
#     isin<-matrix(FALSE,n)
#     for(j in 1:n) {
#       isin[j]<-
#         (length(grep(paste('\\b',R.oo::trim(syms[[i]][j]),'\\b',sep=""),allgenes))>0)
#     }
#     pathways[[i]]<-unique(syms[[i]][isin])
#   }
#   pathways
# }
#
#
# getPathway <- function(sym,allgenes,data){
#   l <- length(sym)
#   x <- NULL
#   isin <- rep(FALSE,l)
#   #---------------------------------------------
#   for(i in 1:l) {
#     ind <- unique(grep(paste('\\b',R.oo::trim(sym[i]),'\\b',sep=""),allgenes))
#     n<-length(ind)
#     if (n>0) {
#       if (n==1)
#         t<-data[ind,]
#       else
#         t<-colMeans(data[ind,])
#       if (var(t)>0) {
#         isin[i]=TRUE;
#         x<-c(x,t)
#       }
#     }
#   }
#
#   #-----------------------------------------------
#   if (is.null(x)) {
#     list(x=NULL,isin=isin)
#   } else {
#     list(x=matrix(x,nrow=ncol(data)),isin=isin)
#   }
# }
#
#
#
# getPathway1 <- function(sym,allgenes,data,threads = 4){
#   l <- length(sym)
#   x <- NULL
#   isin <- rep(FALSE,l)
#   temp.fun <- function(i, sym = sym, allgenes = allgenes, data = data) {
#     ind <- unique(grep(paste('\\b',R.oo::trim(sym[i]),'\\b',sep=""),allgenes))
#     n<-length(ind)
#     if (n>0) {
#       if (n==1)
#         t <- data[ind,]
#       else
#         t <- colMeans(data[ind,])
#       if (var(t)>0) {
#         isin <- TRUE
#         t <- t
#       }else{
#         isin <- FALSE
#         t <- NA
#       }
#       result <- c(list(isin), list(t))
#     }
#   }
#
#   result <- BiocParallel::bplapply(c(1:l), temp.fun,
#                                    BPPARAM =BiocParallel::SnowParam(workers = threads),
#                                    sym = sym, allgenes = allgenes, data = data)
#   result <- result[-which(unlist(lapply(result, is.null)))]
#   isin <- unlist(lapply(result, function(x) x[[1]]))
#   x <- unlist(lapply(result, function(x) x[[2]]))
#   x <- x[!is.na(x)]
#   if (is.null(x)) {
#     list(x=NULL,isin=isin)
#   } else {
#     list(x=matrix(x,nrow=ncol(data)),isin=isin)
#   }
# }
#
#
#
# # system.time(a <- getPathway(sym = pathway.data[[1]][[3]], allgenes = allgenes, data = data))
# # system.time(b <- getPathway1(sym = pathway.data[[1]][[3]], allgenes = allgenes, data = data, threads = 4))
#
#
# setGeneric(name = "scorePathway",
#            def = function(x,m,ranks,
#                           calcerr=FALSE,
#                           thresh = 0.0005,
#                           maxit=200,start,logfile = ""){
#              x<-x[,apply(x,2,sd)>0.001]
#              k<-dim(x)[2]
#              if (k<3) {
#                c <- NULL
#                cat(file=logfile,append=TRUE,'scoring failed (k=',k,').\n')
#              } else {
#                d <- matrix(0,1,m)
#                if (start == "by pca") {
#                  start <- NULL
#                } else if (start == "by ranks") {
#                  start <- aggregate(x, by=list(ranks), FUN=mean)
#                  start <- as.matrix(start[,-1])
#                }
#                c <- princurve::principal.curve(x, ranks, start=start,
#                                                thresh=thresh, maxit=maxit, plot.true=FALSE)
#              }
#              if (!is.null(c)) {
#                d[c$tag[1]]=0
#                for (j in 2:m) {
#                  d[c$tag[j]]<-d[c$tag[j-1]]+dist(c$s[c$tag[(j-1):j],])
#                }
#                d=d/d[c$tag[m]]
#                if (calcerr) {
#                  e<-matrix(0,1,k)
#                  for (i in 1:k) {
#                    e[i]<-mean((c$s[,i]-x[,i])^2)
#                  }
#                } else {
#                  e <- FALSE;
#                }
#                list(score=d,error=e,thecurve=c)
#              } else {
#                cat(file=logfile,append=TRUE,'scoring failed.\n')
#                NULL
#              }
#            })
#
#
#
# setGeneric(name = "samplingsStdev",
#            def = function(m,n,attempts,z,ranks,
#                           samplings,start,logfile = ""){
#              dall<-array(dim=c(attempts,n))
#              skip<-0
#              for(a in 1:attempts) {
#                res<-scorePathway(z[samplings[a,],],m,ranks[samplings[a,]],
#                                  start=start,logfile=logfile)
#                if (!is.null(res)) {
#                  dall[a,samplings[a,]] <- res$score
#                } else {
#                  skip <- skip+1
#                }
#              }
#              if (skip < attempts/2) {
#                mean(apply(dall,2,sd,'na.rm'=TRUE), 'na.rm'=TRUE)
#              } else {
#                Inf
#              }
#            })
#
#
#
#
#
# setGeneric(name = "score_all_pathways_helper",
#            def = function(z, ranks, samplings,
#                           i, attempts, maximize_stability,
#                           logfile = "",start){
#              n<-dim(z)[1]
#              k<-dim(z)[2]
#              m<-dim(samplings)[2]
#              mincheck<-5
#              kmin=max(floor(0.8*k),mincheck+1)
#              mindelta=min(0.009,max(0.002,1.5/k))
#              sig<-matrix(0,1,k)
#              res<-scorePathway(z,n,ranks,calcerr=TRUE,start=start,logfile=logfile)
#              if (is.null(res)) {
#                cat(file=logfile,append=TRUE,'pathway ', i, '> scoring failed 1.\n')
#              } else {
#                sig<-samplingsStdev(m,n,attempts,z,ranks,samplings,start=start)
#                if (sig>10000) {
#                  cat(file=logfile,append=TRUE,
#                      'pathway ', i, '> scoring failed 2 (sig:', sig, ').\n')
#                  res<-NULL
#                } else {
#                  origsig<-sig
#                  cat(file=logfile,append=TRUE,'pathway ', i, '> sig:', sig, '\n')
#                  isin<-1:k
#                  if (maximize_stability) {
#                    testsig<-max(mincheck,floor(0.1*k))
#                    newsig<-rep(0,testsig)
#                    while ((k>=kmin)&(sig>0.05)) {
#                      se<-sort(res$error,index.return=TRUE,decreasing=TRUE)
#                      for (j in 1:testsig) {
#                        newsig[j] <- samplingsStdev(m,n,attempts,z[,-se$ix[j]],
#                                                    ranks,samplings,start=start)
#                      }
#                      wj<-which.min(newsig)
#                      cat(file=logfile,append=TRUE,
#                          'pathway ', i, ' k=', k, '(', ncol(res$thecurve$s),
#                          ') wj=', wj, '>new sig:', newsig[wj])
#                      if (sig-newsig[wj]<mindelta) {
#                        cat(file=logfile,append=TRUE,' x rejected\n')
#                        break
#                      }
#                      cat(file=logfile,append=TRUE,' | accepted!\n')
#                      sig<-newsig[wj]
#                      isin<-isin[-se$ix[wj]];
#                      z<-z[,-se$ix[wj]]
#                      k<-k-1
#                      res<-
#                        scorePathway(z,n,ranks,calcerr=TRUE,
#                                     start=start,logfile=logfile)
#                      if (is.null(res)) {
#                        cat(file=logfile,append=TRUE,
#                            'pathway ', i, '> scoring failed 3.\n')
#                        break;
#                      }
#                    }
#                  }
#                }
#              }
#              if (is.null(res)) {
#                NULL
#              } else {
#                list(score=res$score,thecurve=res$thecurve,z=z,
#                     isin=isin,sig=sig,origsig=origsig,k=k)
#              }
#            })
#
#
#
# setGeneric(name = "single2module",
#            def = function(data, allgenes, syms,
#                           pathwaynames, normals = NULL,
#                           ranks = NULL, attempts = 100,
#                           maximize_stability = TRUE,
#                           logfile = "", samplings=NULL,
#                           min_exp=4,
#                           min_std=0.4){
#              #
#              # cat(file=logfile,append=FALSE,
#              #     'robust_score_bydist. min_exp=',
#              #     min_exp,', min_std=',min_std,'\n')
#              data[data < min_exp] = min_exp;
#              n<-ncol(data)
#              if (is.null(normals)) {
#                normals <- rep(TRUE,n);
#                start <- "by pca";
#              } else {
#                start <- "by ranks"
#              }
#              if (is.null(ranks)) ranks <- !normals;
#              ranks <- rank(ranks)
#              if ((length(normals)!=n)||(length(ranks)!=n)) {
#                stop("invalid dimentions");S
#              }
#              l <- length(syms)
#              nn <- sum(normals)
#              m <- floor(0.8*(n-nn))+nn
#              if (is.null(samplings)) {
#                samplings<-matrix(0,attempts,m)
#                w<-which(!normals)
#                for(a in 1:attempts) {
#                  samplings[a,]<-sort(c(w[sample(n-nn,m-nn)],which(normals)))
#                }
#              }
#              s <- NULL
#              ind <- NULL
#              # s <- foreach::foreach(i = 1:l,
#              #                       .options.multicore = list(preschedule = FALSE),
#              #                       .combine = rbind,
#              #                       .inorder = FALSE,
#              #                       .verbose = FALSE,
#              #                       .errorhandling = 'stop',
#              #                       .multicombine = TRUE) %dopar% {
#              for (i in 1:l) {
#                # if(i == 7)
#                # cat(i); cat(" ")
#                # if using pathifier for large number of pathways,
#                # you might want to use the doMC library to parallelize your code,
#                # in that case replace the above for with:
#                # s <- foreach (i=1:l, .options.multicore=list(preschedule=FALSE),
#                # .combine=rbind, .inorder=FALSE, .verbose=FALSE, .errorhandling='stop',
#                # .multicombine=TRUE) %dopar% {
#                pathway <- syms[[i]]
#                pathwayindata <-
#                  getPathway(pathway, allgenes, data)
#                k1 = sum(pathwayindata$isin)
#                if (k1 < 3) {
#                  si <- NULL
#                  cat(file = logfile,
#                      append = TRUE,
#                      'skipping pathway ',
#                      i,
#                      ' k1=',
#                      k1,
#                      '\n')
#                } else {
#                  x <- pathwayindata$x
#                  pathway <-
#                    pathway[pathwayindata$isin]
#                  xm <- colMeans(x[normals, ])
#                  xs <- apply(x[normals, ], 2, sd)
#                  xs[xs < min_std] = min_std
#
#                  if (0 %in% xs) {
#                    si <- NULL
#                    cat(file = logfile,
#                        append = TRUE,
#                        'skipping pathway ',
#                        i,
#                        ' (0 in xs)\n')
#                  } else {
#                    z <- (x - matrix(rep(xm, each = n), nrow = n)) / (matrix(rep(xs, each =
#                                                                                   n), nrow = n))
#                    t <- prcomp(z)
#                    k2 = max(sum(t$sdev > 1.1), 4)
#                    k2 = min(k2, k1, 0.75 * dim(x)[1], sum(t$sdev >
#                                                             0.25))
#                    if (k2 < 3) {
#                      si <- NULL
#                      cat(file = logfile,
#                          append = TRUE,
#                          'skipping pathway ',
#                          i,
#                          ' k2=',
#                          k2,
#                          '\n')
#                    } else {
#                      pca <- t$x[, 1:k2]
#
#                      check <-
#                        tryCatch({res <-
#                          score_all_pathways_helper(pca,
#                                                    ranks,
#                                                    samplings,
#                                                    i,
#                                                    attempts,
#                                                    maximize_stability,
#                                                    logfile,
#                                                    start = start)},
#                          error = function(e){
#                            NULL
#                          }
#                        )
#
#                      if(is.null(check)){
#                        res <- NULL
#                      }else{
#                        res <-
#                          score_all_pathways_helper(pca,
#                                                    ranks,
#                                                    samplings,
#                                                    i,
#                                                    attempts,
#                                                    maximize_stability,
#                                                    logfile,
#                                                    start = start)
#                      }
#
#
#                      if (is.null(res)) {
#                        si <- NULL
#                        cat(file = logfile,
#                            append = TRUE,
#                            'skipping pathway ',
#                            i,
#                            '\n')
#                      } else {
#                        ind <- c(ind, i)
#                        si <-
#                          list(
#                            res$score,
#                            pathway,
#                            res$sig,
#                            res$origsig,
#                            res$k,
#                            res$thecurve$s,
#                            res$thecurve$tag,
#                            res$z,
#                            res$isin,
#                            xm,
#                            xs,
#                            t$center,
#                            t$rotation,
#                            k2
#                          )
#                      }
#                    }
#                  }
#                }
#                s <- rbind(s,si) # if using doMC foreach above replace this line with simply: si
#                # si
#              }
#              cat(file=logfile,append=TRUE,
#                  length(ind),'pathways processed with start=',
#                  start,'\n')
#              rownames(s)<-pathwaynames[ind]
#              list(scores = s[,1], genesinpathway=s[,2], newmeanstd=s[,3],
#                   origmeanstd=s[,4], pathwaysize=s[,5], curves=s[,6],
#                   curves_order=s[,7], z=s[,8],compin=s[,9],xm=s[,10],
#                   xs=s[,11],center=s[,12],rot=s[,13],pctaken=s[,14],
#                   samplings=samplings,sucess=ind,logfile=logfile)
#
#            })
#
#
#
# # syms <- pathway.data[[1]][1:20]
# # system.time(a <- single2module(data = data, allgenes = allgenes,
# #                                syms = syms,pathwaynames = pathway.data$pathwaynames[1:20],
# #                                normals = normals))
# #
# # system.time(b <- single2module1(data = data, allgenes = allgenes,
# #                                syms = syms,pathwaynames = pathway.data$pathwaynames[1:20],
# #                                normals = normals))
#
#
#
#
# setGeneric(name = "single2module1",
#            def = function(data, allgenes, syms,
#                           pathwaynames, normals = NULL,
#                           ranks = NULL, attempts = 100,
#                           maximize_stability = TRUE,
#                           logfile = "", samplings=NULL,
#                           min_exp=4,
#                           min_std=0.4,
#                           threads = parallel::detectCores()-3){
#              #
#              # cat(file=logfile,append=FALSE,
#              #     'robust_score_bydist. min_exp=',
#              #     min_exp,', min_std=',min_std,'\n')
#              data[data < min_exp] = min_exp;
#              n<-ncol(data)
#              if (is.null(normals)) {
#                normals <- rep(TRUE,n);
#                start <- "by pca";
#              } else {
#                start <- "by ranks"
#              }
#              if (is.null(ranks)) ranks <- !normals;
#              ranks <- rank(ranks)
#              if ((length(normals)!=n)||(length(ranks)!=n)) {
#                stop("invalid dimentions");
#              }
#              l <- length(syms)
#              nn <- sum(normals)
#              m <- floor(0.8*(n-nn))+nn
#              if (is.null(samplings)) {
#                samplings<-matrix(0,attempts,m)
#                w<-which(!normals)
#                for(a in 1:attempts) {
#                  samplings[a,]<-sort(c(w[sample(n-nn,m-nn)],which(normals)))
#                }
#              }
#
#              # s <- NULL
#              # ind <- NULL
#
#
#              temp.fun <- function(index, syms, allgenes, attempts,
#                                   normals, data, start, ranks, maximize_stability,
#                                   min_std, n, samplings){
#                temp.result <- lapply(index, function(i){
#                  pathway <- syms[[i]]
#
#                  pathwayindata <- getPathway(pathway, allgenes, data)
#                  k1 = sum(pathwayindata$isin)
#                  if (k1 < 3) return(NULL)
#
#                  x <- pathwayindata$x
#                  pathway <- pathway[pathwayindata$isin]
#                  xm <- colMeans(x[normals, ])
#                  xs <- apply(x[normals, ], 2, sd)
#                  xs[xs < min_std] = min_std
#                  if (0 %in% xs) return(NULL)
#
#                  z <-
#                    (x - matrix(rep(xm, each = n), nrow = n)) / (matrix(rep(xs, each =
#                                                                              n), nrow = n))
#                  t <- prcomp(z)
#                  k2 = max(sum(t$sdev > 1.1), 4)
#                  k2 = min(k2, k1, 0.75 * dim(x)[1], sum(t$sdev > 0.25))
#                  if (k2 < 3) return(NULL)
#
#                  pca <- t$x[, 1:k2]
#
#                  check <-
#                    tryCatch({res <-
#                      score_all_pathways_helper(pca,
#                                                ranks,
#                                                samplings,
#                                                i,
#                                                attempts,
#                                                maximize_stability,
#                                                logfile,
#                                                start = start)},
#                      error = function(e){
#                        NULL
#                      }
#                    )
#
#                  if(is.null(check)){
#                    res <- NULL
#                  }else{
#                    res <-
#                      score_all_pathways_helper(pca,
#                                                ranks,
#                                                samplings,
#                                                i,
#                                                attempts,
#                                                maximize_stability,
#                                                logfile,
#                                                start = start)
#                  }
#                  # x <- abc
#                  if (is.null(res)) return(NULL)
#
#                  ind <- i
#                  si <- list(res$score, pathway, res$sig, res$origsig,
#                             res$k, res$thecurve$s, res$thecurve$tag, res$z,
#                             res$isin, xm, xs, t$center, t$rotation, k2)
#                  c(list(ind), list(si))
#                })
#                temp.result
#              }
#
#              cl <- snow::makeCluster(spec = threads)
#
#              snow::clusterExport(cl = cl, list = c("getPathway",
#                                                    "score_all_pathways_helper",
#                                                    ".getmeasuredgenesinpathway",
#                                                    "samplingsStdev",
#                                                    "scorePathway"))
#              ichunks <- split(c(1:l), 1:threads)
#
#              result <- snow::clusterApply(cl = cl, ichunks, temp.fun, syms = syms,
#                                           allgenes = allgenes, attempts = attempts,
#                                           normals = normals, data = data, start = start,
#                                           ranks = ranks, maximize_stability = maximize_stability,
#                                           min_std = min_std, n = n,
#                                           samplings = samplings)
#              snow::stopCluster(cl)
#              result1 <- result[[1]]
#              if(length(result) > 1){
#                for(i in 2:length(result)) result1 <- c(result1, result[[i]])
#              }
#              result1 <- result1[order(unlist(ichunks))]
#
#              ind <- lapply(result1, function(x) x[[1]])
#              s <- lapply(result1, function(x) x[[2]])
#              remove.idx <- which(unlist(lapply(ind, is.null)))
#              ind <- ind[-remove.idx]
#              s <- s[-remove.idx]
#
#              ind <- unlist(ind)
#              s <- do.call(rbind, s)
#              rownames(s)<-pathwaynames[ind]
#              list(scores = s[,1], genesinpathway=s[,2], newmeanstd=s[,3],
#                   origmeanstd=s[,4], pathwaysize=s[,5], curves=s[,6],
#                   curves_order=s[,7], z=s[,8],compin=s[,9],xm=s[,10],
#                   xs=s[,11],center=s[,12],rot=s[,13],pctaken=s[,14],
#                   samplings=samplings,sucess=ind,logfile=logfile)
#
#            })
