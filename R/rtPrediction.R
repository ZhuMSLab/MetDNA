# title rtPrediction
# description Predict RTs of kegg metabolites using MS2 matched metabolites.
# author Xiaotao Shen
# \email{shenxt@@sioc.ac.cn}
# param data data From readAnnotation.
# param prefer.adduct The reliable adducts in this LC system. Default is all.
# param threads How many threads do you want to use? Default is the number of
# your PC threads - 3.
# param inHouse.compound.md The molecular descriptors of inhouse compounds.
# param kegg.compound.md The molecular descriptors of kegg compounds.
# return The predicted in house RT and KEGG RT.
# export
setGeneric(name = "rtPrediction",
           def = function(data,
                          prefer.adduct = c("all", "M+Na", "M+H", "M+NH4", "M-H"),
                          threads = 3,
                          inHouse.compound.md,
                          kegg.compound.md,
                          use.default.md = TRUE,
                          column = c("hilic", "rp")
           ){
             #

             prefer.adduct <- match.arg(prefer.adduct)
             column <- match.arg(column)
             tags <- data[[1]]
             sample <- data[[2]]
             sample.int <- apply(sample, 1, median)
             rm(list = "sample")
             gc()

             idx <- which(!is.na(tags$labid))
             tags1 <- tags[idx,]
             sample.int1 <- sample.int[idx]
             rm(list = "sample.int")
             gc()

             labid <- tags1$labid
             labid <- lapply(labid, function(x) {
               strsplit(x, split = ";")[[1]][1]
             })
             labid <- unlist(labid)


             adduct <- tags1$adduct
             adduct <- lapply(adduct, function(x) {
               strsplit(x, split = ";")[[1]][1]
             })
             adduct <- unlist(adduct)

             ##remove multipe peaks matched one metabolite
             dup.id <- unique(labid[duplicated(labid)])
             if(length(dup.id) > 0){
               for(i in 1:length(dup.id)){
                 temp.id <- dup.id[i]
                 temp.idx <- grep(temp.id, labid)
                 temp.int <- sample.int1[temp.idx]
                 # rm(list = c("sample.int1"))
                 # temp.adduct <- adduct[temp.idx]
                 # temp.rt <- tags1$rt[temp.idx]
                 labid[temp.idx[-which.max(temp.int)]]  <- NA
               }
             }

             rt <- tags1$rt
             rm(list = "tags1")
             gc()

             data <- data.frame(labid, rt, adduct, stringsAsFactors = FALSE)
             data <- data[!is.na(data$labid),,drop = FALSE]

             ##filter the metabolite who have no MD in inHouse.compound.md
             temp.idx <- which(data$labid %in% rownames(inHouse.compound.md))

             if(length(temp.idx) == 0) stop("The metabolites from MS2 match have no MD in inHouse.compound.md.\n")
             if(length(temp.idx) < 70){
               prefer.adduct <- 'all'
             }
             data <- data[temp.idx,]


             ##filter data using adduct or not
             adduct.order <- setdiff(names(sort(table(adduct), decreasing = TRUE)), prefer.adduct)
             if (prefer.adduct != 'all') {
               temp.idx <- which(data$adduct %in% prefer.adduct)
               count <- 1
               while(length(temp.idx) < 50 & count < length(adduct.order)){
                 prefer.adduct <- c(prefer.adduct, adduct.order[count])
                 temp.idx <- which(data$adduct %in% prefer.adduct)
                 count <- count + 1
               }
               data <- data[temp.idx,,drop = FALSE]
             }

             if(nrow(data) == 0) stop("No metabolites are identified.")

             cat("There are ", nrow(data),
                 " metabolites are used for RT prediction.\n", sep = "")
             # data <- data[data$adduct == "M+H",]
             #----------------------------------------------------------------
             # x <- table(adduct[which(!is.na(labid))])
             # x <- sort(x, decreasing = T)
             # par(mar = c(5,5,4,2))
             # y <- barplot(x, borde= NA, xlab = "Adduct", ylab = "Peak Number",
             #         cex.lab = 1.8, cex.axis = 1.5, names.arg = "",
             #         col = c("salmon", rep("black", 1),"lightseagreen",
             #                 rep("black", 5),"orchid4", rep("black", 1),rep("black", 8)))
             # par(xpd = T)
             # text(x = y, y = x + 5, labels = names(x), srt = 90)
             # pie(x, border = "white",
             #     col = c("salmon", rep("black", 1),"lightseagreen",
             #             rep("black", 5),"orchid4", rep("black", 1),rep("black", 8)),
             #     radius = 1.3)
             # par(new = T)
             # pie(1, col = "white", radius = 0.9, border = "white", labels = "")
             #----------------------------------------------------------------

             idx <- match(data$labid, rownames(inHouse.compound.md))
             md <- inHouse.compound.md[idx,]
             # rm(list = "inHouse.compound.md")

             ###remove  NA which apper in more than 50% metabolites
             remove.idx1 <-
               which(apply(md, 2, function(x) {sum(is.na(x)/nrow(md))})>0.5)
             md1 <- md[,-remove.idx1]
             rm(list = c("md"))
             gc()
             ##impute NA
             md2 <- t(impute::impute.knn(data = t(md1))[[1]])
             rm(list = "md1")
             gc()

             #remove MD which are same in all metaboites
             remove.idx2 <- which(apply(md2, 2, sd) == 0)
             md3 <- md2[,-remove.idx2]
             rm(list = c("md2"))
             gc()

             ##construct RF model
             train.y <- data$rt
             train.x <- md3
             rm(list = c('data', 'md3'))
             gc()

             if(use.default.md){
               switch(column,
                      "hilic" = {
                        marker.name <- c('XLogP', "tpsaEfficiency", "WTPT.5", "khs.dsCH", "MLogP", "nAcid", "nBase", "BCUTp.1l")
                      },
                      "rp" = {
                        marker.name <- c('XLogP', "WTPT.4", "WTPT.5", "ALogp2", "BCUTp.1l")})
             }else{
               #construct RF 100 times, find the MDs which are always apperar in
               # top 5 as markers
               temp.fun <- function(idx, x, y){
                 suppressMessages(library(randomForest))
                 temp <- idx
                 rf.class <- randomForest::randomForest(x = x, y = y,
                                                        replace = TRUE, importance = TRUE,
                                                        proximity = TRUE)
                 imp <- randomForest::importance(rf.class)
                 rm(list = c("rf.class", "temp"))
                 gc()
                 imp
               }
               cat("\n")
               cat("Find the optimal moleculr descriptor\n")
               imp <- BiocParallel::bplapply(1:100,
                                             temp.fun,
                                             BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                               progressbar = FALSE),
                                             x = train.x,
                                             y = train.y)

               md.name <- list()
               for(i in 1:length(imp)){
                 md.name[[i]] <- names(sort(imp[[i]][,1], decreasing = TRUE)[1:5])
               }

               md.name <- unlist(md.name)
               md.name <- table(md.name)
               md.name <- sort(md.name)
               marker.name <- names(md.name[which(md.name >= 50)])
               rm(list = "imp")
               gc()
               #
               save(marker.name, file = "marker.name", compress = "xz")
               save(md.name, file = "md.name", compress = "xz")
             }




             idx <- match(marker.name, colnames(train.x))

             idx <- idx[!is.na(idx)]
             if(length(idx) == 0){
               stop("Your markers are not in MD data.\n")
             }
             train.x1 <- train.x[,idx]
             rm(list = c("train.x"))
             gc()

             para <- NULL
             ntree1 <- seq(300,1000,by = 200)
             mtry1 <- seq(1,length(marker.name),by = 1)
             for(i in 1:length(ntree1)){
               para <- rbind(para, cbind(ntree1[i],mtry1))
             }
             colnames(para) <- c("ntree", "mtry")
             mse <- NULL
             rf.reg <- list()
             cat("\n")
             cat("Find the optimal parameters\n")
             for(i in 1:nrow(para)){
               cat(i);cat(" ")
               temp.ntree <- para[i,1]
               temp.mtry <- para[i,2]
               rf.reg[[i]] <- randomForest::randomForest(x = train.x1, y = train.y,
                                                         ntree = temp.ntree,mtry = temp.mtry,
                                                         replace = TRUE, importance = TRUE, proximity = TRUE)
               mse[i] <- mean(rf.reg[[i]]$mse)
             }
             cat("\n")
             result <- data.frame(para, mse, stringsAsFactors = FALSE)
             temp.idx <- which.min(result$mse)
             ntree <- result$ntree[temp.idx]
             mtry <- result$mtry[temp.idx]
             ##
             rf.reg <- randomForest::randomForest(x = train.x1, y = train.y,
                                                  ntree = ntree,
                                                  mtry = mtry,
                                                  replace = TRUE,
                                                  importance = TRUE, proximity = TRUE)

             rm(list = c("train.x1"))
             gc()
             ##predict RT in inhouse database
             test.x <- inHouse.compound.md
             rm(list = "inHouse.compound.md")
             gc()
             test.x <- test.x[,match(marker.name, colnames(test.x))]
             test.x <- as.data.frame(test.x)

             inHouse.compound.rt <- rep(NA, nrow(test.x))
             names(inHouse.compound.rt) <- rownames(test.x)
             ##impute NA in test.x
             idx1 <-
               which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
             test.x1 <- test.x[idx1,]
             test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
             temp.rt <- predict(object = rf.reg, newdata = test.x1)
             names(temp.rt) <- rownames(test.x1)
             inHouse.compound.rt[idx1] <- temp.rt
             ##NA give the median RT of all peaks
             inHouse.compound.rt[is.na(inHouse.compound.rt)] <- median(tags$rt)
             attribute <- rep(NA, length(inHouse.compound.rt))
             attribute[idx1] <- "Predicted"
             attribute[is.na(attribute)] <- "Median"
             inHouse.rt <- data.frame(inHouse.compound.rt, attribute, stringsAsFactors = FALSE)
             colnames(inHouse.rt)[1] <- "RT"

             ##cross validation
             ##predict RT in KEGG database
             test.x <- kegg.compound.md
             rm(list = "kegg.compound.md")
             gc()
             test.x <- test.x[,match(marker.name, colnames(test.x))]
             test.x <- as.data.frame(test.x)

             kegg.compound.rt <- rep(NA, nrow(test.x))
             names(kegg.compound.rt) <- rownames(test.x)
             ##impute NA in text.x
             idx1 <-
               which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
             test.x1 <- test.x[idx1,]
             test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
             temp.rt <- predict(object = rf.reg, newdata = test.x1)
             names(temp.rt) <- rownames(test.x1)
             kegg.compound.rt[idx1] <- temp.rt
             ##NA give the median RT of all peaks
             kegg.compound.rt[is.na(kegg.compound.rt)] <- median(tags$rt)
             rm(list = "tags")
             gc()
             attribute <- rep(NA, length(kegg.compound.rt))
             attribute[idx1] <- "Predicted"
             attribute[is.na(attribute)] <- "Median"
             kegg.rt <- data.frame(kegg.compound.rt, attribute, stringsAsFactors = FALSE)
             colnames(kegg.rt)[1] <- "RT"

             rt <- list(inHouse.rt, kegg.rt)
             names(rt) <- c("inHouse.rt", "KEGG.rt")
             rt <- rt
           })




##backup of old version
# setGeneric(name = "rtPrediction",
#            def = function(data,
#                           prefer.adduct = c("all", "M+Na", "M+H", "M+NH4", "M-H"),
#                           threads = 3,
#                           inHouse.compound.md,
#                           kegg.compound.md
#                           ){
#              #
#              prefer.adduct <- match.arg(prefer.adduct)
#              tags <- data[[1]]
#              sample <- data[[2]]
#              sample.int <- apply(sample, 1, median)
#              rm(list = "sample")
#
#              idx <- which(!is.na(tags$labid))
#              tags1 <- tags[idx,]
#              sample.int1 <- sample.int[idx]
#              rm(list = "sample.int")
#
#              labid <- tags1$labid
#              labid <- lapply(labid, function(x) {
#                strsplit(x, split = ";")[[1]][1]
#              })
#              labid <- unlist(labid)
#
#
#              adduct <- tags1$adduct
#              adduct <- lapply(adduct, function(x) {
#                strsplit(x, split = ";")[[1]][1]
#              })
#              adduct <- unlist(adduct)
#
# ##remove multipe peaks matched one metabolite
#              dup.id <- unique(labid[duplicated(labid)])
#              for(i in 1:length(dup.id)){
#                temp.id <- dup.id[i]
#                temp.idx <- grep(temp.id, labid)
#                temp.int <- sample.int1[temp.idx]
#                # rm(list = c("sample.int1"))
#                # temp.adduct <- adduct[temp.idx]
#                # temp.rt <- tags1$rt[temp.idx]
#                labid[temp.idx[-which.max(temp.int)]]  <- NA
#              }
#
#            rt <- tags1$rt
#            rm(list = "tags1")
#
#            data <- data.frame(labid, rt, adduct, stringsAsFactors = FALSE)
#            data <- data[!is.na(data$labid),]
#
#            ##filter data using adduct or not
#            if (prefer.adduct != 'all') {
#              data <- data[data$adduct %in% prefer.adduct,]
#            }
#            if(nrow(data) == 0) stop("No metabolites are identified.")
#
#            cat("There are ", nrow(data),
#                " metabolites are used for RT prediction.\n", sep = "")
#            # data <- data[data$adduct == "M+H",]
#            #----------------------------------------------------------------
#            # x <- table(adduct[which(!is.na(labid))])
#            # x <- sort(x, decreasing = T)
#            # par(mar = c(5,5,4,2))
#            # y <- barplot(x, borde= NA, xlab = "Adduct", ylab = "Peak Number",
#            #         cex.lab = 1.8, cex.axis = 1.5, names.arg = "",
#            #         col = c("salmon", rep("black", 1),"lightseagreen",
#            #                 rep("black", 5),"orchid4", rep("black", 1),rep("black", 8)))
#            # par(xpd = T)
#            # text(x = y, y = x + 5, labels = names(x), srt = 90)
#            # pie(x, border = "white",
#            #     col = c("salmon", rep("black", 1),"lightseagreen",
#            #             rep("black", 5),"orchid4", rep("black", 1),rep("black", 8)),
#            #     radius = 1.3)
#            # par(new = T)
#            # pie(1, col = "white", radius = 0.9, border = "white", labels = "")
#            #----------------------------------------------------------------
#
#            idx <- match(data$labid, rownames(inHouse.compound.md))
#            md <- inHouse.compound.md[idx,]
#            # rm(list = "inHouse.compound.md")
#
#            ###remove  NA which apper in more than 50% metabolites
#            remove.idx1 <-
#              which(apply(md, 2, function(x) {sum(is.na(x)/nrow(md))})>0.5)
#            md1 <- md[,-remove.idx1]
#            rm(list = c("md"))
#            ##impute NA
#            md2 <- t(impute::impute.knn(data = t(md1))[[1]])
#            rm(list = "md1")
#
#            #remove MD which are same in all metaboites
#            remove.idx2 <- which(apply(md2, 2, sd) == 0)
#            md3 <- md2[,-remove.idx2]
#            rm(list = c("md2"))
#
#            ##construct RF model
#            train.y <- data$rt
#            train.x <- md3
#            rm(list = c('data', 'md3'))
#            #construct RF 100 times, find the MDs which are always apperar in
#            # top 5 as markers
#            temp.fun <- function(idx, x, y){
#              suppressMessages(library(randomForest))
#              temp <- idx
#              rf.class <- randomForest::randomForest(x = x, y = y,
#                                                     replace = TRUE, importance = TRUE,
#                                                     proximity = TRUE)
#              imp <- randomForest::importance(rf.class)
#              rm(list = c("rf.class", "temp"))
#              imp
#            }
#            cat("\n")
#            cat("Find the optimal moleculr descriptor\n")
# imp <- BiocParallel::bplapply(1:100,
#                        temp.fun,
#                        BPPARAM = BiocParallel::SnowParam(workers = threads,
#                                                          progressbar = TRUE),
#                        x = train.x,
#                        y = train.y)
#
#            md.name <- list()
#            for(i in 1:length(imp)){
#              md.name[[i]] <- names(sort(imp[[i]][,1], decreasing = TRUE)[1:5])
#            }
#
#            md.name <- unlist(md.name)
#            md.name <- table(md.name)
#            md.name <- sort(md.name)
#            marker.name <- names(md.name[which(md.name >= 50)])
#            rm(list = "imp")
#            #
#            save(marker.name, file = "marker.name")
#            save(md.name, file = "md.name")
#
#            idx <- match(marker.name, colnames(train.x))
#            train.x1 <- train.x[,idx]
#            rm(list = c("train.x"))
#
#            para <- NULL
#            ntree1 <- seq(300,1000,by = 200)
#            mtry1 <- seq(1,length(marker.name),by = 1)
#            for(i in 1:length(ntree1)){
#              para <- rbind(para, cbind(ntree1[i],mtry1))
#            }
#            colnames(para) <- c("ntree", "mtry")
#            mse <- NULL
#            rf.reg <- list()
#            cat("\n")
#            cat("Find the optimal parameters\n")
#            for(i in 1:nrow(para)){
#              cat(i);cat(" ")
#              temp.ntree <- para[i,1]
#              temp.mtry <- para[i,2]
#              rf.reg[[i]] <- randomForest::randomForest(x = train.x1, y = train.y,
#                           ntree = temp.ntree,mtry = temp.mtry,
#                           replace = TRUE, importance = TRUE, proximity = TRUE)
#              mse[i] <- mean(rf.reg[[i]]$mse)
#            }
#            cat("\n")
#            result <- data.frame(para, mse, stringsAsFactors = FALSE)
#            temp.idx <- which.min(result$mse)
#            ntree <- result$ntree[temp.idx]
#            mtry <- result$mtry[temp.idx]
#            ##
#            rf.reg <- randomForest::randomForest(x = train.x1, y = train.y,
#                                                 ntree = ntree,
#                                   mtry = mtry,
#                                   replace = TRUE,
#                                   importance = TRUE, proximity = TRUE)
#
#            rm(list = c("train.x1"))
#            ##predict RT in inhouse database
#            test.x <- inHouse.compound.md
#            rm(list = "inHouse.compound.md")
#            test.x <- test.x[,match(marker.name, colnames(test.x))]
#            test.x <- as.data.frame(test.x)
#
#            inHouse.compound.rt <- rep(NA, nrow(test.x))
#            names(inHouse.compound.rt) <- rownames(test.x)
#            ##impute NA in test.x
#            idx1 <-
#              which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
#            test.x1 <- test.x[idx1,]
#            test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
#            temp.rt <- predict(object = rf.reg, newdata = test.x1)
#            names(temp.rt) <- rownames(test.x1)
#            inHouse.compound.rt[idx1] <- temp.rt
#            ##NA give the median RT of all peaks
#            inHouse.compound.rt[is.na(inHouse.compound.rt)] <- median(tags$rt)
#            attribute <- rep(NA, length(inHouse.compound.rt))
#            attribute[idx1] <- "Predicted"
#            attribute[is.na(attribute)] <- "Median"
#            inHouse.rt <- data.frame(inHouse.compound.rt, attribute, stringsAsFactors = FALSE)
#            colnames(inHouse.rt)[1] <- "RT"
#
#            ##cross validation
#            ##predict RT in KEGG database
#            test.x <- kegg.compound.md
#            rm(list = "kegg.compound.md")
#            test.x <- test.x[,match(marker.name, colnames(test.x))]
#            test.x <- as.data.frame(test.x)
#
#            kegg.compound.rt <- rep(NA, nrow(test.x))
#            names(kegg.compound.rt) <- rownames(test.x)
#            ##impute NA in text.x
#            idx1 <-
#              which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
#            test.x1 <- test.x[idx1,]
#            test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
#            temp.rt <- predict(object = rf.reg, newdata = test.x1)
#            names(temp.rt) <- rownames(test.x1)
#            kegg.compound.rt[idx1] <- temp.rt
#            ##NA give the median RT of all peaks
#            kegg.compound.rt[is.na(kegg.compound.rt)] <- median(tags$rt)
#            rm(list = "tags")
#            attribute <- rep(NA, length(kegg.compound.rt))
#            attribute[idx1] <- "Predicted"
#            attribute[is.na(attribute)] <- "Median"
#            kegg.rt <- data.frame(kegg.compound.rt, attribute, stringsAsFactors = FALSE)
#            colnames(kegg.rt)[1] <- "RT"
#
#            rt <- list(inHouse.rt, kegg.rt)
#            names(rt) <- c("inHouse.rt", "KEGG.rt")
#            rt <- rt
# })


#
# setGeneric(name = "crossValidation",
#            def = function(model, x, y, cross = 7){
#              data <- cbind(y,x)
#              data <- as.data.frame(data)
#              train_control <- caret::trainControl(method="cv", number = cross)
#              grid <- expand.grid(.fL=c(0), .usekernel=c(FALSE))
#              model <- train(y~., data=data,
#                             trControl=train_control, method="rf")
#              return(model)
#            })

#
# load('rt.all')
#
# dim(rt.all)
# colnames(rt.all)
# par(mar = c(5,5,4,2))
# plot(result.rt1[[1]]$RT, result.rt2[[1]]$RT, xlab = "Predicted RT (All adducts)",
#      ylab = "Predicted RT (M+H)", cex.lab = 1.8, cex.axis = 1.5)
# par(xpd = F)
# abline(0,1, lty = 2, col = "tomato", lwd = 1.5)
#
# cor(result.rt1[[1]]$RT, result.rt2[[1]]$RT)
# cor.test(result.rt1[[1]]$RT, result.rt2[[1]]$RT)
# legend("topleft",
#        legend = c("Pearson correlation; 0.86", 'p-value: <0.001'), bty = "n",
#        cex = 1.5)
#
#
# inHouse.rt <- result.rt1[[1]]
# kegg.rt <- result.rt1[[2]]
#
#
# idx1 <- match(rt.all$`in house ID`, rownames(inHouse.rt))
# rt.a <- rt.all$`RT （min）`*60
# rt.b <- inHouse.rt$RT[idx1]
#
# temp1 <- data.frame(rt.a, rt.b, stringsAsFactors = FALSE)
# temp1 <- temp1[!is.na(rt.b),]
#
# relative.error1 <- abs(temp1$rt.a - temp1$rt.b)*100/temp1$rt.a
#
# temp1 <- data.frame(temp1, relative.error1, stringsAsFactors = FALSE)
# plot(temp1$rt.a, temp1$relative.error1, ylim = c(0,1500), pch = 19,
#      xlab = "Experiment RT (second)", ylab = "Relative error (%)",cex.lab = 1.8, cex.axis = 1.5)
# plot(temp1$rt.a, temp1$relative.error1, ylim = c(0,300), pch = 19,
#      xlab = "Experiment RT (second)", ylab = "Relative error (%)",cex.lab = 1.8, cex.axis = 1.5)
# abline(v=300, lty = 2, col = "tomato", lwd = 2)
#
#
# absoult.error1 <- (temp1$rt.b - temp1$rt.a)
# temp1 <- data.frame(temp1, absoult.error1, stringsAsFactors = FALSE)
# plot(temp1$rt.a, abs(temp1$absoult.error1),
#      # ylim = c(-600,600),
#      pch = 19,
#      xlab = "Experiment RT (second)", ylab = "Absoult error (%)",cex.lab = 1.8, cex.axis = 1.5)
# abline(v=300, lty = 2, col = "tomato", lwd = 2)
#
# median.error1 <- NULL
# sd.error1 <- NULL
# for(i in 1:nrow(temp1)){
#   median.error1[i] <- median(abs(temp1$absoult.error1[which(temp1$rt.a>temp1$rt.a[i]-30 & temp1$rt.a<=temp1$rt.a[i]+30)]))
#   sd.error1[i] <- sd(abs(temp1$absoult.error1[which(temp1$rt.a>temp1$rt.a[i]-30 & temp1$rt.a<=temp1$rt.a[i]+30)]))
# }
#
# points(x = sort(temp1$rt.a), y = median.error1[order(temp1$rt.a)], type = "l", col = "tomato", lwd = 2)
#
#
# x1 <- data.frame(temp1$rt.a, median.error1, sd.error1, stringsAsFactors = FALSE)
#
#
# library(ggplot2)
# ggplot(x1, aes(x = temp1.rt.a, y = median.error1)) +
#   geom_ribbon(aes(ymin=median.error1-sd.error1, ymax = median.error1+sd.error1), alpha = 0.3)+
#   geom_line(color = "tomato")
#
#
#
#
#
# inHouse.rt <- result.rt2[[1]]
# kegg.rt <- result.rt2[[2]]
#
#
# idx1 <- match(rt.all$`in house ID`, rownames(inHouse.rt))
# rt.a <- rt.all$`RT （min）`*60
# rt.b <- inHouse.rt$RT[idx1]
#
# temp2 <- data.frame(rt.a, rt.b, stringsAsFactors = FALSE)
# temp2 <- temp2[!is.na(rt.b),]
#
# relative.error2 <- abs(temp2$rt.a - temp2$rt.b)*100/temp2$rt.a
#
# temp2 <- data.frame(temp2, relative.error2, stringsAsFactors = FALSE)
# plot(temp2$rt.a, temp2$relative.error2, ylim = c(0,1500), pch = 19,
#      xlab = "Experiment RT (second)", ylab = "Relative error (%)",cex.lab = 1.8, cex.axis = 1.5)
# plot(temp2$rt.a, temp2$relative.error2, ylim = c(0,300), pch = 19,
#      xlab = "Experiment RT (second)", ylab = "Relative error (%)",cex.lab = 1.8, cex.axis = 1.5)
# abline(v=300, lty = 2, col = "tomato", lwd = 2)
#
#
# absoult.error2 <- (temp2$rt.b - temp2$rt.a)
# temp2 <- data.frame(temp2, absoult.error2, stringsAsFactors = FALSE)
# plot(temp2$rt.a, abs(temp2$absoult.error2),
#      # ylim = c(-600,600),
#      pch = 19,
#      xlab = "Experiment RT (second)", ylab = "Absoult error (%)",cex.lab = 1.8, cex.axis = 1.5)
# abline(v=300, lty = 2, col = "tomato", lwd = 2)
#
# median.error2 <- NULL
# sd.error2 <- NULL
# for(i in 1:nrow(temp2)){
#   median.error2[i] <- median(abs(temp2$absoult.error2[which(temp2$rt.a>temp2$rt.a[i]-30 & temp2$rt.a<=temp2$rt.a[i]+30)]))
#   sd.error2[i] <- sd(abs(temp2$absoult.error2[which(temp2$rt.a>temp2$rt.a[i]-30 & temp2$rt.a<=temp2$rt.a[i]+30)]))
# }
#
# plot(x = sort(temp1$rt.a), y = median.error1[order(temp1$rt.a)], type = "l", col = "black", lwd = 2)
# points(x = sort(temp2$rt.a), y = median.error2[order(temp2$rt.a)], type = "l", col = "tomato", lwd = 2)
#
#
# x2 <- data.frame(temp2$rt.a, median.error2, sd.error2, stringsAsFactors = FALSE)
#
#
# library(ggplot2)
# ggplot(x2, aes(x = temp2.rt.a, y = median.error2)) +
#   geom_ribbon(aes(ymin=median.error2-sd.error2, ymax = median.error2+sd.error2), alpha = 0.2, fill = "salmon")+
#   geom_line(data = x1, aes(x = temp1.rt.a, y = median.error1), color = "lightseagreen") +
#   geom_ribbon(aes(ymin=median.error1-sd.error1, ymax = median.error1+sd.error1), alpha = 0.2, fill = "lightseagreen")+
#   geom_line(color = "salmon") + theme_bw()+
#   theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         axis.title = element_text(),
#         panel.grid.minor = element_blank())+
#   xlab("Experiment RT (second)")+
#   ylab("Absolute error (second)")+
#   scale_color_manual(values = c("All Adduct" = "lightseagreen", "M+H" = "salmon"))
#
#
# #-----------------------------------------------------------------
# median.error1 <- NULL
# sd.error1 <- NULL
# for(i in 1:nrow(temp1)){
#   median.error1[i] <- median(abs(temp1$relative.error1[which(temp1$rt.a>temp1$rt.a[i]-30 & temp1$rt.a<=temp1$rt.a[i]+30)]))
#   sd.error1[i] <- sd(abs(temp1$relative.error1[which(temp1$rt.a>temp1$rt.a[i]-30 & temp1$rt.a<=temp1$rt.a[i]+30)]))
# }
#
#
#
# median.error2 <- NULL
# sd.error2 <- NULL
# for(i in 1:nrow(temp2)){
#   median.error2[i] <- median(abs(temp2$relative.error2[which(temp2$rt.a>temp2$rt.a[i]-30 & temp2$rt.a<=temp2$rt.a[i]+30)]))
#   sd.error2[i] <- sd(abs(temp2$relative.error2[which(temp2$rt.a>temp2$rt.a[i]-30 & temp2$rt.a<=temp2$rt.a[i]+30)]))
# }
#
#
#
# x2_2 <- data.frame(x2$temp2.rt.a, median.error2, x2$sd.error2)
# x1_2 <- data.frame(x1$temp1.rt.a, median.error1, x1$sd.error1)
#
# colnames(x1_2)[1] <- "temp1.rt.a"
# colnames(x2_2)[1] <- "temp2.rt.a"
#
# library(ggplot2)
# ggplot(x2_2, aes(x = temp2.rt.a, y = median.error2)) +
#   geom_ribbon(aes(ymin=median.error2-sd.error2, ymax = median.error2+sd.error2), alpha = 0.2, fill = "salmon")+
#   geom_line(data = x1_2, aes(x = temp1.rt.a, y = median.error1), color = "lightseagreen") +
#   geom_ribbon(aes(ymin=median.error1-sd.error1, ymax = median.error1+sd.error1), alpha = 0.2, fill = "lightseagreen")+
#   geom_line(color = "salmon") + theme_bw()+
#   theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         axis.title = element_text(),
#         panel.grid.minor = element_blank())+
#   xlab("Experiment RT (second)")+
#   ylab("Relative error (%)")+
#   scale_color_manual(values = c("All Adduct" = "lightseagreen", "M+H" = "salmon"))
#
# par(mar = c(5,5,4,2))
# plot(temp2$rt.a, temp2$rt.b, cex.lab = 1.8, cex.axis = 1.5, xlab = "Experiment RT (second)",
#      ylab = "Predicted RT (second)", ylim = c(0, 1200), xlim = c(0,1200))
# abline(0, 1 ,lty = 2, col = "tomato", lwd = 1.5)
# abline(v = 300, lty = 2, col = "tomato", lwd = 1.5)
