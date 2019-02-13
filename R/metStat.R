# dir()
# data.pos <- readr::read_csv("data.pos.csv")
# data.pos <- as.data.frame(data.pos)
#
# data.neg <- readr::read_csv("data.neg.csv")
# data.neg <- as.data.frame(data.neg)
#
# annotation.pos <- readr::read_csv("MRN.annotation.result.pos.csv")
# annotation.pos <- as.data.frame(annotation.pos)
#
# annotation.neg <- readr::read_csv("MRN.annotation.result.neg.csv")
# annotation.neg <- as.data.frame(annotation.neg)
#
# annotation.pos <- removeIsotope(annotation.pos)
# annotation.neg <- removeIsotope(annotation.neg)
#
# sum(data.pos$name == annotation.pos$name)
# sum(data.neg$name == annotation.neg$name)
#
# sample.info <- read.csv("sample.info.csv", stringsAsFactors = FALSE)
#
# colnames(data.pos)
# colnames(annotation.pos)
#
# data.pos <- data.frame(annotation.pos, data.pos[,-c(1:3)], stringsAsFactors = FALSE)
# data.neg <- data.frame(annotation.neg, data.neg[,-c(1:3)], stringsAsFactors = FALSE)
#
# ##remove the NA
# data.pos <- data.pos[which(!is.na(data.pos$ID)),]
# data.neg <- data.neg[which(!is.na(data.neg$ID)),]
#
#
# sample.name <- sample.info$sample.name[sample.info$group %in% c("1011IPS", "1011NPC")]
#
# sample.pos <- data.pos[,match(sample.name, colnames(data.pos))]
# sample.neg <- data.neg[,match(sample.name, colnames(data.neg))]
#
# sample.info <- sample.info[sample.info$group %in% c("1011IPS", "1011NPC"),]
#
# p.pos <- uniTest(sample = sample.pos, sample.info = sample.info,
#                  uni.test = "t", correct = TRUE)
#
# p.neg <- uniTest(sample = sample.neg, sample.info = sample.info,
#                  uni.test = "t", correct = TRUE)
#
# p.index.pos <- which(p.pos < 0.01)
# p.index.neg <- which(p.neg < 0.01)
#
#
# data.pos <- data.pos[p.index.pos,]
# data.neg <- data.neg[p.index.neg,]
#
#
#
# annotation.table <- rbind(annotation.table.pos, annotation.table.neg)
# annotation.table <- annotation.table[,c("name", "ID")]
# colnames(annotation.table) <- c("peak.name", "KEGG.ID")
#
#
#
#
#
#
#
#
