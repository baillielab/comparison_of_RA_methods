# Run BIRRA
# x <- scan("/Users/s1718825/Documents/MAIC/ExistingMethods/1_3_1.txt", what="c", sep="\n")
args=commandArgs(T)
x <- scan(args[1], what="c", sep="\n")
y <- strsplit(x, "[[:space:]]+")
y <- lapply(y,function(x) strsplit(x, "\t"))
rankedy <- list()
i<-1
for (alist in y){
  if (alist[[3]]=="RANKED"){
    rankedy[[i]] <-alist
    i<-i+1
  }
}
y <- rankedy
y <- lapply(y, function(x) x[-1:-4])
y <- lapply(y, function(x) unlist(x))
library(RobustRankAggreg)
r <- -rankMatrix(y, 20000)
source("./runBIRRA.R")
result <- BIRRA(r)
entities <- rownames(r)
rankedEntities = entities[order(result)]
rankedEntities <- entities[order(result)]
rankedEntities <-paste(rankedEntities, collapse='-,' )
print(rankedEntities)
