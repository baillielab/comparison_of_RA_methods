#Boarda's count
library(RobustRankAggreg)
args=commandArgs(T)
Boardamethod = args[2]
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

r = aggregateRanks(glist = y, N = 20000, method = Boardamethod)

rl <- r[1]
rllist <- as.list(rl)
cr <- as.character(rllist[[1]][1:1010])
cr <-paste(cr, collapse='-,' )
print(cr)