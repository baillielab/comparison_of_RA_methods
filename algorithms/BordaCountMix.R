#Boarda's count with ties
#function to build the matrix for the mix of ranked and unranked data
tiedMatrix <-function (glist, isrankedMark, N = NA, full = FALSE)
{
  u = unique(c(glist, recursive = TRUE))
  if (all(is.na(N))) {
    N = length(u)
  }
  if (!full) {
    rmat = matrix(1, nrow = length(u), ncol = length(glist),
                  dimnames = list(u, names(glist)))
    if (length(N) == 1) {
      N = rep(N, ncol(rmat))
    }
  }
  else {
    rmat = matrix(NA, nrow = length(u), ncol = length(glist),
                  dimnames = list(u, names(glist)))
    N = unlist(lapply(glist, length))
  }
  for (i in 1:length(glist)) {
    if(isrankedMark[i]==1){
      rmat[glist[[i]], i] = (1:length(glist[[i]]))/N[i]
    }
    else{
      rmat[glist[[i]], i] = rep(length(glist[[i]])/2, length(glist[[i]]))/N[i]
    }
  }
  return(rmat)
}

library(RobustRankAggreg)
args=commandArgs(T)
Boardamethod = args[2]
x <- scan(args[1], what="c", sep="\n")
#Boardamethod = 'median'
#x <- scan('/Users/s1718825/Documents/MAIC/ExistingMethods/1_3_1.txt', what="c", sep="\n")
y <- strsplit(x, "[[:space:]]+")
y <- lapply(y,function(x) strsplit(x, "\t"))
rankedy <- list()
isrankedMark <- rep(0, length(y))
i<-1
for (alist in y){
  if (alist[[3]]=="RANKED"){
    isrankedMark[i] <- 1
  }
  i<-i+1
}
y <- lapply(y, function(x) x[-1:-4])
y <- lapply(y, function(x) unlist(x))
r = aggregateRanks(glist = y,rmat = tiedMatrix(glist=y, isrankedMark, N = 20000, full = FALSE), N = 20000, method = Boardamethod)
rl <- r[1]
rllist <- as.list(rl)
cr <- as.character(rllist[[1]][1:1010])
cr <-paste(cr, collapse='-,' )
print(cr)
