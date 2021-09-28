#run BiG
source("./BiG_code_platform_changed.R")
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

matrix_transfer<-function (glist, N = NA, full = FALSE)
{
  u = unique(c(glist, recursive = TRUE))
  if (all(is.na(N))) {
    N = length(u)
  }
  rmat = matrix(NA, nrow = length(u), ncol = length(glist),
                dimnames = list(u, names(glist)))
  N = unlist(lapply(glist, length))
  if(!full){
    for (i in 1:length(glist)) {
      rmat[, i] = (1 + length(glist[[i]]))
    }
  }
  for (i in 1:length(glist)) {
    rmat[glist[[i]], i] = (1:length(glist[[i]]))
  }
  return(rmat)
}

if(args[2] == "bottom"){
    r <- matrix_transfer(y, full=FALSE)
}else{
    r <- matrix_transfer(y, full=TRUE)
}

NTlength <- rep(1000,length(y))
for (i in 1:length(y)) {
  NTlength[i] <- length(y[[i]])
}
result <- BiG_diffuse(r=r, n_T=NTlength, n_p1=0, M=2000, burnin=1000, prior="IG")
entities <- rownames(r)
rankedEntities <- entities[order(result,decreasing=TRUE)]
rankedEntities <-paste(rankedEntities, collapse='-,' )
print(rankedEntities)

