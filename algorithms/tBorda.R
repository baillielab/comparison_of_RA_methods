# Run Borda's methods using packages("TopKLists").
# install.packages("TopKLists")
library("TopKLists")
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
outBorda=Borda(y)
if(args[2] == "tMEAN"){
    rllist <- as.list(outBorda$TopK)$mean
}else if(args[2] == "tMED"){
    rllist <- as.list(outBorda$TopK)$median
}else if(args[2] == "tGEO"){
    rllist <- as.list(outBorda$TopK)$geo.mean
}else if(args[2] == "tL2"){
    rllist <- as.list(outBorda$TopK)$l2norm
}else{
    rllist <- "Error: unrecognised name"
}
cr <- as.character(rllist[1:1010])
cr <-paste(cr, collapse='-,' )
print(cr)
