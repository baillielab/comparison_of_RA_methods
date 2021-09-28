# Run MC1-3
# install.packages("TopKLists")
library("TopKLists")
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
outBorda=MC(y)
if(args[2] == "MC1"){
    rllist <- as.list(outBorda$MC1.TopK)
    cr <- as.character(rllist[1:1010])
    cr <-paste(cr, collapse='-,' )
}else if(args[2] == "MC2"){
    rllist <- as.list(outBorda$MC2.TopK)
    cr <- as.character(rllist[1:1010])
    cr <-paste(cr, collapse='-,' )
}else if(args[2] == "MC3"){
    rllist <- as.list(outBorda$MC3.TopK)
    cr <- as.character(rllist[1:1010])
    cr <-paste(cr, collapse='-,' )
}else if(args[2] == "MCall"){
  rllist1 <- as.list(outBorda$MC1.TopK)
  cr1 <- as.character(rllist1[1:1010])
  cr1 <-paste(cr1, collapse='-,' )
  rllist2 <- as.list(outBorda$MC2.TopK)
  cr2 <- as.character(rllist2[1:1010])
  cr2 <-paste(cr2, collapse='-,' )
  rllist3 <- as.list(outBorda$MC3.TopK)
  cr3 <- as.character(rllist3[1:1010])
  cr3 <-paste(cr3, collapse='-,' )
  cr <-paste(cr1,cr2,cr3, sep = '-MCMethodsSeperator-' )
}else{
    rllist <- "Error: unrecognised name"
}

print(cr)
