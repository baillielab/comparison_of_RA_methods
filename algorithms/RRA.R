# Run RRA
#install.packages("RobustRankAggreg")
library(RobustRankAggreg)

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
r = aggregateRanks(glist = y, N = 20000)

rl <- r[1]
rllist <- as.list(rl)
cr <- as.character(rllist[[1]][1:1010])
cr <-paste(cr, collapse='-,' )
print(cr)

#c <- file( "D:/programmer/PythonWorkSpace/crossvalidation_new-4-22/real_data//RRA_result_file.txt", "w" )

#writeLines( cr, c )

#close( c )
