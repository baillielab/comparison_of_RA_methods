#' (Badgeley et al., 2015) provided by author.
BIRRA=function(data, prior=0.05, num.bin=50, num.iter=10, return.all=F, plot.legend=F, grp=NULL, cor.stop=1, ...){
  nr=nrow(data)
  nrp=floor(nrow(data)*prior)
  data=apply(-data,2,rank)/nr
  
  nc=ncol(data)
  TPR=FPR=Bayes.factors=matrix(ncol=nc, nrow=num.bin)
  binned.data=ceiling(data*num.bin)
  
  
  bayes.data=matrix(nrow=nrow(data), ncol=ncol(data))
  
  guess=apply(data,1,mean)
  cprev=0
  #par(mfrow=c(floor(sqrt(num.iter)), ceiling(sqrt(num.iter))), mai=rep(0.7,4))
  for ( iter in 1:num.iter){
    if((cor.stop-cprev)>1e-15){
      guesslast=guess
      oo=order(guess)
      guess[oo[1:nrp]]=1
      guess[oo[(nrp+1):nr]]=0
      
      
      for (i in 1:nc){  
        for (bin in 1:num.bin){
          frac=bin/num.bin
          TPR=sum(guess[binned.data[,i]<=bin])
          FPR=sum((!guess)[binned.data[,i]<=bin])
          
          Bayes.factors[bin,i]=log((TPR+1)/(FPR+1)/(prior/(1-prior)))
          
        }
      }
      
      Bayes.factors=apply(Bayes.factors,2,smooth)
      Bayes.factors=apply(Bayes.factors,2,function(x){rev(cummax(rev(x)))})
      # Plot TPR vs bin for each data set
      # if(is.null(grp)){
      #   matplot(1:num.bin, Bayes.factors, type="l", lwd=2, ...)
      # }
      # else{
      #   matplot(1:num.bin, Bayes.factors, type="l", lwd=2, lty=grp, col=grp)
      # }
      #
      # title(paste("Iteration", iter))
      # if (iter==1&plot.legend){
      #   legend("topright", col=1:5, lty=1:4, legend=colnames(data), lwd=2, ncol=2)
      # }
      for (bin in 1:num.bin){
        oo=order(Bayes.factors[bin,], decreasing=T)
        Bayes.factors[bin, oo[1]]=Bayes.factors[bin, oo[2]]
        
        
      }
      
      for (i in 1:nc){
        
        bayes.data[,i]=Bayes.factors[binned.data[,i],i]
        
      }
      
      
      bb=exp(apply(bayes.data,1, sum))
      f=prior/(1-prior)
      prob=bb*f/(1+bb*f)
      exp=sort(prob, decreasing=F)[nrp]
      
      guess=rank(-apply(bayes.data,1, sum))
      cprev=cor(guess, guesslast)
      if(is.na(cprev)){
        message("correlation with pervious is manually set to be 0 because the standard deviation is zero")
        cprev=0
      }
      else{
        message("correlation with pervious iteration=",cprev)
      }
    }
    else{
      message("Converged");
      break
    }
  }
  if(return.all){
    return(list(result=guess, data=bayes.data, BF=Bayes.factors))
  }
  else{
    guess
  }
}


