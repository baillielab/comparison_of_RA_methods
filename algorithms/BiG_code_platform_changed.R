#' (Li et al., 2018) provided by author.
#' Truncated Gamma distribution
#' 
#' Quantile function and random generation for truncated Gamma distribution with parameters \code{shape} and \code{rate}.
#' 
#' @param p vector of probabilities.
#' @param a vector of lower bounds. These may be \code{-Inf}.
#' @param b vector of upper bounds. These may be \code{Inf}.
#' @param shape,rate shape and rate parameters. Must be positive, \code{rate} strictly.
#' @examples
#' qtruncgamma(0.6,1,2,2,1)
#' rtruncgamma(5,1,2,2,1)
#' @export
qtruncgamma <- function(p,a=-Inf,b=Inf,shape,rate=1)
{stats::qgamma(stats::pgamma(a,shape=shape,rate=rate)+p*(stats::pgamma(b,shape=shape,rate=rate)-stats::pgamma(a,shape=shape,rate=rate)),shape=shape,rate=rate)}
#' @rdname qtruncgamma
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @export
#' 
rtruncgamma <- function(n,a=-Inf,b=Inf,shape,rate=1)
{u = stats::runif(n, min = 0, max = 1)
x = qtruncgamma(a=a,b=b,p=u,shape=shape,rate=rate)
return(x)}

#' Generate initial values for W
#' 
#' Generate initial values for W.
#'
#' @param r matrix that contains the ranked lists to be aggregated. \code{NA} indicates non-inclusion of item.
#' @export
#'
init_W <- function(r){W=base::matrix(NA,base::dim(r)[1],base::dim(r)[2])
rownames(W)=1:(base::dim(r)[1])
for (i in 1:base::dim(r)[2]) {W[base::order(r[,i],na.last=NA),i]=base::sort(stats::rnorm(base::sum(!is.na(r[,i]))),decreasing=TRUE)}
return(W)}

#' BiG with diffuse Inverse Gamma/Uniform prior
#' 
#' BiG implemented with diffuse Inverse Gamma prior or diffuse Uniform prior for the variance/standard deviation parameters.
#'
#' @param r G*S matrix that contains the ranked lists to be aggregated, where G is the total number of items (genes) and S is the total number of ranked lists (studies). \code{NA} indicates non-inclusion of item. Note the matrix needs to be arranged such that all lists that belong to the same platform are next to each other, i.e. the first \code{n_p1} columns are lists from platform 1.
#' @param n_T vector of length \code{S} that contains number of top ranked items in each study.
#' @param n_p1 number of studies belong to platform 1. 
#' @param M number of MCMC iterations.
#' @param burnin number of burn-in iterations.
#' @param prior either \code{"IG"} or \code{"uniform"}
#' @param W G*S matrix that contains initial values for W. Each element of W is the local importance of the corresponding item in the corresponding study, i.e. the latent variable that determines the observed rank.
#' @param mu0 vector of length \code{G} that contains initial values for mu. Each element of mu is global importance of the corresponding item, i.e. the latent variable that determines the true rank. 
#' @param kappa10,kappa20 vectors of length \code{G} that contain initial values for kappa1 and kappa2. Each element of kappa1 and kappa2 is the platform bias for the corresponding item in platform 1 and 2 respectively
#' @param sigma_s0 vector of length \code{S} that contains initial values for the variances of the study bias.
#' @param sigma_p10,sigma_p20 initial values for the variance of the platform bias for platform 1 and platform 2 respectively.
#' @param ds,dp hyperparameter for the prior distributions of variance parameters for study bias and platform bias respectively. Used only when \code{prior="IG"}.
#' @param a,b hyperparameters for the prior distributions of standard deviation parameters. Used only when \code{prior="uniform"}.
#' @examples
#' set.seed(1234)
#' sim = sim_lvm(G=25, S=6, n_p1=3, rho=runif(6,min=0.3,max=0.9), p_p1=0.6, p_p2=0.8, 
#'       lambda=runif(6,min=0.6,max=0.8), n_T=sample(c(5,10,15),6,replace=TRUE))
#' rank(-BiG_diffuse(r=sim$r, n_T=sim$n_T, n_p1=3, M=100, burnin=50, prior="IG"))
#' #rank(-BiG_diffuse(r=sim$r, n_T=sim$n_T, n_p1=3, M=100, burnin=50, prior="uniform"))
#' @export
#' 
BiG_diffuse <- function(r,n_T,n_p1='none',M=20000,burnin,prior,ds=1,dp=1,W=init_W(r),sigma_p10=0.5,sigma_p20=0.5,mu0=numeric(G),
                        kappa10=numeric(G),kappa20=numeric(G),sigma_s0=rep(1,S),a=0.0202,b=98.5025){
  G=dim(r)[1]
  S=dim(r)[2]
  Kappa1=Mu=matrix(0,G,M)
  Sigma_s2=matrix(0,S,M)
  Sigma_p12=numeric(M)
  Mu[,1]=mu0
  Sigma_s2[,1]=sigma_s0
  Kappa1[,1]=kappa10
  Sigma_p12[1]=sigma_p10
  for (m in 2:M)
  {
    for (g in 1:G)
    {
      ## update mu
      c=sum((W[g,is.na(W[g,])==F]-rep(Kappa1[g,m-1],S)[is.na(W[g,])==F])/Sigma_s2[is.na(W[g,])==F,m-1])
      d=sum(1/Sigma_s2[is.na(W[g,])==F,m-1])+1
      Mu[g,m]=stats::rnorm(1,c/d,sqrt(1/d))
      ## update kappa
      c=sum((W[g,is.na(W[g,])==F]-Mu[g,m])/Sigma_s2[is.na(W[g,])==F,m-1])
      d=sum(1/Sigma_s2[is.na(W[g,])==F,m-1])+1/Sigma_p12[m-1]
      Kappa1[g,m]=stats::rnorm(1,c/d,sqrt(1/d))
    }
    ## updata sigma_s2
    for (s in 1:S)
    {
      if (prior=="IG") {Sigma_s2[s,m]=1/stats::rgamma(1,shape=sum(is.na(W[,s])==F)/2+ds,
                                                      rate=sum((W[,s]-Mu[,m]-Kappa1[,m])^2,na.rm=TRUE)/2+ds)}
      else if (prior=="uniform") {Sigma_s2[s,m]=1/rtruncgamma(1,a=a,b=b,shape=(sum(is.na(W[,s])==F)-1)/2,
                                                              rate=sum((W[,s]-Mu[,m]-Kappa1[,m])^2,na.rm=TRUE)/2)}
    }
    #update sigma_p2
    P1=as.numeric(base::names(base::apply(is.na(W),1,prod)==0))
    if (prior=="IG") {Sigma_p12[m]=1/stats::rgamma(1,shape=length(P1)/2+dp, rate=sum(Kappa1[P1,m]^2,na.rm=TRUE)/2+dp)}
    # Sigma_p22[m]=1/stats::rgamma(1,shape=length(P2)/2+dp,rate=sum(Kappa2[P2,m]^2,na.rm=TRUE)/2+dp)}
    else if (prior=="uniform") {Sigma_p12[m]=1/rtruncgamma(1,a=a,b=b,shape=(length(P1)-1)/2, rate=sum(Kappa1[P1,m]^2,na.rm=TRUE)/2)}
    # Sigma_p22[m]=1/rtruncgamma(1,a=a,b=b,shape=(length(P2)-1)/2, rate=sum(Kappa2[P2,m]^2,na.rm=TRUE)/2)}
    #update W
    for (l in 1:10)
    {
      for (s in 1:(S))
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {
            if(n_T[s]<length(W[is.na(W[,s])==F,s])){
              aa=max(W[which(r[,s]==(n_T[s]+1)),s])
            }
            else{
              aa=-Inf
            }
            bb=W[which(r[,s]==(n_T[s]-1)),s]
          } 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=truncnorm::rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+Kappa1[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        if(n_T[s]<length(W[is.na(W[,s])==F,s])){
          W[which(r[,s]==(n_T[s]+1)),s]=truncnorm::rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=TRUE),b=W[which(r[,s]==n_T[s]),s],
                                                              mean=Mu[which(r[,s]==n_T[s]+1),m]+Kappa1[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
        }
      }
    }
  }
  post.mean.mu=base::apply(Mu[,ceiling(burnin):M],1,mean)
  return(post.mean.mu)
}

#' BiG with half-t prior through a data augmentation approach
#' 
#' BiG implemented with half-t prior for the standard deviation parameters of the platform bias and diffuse uniform prior for the variance parameters of the study bias.
#'
#' @inheritParams BiG_diffuse
#' @param xi10,xi20 vectors of length \code{G} that contain initial values for xi1 and xi2.
#' @examples 
#' set.seed(1234)
#' sim = sim_lvm(G=25, S=6, n_p1=3, rho=runif(6,min=0.3,max=0.9), p_p1=0.6, p_p2=0.8, 
#'       lambda=runif(6,min=0.6,max=0.8), n_T=sample(c(5,10,15),6,replace=TRUE))
#' BiG_DA(r=sim$r, n_T=sim$n_T, n_p1=3, M=100, burnin=50)
#' @export
BiG_DA <- function(r,n_T,n_p1,M=20000,burnin,a=0.0202,b=98.5025,dp=1,W=init_W(r),sigma_p10=0.5,sigma_p20=0.5,mu0=numeric(G),
                   xi10=numeric(G),xi20=numeric(G),sigma_s0=rep(1,S)){
  G=dim(r)[1]
  S=dim(r)[2]
  xi2=xi1=Mu=matrix(0,G,M)
  Sigma_s2=matrix(0,S,M)
  alpha1=alpha2=Sigma_p12=Sigma_p22=numeric(M)
  Mu[,1]=mu0
  Sigma_s2[,1]=sigma_s0
  xi1[,1]=xi10
  xi2[,1]=xi20
  alpha1[1]=alpha2[1]=1
  Sigma_p12[1]=sigma_p10
  Sigma_p22[1]=sigma_p20
  for (m in 2:M)
  {
    for (g in 1:G)
    {
      ## update mu
      c=sum((W[g,is.na(W[g,])==F]-c(rep(alpha1[m-1]*xi1[g,m-1],n_p1),rep(alpha2[m-1]*xi2[g,m-1],n_p1))[is.na(W[g,])==F])/Sigma_s2[is.na(W[g,])==F,m-1])
      d=sum(1/Sigma_s2[is.na(W[g,])==F,m-1])+1
      Mu[g,m]=stats::rnorm(1,c/d,sqrt(1/d))
      ## update xi
      c=alpha1[m-1]*sum((W[g,c(is.na(W[g,1:(n_p1)])==F,rep(F,n_p1))]-Mu[g,m])/Sigma_s2[c(is.na(W[g,1:(n_p1)])==F,rep(F,n_p1)),m-1])
      d=alpha1[m-1]^2*sum(1/Sigma_s2[c(is.na(W[g,1:(n_p1)])==F,rep(F,n_p1)),m-1])+1/Sigma_p12[m-1]
      xi1[g,m]=stats::rnorm(1,c/d,sqrt(1/d))
      c=alpha2[m-1]*sum((W[g,c(rep(F,n_p1),is.na(W[g,(n_p1+1):S])==F)]-Mu[g,m])/Sigma_s2[c(rep(F,n_p1),is.na(W[g,(n_p1+1):S])==F),m-1])
      d=alpha2[m-1]^2*sum(1/Sigma_s2[c(rep(F,n_p1),is.na(W[g,(n_p1+1):S])==F),m-1])+1/Sigma_p22[m-1]
      xi2[g,m]=stats::rnorm(1,c/d,sqrt(1/d))
    }
    ##update alpha
    P1=as.numeric(base::names(base::apply(is.na(W[,1:(n_p1)]),1,prod)==0))
    P2=as.numeric(base::names(base::apply(is.na(W[,(n_p1+1):S]),1,prod)==0))
    c=sum(t((W[P1,1:(n_p1)]-Mu[P1,m])*xi1[P1,m])/Sigma_s2[1:(n_p1),m-1],na.rm=TRUE)
    d=sum(xi1[P1,m]^2)*sum(1/Sigma_s2[1:(n_p1),m-1])+1
    alpha1[m]=stats::rnorm(1,c/d,sqrt(1/d))
    c=sum(t((W[P2,(n_p1+1):S]-Mu[g,m])*xi2[P2,m])/Sigma_s2[(n_p1+1):S,m-1],na.rm=TRUE)
    d=sum(xi2[P2,m]^2)*sum(1/Sigma_s2[(n_p1+1):S,m-1])+1
    alpha2[m]=stats::rnorm(1,c/d,sqrt(1/d))
    #update sigma_p2
    Sigma_p12[m]=1/stats::rgamma(1,shape=length(P1)/2+1,
                                 rate=sum(xi1[P1,m]^2,na.rm=TRUE)/2+1)
    Sigma_p22[m]=1/stats::rgamma(1,shape=length(P2)/2+1,
                                 rate=sum(xi2[P2,m]^2,na.rm=TRUE)/2+1)
    ##update sigma_s2
    for (s in 1:(n_p1))
    {
      Sigma_s2[s,m]=1/rtruncgamma(1,a=a,b=b,shape=(sum(is.na(W[,s])==F)-1)/2,
                                  rate=sum((W[,s]-Mu[,m]-alpha1[m]*xi1[,m])^2,na.rm=TRUE)/2)
    }
    for (s in (n_p1+1):S)
    {
      Sigma_s2[s,m]=1/rtruncgamma(1,a=a,b=b,shape=(sum(is.na(W[,s])==F)-1)/2,
                                  rate=sum((W[,s]-Mu[,m]-alpha2[m]*xi2[,m])^2,na.rm=TRUE)/2)
    }
    
    #update W
    for (l in 1:10)
    {
      for (s in 1:(n_p1))
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=truncnorm::rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+alpha1[m]*xi1[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=truncnorm::rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=TRUE),b=W[which(r[,s]==n_T[s]),s],
                                                            mean=Mu[which(r[,s]==n_T[s]+1),m]+alpha1[m]*xi1[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
      for (s in (n_p1+1):S)
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=truncnorm::rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+alpha2[m]*xi2[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=truncnorm::rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=TRUE),b=W[which(r[,s]==n_T[s]),s],
                                                            mean=Mu[which(r[,s]==n_T[s]+1),m]+alpha2[m]*xi2[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
    }
  }
  post.mean.mu=base::apply(Mu[,ceiling(burnin):M],1,mean)
  return(post.mean.mu)
}


#' Simulate rank data from latent variable model
#' 
#' Simulate rank data from latent variable model
#'
#' @inheritParams BiG_diffuse
#' @param G total number of genes involved in all of the studies.
#' @param S number of studies (ranked lists) to be aggregated.
#' @param rho correlation between local importance (w) and global importance (mu) for each study, which determines the total variance of w.
#' @param p_p1,p_p2 percentage of total variance of w contributed by platform variance from platform 1 and platform 2 respectively for the study with the lowest total variance. 
#' @param lambda inclusion rate for each study.
#' @examples 
#' set.seed(1234)
#' sim_lvm(G=25, S=6, n_p1=3, rho=runif(6,min=0.3,max=0.9), p_p1=0.6, p_p2=0.8, 
#'         lambda=runif(6,min=0.6,max=0.8), n_T=sample(c(5,10,15),6,replace=TRUE))
#' @export
#'
sim_lvm <- function(G,S,n_p1,rho,p_p1,p_p2,lambda,n_T){
  sigma_p2=c((max(rho)^(-2)-1)*p_p1,(max(rho)^(-2)-1)*p_p2) #platform variance 
  sigma_s2=rho^(-2)-1-rep(sigma_p2,c(n_p1,S-n_p1)) #study variance
  mu_g=stats::rnorm(G)
  names(mu_g)=1:G
  true_rank=rank(-mu_g)
  kappa=cbind(stats::rnorm(G,mean=0,sd=sqrt(sigma_p2[1])),stats::rnorm(G,mean=0,sd=sqrt(sigma_p2[2])))
  ## simulate data for each study ##
  ss=1 ## initialize ss to gurantee each gene is included in at least one study
  while(ss>0)
  {
    w=r=r.k=matrix(NA,G,S)
    rownames(w)=rownames(r)=rownames(r.k)=1:G
    colnames(r)=colnames(r.k)=1:S
    for (s in 1:(n_p1))
    {
      index=which(stats::rbinom(G,1,lambda[s])==1)
      n_T[s]=min(n_T[s],length(index))
      w[index,s]=mu_g[index]+stats::rnorm(length(index),mean=0,sd=sqrt(sigma_s2[s]))+kappa[index,1]
      r[index,s]=r.k[index,s]=rank(-w[index,s])
      r[is.na(r[,s])==F&r[,s]>n_T[s],s]=n_T[s]+1
      r.k[is.na(r[,s])==F&r[,s]>n_T[s],s]=NA
    }
    for (s in (n_p1+1):S)
    {
      index=which(stats::rbinom(G,1,lambda[s])==1)
      n_T[s]=min(n_T[s],length(index))
      w[index,s]=mu_g[index]+stats::rnorm(length(index),mean=0,sd=sqrt(sigma_s2[s]))+kappa[index,2]
      r[index,s]=r.k[index,s]=rank(-w[index,s])
      r[is.na(r[,s])==F&r[,s]>n_T[s],s]=n_T[s]+1
      r.k[is.na(r[,s])==F&r[,s]>n_T[s],s]=NA
    }
    ss=sum(apply(is.na(w),1,sum)==S)
  }  
  return(list(r=r,true_rank=true_rank,n_T=n_T))
}

