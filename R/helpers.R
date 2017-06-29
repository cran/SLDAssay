# Helper functions

# Calculate the log likelihood of a result given the log of IUPM
like <- function(logiupm,pos,replicates,dilutions){
  iupm <- exp(logiupm)
  ans <- 0
  nlevels <- length(replicates)
  for (ii in 1:nlevels){
    xx <- iupm*dilutions[ii]/1e6
    ans <- ans + pos[ii]*log(1-exp(-xx)) - (replicates[ii]-pos[ii])*xx
  }
  -ans
}

# For Gart: Log Likelihood, Fisher's information, first derivative of Fisher's information, and
# expectation of third derivative of log likelihood

like.gart<-function(logiupm,pos,replicates,dilutions){
  iupm<-exp(logiupm)
  ans <- 0
  nlevels <- length(replicates)
  p <- iupm/1e6
  for (ii in 1:nlevels){
    u<-dilutions[ii]
    ans<- ans+ pos[ii]*log((1-(1-p)^u)/(1-p)^u) + (replicates[ii]*u)*log((1-p))
  }
  -ans
}

fisher.der.gart<-function(p,pos,replicates,dilutions){
  ans<-0
  nlevels<-length(replicates)
  q<-1-p
  for(ii in 1:nlevels){
    k<-dilutions[ii]
    n<-replicates[ii]
    ans<-ans+ (-k^2*n*q^(k-3))*(k-2+2*q^k)/(1-q^k)^2
  }
  return(ans)
}

fisher.gart<-function(p,pos,replicates,dilutions){
  ans<-0
  nlevels<-length(replicates)
  q<-1-p
  for(ii in 1:nlevels){
    k<-dilutions[ii]
    n<-replicates[ii]
    ans<-ans+ (k^2*n*q^(k-2))/(1-q^k)
  }
  return(ans)
}

third.der.gart<-function(p,pos,replicates,dilutions){
  ans<-0
  nlevels<-length(replicates)
  q<-1-p
  for(ii in 1:nlevels){
    k<-dilutions[ii]
    n<-replicates[ii]
    ans<-ans+ (-k^2*n*q^(k-3))*(3-3*q^k-k-k*q^k)/(1-q^k)^2
  }
  return(ans)
}

B.gart<-function(p,pos,replicates,dilutions){
  dIdP<-fisher.der.gart(p,pos,replicates,dilutions)
  Ed3dp<-third.der.gart(p,pos,replicates,dilutions)
  I<-fisher.gart(p,pos,replicates,dilutions)
  ans<- -(2*dIdP+Ed3dp)/(2*I^2)
  return(ans)
}

# For exact CI - produce outcome of positive wells from an integer (change of base)
convert <- function(rows,replicates){
  # rows are the randomly selected row indices where one row is a dilution assay outcome
  nlevels <- length(replicates) #number of dilution levels
  nrows <- length(rows) #number of mc samples
  output <- matrix(0,nrows,nlevels) #initialize matrix where each row is one assay outcome
  newrows <- rows
  for (ii in 1:nlevels){
    # Ilana changed this
    # original: output[,nlevels-ii+1] <- newrows%%(replicates[nlevels-ii+1]+1)
    output[,ii] <- newrows%%(replicates[ii]+1)
    
    # Ilana changed this
    # orig: newrows <- (newrows-output[,nlevels-ii+1])/(replicates[nlevels-ii+1]+1)
    newrows <- (newrows-output[,ii])/(replicates[ii]+1)
  }
  return(output)
}

iconvert<-function(crow,replicates){
  #browser()
  nlevels<-length(replicates)
  nrow1<-0
  for(ii in 1:nlevels){
    nrow0<-(nrow1*(replicates[nlevels-ii+1]+1))+crow[nlevels-ii+1]
    #print(nrow0)
    nrow1<-nrow0
  }
  return(nrow1)
}

# Calculate the likelihood of a result given IUPM
like.long <-  function(iupm,pos,replicates,dilutions){
  # iupm can be scalar or vector
  # pos needs to be a matrix
  ans <- 1
  nlevels <- length(replicates)
  for (ii in 1:nlevels){
    xx <- iupm*dilutions[ii]/10^6
    ans <- ans * dbinom(pos[,ii],replicates[ii],1-exp(-xx))
  }
  ans
}


# LRT p-value - for confidence intervals
# pobsj - probability of observed outcome at IUPM tauj, num
# pobsmle - probability of observed outcome at MLE (tau hat), denom
# ratio2 - ratio all others being compared to
# pijs - nums for relative likelihoods we are going to compare to ratio2
# piis - denoms for relative likelihoods we are going to compare to ratio2
# mc- number of MC simulations
getp <- function(tauj,pos,replicates,dilutions,mle,mles,crows,mc){
  if (tauj==mle){
    return(1)
  }
  ans <- 0
  posvec <- matrix(pos,1,length(dilutions))
  pobsj <- 	like.long(tauj,	pos=posvec,replicates=replicates,dilutions=dilutions)
  pobsmle <- 	like.long(mle, pos=posvec,replicates=replicates,dilutions=dilutions)
  ratio2 <- pobsj/pobsmle
  pijs <- like.long(tauj,	pos=crows,replicates=replicates,dilutions=dilutions)
  piis <- like.long(mles,	pos=crows,replicates=replicates,dilutions=dilutions)
  ratio1 <- pijs/piis
  index <- ratio1 <= ratio2
  # It is possible for piis to be 0, thus causing ratio1 to be NA. The equality
  # will not hold in these cases, so exclude NA's from the sum
  ans <- sum(pijs[index], na.rm=T) #sum probabilities where ratio less than observed obs's ratio
  # slower version using loop
  # for (ii in 1:mc){
  #	crowsii <- matrix(crows[ii,],1,length(dilutions))
  # 	pij <- like.long(tauj,		pos=crowsii,replicates=replicates,dilutions=dilutions)
  #	pii <- like.long(mles[ii],	pos=crowsii,replicates=replicates,dilutions=dilutions)
  #	ratio1 <- pijs[ii]/pii
  #	if (ratio1<=ratio2) ans <- ans + pijs[ii]
  #	}
  weight <- (prod(replicates+1))/(mc) # weight for monte carlo, add 1 cuz need to count 0 case as option
  sumweight<-ans*weight
  return(sumweight)
}

# CMPfun - Speeds up code using just in time compiler
getpc <- cmpfun(getp)

# LRT p-value minus alpha (find root of this function for CI)
getp.alpha <- function(tauj,pos,replicates,dilutions,mle,mles,crows,mc,alpha){
  getpc(tauj,pos,replicates,dilutions,mle,mles,crows,mc)-alpha}

getp.alpha.log<-function(logtauj,pos,replicates,dilutions,mle,mles,crows,mc,alpha){
  tauj<-exp(logtauj)
  getpc(tauj,pos,replicates,dilutions,mle,mles,crows,mc)-alpha}

getp.alphac <- cmpfun(getp.alpha)
getp.alphac.log<-cmpfun(getp.alpha.log)

# GOF p-value - for PGOF
gofp <- function(pos,replicates,dilutions,mle,mles,crows,mc){
  ans <- 0
  posvec <- matrix(pos,1,length(dilutions))
  pobsmle <- 	like.long(mle,	pos=posvec,replicates=replicates,dilutions=dilutions)
  pimle <- 	like.long(mle,	pos=crows,replicates=replicates,dilutions=dilutions)
  index <- pimle <= pobsmle
  ans <- sum(pimle[index])
  # slower version using loop
  # for (ii in 1:mc){
  #	crowsii <- matrix(crows[ii,],1,length(dilutions))
  #	pimle <- like.long(mle,		pos=crowsii,replicates=replicates,dilutions=dilutions)
  #	if (pimle<=pobsmle) ans <- ans + pimle
  #	}
  weight <- (prod(replicates+1))/(mc)
  weightsum<-ans*weight
  return(weightsum)
}
