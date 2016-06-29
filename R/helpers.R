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

# For exact CI - produce outcome of positive wells from an integer (change of base)
convert <- function(rows,replicates){
  nlevels <- length(replicates)
  nrows <- length(rows)
  output <- matrix(0,nrows,nlevels)
  newrows <- rows
  for (ii in 1:nlevels){
    output[,nlevels-ii+1] <- newrows%%(replicates[nlevels-ii+1]+1)
    newrows <- (newrows-output[,nlevels-ii+1])/(replicates[nlevels-ii+1]+1)
  }
  output
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
  ans <- 0
  posvec <- matrix(pos,1,length(dilutions))
  pobsj <- 	like.long(tauj,		pos=posvec,replicates=replicates,dilutions=dilutions)
  pobsmle <- 	like.long(mle,		pos=posvec,replicates=replicates,dilutions=dilutions)
  ratio2 <- pobsj/pobsmle
  pijs <- like.long(tauj,		pos=crows,replicates=replicates,dilutions=dilutions)
  piis <- like.long(mles,		pos=crows,replicates=replicates,dilutions=dilutions)
  ratio1 <- pijs/piis
  index <- ratio1 <= ratio2
  ans <- sum(pijs[index]) #sum probabilities where ratio less than observed obs's ratio
  # slower version using loop
  # for (ii in 1:mc){
  #	crowsii <- matrix(crows[ii,],1,length(dilutions))
  # 	pij <- like.long(tauj,		pos=crowsii,replicates=replicates,dilutions=dilutions)
  #	pii <- like.long(mles[ii],	pos=crowsii,replicates=replicates,dilutions=dilutions)
  #	ratio1 <- pijs[ii]/pii
  #	if (ratio1<=ratio2) ans <- ans + pijs[ii]
  #	}
  weight <- prod(replicates+1)/mc # weight for monte carlo, add 1 cuz need to count 0 case as option
  ans*weight
}

# CMPfun - Speeds up code using just in time compiler
getpc <- cmpfun(getp)

# LRT p-value minus alpha (find root of this function for CI)
getp.alpha <- function(tauj,pos,replicates,dilutions,mle,mles,crows,mc,alpha){
  getpc(tauj,pos,replicates,dilutions,mle,mles,crows,mc)-alpha}

getp.alphac <- cmpfun(getp.alpha)

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
  weight <- prod(replicates+1)/mc
  ans*weight
}
