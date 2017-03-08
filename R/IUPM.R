#' IUPM, PGOF, and CI
#'
#' Calculates the maximum likelihood estimate of infectious units per million (IUPM)
#' from a single serial limiting dilution (SLD) assay. Also calculates corresponding
#' exact and asymptotic confidence intervals, and a goodness-of-fit p-value. While
#' this package was developed with the purpose of estimating IUPM, it is applicable to SLD
#' assays in general.
#'
#' @importFrom stats dbinom optim qnorm uniroot pchisq
#'
#' @param pos  Vector of number of positive wells at each dilution level (outcome of SLD Assay)
#' @param replicates  Vector of number of replicates at each dilution level
#' @param dilutions  Vector of number of cells per well at each dilution level
#' @param monte  Number of Monte Carlo samples. Default is exact (no MC sampling), unless
#' more than 15,000 possible positive well outcomes exist, in which case 15,000 MC samples are taken.
#'  Use monte=F for exact computation.
#' @param conf.level Confidence level of the interval.
#'
#' @return \item{IUPM_MLE}{Maximum likelihood estimate of IUPM for the given outcome vector.}
#' @return \item{Exact_PGOF}{P value for goodness of fit. PGOF is the probability of an
#' experimental result as rare as or rarer than that obtained, assuming that the model is
#' correct. Low values of PGOF, (e.g. PGOF < 0.01), indicate rare or implausible experimental
#' results. Samples with a very low PGOF might be considered for retesting.}
#' @return \item{Asymp_PGOF}{P value calculated using an asymptotic Chi-Squared distribution with D-1 degrees of freedom,
#' where D is the number of dilution levels in an SLD assay.}
#' @return \item{Exact_CI}{Exact confidence interval, computed from the likelihood ratio test (recommended)}
#' @return \item{Asymp_CI}{Wald asymptotic confidence interval, based on the normal approximation to the binomial distribution.}
#'
#' @references Myers, L. E., McQuay, L. J., & Hollinger, F. B. (1994). Dilution assay statistics.
#'             Journal of Clinical Microbiology, 32(3), 732-739. DOI:10.1.1.116.1568
#'
#' @examples
#' # Duplicates row 4 of Table 4 from Myers, et. al.
#' # Myers et. al. divides IUPM space into discrete values. This package searches
#' # entire parameter space, yielding a slightly different and more accurate MLE.
#' row4 <- get.iupm(pos=c(2,1,0,0,0,0),  # Number of positive wells per dilution level
#'                  replicates=rep(2,6), # Number of replicates per dilution level
#'                  dilutions=c(1e6,2e5,4e4,8e3,1600,320), # Cells per dilution level
#'                  conf.level=0.95   # Significance level
#'                  )
#'
#' # Duplicates row 21 of Table 4 from Myers, et. al.
#' # Low PGOF example
#' # Myers et. al. divides IUPM space into discrete values. This package searches
#' # entire parameter space, yielding a slightly different and more accurate MLE.
#' row21 <- get.iupm(pos=c(2,2,2,0,1,0),
#'                  replicates=rep(2,6),
#'                  dilutions=c(1e6,2e5,4e4,8e3,1600,320),
#'                  conf.level=0.95 )
#'
#' # Monte Carlo example
#' # 67,081 total possible positive well outcomes, therefore
#' # Monte Carlo sampling is used to reduce computation time.
#' MC.example <- get.iupm(pos=c(30,9,1,0),
#'                        replicates=c(36,36,6,6),
#'                        dilutions=c(2.5e6,5e5,1e5,2.5e4),
#'                        conf.level=0.95,
#'                        monte = 5000 )

#' @import compiler
#'
#' @export

# Calculate IUPM, PGOF, and CIs (main function)
# Exact CI is default, only use Monte Carlo Sims if total rows is greater than 15000

get.iupm <- function(pos,replicates,dilutions,monte=15000,conf.level=0.95){

  alpha <- 1-conf.level
  # Calculate MLE
  opt <- optim(-6,like,pos=pos,replicates=replicates,dilutions=dilutions,lower=-18,upper=14,method="Brent",hessian=T)
  mle <- exp(opt$par)
  se <- sqrt(1/opt$hessian)

  # Upper and lower asymptotic confidence intervals
  lower.asy <- exp(opt$par-qnorm(1-alpha/2)*se)
  upper.asy <- exp(opt$par+qnorm(1-alpha/2)*se)

  # Asymptotic PGOF
  asy.pvalue<-function(n,repsize,numsuccess,f)
  {
    chisq=0
    for(ii in 1:length(n))
    {
      E.i=repsize[ii]*(1-exp(-1*f*n[ii]))
      O.i=numsuccess[ii]
      val1=((E.i-O.i)^2)/(E.i)
      val2=((E.i-O.i)^2)/(repsize[ii]-E.i)
      chisq=chisq+val1+val2
    }
    df=length(n)-1
    pval=1-pchisq(chisq,df)
    return(pval)
  }
  pgof.asy<- asy.pvalue(n=dilutions,repsize=replicates,numsuccess=pos,f=mle/10^6)

  # Define number of Monte Carlo samples
  if (monte==F) mc <- prod(replicates+1)
  else mc <- min(monte,prod(replicates+1))

  # convert - randomly generate outcomes of positive wells "rows" times, without replacement
  rows <- sample(1:prod(replicates+1),mc,replace=F)
  crows <- convert(rows,replicates)

  # MLEs for all outcomes
  mles <- matrix(0,mc,1)
  for (ii in 1:mc){
    opt <- optim(-6,like,pos=crows[ii,],replicates=replicates,dilutions=dilutions,lower=-18,upper=14,method="Brent")
    mles[ii] <- exp(opt$par)
  }

  # PGOF
  gof.p <- gofp(pos,replicates,dilutions,mle,mles,crows,mc)

  # Exact Confidence Intervals

  # If p-value at 0 is not less than 0.05, set lower exact interval to 0
  if (getp.alphac(0,pos,replicates,dilutions,mle,mles,crows,mc,alpha)>0){
    lower <- 0
  }else{
    lower <- uniroot(getp.alphac,c(0,mle),pos,replicates,dilutions,mle,mles,crows,mc,alpha,extendInt="up")$root
  }

  # If p-value at 1e6 is not less than 0.05, set upper exact interval to infinity
  # (We cannot say anything about the upper bound in this case)
  if (getp.alphac(1e6,pos,replicates,dilutions,mle,mles,crows,mc,alpha)>0){
    upper <- Inf
    mle <- Inf
    lower.asy= NA
    upper.asy= NA
  }else{
    upper <- uniroot(getp.alphac,c(mle,1e6),pos,replicates,dilutions,mle,mles,crows,mc,alpha,extendInt="down")$root
  }

  # If Monte Carlo sampling is used, print sampling fraction
  if (mc<prod(replicates+1)){
    sampfrac=mc/prod(replicates+1)
    print(paste("NOTE: MC approximation used. Sampling fraction is",round(sampfrac,4)))
  }

  # If all wells are positive, print warning
  if (identical(pos,replicates)){
    print("WARNING: All wells are positive, therefore MLE does not exist")
  }

  asy=c(lower.asy,upper.asy)
  exact=c(lower,upper)
  return(list("IUPM_MLE"=mle,"Exact_PGOF"=gof.p,"Exact_CI"=exact, "Asymp_PGOF"=pgof.asy,"Asymp_CI"=asy))
}
