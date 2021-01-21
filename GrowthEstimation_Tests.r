#///////////////////////////////////////////////////////////////////////////////
#### GrowthEstimation v0.3.1 ####
#///////////////////////////////////////////////////////////////////////////////


#/////////////////////////////////////////////////////////////////////////
#### LRT_2pop_fa65: LRT for two populations comparison, Fabens (1965) ####
#/////////////////////////////////////////////////////////////////////////

# LRT = likelihood ratio test

# require(TMB)
# # compile("growthestimation/FabensTwoPop.cpp")
# dyn.load(dynlib("growthestimation/FabensTwoPop"))


LRT_2pop_fa65 <- function(par=NULL,alpha=0.05,
                          L1.pop1,L2.pop1,T1.pop1,T2.pop1,
                          L1.pop2,L2.pop2,T1.pop2,T2.pop2){
  
  ### setup
  if (is.null(par)){
    par.ini <- c(0,0) # bad starting values
  } else {
    par.ini <- c(par[1],par[2]) # keep only first two for Linf and K in pop1
  }
  
  datalist <- list(
    'L1_1'=L1.pop1,
    'L2_1'=L2.pop1,
    'T1_1'=T1.pop1,
    'T2_1'=T2.pop1,
    'L1_2'=L1.pop2,
    'L2_2'=L2.pop2,
    'T1_2'=T1.pop2,
    'T2_2'=T2.pop2
  )
  parlist <- list(
    'logLinf'=log(par.ini[1]),
    'logK'=log(par.ini[2]),
    'logsigma1'=0, # pop1
    'deltaLinf'=0,
    'deltaK'=0,
    'logsigma2'=0 # pop2
  )
  
  ### fit1: unconstrained, full model where all param estimated
  obj1 <- MakeADFun(data=datalist,
                    parameters=parlist,
                    DLL="FabensTwoPop",silent=T)

  opt1 <- nlminb(start=obj1$par,obj=obj1$fn,gr=obj1$gr,
                 control=list(eval.max=2000,iter.max=2000))
  # opt1$mess # ok
  
  summ.rep1 <- summary(sdreport(obj1))
  
  
  ### fit0: constrained under H0: (deltaLinf,deltaK) = 0
  obj0 <- MakeADFun(data=datalist,
                    parameters=parlist,
                    map=list('deltaLinf'=factor(NA), # fixed at 0
                             'deltaK'=factor(NA)), # fixed at 0
                    DLL="FabensTwoPop",silent=T)
  
  opt0 <- nlminb(start=obj0$par,obj=obj0$fn,gr=obj0$gr,
                 control=list(eval.max=2000,iter.max=2000))
  # opt0$mess # ok
  
  summ.rep0 <- summary(sdreport(obj0))
  
  
  ### compute LRT test stat, chi2 critical value, p-value
  lrt <- 2*(opt0$obj-opt1$obj) # LRT stat value, compare to chi^2 with df=2
  chi2quant <- qchisq(p=alpha,df=2,lower.tail=F)
  pval <- pchisq(q=lrt,df=2,lower.tail=F)
  

  ### output
  
  par.pop1 <- summ.rep1[c('Linf','K','sigma1'),1] # estimates only
  par.pop2 <- c(summ.rep1[c('Linf','K'),1]+summ.rep1[c('deltaLinf','deltaK'),1],
                summ.rep1['sigma2',1])# estimates only
  names(par.pop1) <- c('Linf','K','sigma')
  names(par.pop2) <- c('Linf','K','sigma')
  
  res <- list('lrt.teststat'=lrt,
              'chi2.critval'=chi2quant,
              'pvalue'=pval,
              'loglik.full'=-opt1$obj,
              'loglik.constr'=-opt0$obj,
              'est.pop1'=par.pop1,
              'est.pop2'=par.pop2)
  # ^ simple output for now, no se
  
  return(res)
}




# END GrowthEstimation_Tests
