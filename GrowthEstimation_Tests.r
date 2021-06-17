#///////////////////////////////////////////////////////////////////////////////
#### GrowthEstimation v0.4.3 ####
#///////////////////////////////////////////////////////////////////////////////


#/////////////////////////////////////////////////////////////////////////////////
#### BF_2pop_Bfa65: BF for two populations comparison, Bayesian Fabens (1965) ####
#/////////////////////////////////////////////////////////////////////////////////

# BF = Bayes factor
# M0 = constrained model, M1 = full model

# require(TMB)
# compile("growthestimation/FabensTwoPopBayesian_M0.cpp")
# compile("growthestimation/FabensTwoPopBayesian_M1.cpp")
# dyn.load(dynlib("growthestimation/FabensTwoPopBayesian_M0")) 
# dyn.load(dynlib("growthestimation/FabensTwoPopBayesian_M1"))

BF_2pop_Bfa65 <- function(par=NULL,
                          L1.pop1,L2.pop1,deltaT.pop1,
                          L1.pop2,L2.pop2,deltaT.pop2,
                          priordist.M0,priordist.M1.pop1,priordist.M1.pop2,
                          hyperpar.M0,hyperpar.M1.pop1,hyperpar.M1.pop2){
  
  # priordist: 0=no prior, 1=unif, 2=Gaussian, 3=lognormal
  # hyperpar: list of 3 elements (Linf,K,sigma), each a vector of length 2
  
  ### setup
  if (is.null(par)){
    par.ini <- c(1,0.5) # bad starting values
    warning('No par supplied, using default ini (Linf,K) = c(1,0.5).')
  } else {
    par.ini <- c(par[1],par[2]) # keep only first two for Linf and K
  }
  
  ### M0: constrained model, same (Linf,K) but sigma1 != sigma2
  datalist.M0 <- list(
    'L1_1'=L1.pop1,
    'L2_1'=L2.pop1,
    'deltaT_1'=deltaT.pop1,
    'L1_2'=L1.pop2,
    'L2_2'=L2.pop2,
    'deltaT_2'=deltaT.pop2,
    'hp_Linf'=hyperpar.M0[[1]],
    'hp_K'=hyperpar.M0[[2]],
    'hp_sigma1'=hyperpar.M0[[3]], # same prior for both sigma, for simplicity
    'hp_sigma2'=hyperpar.M0[[3]], # same prior for both sigma, for simplicity
    'priordist'=priordist.M0
    # ^ 0=no prior, 1=unif, 2=Gaussian, 3=lognormal
  )
  parlist.M0 <- list(
    'logLinf'=log(par.ini[1]),
    'logK'=log(par.ini[2]),
    'logsigma1'=0, # default value
    'logsigma2'=0  # default value
  )
  
  obj.M0 <- MakeADFun(data=datalist.M0, parameters=parlist.M0,
                      DLL="FabensTwoPopBayesian_M0",silent=T) # no random effect
  opt.M0 <- nlminb(start=obj.M0$par, obj=obj.M0$fn, gr=obj.M0$gr,
                   control=list(eval.max=5000,iter.max=5000))
  # unlist(obj.M0$rep()[c('Linf','K','sigma1','sigma2')])
  # ^ TMB est = posterior modes, used as ini for Laplace approx of int
  
  
  ### M1: full model, different (Linf,K,sigma) between two pop
  datalist.M1 <- list(
    'L1_1'=L1.pop1,
    'L2_1'=L2.pop1,
    'deltaT_1'=deltaT.pop1,
    'L1_2'=L1.pop2,
    'L2_2'=L2.pop2,
    'deltaT_2'=deltaT.pop2,
    'hp_Linf1'=hyperpar.M1.pop1[[1]],
    'hp_K1'=hyperpar.M1.pop1[[2]],
    'hp_sigma1'=hyperpar.M1.pop1[[3]], # allowed to be different here
    'hp_Linf2'=hyperpar.M1.pop2[[1]],
    'hp_K2'=hyperpar.M1.pop2[[2]],
    'hp_sigma2'=hyperpar.M1.pop2[[3]], # allowed to be different here
    'priordist1'=priordist.M1.pop1,
    'priordist2'=priordist.M1.pop2
    # ^ 0=no prior, 1=unif, 2=Gaussian, 3=lognormal
  )
  parlist.M1 <- list(
    'logLinf1'=log(par.ini[1]),
    'logK1'=log(par.ini[2]),
    'logsigma1'=0, # default value
    'logLinf2'=log(par.ini[1]),
    'logK2'=log(par.ini[2]),
    'logsigma2'=0  # default value
  ) # same par.ini for both pop
  
  obj.M1 <- MakeADFun(data=datalist.M1, parameters=parlist.M1,
                      DLL="FabensTwoPopBayesian_M1",silent=T) # no random effect
  opt.M1 <- nlminb(start=obj.M1$par, obj=obj.M1$fn, gr=obj.M1$gr,
                   control=list(eval.max=5000,iter.max=5000))
  # unlist(obj.M1$rep()[c('Linf1','K1','sigma1','Linf2','K2','sigma2')])
  # ^ TMB est = posterior modes, used as ini for Laplace approx of int
  
  
  ### compute BF from Laplace approx of marginal pdf
  obj.M0.rand <- MakeADFun(data=datalist.M0,
                           parameters=as.list(opt.M0$par),
                           random=names(parlist.M0), # all integrated out
                           DLL="FabensTwoPopBayesian_M0",silent=T)
  # exp(-obj.M0.rand$fn(obj.M0.rand$par)) # marginal pdf of data|M0
  
  obj.M1.rand <- MakeADFun(data=datalist.M1,
                           parameters=as.list(opt.M1$par),
                           random=names(parlist.M1), # all integrated out
                           DLL="FabensTwoPopBayesian_M1",silent=T)
  # exp(-obj.M1.rand$fn(obj.M1.rand$par)) # marginal pdf of data|M1
  
  BF <- as.numeric(exp(-obj.M1.rand$fn(obj.M1.rand$par) +
                         + obj.M0.rand$fn(obj.M0.rand$par)))
  # ^ Laplace-approximated BF = P(data|M1)/P(data|M0), where both numerator and
  #   denominator integrals are approximated by Laplace's method (TMB's obj$fn) 
  
  # BF interpretation table from Kass and Raftery (1995, p. 777, JASA), adapted
  # from Jeffreys (1961) Theory of Probability, 3rd Ed., Appendix B:
  # 1 - 3.2  | Not worth more than a bare mention
  # 3.2 - 10 | Substantial
  # 10 - 100 | Strong
  # >100     | Decisive
  
  
  ### output
  par.M0 <- unlist(obj.M0$rep()[c('Linf','K','sigma1','sigma2')]) # est only
  par.M1 <- unlist(obj.M1$rep()[c('Linf1','K1','sigma1',
                                    'Linf2','K2','sigma2')]) # est only
  
  res <- list('BF'=BF,
              'marg.logpdf.M0'=as.numeric(-obj.M0.rand$fn(obj.M0.rand$par)),
              'marg.logpdf.M1'=as.numeric(-obj.M1.rand$fn(obj.M1.rand$par)),
              'estTMB.M0'=par.M0,
              'estTMB.M1'=par.M1)
  # ^ simple output for now, no se
  
  return(res)
}




#/////////////////////////////////////////////////////////////////////////
#### LRT_2pop_fa65: LRT for two populations comparison, Fabens (1965) ####
#/////////////////////////////////////////////////////////////////////////

# LRT = likelihood ratio test

# require(TMB)
# # compile("growthestimation/FabensTwoPop.cpp")
# dyn.load(dynlib("growthestimation/FabensTwoPop"))


LRT_2pop_fa65 <- function(par=NULL,alpha=0.05,
                          L1.pop1,L2.pop1,deltaT.pop1,
                          # T1.pop1,T2.pop1,T1.pop2,T2.pop2 # as of v0.4.3
                          L1.pop2,L2.pop2,deltaT.pop2){
  
  ### setup
  if (is.null(par)){
    par.ini <- c(1,0.5) # bad starting values
    warning('No par supplied, using default ini (Linf,K) = c(1,0.5).')
  } else {
    par.ini <- c(par[1],par[2]) # keep only first two for Linf and K in pop1
  }
  
  datalist <- list(
    'L1_1'=L1.pop1,
    'L2_1'=L2.pop1,
    # 'T1_1'=T1.pop1,       # as of v0.4.3
    # 'T2_1'=T2.pop1,       # as of v0.4.3
    'deltaT_1'=deltaT.pop1, # as of v0.4.3
    'L1_2'=L1.pop2,
    'L2_2'=L2.pop2,
    # 'T1_2'=T1.pop2,       # as of v0.4.3
    # 'T2_2'=T2.pop2        # as of v0.4.3
    'deltaT_2'=deltaT.pop2 # as of v0.4.3
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
  
  # summ.rep1 <- summary(sdreport(obj1)) # as of v0.4.3: don't need se
  summ.rep1 <- obj1$rep()
  
  ### fit0: constrained under H0: (deltaLinf,deltaK) = (0,0)
  obj0 <- MakeADFun(data=datalist,
                    parameters=parlist,
                    map=list('deltaLinf'=factor(NA), # fixed at 0
                             'deltaK'=factor(NA)), # fixed at 0
                    DLL="FabensTwoPop",silent=T)
  
  opt0 <- nlminb(start=obj0$par,obj=obj0$fn,gr=obj0$gr,
                 control=list(eval.max=2000,iter.max=2000))
  # opt0$mess # ok
  
  # summ.rep0 <- summary(sdreport(obj0)) # as of v0.4.3: don't need se
  summ.rep0 <- obj0$rep()
  
  
  ### compute LRT test stat, chi2 critical value, p-value
  lrt <- 2*(opt0$obj-opt1$obj) # LRT stat value, compare to chi^2 with df=2
  chi2quant <- qchisq(p=alpha,df=2,lower.tail=F)
  pval <- pchisq(q=lrt,df=2,lower.tail=F)
  

  ### output
  par.pop1 <- unlist(summ.rep1[c('Linf','K','sigma1')]) # estimates only
  par.pop2 <- c(unlist(summ.rep1[c('Linf','K')]) +
                  + unlist(summ.rep1[c('deltaLinf','deltaK')]),
                unlist(summ.rep1['sigma2'])) # estimates only
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
