#///////////////////////////////////////////////////////////////////////////////
#### GrowthEstimation v0.3 ####
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#### GrowthPriors: uniform prior lower/upper bounds for Linf and K  ####
#///////////////////////////////////////////////////////////////////////////////


GrowthPriors <- function(Lmax=133, species="Mustelus asterias",
                         category="BodyShape", LowP=0.8, UpP=1.2, LQ=0, UQ=1){
  # Documentation v0.2.1
  #' @title Calculates priors for the von Bertalanffy growth parameters
  # Linf and k using the growth performance index and maximum size Lmax
  #' @description
  #' TODO
  #' @param Lmax positive scalar, maximum (observed) size of the species
  #' @param species character string, target species
  #' @param category character string, criterion for selecting other species in
  #' fishbase, choice among "BodyShape", "Habitat", and "Environment"
  #' @param LowP positive scalar, proportion of maximum length that defines the
  #' uniform prior lower bound for the asymptotic size Linf
  #' @param UpP positive scalar, proportion of maximum length that defines the
  #' uniform prior upper bound for the asymptotic size Linf
  #' @param LQ scalar within (0,1), probability for lower quantile of growth
  #' performance index phi; if 0 then the minimum phi is used
  #' @param UQ scalar within (0,1), probability for upper quantile of growth
  #' performance index phi; if 1 then the maximum phi is used
  
  require(rfishbase) # library(rfishbase)
  
  ### Get species information from fishbase
  SpeciesDF <- species((species_list=c(species_list(Species=species))))
  
  # SIMILAR HABITAT:
  # can be: #levels(as.factor(AllSpeciesDF$DemersPelag))
  # "bathydemersal", "bathypelagic", "benthopelagic", "demersal", "pelagic",
  # "pelagic-neritic", "pelagic-oceanic", or "reef-associated"
  
  ### Extract species habitat
  Habit <- levels(as.factor(SpeciesDF$DemersPelag))
  
  ### Find all species with same habitat in fishbase
  AllSpeciesDF <- data.frame(species())
  AllSpeciesDF$Exclude <- match(AllSpeciesDF$DemersPelag, c(Habit))
  HABITAT <- subset(AllSpeciesDF, is.na(Exclude)=="FALSE")
  
  
  # SIMILAR BODY SHAPE:
  # can be: #levels(as.factor(AllSpeciesDF$BodyShapeI))
  # "eel-like", "elongated", "Elongated", "fusiform / normal", "other",
  # "other (see remarks)", or "short and / or deep"
  
  ### Extract species body shape
  Body <- levels(as.factor(SpeciesDF$BodyShapeI))
  
  ### Find all species with same body shape in fishbase
  AllSpeciesDF$Exclude <- match(AllSpeciesDF$BodyShapeI, c(Body))
  BODY <- subset(AllSpeciesDF, is.na(Exclude)=="FALSE")
  
  
  # SIMILAR ENVIRONMENT:
  # can be: levels(as.factor(allstocksdf$EnvTemp))
  # "boreal", "deep-water", "high altitude", "polar", "subtropical",
  # "temperate", or "tropical"
  
  ### Extract species environment
  stocksdf <- data.frame(stocks(species_list=c(species_list(Species=species))))
  environ <- levels(as.factor(stocksdf$EnvTemp))
  
  ### Find all species with same environment in fishbase
  allstocksdf <- data.frame(stocks())
  allstocksdf$Exclude <- match(allstocksdf$EnvTemp, c(environ))
  
  ENV <- subset(allstocksdf, is.na(Exclude)=="FALSE")
  
  
  ### Match all species (SpecCodes)
  
  habitat <- levels(as.factor(HABITAT$SpecCode))
  bodyshape <- levels(as.factor(BODY$SpecCode))
  environ <- levels(as.factor(ENV$SpecCode))
  
  SPECCODELIST1 <- bodyshape
  SPECCODELIST2 <- intersect(bodyshape, habitat)
  SPECCODELIST3 <- intersect(bodyshape, environ)
  SPECCODELIST4 <- intersect(habitat, environ)
  
  if (category=="BodyShape"){
    SPECCODELIST <- SPECCODELIST1
  } else if (category=="Habitat"){
    SPECCODELIST <- SPECCODELIST2 
  } else if (category=="Environment"){
    SPECCODELIST <- SPECCODELIST3
  } else { # then use interesection of habitat and environ as default
    SPECCODELIST <- SPECCODELIST4
  }
  
  
  ### Get growth parameters for all species
  popdf <- data.frame(popgrowth())
  
  ### only take total length measurements for comparibility
  dat <- subset(popdf, !is.na(TLinfinity))
  
  ### make sure von Bertalanffy K is available for all of them
  dat <- subset(dat, !is.na(K))
  
  ### exclude growth experiments in captivity
  dat <- subset(dat, GrowthEnviron!="captivity")
  
  ### exclude doubtful estimates
  dat <- subset(dat, Auxim!="doubtful")
  dat <- droplevels(dat)
  
  ### select only species of interest (as defined in Category)
  dat$Exclude <- match(dat$SpecCode, SPECCODELIST)
  POP <- subset(dat, !is.na(Exclude))
  
  ### "sample" size
  nsize <- length(POP$TLinfinity)
  
  ### prior for Linf from Lmax 
  # Linf <- 10^(0.044+0.9841*log10(Lmax)) # Froese and Binohlan (2000)
  Linf <- Lmax/0.99 # simpler assumption, anyway close to Froese and Binohlan (2000)
  
  # calculate prior boundaries for Linf (uniform prior)
  # take priors for Linf as X% of Lmax
  LowLinf <- LowP*Lmax # uniform prior lower bound
  UpLinf <- UpP*Lmax # uniform prior upper bound
  
  
  ### Growth performance index phi according to Pauly (slope of -2)
  POP$phi <- log10(POP$K)+2*log10(POP$TLinfinity)
  
  # get upper and lower quantiles and median
  MedPhi <- as.numeric(quantile(POP$phi, 0.5)) # median
  LowPhi <- as.numeric(quantile(POP$phi, LQ))
  UpPhi <- as.numeric(quantile(POP$phi, UQ))
  
  
  ### Corresponding K with Linf set to its center value and 3 phi quantiles
  K <- 10^(MedPhi-2*log10(Linf))
  LowK <- 10^(LowPhi-2*log10(Linf)) 
  UpK <- 10^(UpPhi-2*log10(Linf)) 
  
  
  ### Output
  df <- data.frame('species'=species, 
                   'category'=category,
                   'Lmax'=Lmax,
                   'n'=nsize,
                   'Linf'=Linf,
                   'LowLinf'=LowLinf,
                   'UpLinf'=UpLinf,
                   'K'=K,
                   'LowK'=LowK,
                   'UpK'=UpK)
  
  return(df)
}



#///////////////////////////////////////////////////////////////////////////////
#### Bla02: Bayesian Laslett, Eveson and Polacheck (2002) ####
#///////////////////////////////////////////////////////////////////////////////

Bla02 <- function(par,L1,L2,deltaT,hyperpar=NULL,meth='nlminb',compute.se=T,
                  output.post.draws=F,lb.sd=NULL,
                  mcmc.control=list('nchains'=3,'iter'=5000,'warmup'=4000)){
  # uniform priors with user-supplied hyperparam
  # hyperpar=list(lbubmuinf,lbubK) for consisetncy across methods, while bounds
  # for lbubsigmainf, lbubmuA, lbubsigmaA and lbubsigmaeps are hard-coded below.
  # lb.sd = lower bound for sigmainf and sigmaeps, in TMB optimization only.
  
  ### setup
  if (is.null(hyperpar)){
    stop('hyperpar must be supplied as list(c(lb,ub),c(lb,ub)) respectively for
         the mean of Linf and for K.')
  }
  if (!compute.se){warning('se will be computed anyway.')}
  
  n <- length(L1)
  parlist <- list('logmuinf'=log(par[1]),'logsigmainf'=0,
                  'logK'=log(par[2]),
                  'logmuA'=1,'logsigmaA'=0,
                  'logsigmaeps'=0,
                  'logLinf'=rep(log(par[1]),n),
                  'logA'=rep(1,n))
  length.theta <- 6
  datalist <- list('L1'=L1,'L2'=L2,
                   # 'T1'=T1,'T2'=T2,
                   'deltaT'=deltaT, # as of v0.3: directly supply deltaT
                   'lbubmuinf'=hyperpar[[1]], # user-supplied
                   'lbubK'=hyperpar[[2]], # user-supplied
                   'lbubsigmainf'=c(1e-10,1e4),
                   'lbubmuA'=c(1e-10,1000),
                   'lbubsigmaA'=c(1e-10,1e4),
                   'lbubsigmaeps'=c(1e-10,1e4),
                   'enablepriors'=1)
  obj <- MakeADFun(data=datalist,parameters=parlist,
                   random=c('logLinf','logA'),
                   DLL="Laslett",silent=T)
  
  ### TMB optimization and sdreport
  if (meth=='nlminb'){
    if (is.null(lb.sd)){
      lb.par <- rep(-Inf,length.theta) # no bounds
    } else {
      lb.par <- c(-Inf,log(lb.sd),-Inf,-Inf,-Inf,log(lb.sd))
      # ^ helps if data too close to pure vB, otherwise Hessian not pos def
    }
    
    opt <- nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                  control=list(eval.max=5000,iter.max=5000),lower=lb.par)
    rep <- sdreport(obj)
    summary.rep <- summary(rep)
    theta.tmb <- summary.rep[c('muinf','K'),1]
    se.theta.tmb <- summary.rep[c('muinf','K'),2]
    # val.tmb <- obj$fn()
  } else {stop('Only meth="nlminb" is allowed for now.')}
  
  ### MCMC based on tmbstan::tmbstan, by default uses NUTS
  mcmc.obj <- tmbstan(obj=obj,
                      lower=rep(-Inf,length(unlist(parlist))),
                      upper=rep(Inf,length(unlist(parlist))),
                      silent=T,laplace=F,
                      chains=mcmc.control$nchains,
                      warmup=mcmc.control$warmup,
                      iter=mcmc.control$iter,
                      # init='random'
                      init='last.par.best' # start from MLE above
  )
  
  # traceplot(mcmc.obj, pars=c('logmuinf','logK'), inc_warmup=TRUE) # check conv
  # pairs(mcmc.obj, pars=c('logmuinf','logK')) # post dist and scatterplots
  # # ^ Linf and K typically correlate a lot (negatively), but not so linearly
  
  # extract MCMC post draws for derived quantities specified in obj's REPORT
  mcmc.post <- as.matrix(mcmc.obj)
  mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=2) # only Linf and K
  for (i in 1:nrow(mcmc.post)){
    mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[c('muinf','K')])
  }
  # colMeans(mcmc.est) # post means for Linf and K
  
  
  ### output
  res <- list('par'=c(mean(mcmc.est[,1]),mean(mcmc.est[,2])),
              'par.TMB'=theta.tmb,
              'se'=c(sqrt(var(mcmc.est[,1])),  # /length(mcmc.est[,1])
                     sqrt(var(mcmc.est[,2]))), # /length(mcmc.est[,2])
              'se.TMB'=se.theta.tmb#,
              # 'other.par'=theta[-c(1,3)],'other.par.se'=se.theta[-c(1,3)],
  )
  names(res$par) <- c('Linf','K')
  names(res$par.TMB) <- c('Linf','K')
  names(res$se) <- c('Linf','K')
  names(res$se.TMB) <- c('Linf','K')
  
  if (output.post.draws){
    # if (onlyTMB){
    #   warning('Cannot output posterior draws if onlyTMB=TRUE.')
    # } else {
    res$post.draws <- list('Linf'=mcmc.est[,1],'K'=mcmc.est[,2])
    # res$post.draws <- as.list(mcmc.post)
    # }
  }
  
  # if (!onlyTMB){ # equal-tailed 95% credible intervals based on MCMC draws
  res$cred.int <- list('Linf'=quantile(mcmc.est[,1],probs=c(0.025,0.975)),
                       'K'=quantile(mcmc.est[,2],probs=c(0.025,0.975)))
  # }
  
  res$value.TMB <- opt$obj
  res$conv.TMB <- opt$mess
  
  return(res)
}




#///////////////////////////////////////////////////////////////////////////////
#### Bfa65: Bayesian Fabens (1965) ####
#///////////////////////////////////////////////////////////////////////////////


Bfa65 <- function(par,L1,L2,deltaT,
                  priordist.Linf='uniform',
                  priordist.K='uniform',
                  priordist.sigma='uniform',
                  hyperpar=NULL,
                  meth='nlminb',compute.se=T,
                  onlyTMB=F,output.post.draws=F,
                  mcmc.control=list('nchains'=3,'iter'=5000,'warmup'=4000)){
  
  ### setup priors
  
  priordist.code <- integer(3) # codes for Linf, K, and sigma
  
  # prior for Linf
  
  if (priordist.Linf=='uniform'){
    priordist.code[1] <- 1L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.Linf="uniform" the 1st one is a vector of ',
           'length 2 of lower and upper bounds.')
    }
  } else if (priordist.Linf=='normal' | priordist.Linf=='Gaussian'){
    priordist.code[1] <- 2L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.Linf="normal" the 1st one is a vector of ',
           'length 2 of mean and sd.')
    }
  } else if (priordist.Linf=='lognormal'){
    priordist.code[1] <- 3L
    hyperpar[[1]] <- c(2*log(hyperpar[[1]][1]) +
                         - log(hyperpar[[1]][2]^2+hyperpar[[1]][1]^2)/2,
                       sqrt(-2*log(hyperpar[[1]][1]) +
                              + log(hyperpar[[1]][2]^2+hyperpar[[1]][1]^2)))
    # ^ mean and sd on log scale, user supplies mean and sd on exp scale
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.Linf="lognormal" the 1st one is a vector of ',
           'length 2 of mean and sd on the original (exponential) scale.')
    }
  } else { # then no prior specified, straight frequentist Fabens, no MCMC
    priordist.code[1] <- 0L
  }
  
  # prior for K
  
  if (priordist.K=='uniform'){
    priordist.code[2] <- 1L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.K="uniform" the 2nd one is a vector of ',
           'length 2 of lower and upper bounds.')
    }
  } else if (priordist.K=='normal' | priordist.K=='Gaussian'){
    priordist.code[2] <- 2L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.K="normal" the 2nd one is a vector of ',
           'length 2 of mean and sd.')
    }
  } else if (priordist.K=='lognormal'){
    priordist.code[2] <- 3L
    hyperpar[[2]] <- c(2*log(hyperpar[[2]][1]) +
                         - log(hyperpar[[2]][2]^2+hyperpar[[2]][1]^2)/2,
                       sqrt(-2*log(hyperpar[[2]][1]) +
                              + log(hyperpar[[2]][2]^2+hyperpar[[2]][1]^2)))
    # ^ mean and sd on log scale, user supplies mean and sd on exp scale
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.K="lognormal" the 2nd one is a vector of ',
           'length 2 of mean and sd on the original (exponential) scale.')
    }
  } else { # then no prior specified, straight frequentist Fabens, no MCMC
    priordist.code[2] <- 0L
  }
  
  # prior for sigma
  
  if (priordist.sigma=='uniform'){
    priordist.code[3] <- 1L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.sigma="uniform" the 3rd one is a vector of ',
           'length 2 of lower and upper bounds.')
    }
  } else if (priordist.sigma=='normal' | priordist.sigma=='Gaussian'){
    priordist.code[3] <- 2L
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.sigma="normal" the 3rd one is a vector of ',
           'length 2 of mean and sd.')
    }
  } else if (priordist.sigma=='lognormal'){
    priordist.code[3] <- 3L
    hyperpar[[3]] <- c(2*log(hyperpar[[3]][1]) +
                         - log(hyperpar[[3]][2]^2+hyperpar[[3]][1]^2)/2,
                       sqrt(-2*log(hyperpar[[3]][1]) +
                              + log(hyperpar[[3]][2]^2+hyperpar[[3]][1]^2)))
    # ^ mean and sd on log scale, user supplies mean and sd on exp scale
    if (!is.list(hyperpar)){ # is.null(hyperpar) # NULL is not a list
      stop('hyperpar must be supplied as a list of 3 elements (Linf, K, sigma)',
           ', where for priordist.sigma="lognormal" the 3rd one is a vector of ',
           'length 2 of mean and sd on the original (exponential) scale.')
    }
  } else { # then no prior specified, straight frequentist Fabens, no MCMC
    priordist.code[3] <- 0L
  }
  
  if (any(priordist.code==0L) & !all(priordist.code==0L)){
    warning('At least one prior distribution specified, but not for all three ',
            'Linf, K, and sigma, so reverting to TMB estimation only and no MCMC.')
  }
  
  ### setup
  if (meth!='nlminb'){warning('Only meth="nlminb" supported for now.')}
  # if (!compute.se){warning('se will be computed anyway.')}
  
  # n <- length(L1)
  
  datalist <- list('L1'=L1,'L2'=L2,
                   # 'T1'=T1,'T2'=T2,
                   'deltaT'=deltaT, # as of v0.3: directly supply deltaT
                   'hp_Linf'=hyperpar[[1]], # user-supplied
                   'hp_K'=hyperpar[[2]], # user-supplied
                   'hp_sigma'=hyperpar[[3]], # user-supplied
                   'priordist'=priordist.code) # as of v0.2.1: vector dim 3
  parlist <- list('logLinf'=log(par[1]),
                  'logK'=log(par[2]),
                  'logsigma'=0)
  
  ### TMB estimation (not necessary but good benchmark)
  obj <- MakeADFun(data=datalist,parameters=parlist,
                   DLL="FabensBayesian",silent=T)
  opt <- nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                control=list(eval.max=5000,iter.max=5000))
  theta.tmb <- exp(opt$par[1:2]) # naive estimates without sdreport()
  
  ### MCMC based on tmbstan::tmbstan, by default uses NUTS
  if (!onlyTMB & all(priordist.code!=0)){
    mcmc.obj <- tmbstan(obj=obj,lower=rep(-Inf,3),upper=rep(Inf,3),
                        silent=T,laplace=F,
                        chains=mcmc.control$nchains,
                        warmup=mcmc.control$warmup,
                        iter=mcmc.control$iter,
                        # init='random'
                        init='last.par.best' # start from MLE above
    )
    
    # traceplot(mcmc.obj, pars=names(obj$par), inc_warmup=TRUE) # check conv
    # pairs(mcmc.obj, pars=names(obj$par)) # post dist and scatterplots
    # # ^ Linf and K typically correlate a lot (negatively)
    
    # extract MCMC post draws for derived quantities specified in obj's REPORT
    mcmc.post <- as.matrix(mcmc.obj)
    mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=2) # only Linf and K
    for (i in 1:nrow(mcmc.post)){
      mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[c('Linf','K')])
    }
    # colMeans(mcmc.est) # post means for Linf and K
    
    res <- list('par'=c(mean(mcmc.est[,1]),mean(mcmc.est[,2])),
                'par.TMB'=theta.tmb)
    # ^ MCMC point estimates are posterior means
    names(res$par) <- c('Linf','K')
    names(res$par.TMB) <- c('Linf','K')
  } else {
    res <- list('par'=c(NA,NA),'par.TMB'=theta.tmb)
    names(res$par.TMB) <- c('Linf','K')
  }
  
  ### optional: compute standard errors
  if (compute.se){
    if (!onlyTMB & all(priordist.code!=0)){
      res$se <- c(sqrt(var(mcmc.est[,1])), # /length(mcmc.est[,1])
                  sqrt(var(mcmc.est[,2]))) # /length(mcmc.est[,2])
      # ^ posterior naive se
      names(res$se) <- c('Linf','K')
    } else {
      res$se <- c(NA,NA)
    }
    rep <- sdreport(obj)
    res$se.TMB <- summary(rep)[c('Linf','K'),2] # delta method std errors
    names(res$se.TMB) <- c('Linf','K')
  } else {
    res$se.TMB <- c(NA,NA)
    if (!onlyTMB){
      res$se <- c(NA,NA)
    }
  }
  
  ### output
  if (output.post.draws){
    if (onlyTMB | any(priordist.code==0)){
      warning('No posterior draws if onlyTMB=TRUE or not all priordist set.')
    } else {
      # res$post.draws <- list('Linf'=mcmc.Linf,'K'=mcmc.K)
      # res$post.draws <- as.list(mcmc.post)
      res$post.draws <- list('Linf'=mcmc.est[,1],
                             'K'=mcmc.est[,2])
    }
  }
  
  if (!onlyTMB & all(priordist.code!=0)){
    res$cred.int <- list('Linf'=quantile(mcmc.est[,1],probs=c(0.025,0.975)),
                         'K'=quantile(mcmc.est[,2],probs=c(0.025,0.975)))
    # ^ equal-tailed 95% credible intervals based on MCMC draws
  }
  
  res$value.TMB <- opt$obj
  res$conv.TMB <- opt$mess
  
  return(res)
}




#///////////////////////////////////////////////////////////////////////////////
#### zh09: Zhang, Lessard and Campbell (2009), Bayesian ####
#///////////////////////////////////////////////////////////////////////////////

# Notes:
# * Estimation difficult, need many different random starting values.
# * Model could be overparameterized, 3 random effects seems a bit overkill.
# * Many parameters badly estimated: Linf overestimated, sdLinf underestimated,
#   muA highly overestimated and sdeps overestimated.
# * No noticeable impact of original priors, all fairly flat, but necessary to
#   help numerical stability in estimation (rather than box constraints).

zh09 <- function(par,L1,L2,deltaT,hyperpar=NULL,meth='nlminb',compute.se=T,
                 rand.ini=list(perform=F,'n'=10,'radius'=0.1),
                 onlyTMB=F,output.post.draws=F,lb.sd=NULL,enablepriors=1,
                 mcmc.control=list('nchains'=3,'iter'=5000,'warmup'=4000)){
  # uniform priors with user-supplied hyperparam
  # hyperpar=list(lbubmuinf,lbubK) for consistency across methods, while bounds
  # for sigma are hard-coded below.
  # lb.sd = lower bound for sigmainf, sigmaK and sigmaeps, in TMB optim only
  
  ### setup
  if (is.null(hyperpar)){
    stop('hyperpar must be supplied as list(c(lb,ub),c(lb,ub)) respectively for
         Linf and K.')
  }
  if (!compute.se){warning('se will be computed anyway.')}
  
  length.theta <- 7
  n <- length(L1)
  rad <- rand.ini$radius
  datalist <- list('L1'=L1,'L2'=L2,
                   # 'T1'=T1,'T2'=T2,
                   'deltaT'=deltaT, # as of v0.3: directly supply deltaT
                   'lbubmuinf'=hyperpar[[1]], # user-supplied
                   'lbubK'=hyperpar[[2]], # user-supplied
                   'lbubgammashapeA'=c(1e-3,1000),
                   'lbubgammarateA'=c(1e-3,1000),
                   'lbubsigmainf'=c(1e-10,1e4),
                   'lbubsigmaK'=c(1e-10,1e4),
                   'lbubsigmaeps'=c(1e-10,1e4),
                   'enablepriors'=enablepriors)
  
  ### TMB optimization and sdreport
  if (meth!='nlminb'){stop('Only meth="nlminb" is allowed for now.')}
  
  if (is.null(lb.sd)){
    lb.par <- rep(-Inf,length.theta) # no bounds
  } else {
    lb.par <- c(-Inf,log(lb.sd),-Inf,log(lb.sd),-Inf,-Inf,log(lb.sd))
    # ^ helps if data too close to pure vB, otherwise Hessian not pos def
  }
  
  if (rand.ini$perform){ # try multiple sets of random initial values
    it.ini <- 1
    obj <- NULL
    opt <- NULL
    while (it.ini <= rand.ini$n){
      parlist <- list('logmuinf'=log(par[1])+runif(1,-rad,rad),
                      'logsigmainf'=runif(1,-rad,rad),
                      'logmuK'=log(par[2])+runif(1,-rad,rad),
                      'logsigmaK'=runif(1,-rad,rad),
                      'logmuA'=1+runif(1,-rad,rad),
                      'logsigmaA'=runif(1,-rad,rad),
                      'logsigmaeps'=runif(1,-rad,rad),
                      'logLinf'=rep(log(par[1]),n)+runif(n,-rad,rad),
                      'logK'=rep(log(par[2]),n)+runif(n,-rad,rad),
                      'logA'=rep(1,n)+runif(n,-rad,rad)) # jittering
      obj0 <- MakeADFun(data=datalist,parameters=parlist,
                        random=c('logLinf','logK','logA'),
                        DLL="Zhang",silent=T)
      opt0 <- try(nlminb(start=obj0$par,obj=obj0$fn,gr=obj0$gr,
                         control=list(eval.max=5000,iter.max=5000),
                         lower=lb.par),T)
      if (!is(opt0,'try-error')){
        if (opt0$conv%in%c(0,1)){ # opt0$conv==0 too strict
          if (is.null(obj)){ # then no good conv so far, keep this one
            obj <- obj0
            opt <- opt0
          } else if (opt0$obj < opt$obj){ # then compare to best conv so far
            obj <- obj0
            opt <- opt0
          }
        }
      }
      it.ini <- it.ini+1
    }
    if (is.null(obj)){
      stop('Fitting did not converge, try more random initial values.')
    } else {
      rep <- sdreport(obj)
      summary.rep <- summary(rep)
      theta.tmb <- summary.rep[c('muinf','K'),1]
      se.theta.tmb <- summary.rep[c('muinf','K'),2]
      val.tmb <- opt$obj # obj$fn()
    }
  } else { # single set of non-random initial values
    parlist <- list('logmuinf'=log(par[1]),
                    'logsigmainf'=0,
                    'logmuK'=log(par[2]),
                    'logsigmaK'=0,
                    'logmuA'=1,
                    'logsigmaA'=0,
                    'logsigmaeps'=0,
                    'logLinf'=rep(log(par[1]),n),
                    'logK'=rep(log(par[2]),n),
                    'logA'=rep(1,n)) # no jittering
    obj <- MakeADFun(data=datalist,parameters=parlist,
                     random=c('logLinf','logK','logA'),
                     DLL="Zhang",silent=T)
    opt <- try(nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                      control=list(eval.max=5000,iter.max=5000),
                      lower=lb.par),T)
    if (is(opt,'try-error')){
      stop('Optimization failed, try random initial values.')
    }
    rep <- sdreport(obj)
    summary.rep <- summary(rep)
    theta.tmb <- summary.rep[c('muinf','K'),1]
    se.theta.tmb <- summary.rep[c('muinf','K'),2]
    # val.tmb <- opt$obj # obj$fn()
    # mess.tmb <- opt$mess
  }
  
  ### MCMC based on tmbstan::tmbstan, by default uses NUTS
  if (!onlyTMB){
    mcmc.obj <- tmbstan(obj=obj,
                        lower=rep(-Inf,length(unlist(parlist))),
                        upper=rep(Inf,length(unlist(parlist))),
                        silent=T,laplace=F,
                        chains=mcmc.control$nchains,
                        warmup=mcmc.control$warmup,
                        iter=mcmc.control$iter,
                        # init='random'
                        init='last.par.best' # start from MLE above
    )
    
    # traceplot(mcmc.obj, pars=c('logmuinf','logmuK'), inc_warmup=TRUE) # check conv
    # pairs(mcmc.obj, pars=c('logmuinf','logmuK')) # post dist and scatterplots
    
    # extract MCMC post draws for derived quantities specified in obj's REPORT
    mcmc.post <- as.matrix(mcmc.obj)
    mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=2) # only Linf and K
    for (i in 1:nrow(mcmc.post)){
      mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[c('muinf','muK')])
    }
    # colMeans(mcmc.est) # post means for (expectation of) Linf and K

    res <- list('par'=c(mean(mcmc.est[,1]),mean(mcmc.est[,2])),'par.TMB'=theta.tmb,
                'se'=c(sqrt(var(mcmc.est[,1])),  # /length(mcmc.est[,1])
                       sqrt(var(mcmc.est[,2]))), # /length(mcmc.est[,2])
                'se.TMB'=se.theta.tmb)
    names(res$par) <- c('Linf','K')
    names(res$par.TMB) <- c('Linf','K')
    names(res$se) <- c('Linf','K')
    names(res$se.TMB) <- c('Linf','K')
  } else {
    res <- list('par'=c(NA,NA),'par.TMB'=theta.tmb,
                'se'=c(NA,NA),'se.TMB'=se.theta.tmb)
    names(res$par.TMB) <- c('Linf','K')
    names(res$se.TMB) <- c('Linf','K')
  }
  
  ### output
  if (output.post.draws){
    if (onlyTMB){
      warning('Cannot output posterior draws if onlyTMB=TRUE.')
    } else {
      res$post.draws <- list('Linf'=mcmc.est[,1],'K'=mcmc.est[,2])
      # res$post.draws <- as.list(mcmc.post)
    }
  }
  if (!onlyTMB){ # equal-tailed 95% credible intervals based on MCMC draws
    res$cred.int <- list('Linf'=quantile(mcmc.est[,1],probs=c(0.025,0.975)),
                         'K'=quantile(mcmc.est[,2],probs=c(0.025,0.975)))
  }
  
  res$value.TMB <- opt$obj
  res$conv.TMB <- opt$mess
  
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### la02: Laslett, Eveson and Polacheck (2002) ####
#///////////////////////////////////////////////////////////////////////////////

# Notes:
# * Mean of age at capture likely overestimated, with sd underestimated.
# * Already with n=100 all estimates seem symmetric.
# * Naive CIs on original scale cover too much.

la02 <- function(par,L1,L2,deltaT,hyperpar=NULL,meth='nlminb',
                 compute.se=T,enable.priors=F,lb.sd=NULL){
  # uniform priors with user-supplied hyperparam
  # hyperpar=list(lbubmuinf,lbubK) for consistency across methods, while bounds
  # for lbubsigmainf, lbubmuA, lbubsigmaA and lbubsigmaeps are hard-coded below.
  # lb.sd = lower bound for sigmainf and sigmaeps, in TMB optimization only
  
  ### setup
  if (is.null(hyperpar)){
    if (enable.priors){ # then hyperpar must be supplied
      stop('hyperpar must be supplied as list(c(lb,ub),c(lb,ub)) respectively for
         the mean of Linf and for K.')
    } else { # then set arbitrary values for hyperpar, ignored in TMB
      hyperpar <- list(c(0,1),c(0,1))
    }
  }
  if (!compute.se){warning('se will be computed anyway.')}
  
  n <- length(L1)
  parlist <- list('logmuinf'=log(par[1]),'logsigmainf'=0,
                  'logK'=log(par[2]),
                  'logmuA'=1,'logsigmaA'=0,
                  'logsigmaeps'=0,
                  'logLinf'=rep(log(par[1]),n),
                  'logA'=rep(1,n))
  length.theta <- 6
  datalist <- list('L1'=L1,'L2'=L2,
                   # 'T1'=T1,'T2'=T2,
                   'deltaT'=deltaT, # as of v0.3: directly supply deltaT
                   'lbubmuinf'=hyperpar[[1]], # user-supplied
                   'lbubK'=hyperpar[[2]], # user-supplied
                   'lbubsigmainf'=c(1e-10,1e4),
                   'lbubmuA'=c(1e-10,1000),
                   'lbubsigmaA'=c(1e-10,1e4),
                   'lbubsigmaeps'=c(1e-10,1e4),
                   'enablepriors'=as.integer(enable.priors))
  obj <- MakeADFun(data=datalist,parameters=parlist,
                   random=c('logLinf','logA'),
                   DLL="Laslett",silent=T)
  
  ### optimization and sdreport
  if (meth=='nlminb'){
    if (is.null(lb.sd)){
      lb.par <- rep(-Inf,length.theta) # no bounds
    } else {
      lb.par <- c(-Inf,log(lb.sd),-Inf,-Inf,-Inf,log(lb.sd))
      # ^ helps if data too close to pure vB, otherwise Hessian not pos def
    }
    opt <- nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                  control=list(eval.max=5000,iter.max=5000),lower=lb.par)
    rep <- sdreport(obj)
    summary.rep <- summary(rep)
  } else {stop('Only meth="nlminb" is allowed for now.')}
  
  ### output
  theta <- summary.rep[(length.theta+2*n+1):(2*length.theta+2*n),1]
  se.theta <- summary.rep[(length.theta+2*n+1):(2*length.theta+2*n),2]
  res <- list('par'=theta[c(1,3)],'se'=se.theta[c(1,3)],
              'other.par'=theta[-c(1,3)],'other.par.se'=se.theta[-c(1,3)],
              'random.effects'=list(
                'Linf'=summary.rep[dimnames(summary.rep)[[1]]=='Linf',1],
                'se.Linf'=summary.rep[dimnames(summary.rep)[[1]]=='Linf',2],
                'A'=summary.rep[dimnames(summary.rep)[[1]]=='A',1],
                'se.A'=summary.rep[dimnames(summary.rep)[[1]]=='A',2]),
              'value'=obj$fn(),'conv'=opt$mess)
  names(res$par) <- c('Linf','K')
  names(res$se) <- c('Linf','K')
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### ja91: James (1991) ####
#///////////////////////////////////////////////////////////////////////////////

# Notes:
# * This is the WLS from Section 2.2, presented as improved version of fa65.
# * Original specification of Y_{1i} and Y_{2i} with indep epsi_1 and epsi_2
#   leads to eta_i = epsi_{2i} - epsi_{1i}*exp(-K*deltaT), i.e. does not
#   technically depend on Linf. So corresponding Var[eta_i] = sigma2*wi may not
#   match a DGP where the original specification is satisfied.
# * Consequence: wasn't able to find correct DGP to verify hand deriv with
#   expectations.
# * Consequence: computed s.e. may corespond to nothing. They generally
#   overestimate the empirical s.e. under best guess of DGP (which is: generate
#   Y1 and Y2 separately with iid additive Gaussian mean 0 and var sigma2).
# * With large sigma2 (>0.4) and large K (>0.3): estimates of K highly
#   right-skewed and overestimated s.e. of Linf but dramatically underestimated
#   s.e. of K, even with large n. 

ja91 <- function(par,L1,L2,deltaT,meth='nlminb',compute.se=T){
  ### objective function
  rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    expKdt <- exp(-K*deltaT)
    RSS <- sum((L2-L1-(Linf-L1)*(1-expKdt))^2/(1+expKdt^2))
    return(RSS) # weighted version of fa65
  }
  gr.rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    expKdt <- exp(-K*deltaT)
    exp2Kdt <- expKdt^2 #  exp(-2*K*deltaT)
    wi <- 1+exp2Kdt # weight in James' WLS
    resid <- L2-L1-(Linf-L1)*(1-expKdt) # eta_i in notation of James (1991)
    psi1 <- sum(resid*(1-expKdt)/wi)
    psi2 <- sum(resid*(Linf-L1)*deltaT*expKdt/wi-resid^2*exp2Kdt*deltaT/wi^2)
    return(-2*c(psi1,psi2))
  }
  ### setup
  # deltaT <- as.numeric(T2-T1) # asof v0.3: directly supply deltaT
  ### optimization
  if (meth=='nlminb'){
    opt <- nlminb(start=par,objective=rss,gradient=gr.rss,
                  L1=L1,L2=L2,deltaT=deltaT)
    res <- list('par'=opt$par,'value'=opt$obj,'conv'=opt$mess)
  } else {
    opt <- optim(par=par,fn=rss,gr=gr.rss,
                 L1=L1,L2=L2,deltaT=deltaT,method=meth)
    res <- list('par'=opt$par,'value'=opt$val,'conv'=opt$mess)
  }
  ### optional: compute standard errors
  if (compute.se){
    # setup
    n <- length(L1)
    Linf <- opt$par[1] # plug-in estimated values
    K <- opt$par[2] # plug-in estimated values
    expKdt <- exp(-K*deltaT)
    exp2Kdt <- expKdt^2 # exp(-2*K*deltaT)
    wi <- 1+exp2Kdt # weight in James' WLS
    # estimate residual variance under James' model (corrected WLS)
    y <- L2-L1 # deltaL
    sigma2 <- sum((y-(Linf-L1)*(1-exp(-K*deltaT)))^2/wi)/(n-2) # bottom p. 1521
    # compute M and Q matrices (need full sandwich because M22!=Q22)
    M11 <- sum((1-expKdt)^2/wi)
    M12 <- sum((Linf-L1)*deltaT*expKdt*(1-expKdt)/wi)
    M22 <- sum(deltaT^2*exp2Kdt*((Linf-L1)^2/wi+2*sigma2*(1-exp2Kdt)/wi^2))
    # Mmat <- -2*cbind(c(M11,M12),c(M12,M22))#/n # not required
    Minv <- cbind(c(M22,-M12),c(-M12,M11))/(M11*M22-M12^2)#/(-2)*n
    # Q11 <- sum((1-expKdt)^2/wi) # identical to M11
    # Q12 <- sum((Linf-L1)*deltaT*expKdt*(1-expKdt)/wi) # identical to M12
    Q22 <- sum(deltaT^2*exp2Kdt*((Linf-L1)^2/wi+3*sigma2*exp2Kdt/wi^2))
    Qmat <- sigma2*cbind(c(M11,M12),c(M12,Q22))#*4/n
    # compute var-cov sandwich matrix (based on expectations)
    varcov <- Minv%*%Qmat%*%Minv#/n
    res$se <- sqrt(diag(varcov))
  } else {
    res$se <- c(NA,NA)
  }
  ### output
  names(res$par) <- c('Linf','K')
  names(res$se) <- c('Linf','K')
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### fr88: Francis (1988) ####
#///////////////////////////////////////////////////////////////////////////////

# Notes:
# * Wasn't able to find a DGP that leads to unbiased estimates. Both Linf and K
#   are biased, even asymptotically: Linf overestimated and K underestimated.
# * Both estimates look normally distributed empirically, already with n=500.
# * Sampling variance of estimates decreases with nu decreasing, but not bias.
# * Problem of varcovar function not using asymptotic sandwich, but relying on
#   numerical Hessian and delta method for (Linf,K). Since estimates not
#   consistent, delta method is probably returning garbage.

fr88 <- function(par,L1,L2,deltaT,meth='nlminb',compute.se=T,tol.sd=1e-4){
  # par argument not necessary, but for consistency with other methods
  ### objective function
  nll <- function(par,alpha,beta,deltaT,L1,L2,tol.sd){ # par=c(ga,gb,nu)
    ga <- par[1]
    gb <- par[2]
    nu <- par[3] # sd = nu*mu, Francis' eq. (5)
    deltaL <- L2-L1
    mu <- ((beta*ga-alpha*gb)/(ga-gb)-L1)*(1-(1+(ga-gb)/(alpha-beta))^deltaT)
    std.dev <- nu*ifelse(mu>0, mu, tol.sd) # Francis' eq. (5)
    pllvec <- dnorm(deltaL, mean=mu, sd=std.dev, log=T)
    return(-sum(pllvec)) # negloglik
  }
  ### setup
  alpha <- as.numeric(quantile(L1,0.1))
  beta <- as.numeric(quantile(L1,0.9))
  deltaL <- as.numeric(L2-L1)
  # deltaT <- as.numeric(T2-T1) # as of v0.3: directly supply deltaT
  ### starting values, from fishmethods::grotag
  subs.a <- which(L1>(0.8*alpha) & L1<(1.2*alpha) & deltaT>median(deltaT))
  if (length(subs.a)==0){ # empty subset
    subs.a <- which(L1>(0.7*alpha) & L1<(1.3*alpha) & deltaT>median(deltaT))
  } # ^ ad hoc but should work on most cases, NBG as long as non-empty subset
  subs.b <- which(L1>(0.8*beta) & L1<(1.2*beta) & deltaT>median(deltaT))
  if (length(subs.b)==0){ # empty subset
    subs.b <- which(L1>(0.5*beta) & L1<(1.3*beta) & deltaT>median(deltaT))
  } # ^ ad hoc but should work on most cases, NBG as long as non-empty subset
  ga.ini <- mean(deltaL[subs.a]/deltaT[subs.a])
  gb.ini <- mean(deltaL[subs.b]/deltaT[subs.b])
  par.ini <- c(ga.ini,gb.ini,1) # hardcoded starting value for nu
  ### lower and upper bounds, from fishmethods::grotag
  lb <- c(0.5*ga.ini, 0.5*gb.ini, 1e-6) # hardcoded lb for nu
  ub <- c(1.5*ga.ini, 1.5*gb.ini, 100) # hardcoded ub for nu
  ### optimization
  if (meth=='nlminb'){
    opt <- nlminb(start=par.ini,objective=nll,lower=lb,upper=ub,
                  alpha=alpha,beta=beta,L1=L1,L2=L2,deltaT=deltaT,tol.sd=tol.sd)
    opt$value <- opt$obj
  } else if (meth=='L-BFGS-B'){
    opt <- optim(par=par.ini,fn=nll,lower=lb,upper=ub,
                 method='L-BFGS-B',hessian=F, # only compute hess if compute.se=T
                 alpha=alpha,beta=beta,L1=L1,L2=L2,deltaT=deltaT,tol.sd=tol.sd)
  } else {stop('Only "nlminb" and "L-BFGS-B" are allowed for meth.')}
  ### optional: compute standard errors
  if (compute.se){
    # suppressWarnings(hess <- try(optimHess(par=opt$par,fn=nll,
    #                                        alpha=alpha,beta=beta,L1=L1,L2=L2,
    #                                        deltaT=deltaT,tol.sd=tol.sd),T))
    # if (is(hess,'try-error')){ # likely failed because of nu hitting lb
    hess <- optimHess(par=opt$par[1:2],
                      fn=function(ga.gb,nu,alpha,beta,deltaT,L1,L2,tol.sd){
                        nll(c(ga.gb,nu),alpha,beta,deltaT,L1,L2,tol.sd)},
                      nu=opt$par[3],alpha=alpha,beta=beta,L1=L1,L2=L2,
                      deltaT=deltaT,tol.sd=tol.sd)
    # } # ^ Hessian wrt ga and gb only, nu held fixed, more stable
    ga <- opt$par[1]
    gb <- opt$par[2]
    jac.g <- cbind(c(gb*(alpha-beta)/(ga-gb)^2,   # jacobian of transformation
                     -1/(ga-gb+alpha-beta)),      # from (ga,gb) to (Linf,K) 
                   c(-ga*(alpha-beta)/(ga-gb)^2,
                     1/(ga-gb+alpha-beta)))
    # lazy: numerical Hessian instead of analytical Fisher info of (ga,gb,nu)
    # (same laziness as in fishmethods::grotag)
    varcov <- t(jac.g)%*%solve(hess[1:2,1:2])%*%jac.g # delta method
    # ^ nu ignored in varcov, because Hessian rarely pos def with it
    se <- sqrt(diag(varcov))
  } else {
    hess <- NA
    se <- c(NA,NA)
  }
  ### output
  Linf <- (beta*opt$par[1]-alpha*opt$par[2])/(opt$par[1]-opt$par[2])
  K <- -log(1+(opt$par[1]-opt$par[2])/(alpha-beta))
  names(opt$par) <- c('g.alpha','g.beta','nu')
  res <- list('par'=c(Linf,K),'value'=opt$value,'conv'=opt$mess,'se'=se,
              'original.par'=opt$par,'hessian'=hess,
              'alpha'=alpha,'beta'=beta)
  names(res$par) <- c('Linf','K')
  names(res$se) <- c('Linf','K')
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### fa65: Fabens (1965) ####
#///////////////////////////////////////////////////////////////////////////////

# Notes:
# * Empirical variance of estimator converges to sandwich with expected values
# * Linf and K seem to converge at about same speed (bias and as. norm.)
# * n=500 is already good enough for getting variances in right ballpark.

fa65 <- function(par,L1,L2,deltaT,meth='nlminb',compute.se=T){
  ### objective function and gradient wrt (Linf,K)
  rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    dL.fitted <- (Linf-L1)*(1-exp(-K*deltaT))
    dL <- L2-L1
    RSS <- sum((dL-dL.fitted)^2)
    return(RSS)
  }
  gr.rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    dL.fitted <- (Linf-L1)*(1-exp(-K*deltaT))
    dL <- L2-L1
    resid <- dL-dL.fitted
    expKdt <- exp(-K*deltaT)
    gr <- -2*c(sum(resid*(1-expKdt)),sum(resid*(Linf-L1)*deltaT*expKdt))
    return(gr)
  }
  ### setup
  # deltaT <- as.numeric(T2-T1) # as of v0.3: directly supply deltaT
  ### optimization
  if (meth=='nlminb'){
    opt <- nlminb(start=par,objective=rss,gradient=gr.rss,L1=L1,L2=L2,deltaT=deltaT)
    res <- list('par'=opt$par,'value'=opt$obj,'conv'=opt$mess)
  } else {
    opt <- optim(par=par,fn=rss,gr=gr.rss,L1=L1,L2=L2,deltaT=deltaT,method=meth)
    res <- list('par'=opt$par,'value'=opt$val,'conv'=opt$mess)
  }
  ### optional: compute standard errors
  if (compute.se){
    # setup
    n <- length(L1)
    Linf <- opt$par[1] # plug-in estimated values
    K <- opt$par[2] # plug-in estimated values
    # estimate residual variance under Fabens' model (implicit Gaussian error)
    y <- L2-L1 # deltaL
    sigma2 <- sum((y-(Linf-L1)*(1-exp(-K*deltaT)))^2)/(n-2) # residual sum of squares
    # compute M=Q
    expKdt <- exp(-K*deltaT)
    M11 <- sum((1-expKdt)^2)
    M12 <- sum(deltaT*expKdt*(Linf-L1)*(1-expKdt))
    M22 <- sum((deltaT*expKdt*(Linf-L1))^2)
    Minv <- cbind(c(M22,-M12),c(-M12,M11))/(M11*M22-M12^2)
    # compute var-cov sandwich matrix (based on expectations)
    varcov <- sigma2*Minv # (2x2) variance-covariance matrix
    res$se <- sqrt(diag(varcov)) # standard errors
  } else {
    res$se <- c(NA,NA)
  }
  ### output
  names(res$par) <- c('Linf','K')
  names(res$se) <- c('Linf','K')
  return(res)
}


#///////////////////////////////////////////////////////////////////////////////
#### gh59: Gulland and Holt (1959) ####
#///////////////////////////////////////////////////////////////////////////////

# Notes:
# * Emp var of estimator converges to sandwich with expected values.
# * Emp var of K converges much more quickly than Linf, maybe because
#   of bias that vanishes only with large sample sizes (>5000).
# * Consequence: unless n is really large, sandwich underestimates greatly emp
#   variance of Linf (let alone MSE with non-zero bias)
# * Related: estimated Linf is quite positively skewed unless n very large.

gh59 <- function(par,L1,L2,deltaT,meth='nlminb',compute.se=T){
  ### objective function and gradient wrt (Linf,K)
  rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    Lmean <- (L2+L1)/2
    dLdt.fitted <- K*(Linf-Lmean)
    dLdt <- (L2-L1)/deltaT
    RSS <- sum((dLdt-dLdt.fitted)^2)
    return(RSS)
  }
  gr.rss <- function(par,L1,L2,deltaT){
    Linf <- par[1]
    K <- par[2]
    Lmean <- (L2+L1)/2
    dLdt.fitted <- K*(Linf-Lmean)
    dLdt <- (L2-L1)/deltaT
    resid <- dLdt-dLdt.fitted
    gr <- -2*c(sum(resid*K),sum(resid*(Linf-Lmean)))
    return(gr)
  }
  ### setup
  # deltaT <- as.numeric(T2-T1) # as of v0.3: directly supply deltaT
  ### optimization
  if (meth=='nlminb'){
    opt <- nlminb(start=par,objective=rss,gradient=gr.rss,L1=L1,L2=L2,deltaT=deltaT)
    res <- list('par'=opt$par,'value'=opt$obj,'conv'=opt$mess)
  } else {
    opt <- optim(par=par,fn=rss,gr=gr.rss,L1=L1,L2=L2,deltaT=deltaT,method=meth)
    res <- list('par'=opt$par,'value'=opt$val,'conv'=opt$mess)
  }
  ### optional: compute standard errors
  if (compute.se){
    # setup
    n <- length(L1)
    Linf <- opt$par[1] # plug-in estimated values
    K <- opt$par[2] # plug-in estimated values
    # estimate residual variance under G&H's linear model
    # y <- (L2-L1)/(T2-T1) # response in G&H
    y <- (L2-L1)/deltaT # response in G&H
    x <- (L2+L1)/2 # fixed covariate in G&H
    sigma2 <- sum((y-K*Linf+K*x)^2)/(n-2) # residual sum of squares
    # compute M=Q
    M11 <- K^2
    M12 <- K*sum(Linf-x)/n
    M22 <- sum((Linf-x)^2)/n
    # compute var-cov sandwich matrix (based on expectations)
    Minv <- cbind(c(M22,-M12),c(-M12,M11))/(M11*M22-M12^2)
    varcov <- sigma2*Minv/n # (2x2) variance-covariance matrix
    res$se <- sqrt(diag(varcov)) # standard errors
  } else {
    res$se <- c(NA,NA)
  }
  ### output
  names(res$par) <- c('Linf','K')
  names(res$se) <- c('Linf','K')
  return(res)
}

# END GrowthEstimation_Methods
