#///////////////////////////////////////////////////////////////////////////////
#### GrowthEstimation v0.4.5 ####
#///////////////////////////////////////////////////////////////////////////////


#///////////////////////////////////////////////////////////////////////////////
#### hp.normal and hp.lognormal: prior hyperparameters for Linf and K ####
#///////////////////////////////////////////////////////////////////////////////

hp.normal <- function(median, upper.bound, plot=T){
  sd <- (upper.bound-median)/qnorm(0.99)
  # ub=0.99 quantile, mean=median for Gaussian dist
  
  if (plot){
    xgrid <- seq(qnorm(1e-3,mean=median,sd=sd),
                 qnorm(1-1e-3,mean=median,sd=sd),length.out=200)
    plot(xgrid,dnorm(xgrid,mean=median,sd=sd),type='l',
         xlab='',ylab='Normal pdf')
    abline(v=median,col='red',lty=2)
    text(x=median,y=0,labels=paste0('median = ',median),
         pos=4,col='red',cex=0.6)
    abline(v=qnorm(c(0.99),mean=median,sd=sd),col='red',lty=2)
    text(x=upper.bound,y=1/sqrt(2*pi)/sd,
         labels=paste0('upper bound = 0.99 quantile\n= ',upper.bound),
         pos=2,col='red',cex=0.6)
    abline(v=qnorm(0.01,mean=median,sd=sd),lty=2)
    text(x=qnorm(0.01,mean=median,sd=sd),y=1/sqrt(2*pi)/sd,
         labels=paste0('0.01 quantile\n~= ',
                       round(qnorm(0.01,mean=median,sd=sd),2)),
         pos=4,cex=0.6)
  }
  
  return(c('mean'=median,'sd'=sd))
}

hp.lognormal <- function(median, upper.bound, plot=T, interval.sdlog=c(1e-5,5)){
  meanlog <- log(median) # lognormal median = exp(meanlog)
  
  diffq <- function(sdlog,q,meanlog,prob=0.99){
    qlnorm(p=prob,meanlog=meanlog,sdlog=sdlog)-q
  }
  
  sdlog <- uniroot(f=diffq,interval=interval.sdlog,
                   q=upper.bound,meanlog=meanlog)$root
  
  if (plot){
    xgrid <- seq(qlnorm(1e-3,meanlog=meanlog,sdlog=sdlog),
                 qlnorm(1-1e-3,meanlog=meanlog,sdlog=sdlog),length.out=200)
    plot(xgrid,dlnorm(xgrid,meanlog=meanlog,sdlog=sdlog),type='l',
         xlab='',ylab='Lognormal pdf')
    abline(v=median,col='red',lty=2)
    text(x=median,y=0,labels=paste0('median = ',median),
         pos=4,col='red',cex=0.6)
    abline(v=qlnorm(0.99,meanlog=meanlog,sdlog=sdlog),col='red',lty=2)
    text(x=upper.bound,y=1/sqrt(2*pi)/sdlog*exp(sdlog^2/2 - meanlog),
         labels=paste0('upper bound = 0.99 quantile\n= ',upper.bound),
         pos=2,col='red',cex=0.6)
    abline(v=qlnorm(0.01,meanlog=meanlog,sdlog=sdlog),lty=2)
    text(x=qlnorm(0.01,meanlog=meanlog,sdlog=sdlog),
         y=1/sqrt(2*pi)/sdlog*exp(sdlog^2/2 - meanlog),
         labels=paste0('0.01 quantile\n~= ',
                       round(qlnorm(0.01,meanlog=meanlog,sdlog=sdlog),2)),
         pos=4,cex=0.6)
  }
  
  return(c('mean'=exp(meanlog+sdlog^2/2), # expectation ori scale
           'sd'=sqrt((exp(sdlog^2)-1)*exp(2*meanlog+sdlog^2)))) # sd ori scale
}



#///////////////////////////////////////////////////////////////////////////////
#### GrowthPriors: prior hyperparameters for Linf and K from fishbase ####
#///////////////////////////////////////////////////////////////////////////////


GrowthPriors <- function(Lmax=133, species="Mustelus asterias",
                         category="BodyShape", LowP=0.8, UpP=1.2,
                         lnorm.coef=c(1.2,3), LQ=0, UQ=1){
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
  #' @param lnorm.coef vector of 2 positive scalars defining narrow and wide
  #' lognormal priors for Linf; first element * median is set as 0.99 quantile
  #' to compute the narrow lognormal sdlog hyperparameter, while the second
  #' element * median set as 0.99 quantile defines the wide lognormal.
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
  
  ### uniform prior bounds for Linf from Lmax 
  # Linf <- 10^(0.044+0.9841*log10(Lmax)) # Froese and Binohlan (2000)
  Linf <- Lmax/0.99 # simpler assumption, anyway close to Froese and Binohlan (2000)
  
  # calculate prior boundaries for Linf (uniform prior)
  # take priors for Linf as X% of Lmax
  # LowLinf <- LowP*Lmax # uniform prior lower bound
  # UpLinf <- UpP*Lmax # uniform prior upper bound
  LowLinf <- LowP*Linf # uniform prior lower bound
  UpLinf <- UpP*Linf # uniform prior upper bound
  # ^ as of v0.3: center on Linf rather than Lmax
  
  
  ### Gaussian prior mean and sd for Linf
  meanLinf <- Linf # (LowLinf+UpLinf)/2 # center of uniform prior
  sdLinf <- uniroot(f=function(sdLinf){meanLinf+3*sdLinf-UpLinf},
                    interval=c(0,1e5))$root
  # ^ find sd such that uniform upper bound matches mean+3*sd. So prob under
  #   Gaussian within uniform bounds is pnorm(3)-pnorm(-3) = 99.7% roughly
  
  
  ### lognormal mean and sd (original scale) for Linf
  # only requires Lmax and coef 1.2 (narrow) and 3 (wide) for finding sdlog
  diffq <- function(sdlog,q,Lmax,prob=0.99){
    meanlog <- log(Lmax/0.99) # median = Lmax/0.99
    return(qlnorm(p=prob,meanlog=meanlog,sdlog=sdlog)-q)
  }
  meanlog1 <- log(Lmax/0.99)
  sdlog1 <- uniroot(f=diffq, interval=c(1e-5,3),
                    q=Lmax/0.99*lnorm.coef[1], Lmax=Lmax)$root
  meanlog2 <- log(Lmax/0.99) # same median as narrow
  sdlog2 <- uniroot(f=diffq, interval=c(1e-5,3),
                    q=Lmax/0.99*lnorm.coef[2], Lmax=Lmax)$root
  
  logn.mean1 <- exp(meanlog1+sdlog1^2/2) # expectation on original scale
  # logn.sd1 <- sqrt(exp(sdlog1^2-1)*exp(2*meanlog1+sdlog1^2)) # sd original scale <= WRONG
  logn.sd1 <- sqrt((exp(sdlog1^2)-1)*exp(2*meanlog1+sdlog1^2)) # sd original scale <= CORRECT

  logn.mean2 <- exp(meanlog2+sdlog2^2/2) # expectation on original scale
  # logn.sd2 <- sqrt(exp(sdlog2^2-1)*exp(2*meanlog2+sdlog2^2)) # sd original scale <= WRONG
  logn.sd2 <- sqrt((exp(sdlog2^2)-1)*exp(2*meanlog2+sdlog2^2)) # sd original scale <= CORRECT
  
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
                   'Linf'=Linf,       # for uniform prior on Linf
                   'LowLinf'=LowLinf, # for uniform prior on Linf
                   'UpLinf'=UpLinf,   # for uniform prior on Linf
                   'meanLinf'=meanLinf, # for Gaussian prior on Linf
                   'sdLinf'=sdLinf,     # for Gaussian prior on Linf
                   'lnorm.meanLinf.narrow'=logn.mean1, # narrow lognormal Linf
                   'lnorm.sdLinf.narrow'=logn.sd1,     # narrow lognormal Linf
                   'lnorm.meanLinf.wide'=logn.mean2, # wide lognormal Linf
                   'lnorm.sdLinf.wide'=logn.sd2,     # wide lognormal Linf
                   'K'=K,       # uniform prior on K
                   'LowK'=LowK, # uniform prior on K
                   'UpK'=UpK)   # uniform prior on K
  
  return(df)
}







#///////////////////////////////////////////////////////////////////////////////
#### Bfa65: Bayesian Fabens (1965) ####
#///////////////////////////////////////////////////////////////////////////////

Bfa65 <- function(par,L1,L2,deltaT,
                  priordist.Linf='lognormal',
                  priordist.K='uniform',
                  priordist.sigma='uniform',
                  hyperpar=NULL,
                  meth='nlminb',compute.se=T,
                  onlyTMB=F,output.post.draws=F,
                  mcmc.control=list('nchains'=5,'iter'=20000,'warmup'=10000,
                                    'adapt_delta'=0.8,'max_treedepth'=20)){
  
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
                        control=list(adapt_delta=mcmc.control$adapt_delta,
                                     max_treedepth=mcmc.control$max_treedepth
                                     # adapt_delta=0.99, # def=0.8, larger=safer
                                     # max_treedepth=15 # def=10, larger helps
                        ),
                        # init='random'
                        init='last.par.best' # start from MLE above
    )
    
    # traceplot(mcmc.obj, pars=names(obj$par), inc_warmup=TRUE) # check conv
    # pairs(mcmc.obj, pars=names(obj$par)) # post dist and scatterplots
    # # ^ Linf and K typically correlate a lot (negatively)
    
    # extract MCMC post draws for derived quantities specified in obj's REPORT
    mcmc.post <- as.matrix(mcmc.obj)
    mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=3) # Linf, K, sigma
    for (i in 1:nrow(mcmc.post)){
      mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
        c('Linf','K','sigma')]) # need all param for DIC
    }
    # colMeans(mcmc.est[,1:2]) # post means for Linf and K
    
    res <- list('par'=c(median(mcmc.est[,1]),median(mcmc.est[,2])),
                # 'par'=c(mean(mcmc.est[,1]),mean(mcmc.est[,2])),
                # ^ post median better, dist can be very skewed with small n
                'par.TMB'=theta.tmb)
    # ^ MCMC point estimates are posterior median as of v0.3.1
    names(res$par) <- c('Linf','K')
    names(res$par.TMB) <- c('Linf','K')
  } else {
    res <- list('par'=c(NA,NA),'par.TMB'=theta.tmb)
    names(res$par.TMB) <- c('Linf','K')
  }
  
  ### optional: compute standard errors
  if (compute.se){
    if (!onlyTMB & all(priordist.code!=0)){
      # res$se <- c(sqrt(var(mcmc.est[,1])),sqrt(var(mcmc.est[,2])))
      # ^ posterior naive se numerically unstable if low n
      res$se <- c(median(abs(mcmc.est[,1]-res$par[1])),
                  median(abs(mcmc.est[,2]-res$par[2]))) # MADAM
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
      res$post.draws <- list('Linf'=mcmc.est[,1],
                             'K'=mcmc.est[,2],
                             'sigma'=mcmc.est[,3])
    }
  }
  
  if (!onlyTMB & all(priordist.code!=0)){
    res$cred.int <- list('Linf'=quantile(mcmc.est[,1],probs=c(0.025,0.975)),
                         'K'=quantile(mcmc.est[,2],probs=c(0.025,0.975)),
                         'sigma'=quantile(mcmc.est[,3],probs=c(0.025,0.975))
    )
    # ^ equal-tailed 95% credible intervals based on MCMC draws
    
    # DIC, both p_D and p_V versions
    loglikvec <- function(par,L1,L2,deltaT,log=T){
      Linf <- par[1]
      K <- par[2]
      sigma <- par[3]
      meanL2 <- Linf-(Linf-L1)*exp(-K*deltaT)
      return(dnorm(x=L2, mean=meanL2, sd=sigma, log=log))
    }
    loglik <- function(par,L1,L2,deltaT){
      return(sum(loglikvec(par,L1,L2,deltaT)))
    } # common wrapper for sum of log-likelihoods over obs
    
    postdev <- -2*apply(X=mcmc.est,MARGIN=1,FUN=loglik,L1=L1,L2=L2,deltaT=deltaT)
    # ^ post draws of deviance = -2*loglik
    ll.mean <- loglik(par=colMeans(mcmc.est),L1=L1,L2=L2,deltaT=deltaT)
    
    dic1 <- 2*(mean(postdev)+ll.mean)
    # ^ original DIC def with p_D, Spiegelhalter et al. (2002, JRSSB)
    dic2 <- mean((postdev-mean(postdev))^2)-2*ll.mean
    # ^ alt DIC def with p_V = 0.5*post var of deviance, Gelman et al. (2014)
    res$DIC <- c(dic1,dic2)
    names(res$DIC) <- c('pD','pV')
    
    # WAIC with variance-based "bias correction" term
    llmat <- apply(X=mcmc.est, MARGIN=1, FUN=loglikvec,
                   L1=L1, L2=L2, deltaT=deltaT)
    # ^ n rows, iter cols
    pWAIC2 <- mean(apply(X=llmat,MARGIN=1,FUN=var))
    # ^ sum over obs of emp var over MCMC draws, eq. (7.12) in Gelman et al. (2014)
    
    likmat <- apply(X=mcmc.est, MARGIN=1, FUN=loglikvec,
                    L1=L1, L2=L2, deltaT=deltaT, log=F) # post density
    lppd <- sum(log(rowMeans(likmat)))
    # ^ log pointwise pred density, eq. (7.5) in Gelman et al. (2014)
    
    res$WAIC <- -2*(lppd-pWAIC2)
  } else {
    res$cred.int <- NA
    res$DIC <- c(NA,NA)
    res$WAIC <- NA
  }
  
  res$AIC <- 2*opt$obj+2*length(opt$par) # output AIC anyway
  res$value.TMB <- opt$obj
  res$conv.TMB <- opt$mess
  
  return(res)
}







#///////////////////////////////////////////////////////////////////////////////
#### ja91: James (1991) ####
#///////////////////////////////////////////////////////////////////////////////

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
#### Bfr88: Bayesian Francis (1988), "our take" ####
#///////////////////////////////////////////////////////////////////////////////

Bfr88 <- function(par,L1,L2,deltaT,
                  priordist.Linf='lognormal',
                  priordist.K='uniform',
                  priordist.sigma='uniform',
                  hyperpar,
                  varfunc='constant',enablepriorsd=1L,
                  meth='nlminb',compute.se=T,
                  onlyTMB=F,output.post.draws=F,
                  mcmc.control=list('nchains'=5,'iter'=20000,'warmup'=10000,
                                    'adapt_delta'=0.8,'max_treedepth'=20)){
  
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
  
  ### setup varfunc code
  
  if (varfunc=='constant'){
    varfunc.code <- 0L
  } else if (varfunc=='prop.dL'){
    varfunc.code <- 1L
  } else if (varfunc=='prop.dT'){
    varfunc.code <- 2L
  } else if (varfunc=='prop.L2'){
    varfunc.code <- 3L
  } else if (varfunc=='exp.dL'){
    varfunc.code <- 4L
  } else if (varfunc=='exp.L2'){
    varfunc.code <- 5L
  } else if (varfunc=='pow.dL'){
    varfunc.code <- 6L
  } else if (varfunc=='pow.L2'){
    varfunc.code <- 7L
  } else {
    stop('varfunc must be one of: "constant", "prop.dL", "prop.dT", "prop.L2",',
         ' "exp.dL", "exp.L2", "pow.dL", or "pow.L2".')
  }
  
  ### setup
  if (meth!='nlminb'){warning('Only meth="nlminb" supported for now.')}
  # if (!compute.se){warning('se will be computed anyway.')}
  
  # n <- length(L1)
  
  datalist <- list('L1'=L1,'L2'=L2,
                   'deltaT'=deltaT, # as of v0.3: directly supply deltaT
                   'hp_Linf'=hyperpar[[1]],  # user-supplied
                   'hp_K'=hyperpar[[2]],     # user-supplied
                   'hp_sigma'=hyperpar[[3]], # user-supplied
                   'priordist'=priordist.code, # as of v0.2.1: vector dim 3
                   'priorsd'=as.integer(enablepriorsd), # 0=disabled, 1=enabled
                   'varfunc'=varfunc.code
  )
  parlist <- list('logLinf'=log(par[1]),
                  'logK'=log(par[2]),
                  'logsigma'=0, # "intercept" in error sd
                  'logsda'=-1,'logsdb'=-1 # multiplicative terms, smallish
  )
  
  ### TMB estimation (good for ini MCMC and as benchmark)
  
  if (varfunc.code==0L){ # constant
    obj <- MakeADFun(data=datalist,parameters=parlist,
                     map=list('logsda'=factor(NA),'logsdb'=factor(NA)),
                     # ^ do not estimate sda and sdb
                     DLL="FrancisBayesian",silent=T)
  }  else if (varfunc.code%in%c(1L,2L,3L)){ # prop.dL, prop.dT, prop.L2
    obj <- MakeADFun(data=datalist,parameters=parlist,
                     map=list('logsdb'=factor(NA)),
                     # ^ do not estimate sdb
                     DLL="FrancisBayesian",silent=T)
  } else { # exp.dL, exp.L2, pow.dL, pow.L2
    obj <- MakeADFun(data=datalist,parameters=parlist,
                     DLL="FrancisBayesian",silent=T)
  }
  
  opt <- nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                control=list(eval.max=5000,iter.max=5000))
  theta.tmb <- exp(opt$par[1:2]) # naive estimates without sdreport()
  
  ### MCMC based on tmbstan::tmbstan, by default uses NUTS
  if (!onlyTMB & all(priordist.code!=0)){
    mcmc.obj <- tmbstan(obj=obj,
                        lower=rep(-Inf,length(opt$par)),
                        upper=rep(Inf,length(opt$par)),
                        silent=T,laplace=F,
                        chains=mcmc.control$nchains,
                        warmup=mcmc.control$warmup,
                        iter=mcmc.control$iter,
                        control=list(adapt_delta=mcmc.control$adapt_delta,
                                     max_treedepth=mcmc.control$max_treedepth
                        ),
                        # control=list(adapt_delta=0.9, # def=0.8, larger=safer
                        #              max_treedepth=20), # def=10, larger helps
                        # init='random'
                        init='last.par.best' # start from MLE above
    )
    
    # traceplot(mcmc.obj, pars=names(obj$par), inc_warmup=TRUE) # check conv
    # pairs(mcmc.obj, pars=names(obj$par)) # post dist and scatterplots
    # # ^ Linf and K typically correlate a lot (negatively)
    
    # extract MCMC post draws for derived quantities specified in obj's REPORT
    mcmc.post <- as.matrix(mcmc.obj)
    
    if (varfunc=='constant'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=3) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        meanL2 <- Linf-(Linf-L1)*exp(-K*deltaT)
        return(dnorm(x=L2, mean=meanL2, sd=sigma, log=log))
      }
    } else if (varfunc=='prop.dL'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=4) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma','sda')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        sda <- par[4]
        meandL <- (Linf-L1)*(1-exp(-K*deltaT))
        return(dnorm(x=L2-L1, mean=meandL, sd=sigma+sda*meandL, log=log))
      }
    } else if (varfunc=='prop.dT'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=4) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma','sda')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        sda <- par[4]
        meandL <- (Linf-L1)*(1-exp(-K*deltaT))
        return(dnorm(x=L2-L1, mean=meandL, sd=sigma+sda*deltaT, log=log))
      }
    } else if (varfunc=='prop.L2'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=4) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma','sda')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        sda <- par[4]
        meanL2 <- Linf-(Linf-L1)*exp(-K*deltaT)
        return(dnorm(x=L2, mean=meanL2, sd=sigma+sda*meanL2, log=log))
      }
    } else if (varfunc=='exp.dL'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=5) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma','sda','sdb')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        sda <- par[4]
        sdb <- par[5]
        meandL <- (Linf-L1)*(1-exp(-K*deltaT))
        sd_expdL <- sigma + sdb*(1-exp(-sda*meandL))
        return(dnorm(x=L2-L1, mean=meandL, sd=sd_expdL, log=log))
      }
    } else if (varfunc=='exp.L2'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=5) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma','sda','sdb')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        sda <- par[4]
        sdb <- par[5]
        meanL2 <- Linf-(Linf-L1)*exp(-K*deltaT)
        sd_expL2 <- sigma + sdb*(1-exp(-sda*meanL2))
        return(dnorm(x=L2, mean=meanL2, sd=sd_expL2, log=log))
      }
    } else if (varfunc=='pow.dL'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=5) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma','sda','sdb')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        sda <- par[4]
        sdb <- par[5]
        meandL <- (Linf-L1)*(1-exp(-K*deltaT))
        sd_powdL <- sigma + sda*meandL^sdb
        return(dnorm(x=L2-L1, mean=meandL, sd=sd_powdL, log=log))
      }
    } else if (varfunc=='pow.L2'){
      mcmc.est <- matrix(NA_real_,nrow=nrow(mcmc.post),ncol=5) # all param
      for (i in 1:nrow(mcmc.post)){
        mcmc.est[i,] <- unlist(obj$report(mcmc.post[i,-ncol(mcmc.post)])[
          c('Linf','K','sigma','sda','sdb')])
      }
      
      loglikvec <- function(par,L1,L2,deltaT,log=T){
        Linf <- par[1]
        K <- par[2]
        sigma <- par[3]
        sda <- par[4]
        sdb <- par[5]
        meanL2 <- Linf-(Linf-L1)*exp(-K*deltaT)
        sd_powL2 <- sigma + sda*meanL2^sdb
        return(dnorm(x=L2, mean=meanL2, sd=sd_powL2, log=log))
      }
    }
    loglik <- function(par,L1,L2,deltaT){
      return(sum(loglikvec(par,L1,L2,deltaT)))
    } # common wrapper for sum of log-likelihoods over obs
    # colMeans(mcmc.est[,1:2]) # post means for Linf and K
    
    res <- list('par'=c(median(mcmc.est[,1]),median(mcmc.est[,2])),
                # 'par'=c(mean(mcmc.est[,1]),mean(mcmc.est[,2])),
                # ^ post median better, dist can be very skewed with small n
                'par.TMB'=theta.tmb)
    # ^ MCMC point estimates are posterior median as of v0.3.1
    names(res$par) <- c('Linf','K')
  } else { # then no MCMC and only report TMB (frequentist) estimate
    res <- list('par'=c(NA,NA),'par.TMB'=theta.tmb)
  }
  names(res$par.TMB) <- c('Linf','K')
  
  ### optional: compute standard errors
  if (compute.se){
    if (!onlyTMB & all(priordist.code!=0)){
      # res$se <- c(sqrt(var(mcmc.est[,1])),sqrt(var(mcmc.est[,2])))
      # ^ posterior naive se numerically unstable if low n
      res$se <- c(median(abs(mcmc.est[,1]-res$par[1])),
                  median(abs(mcmc.est[,2]-res$par[2]))) # MADAM
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
  
  ### output, AIC/DIC/WAIC
  if (output.post.draws){
    if (onlyTMB | any(priordist.code==0)){
      warning('No posterior draws if onlyTMB=TRUE or not all priordist set.')
    } else {
      res$post.draws <- list('Linf'=mcmc.est[,1],
                             'K'=mcmc.est[,2],
                             'sigma'=mcmc.est[,3])
    }
  }
  
  if (!onlyTMB & all(priordist.code!=0)){
    res$cred.int <- list('Linf'=quantile(mcmc.est[,1],probs=c(0.025,0.975)),
                         'K'=quantile(mcmc.est[,2],probs=c(0.025,0.975)),
                         'sigma'=quantile(mcmc.est[,3],probs=c(0.025,0.975))
    )
    # ^ equal-tailed 95% credible intervals based on MCMC draws
    
    # DIC, both p_D and p_V versions
    postdev <- -2*apply(X=mcmc.est, MARGIN=1, FUN=loglik,
                        L1=L1, L2=L2, deltaT=deltaT)
    ll.mean <- loglik(par=colMeans(mcmc.est),L1=L1,L2=L2,deltaT=deltaT)
    
    dic1 <- 2*(mean(postdev)+ll.mean)
    # ^ original DIC def with p_D, Spiegelhalter et al. (2002, JRSSB)
    dic2 <- mean((postdev-mean(postdev))^2)-2*ll.mean
    # ^ alt DIC def with p_V = 0.5*post var of deviance, Gelman et al. (2014)
    res$DIC <- c(dic1,dic2)
    names(res$DIC) <- c('pD','pV')
    
    # WAIC with variance-based "bias correction" term
    llmat <- apply(X=mcmc.est, MARGIN=1, FUN=loglikvec,
                   L1=L1, L2=L2, deltaT=deltaT)
    # ^ n rows, iter cols
    pWAIC2 <- mean(apply(X=llmat,MARGIN=1,FUN=var))
    # ^ sum over obs of emp var over MCMC draws, eq. (7.12) in Gelman et al. (2014)
    
    likmat <- apply(X=mcmc.est, MARGIN=1, FUN=loglikvec,
                    L1=L1, L2=L2, deltaT=deltaT, log=F) # post density
    lppd <- sum(log(rowMeans(likmat)))
    # ^ log pointwise pred density, eq. (7.5) in Gelman et al. (2014)
    
    res$WAIC <- -2*(lppd-pWAIC2)
  } else {
    res$cred.int <- NA
    res$DIC <- c(NA,NA)
    res$WAIC <- NA
  }
  
  res$AIC <- 2*opt$obj+2*length(opt$par) # output AIC anyway
  res$value.TMB <- opt$obj
  res$conv.TMB <- opt$mess
  
  return(res)
}



### wrapper for Bfr88 with model selection by min WAIC/AIC
Bfr88.minIC <- function(par,L1,L2,deltaT,
                        priordist.Linf='lognormal',
                        priordist.K='uniform',
                        priordist.sigma='uniform',
                        hyperpar,
                        enablepriorsd=1L,
                        varfunc='all',
                        meth='nlminb',compute.se=T,
                        onlyTMB=F,output.post.draws=F,
                        mcmc.control=list('nchains'=5,
                                          'iter'=20000,'warmup'=10000,
                                          'adapt_delta'=0.8,'max_treedepth'=20)){
  # could be improved by setting compute.se=F for all candidate models and
  # separate se comp without estim
  
  # except "constant" (= Fabens), all other variance functions are non-nested,
  # so simply fit all and keep on with smallest WAIC
  if (any(varfunc=='all')){
    varfunc.charvec <- c("constant","prop.dL","prop.dT","prop.L2",
                         "exp.dL","exp.L2","pow.dL","pow.L2") # test all
  } else { # then assumed subset of variance functions supplied
    varfunc.charvec <- varfunc
  }
  fitlist <- vector('list',length(varfunc.charvec))
  
  # wallclock <- proc.time()[3]
  for (j in 1:length(varfunc.charvec)){
    system.time(fitlist[[j]] <- Bfr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                                      priordist.Linf=priordist.Linf,
                                      priordist.K=priordist.K,
                                      priordist.sigma=priordist.sigma,
                                      hyperpar=hyperpar,
                                      varfunc=varfunc.charvec[j], 
                                      enablepriorsd=enablepriorsd,meth=meth,
                                      compute.se=compute.se,onlyTMB=onlyTMB,
                                      output.post.draws=output.post.draws,
                                      mcmc.control=mcmc.control))
  }
  # print(proc.time()[3]-wallclock)
  
  # rbind(sapply(fitlist,'[[','DIC'),
  #       sapply(fitlist,'[[','WAIC'),
  #       sapply(fitlist,'[[','AIC')) # meaningful only if frequentist est
  
  if (!all(c(priordist.Linf,priordist.K,priordist.sigma)%in%
           c('uniform','normal','Gaussian','lognormal')) |
      !as.logical(enablepriorsd) | onlyTMB){
    # then consider est as frequentist and use AIC for model selection
    ind.best <- which.min(sapply(fitlist,'[[','AIC')) # AIC
    message('Best (smallest AIC) Bfr88 model: varfunc = "',
            varfunc.charvec[ind.best],'"')
  } else { # then legit Bayesian est, use DIC computed from MCMC draws
    # ind.best <- which.min(sapply(fitlist,function(x){x$DIC[2]})) # DIC with pV
    ind.best <- which.min(sapply(fitlist,'[[','WAIC')) # WAIC
    message('Best (smallest WAIC) Bfr88 model: varfunc = "',
            varfunc.charvec[ind.best],'"')
  }
  
  return(c(fitlist[[ind.best]],'best.model'=varfunc.charvec[ind.best]))
}







#///////////////////////////////////////////////////////////////////////////////
#### fr88: "standard" Francis (1988) full, with min AIC for best sub-model ####
#///////////////////////////////////////////////////////////////////////////////

fr88 <- function(par,L1,L2,deltaT,
                 par.nu=NULL,par.m=NULL,par.p=NULL,
                 meth='nlminb',try.good.ini=T,compute.se=T,tol.sd=1e-5){
  # par=(Linf,K) argument now necessary, used as first ini to try and only one
  # if try.good.ini=F
  
  # par.nu, par.m ,par.p determine any additional "component" of the model. By
  # default, leaving all three NULL fits a Fabens model (with error sd s) with
  # Francis' reparametrization in terms of g_alpha and g_beta. Providing a
  # numerical value for any of the three enables that "component" and uses the 
  # suplied value as ini.
  
  if (meth!='nlminb'){stop('Only "nlminb" allowed for meth for now.')}
  
  ### setup
  alpha <- as.numeric(quantile(L1,0.1))
  beta <- as.numeric(quantile(L1,0.9))
  
  par.ini <- c((par[1]-alpha)*(1-exp(-par[2])), # g.alpha
               (par[1]-beta)*(1-exp(-par[2])),  # g.beta
               0)                               # log(nu)
  # ^ transform (Linf,K) into (ga,gb), hardcoded ini value for s on log scale,
  #   and then append below extra ini among c(log(nu), m, logit(p))
  
  ### objective function: enable components according to par.nu, par.m, par.p
  if (is.null(par.nu)){ # then no additive variance term nu*E[deltaL]
    if (is.null(par.m)){ # then no "measurement error" additive term m on mean
      if (is.null(par.p)){ # then no mixture
        # model: ga, gb, s
        names.par <- c('g_alpha','g_beta','log(s)')
        
        untransfo <- function(x){exp(x[3])}
        names.otherpar <- 's'
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s))
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          
          nll <- -sum(dnorm(x=L2-L1, mean=mean.dL, sd=s, log=T))
          
          return(nll) # negloglik
        }
      } else { # mixture model with "contamination" uniform within deltaL range
        # model: ga, gb, s, p
        par.ini <- c(par.ini, log(par.p/(1-par.p))) # logit transfo
        
        names.par <- c('g_alpha','g_beta','log(s)','logit(p)')
        
        untransfo <- function(x){c(exp(x[3]), 1/(1+exp(-x[4])))}
        names.otherpar <- c('s','p')
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s), logit(p))
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          p <- 1/(1+exp(-par[4])) # inv logit transfo
          
          dL <- L2-L1
          range.dL <- max(dL)-min(dL) # R in Francis' eq. (9)
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          
          pllvec <- log((1-p)*dnorm(x=dL,mean=mean.dL,sd=s,log=F) + p/range.dL)
          # ^ mixture: (1-p) Gaussian + p uniform [min deltaL, max deltaL]
          
          return(-sum(pllvec)) # negloglik
        }
      }
    } else { # "measurement error" additive term m on scale of E[deltaL]
      if (is.null(par.p)){ # then no mixture
        # model: ga, gb, s, m
        par.ini <- c(par.ini, par.m) # no transfo
        
        names.par <- c('g_alpha','g_beta','log(s)','m')
        
        untransfo <- function(x){c(exp(x[3]), x[4])}
        names.otherpar <- c('s','m')
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s), m)
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          m <- par[4] # "bias" in E[deltaL], Francis' lambda_i
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          
          nll <- -sum(dnorm(x=L2-L1, mean=mean.dL+m, sd=s, log=T))
          
          return(nll) # negloglik
        }
      } else { # fit a mixture model with "contamination" uniform within [minL,maxL]
        # model: ga, gb, s, m, p
        par.ini <- c(par.ini, par.m, log(par.p/(1-par.p)))
        # ^ no transfo for m, logit transfo for p
        
        names.par <- c('g_alpha','g_beta','log(s)','m','logit(p)')
        
        untransfo <- function(x){c(exp(x[3]), x[4], 1/(1+exp(-x[5])))}
        names.otherpar <- c('s','m','p')
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s), m, logit(p))
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          m <- par[4] # "bias" in E[deltaL], Francis' lambda_i
          p <- 1/(1+exp(-par[5])) # inv logit
          
          dL <- L2-L1
          range.dL <- max(dL)-min(dL) # R in Francis' eq. (9)
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          
          pllvec <- log((1-p)*dnorm(x=dL,mean=mean.dL+m,sd=s,log=F) + p/range.dL)
          # ^ mixture: (1-p) Gaussian + p uniform [min deltaL, max deltaL]
          
          return(-sum(pllvec)) # negloglik
        }
      }
    }
  } else { # additive variance term prop to E[deltaL]
    if (is.null(par.m)){ # then no "measurement error" additive term m on mean
      if (is.null(par.p)){ # then no mixture
        # model: ga, gb, s, nu
        par.ini <- c(par.ini, log(par.nu)) # log transfo
        
        names.par <- c('g_alpha','g_beta','log(s)','log(nu)')
        
        untransfo <- function(x){c(exp(x[3]), exp(x[4]))}
        names.otherpar <- c('s','nu')
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s), log(nu))
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          nu <- exp(par[4]) # "growth variability", heteroscedatic errors
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          posmean.dL <- ifelse(mean.dL>0, mean.dL, tol.sd) # for neg growth
          
          nll <- -sum(dnorm(x=L2-L1, mean=mean.dL, sd=s+nu*posmean.dL, log=T))
          # ^ affine variance function, Francis' eq. (5)
          
          return(nll) # negloglik
        }
      } else { # fit a mixture model with "contamination" uniform within [minL,maxL]
        # model: ga, gb, s, nu, p
        par.ini <- c(par.ini, log(par.nu), log(par.p/(1-par.p)))
        # ^ log transfo for nu, logit transfo for p
        
        names.par <- c('g_alpha','g_beta','log(s)','log(nu)','logit(p)')
        
        untransfo <- function(x){c(exp(x[3]), exp(x[4]), 1/(1+exp(-x[5])))}
        names.otherpar <- c('s','nu','p')
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s), log(nu), logit(p))
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          nu <- exp(par[4]) # "growth variability", heteroscedatic errors
          p <- 1/(1+exp(-par[5])) # inv logit
          
          dL <- L2-L1
          range.dL <- max(dL)-min(dL) # R in Francis' eq. (9)
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          posmean.dL <- ifelse(mean.dL>0, mean.dL, tol.sd) # for neg growth
          
          pllvec <- log((1-p)*dnorm(x=dL,mean=mean.dL,sd=s+nu*posmean.dL,log=F) +
                          + p/range.dL)
          # ^ mixture: (1-p) Gaussian + p uniform [min deltaL, max deltaL]
          
          return(-sum(pllvec)) # negloglik
        }
      }
    } else { # "measurement error" additive term m on scale of E[deltaL]
      if (is.null(par.p)){ # then no mixture
        # model: ga, gb, s, nu, m
        par.ini <- c(par.ini, log(par.nu), par.m)
        # ^ log transfo for nu, no transfo for m
        
        names.par <- c('g_alpha','g_beta','log(s)','log(nu)','m')
        
        untransfo <- function(x){c(exp(x[3]), exp(x[4]), x[5])}
        names.otherpar <- c('s','nu','m')
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s), log(nu), m)
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          nu <- exp(par[4]) # "growth variability", heteroscedatic errors
          m <- par[5] # "bias" in E[deltaL], Francis' lambda_i
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          posmean.dL <- ifelse(mean.dL>0, mean.dL, tol.sd) # for neg growth
          
          nll <- -sum(dnorm(x=L2-L1, mean=mean.dL+m, sd=s+nu*posmean.dL, log=T))
          # ^ affine variance function, Francis' eq. (5)
          
          return(nll) # negloglik
        }
      } else { # fit a mixture model with "contamination" uniform within [minL,maxL]
        # full model: ga, gb, s, nu, m, p
        par.ini <- c(par.ini, log(par.nu), par.m, log(par.p/(1-par.p)))
        # ^ log transfo for nu, no transfo for m, logit transfo for p
        
        names.par <- c('g_alpha','g_beta','log(s)','log(nu)','m','logit(p)')
        
        untransfo <- function(x){c(exp(x[3]), exp(x[4]), x[5], 1/(1+exp(-x[6])))}
        names.otherpar <- c('s','nu','m','p')
        
        nll <- function(par,alpha,beta,L1,L2,deltaT,tol.sd){
          # par=c(ga, gb, log(s), log(nu), m, logit(p))
          ga <- par[1] # > 0 but generally far away enough from 0
          gb <- par[2] # > 0 but generally far away enough from 0
          s <- exp(par[3]) # est on log scale because > 0
          
          nu <- exp(par[4]) # "growth variability", heteroscedatic errors
          m <- par[5] # "bias" in E[deltaL], Francis' lambda_i
          p <- 1/(1+exp(-par[6])) # inv logit
          
          dL <- L2-L1
          range.dL <- max(dL)-min(dL) # R in Francis' eq. (9)
          
          Linf <- (beta*ga-alpha*gb)/(ga-gb)
          expmK <- 1 + (ga-gb)/(alpha-beta) # exp(-K)
          mean.dL <- (Linf-L1)*(1-expmK^deltaT) # Francis' eq. (2)
          posmean.dL <- ifelse(mean.dL>0, mean.dL, tol.sd) # for neg growth
          
          pllvec <- log((1-p)*dnorm(x=dL,mean=mean.dL+m,sd=s+nu*posmean.dL,log=F) +
                          + p/range.dL)
          # ^ mixture: (1-p) Gaussian + p uniform [min deltaL, max deltaL]
          
          return(-sum(pllvec)) # negloglik
        }
      }
    }
  }
  
  
  ### optimization with supplied par ini
  opt <- try(nlminb(start=par.ini,objective=nll,
                    alpha=alpha,beta=beta,
                    L1=L1,L2=L2,deltaT=deltaT,
                    tol.sd=tol.sd),T)
  if (class(opt)=='try-error'){
    # opt$value <- Inf # so comparison with other loglik succeeds below
    warning('Supplied ini failed in nlminb, trying other ini.')
  } else {
    opt$value <- opt$obj
  }
  
  
  ### optional: try many different starting values, opt actually not convex
  if (try.good.ini | class(opt)=='try-error'){
    # then use my better ini, no subsets
    # ini for ga and gb based on Linf = max obs length at recap, and K = vB eq
    # for obs lengths with deltaT closest to 1 and using Linf = max obs length
    
    ind.dt1 <- which.min(abs(deltaT-1)) # obs with deltaT closest to 1
    ga.ini <- (max(L2)-alpha)*(L2[ind.dt1]-L1[ind.dt1])/(max(L2)-L1[ind.dt1])
    gb.ini <- (max(L2)-beta)*(L2[ind.dt1]-L1[ind.dt1])/(max(L2)-L1[ind.dt1])
    par.ini1 <- par.ini
    par.ini1[1:2] <- c(ga.ini, gb.ini)
    # ^ change only ini for g_alpha and g_beta, keep others as supplied
    
    opt1 <- nlminb(start=par.ini1,objective=nll,
                   alpha=alpha,beta=beta,
                   L1=L1,L2=L2,deltaT=deltaT,
                   tol.sd=tol.sd)
    opt1$value <- opt1$obj
    
    if (class(opt)=='try-error'){ # then replace
      opt <- opt1
    } else { # then need improvement in loglik to replace
      if (opt1$val < opt$val){ # then improves over first opt with supplied ini
        opt <- opt1
      }
    }
  } # otherwise only fit with supplied par as ini
  
  
  ### optional: compute standard errors
  if (compute.se){
    other.par <- opt$par[3:length(par.ini)] # all except g_alpha and g_beta
    hess <- optimHess(par=opt$par[1:2],
                      fn=function(ga.gb, other.par, alpha, beta,
                                  L1, L2, deltaT, tol.sd){
                        nll(c(ga.gb,other.par),alpha,beta,L1,L2,deltaT,tol.sd)
                      },other.par=other.par, alpha=alpha, beta=beta,
                      L1=L1, L2=L2, deltaT=deltaT, tol.sd=tol.sd)
    # ^ Hessian wrt ga and gb only, all other pars held fixed, more stable
    ga <- opt$par[1]
    gb <- opt$par[2]
    jac.g <- cbind(c(gb*(alpha-beta)/(ga-gb)^2,   # jacobian of transformation
                     -1/(ga-gb+alpha-beta)),      # from (ga,gb) to (Linf,K) 
                   c(-ga*(alpha-beta)/(ga-gb)^2,
                     1/(ga-gb+alpha-beta)))
    # lazy: numerical Hessian instead of analytical Fisher info of (ga,gb,...)
    # (same laziness as in fishmethods::grotag)
    varcov <- t(jac.g)%*%solve(hess[1:2,1:2])%*%jac.g # delta method
    # ^ other pars ignored in varcov, because Hessian rarely pos def with them
    se <- sqrt(diag(varcov))
  } else {
    hess <- NA
    se <- c(NA,NA)
  }
  
  ### output
  Linf <- (beta*opt$par[1] - alpha*opt$par[2])/(opt$par[1]-opt$par[2])
  K <- -log(1+(opt$par[1]-opt$par[2])/(alpha-beta))
  
  names(opt$par) <- names.par # all transformed param in optim
  
  otherpar <- untransfo(opt$par) # all other par except ga and gb
  names(otherpar) <- names.otherpar
  
  res <- list('par'=c(Linf,K),'se'=se,
              'other.par'=otherpar, # no se for other par for now
              'value'=opt$val,'AIC'=2*opt$val+2*length(par.ini),
              'original.par'=opt$par,'hessian'=hess,
              'alpha'=alpha,'beta'=beta,
              'conv'=opt$mess)
  names(res$par) <- c('Linf','K')
  names(res$se) <- c('Linf','K')
  
  return(res)
}


### wrapper for fr88 with model selection by min AIC
fr88.minAIC <- function(par,L1,L2,deltaT,
                        par.nu,par.m,par.p,
                        meth='nlminb',try.good.ini=T,
                        compute.se=T,tol.sd=1e-5){
  # could be improved by setting compute.se=F for all candidate models and
  # separate se comp without estim
  
  # Step 1: m1 = equiv to Fabens, baseline model
  fit <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
              par.nu=NULL,par.m=NULL,par.p=NULL,
              meth=meth,try.good.ini=try.good.ini,
              compute.se=compute.se,tol.sd=tol.sd) # equiv to Fabens
  # fit <- m1
  
  # Step 2: m2 = m1 + 1 param among 3 possible (nu, m, and p)
  m2 <- vector('list',3)
  m2[[1]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                  par.nu=par.nu,par.m=NULL,par.p=NULL,
                  meth=meth,try.good.ini=try.good.ini,
                  compute.se=compute.se,tol.sd=tol.sd)
  m2[[2]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                  par.nu=NULL,par.m=par.m,par.p=NULL,
                  meth=meth,try.good.ini=try.good.ini,
                  compute.se=compute.se,tol.sd=tol.sd)
  m2[[3]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                  par.nu=NULL,par.m=NULL,par.p=par.p,
                  meth=meth,try.good.ini=try.good.ini,
                  compute.se=compute.se,tol.sd=tol.sd)
  
  if (min(sapply(m2,'[[','AIC')) < fit$AIC){
    bestm2 <- which.min(sapply(m2,'[[','AIC'))
    fit <- m2[[bestm2]]
    # message('m2 improves over m1')
    
    # Step 3: m3 = m2 + 1 param among 2 remaining
    if (bestm2==1){
      # nu in m2, choice among m and p for m3
      m3 <- vector('list',2)
      m3[[1]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                      par.nu=par.nu,par.m=par.m,par.p=NULL,
                      meth=meth,try.good.ini=try.good.ini,
                      compute.se=compute.se,tol.sd=tol.sd)
      m3[[2]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                      par.nu=par.nu,par.m=NULL,par.p=par.p,
                      meth=meth,try.good.ini=try.good.ini,
                      compute.se=compute.se,tol.sd=tol.sd)
    } else if (bestm2==2){
      # m in m2, choice among nu and p for m3
      m3 <- vector('list',2)
      m3[[1]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                      par.nu=par.nu,par.m=par.m,par.p=NULL,
                      meth=meth,try.good.ini=try.good.ini,
                      compute.se=compute.se,tol.sd=tol.sd)
      m3[[2]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                      par.nu=NULL,par.m=par.m,par.p=par.p,
                      meth=meth,try.good.ini=try.good.ini,
                      compute.se=compute.se,tol.sd=tol.sd)
    } else {
      # p in m2, choice among nu and m for m3
      m3 <- vector('list',2)
      m3[[1]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                      par.nu=par.nu,par.m=NULL,par.p=par.p,
                      meth=meth,try.good.ini=try.good.ini,
                      compute.se=compute.se,tol.sd=tol.sd)
      m3[[2]] <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                      par.nu=NULL,par.m=par.m,par.p=par.p,
                      meth=meth,try.good.ini=try.good.ini,
                      compute.se=compute.se,tol.sd=tol.sd)
    }
    
    if (min(sapply(m3,'[[','AIC')) < fit$AIC){
      bestm3 <- which.min(sapply(m3,'[[','AIC'))
      fit <- m3[[bestm3]]
      # message('m3 improves over m2')
      
      # Step 4: m4 = full model
      m4 <- fr88(par=par,L1=L1,L2=L2,deltaT=deltaT,
                 par.nu=par.nu,par.m=par.m,par.p=par.p,
                 meth=meth,try.good.ini=try.good.ini,
                 compute.se=compute.se,tol.sd=tol.sd)
      
      if (m4$AIC < fit$AIC){
        fit <- m4
        # message('m4 improves over m3')
        bestmodel.charstr <- '(ga,gb,s) + (nu,m,p)'
        message('Best (smallest AIC) fr88 model: full with ',bestmodel.charstr)
      } else {
        # message('stick to m3')
        bestmodel.charstr <- paste0('(ga,gb,s) + (',c('nu','m','p')[bestm2],',',
                                    c('nu','m','p')[-bestm2][bestm3],')')
        message('Best (smallest AIC) fr88 model: ',bestmodel.charstr)
      }
    } else {
      # message('stick to m2')
      bestmodel.charstr <- paste0('(ga,gb,s) + (',c('nu','m','p')[bestm2],')')
      message('Best (smallest AIC) fr88 model: ',bestmodel.charstr)
    }
  } else {
    # message('stick to m1')
    bestmodel.charstr <- '(ga,gb,s)'
    message('Best (smallest AIC) fr88 model: simplest with only ',bestmodel.charstr)
  }
  return(c(fit,'best.model'=bestmodel.charstr))
}







#///////////////////////////////////////////////////////////////////////////////
#### fa65: Fabens (1965) ####
#///////////////////////////////////////////////////////////////////////////////

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
