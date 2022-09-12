########################################################################
# Bayesian implementation of Fabens method (BFa) 
########################################################################
 
########################################################################
# For your species modify code as explained in the Supplement and below 
# from rows 11-85
########################################################################
 
### Set Working Directory (folder where all files are saved):
setwd("D:/GrowthEstimation") # this needs to be changed to your folder path!

### Before first use only make sure Rtools is installed:
# https://cran.r-project.org/bin/windows/Rtools/
# as well as the packages TMB and tmbstan (see Supplementary Information for details)

############################################
require(TMB) 
require(tmbstan)
source('GrowthEstimation_CapRecapSim.r')
source('GrowthEstimation_Methods.r')  
############################################
compile("FabensBayesian.cpp") 
dyn.load(dynlib("FabensBayesian"))

### Number of iterations, warm-up and chains used for Bayesian method:
ITER<-10000 # Recommended are at least 10000 iterations, more iterations are better but will take longer 
WARM<-5000 # Recommended is a warm-up of at least 5000
CHAINS<-3 # Recommended are at least 3 chains, more chains are better but will take more time

### Y = data and outlier removal as described in the manuscript is applied
### N = No data point is removed (note that also negative growth is not removed then)
Exclude<-"N" # alternative is Exclude<-"Y"

### Min number of days between a mark and a recapture event (minimum deltaT)
MinDays<-90 #if possible times at liberty should be larger 365 days minimum, given that growth is often seasonal

### Define prior distribution for the parameters
### Choices are: uniform, normal or lognormal
LinfPriorDist<-'lognormal' # uniform # normal # lognormal
KPriorDist<-'uniform' # uniform # normal # lognormal
SigmaPriorDist<-'uniform' # uniform # normal # lognormal

### Define the maximum length (your best guess of the average maximum length in the population, Linf)
Lmax<-133
### Define an upper value for your best guess on Linf (this value needs to be at least slightly larger than Lmax)
UpperLmax<-200

### Plot the prior distribution and generate the input mean and sd
LinfPrior<-hp.lognormal(Lmax, UpperLmax, plot=T) # for a normal distributed prior on Linf use LinfPrior<-hp.normal(Lmax, UpperLmax, plot=T)

### Specify the priors here, when prior distribution is uniform specify min and max,
### when otherwise specify mean and sd. 
LinfPr1<-as.numeric(LinfPrior[1]) #Linf prior, if uniform than this is the minimum bound, else it is the mean
LinfPr2<-as.numeric(LinfPrior[2]) #Linf prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
KPr1<-10^-10 #k prior, if uniform than this is the minimum bound, else it is the mean
KPr2<-100 #k prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)
SPr1<-10^-10 #sigma prior, if uniform than this is the minimum bound, else it is the mean
SPr2<-10^3 #sigma prior, if uniform than this is the maximum bound, else it is the standard deviation (sd)

priorlist<-list(c(LinfPr1,LinfPr2),
                c(KPr1,KPr2),
                c(SPr1, SPr2))

### Load the tagging data from the Excel CSV file:
rawdata<- "ExampleMarkRecapture.csv" # rename according to your file name
rawdata<- read.csv(rawdata, header=T, sep=",")
rawdata 
# make sure the data in the console looks like
#> rawdata
#L1        L2    deltaT
#1   72.26157  90.28036  1.281660
#2   81.46210  88.23489  3.656536
#3   91.63498 112.53977  6.822329
#4  112.46007 114.14646  6.438745
#5   60.31395  97.22382  6.106889
#6  106.73469 112.84647  5.146008
#7   60.23654  95.57799  5.695239
#8   93.64087 122.98888  9.128786
#9   83.24406 116.28717 13.431396
#10  83.11262 108.11820  3.111449
#If not, try a different separator, such as a ";" instead of a ",": 
#rawdata<- "ExampleMarkRecapture.csv" 
#rawdata<- read.csv(rawdata, header=T, sep=";")

########################################################################################
### You can now highlight the code from here until the end of the script and click 'Run'
########################################################################################

### Delete observations with short times at liberty
ndays <- 365.2425 # nb of days in 1 year, for conversion
rawdata2<-subset(rawdata, deltaT>(MinDays/ndays))
rawdata2<-droplevels(rawdata2)

### Calculate growth rate per year
rawdata2$Growth_rate_cm_Year<-(rawdata2$L2-rawdata2$L1)/rawdata2$deltaT

# remove outliers
#windows()
boxplot(rawdata2$Growth_rate_cm_Year,horizontal = TRUE, xlab="Growth per year in cm")

outs<-boxplot.stats(rawdata2$Growth_rate_cm_Year)
RemovedN<-length(outs$out)
outlier<-outs$out

rawdata2$MATCH<-match(rawdata2$Growth_rate_cm_Year, outlier)
rawdata3<-subset(rawdata2, is.na(MATCH)=="TRUE")
rawdata3<-droplevels(rawdata3)

### Exclude negative growth rates
rawdata3$Length_change<-rawdata3$L2-rawdata3$L1
data<-subset(rawdata3, Length_change>=0)
data<-droplevels(data)

### Use original dataset (no data points excluded) or approach as described in the manuscript (Exclude = Y):
invisible(ifelse(Exclude=="Y", data<-data,
                 data<-rawdata2))

### Transform dataframe to list
datalist <- as.list(as.data.frame(data))

### Generate starting values
Linfstart<-mean(as.numeric(Lmax)/0.99)
Kstart <- as.numeric(median(-log((Linfstart-data$L2)/(Linfstart-data$L1))/data$deltaT, na.rm=T))
starts<-c(Linfstart, Kstart)

### Run BFa:
Bfa<-Bfa65(par=starts,L1=datalist$L1,L2=datalist$L2,
           deltaT=datalist$deltaT,
           priordist.Linf=LinfPriorDist,
           priordist.K=KPriorDist,
           priordist.sigma=SigmaPriorDist,
           hyperpar=priorlist,
           meth='nlminb',compute.se=T,
           onlyTMB=F,output.post.draws=F,
           mcmc.control=list('nchains'=CHAINS,'iter'=ITER,'warmup'=WARM,
                             'adapt_delta'=.8,'max_treedepth'=20))

# Initial sample size before data point removal (if selected)
InitialN<-length(rawdata$L1)
InitialN 
# Get final sample after data pint removal
FinalN<-length(data$L1)
FinalN
# Get priorlist 
# [[1]]== Linf min max or mean and sd for uniform or normal/lognormal priors respectively
# [[2]]== K min max or mean and sd for uniform or normal/lognormal priors respectively
# [[3]]== Sigma min max or mean and sd for uniform or normal/lognormal priors respectively
priorlist
# Get parameter estimates
Bfa$par
# Get parameter uncertainty (95% credible intervals)
Bfa$cred.int