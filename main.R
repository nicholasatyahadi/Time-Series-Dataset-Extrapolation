# library(tidyverse)
# library(tseries)
# library(TSA)
# library(extraDistr)
# library(fitdistrplus)
# library(actuar)
# setwd("~/R/Projects/TS Simulation")

stat_rep <- function(data){
  fields <- c("Min","Max","Mean","SD","Q1","Q2","Q3")
  res <- c(min(data),
           max(data),
           mean(data),
           sd(data),
           as.numeric(quantile(data,c(0.25,0.5,0.75))))
  resdf<-data.frame(Fields=fields, Val = res)
  return(resdf)
}

#import the dataset
sunspots <- read.table("sunspots.dat",header=T, sep=",")
income <- read.csv("s1902_sum_ann.csv",header = T, sep = ";")

#making sure both are dataframes, and the the fields in both datasets are appropiate
class(sunspots);class(income)
str(sunspots);str(income)

#extract sunspots dataset to a smaller number of data
nrow(sunspots) #with 2820 rows, we will use the last 3 years of observations as simplification
#we cannot use random sampling since this is a time series dataset
#which are time-affected dataset, so it's best to extract like so
n<-36
ex_sunspots <- sunspots[c((nrow(sunspots)-n+1):nrow(sunspots)),]

#Univariate dataset with MOE - US Mean Income Dataset
income <- income %>% mutate(min = mean-moe, max = mean+moe)

#We will randomize the data per year, to simulate a fluctuating income
#in this example, we will use uniform(0,1) and lognormal distribution for this example
#for lognormal distribution, we will fit the data to at least get a better view of it
inc_lnorm <- fitdist(income[,2],"lnorm",method = "mle");gofstat(inc_lnorm);summary(inc_lnorm)
incparam <- as.vector(unlist(inc_lnorm[1]))
income <- income %>% mutate(plmin = plnorm(min,meanlog=incparam[1],sdlog=incparam[2]),
                  plmax = plnorm(max,meanlog=incparam[1],sdlog=incparam[2])) #pikirin lagi gini gimana

#say that we want to make monthly values
u_income <- sapply(c(0:(nrow(income)*12-1)),function(x){
  runif(1)*(income$max[floor(x)/12+1]-income$min[floor(x)/12+1])+income$min[floor(x)/12+1]
})

#turn the dataset into time series dataset
ts.income <- ts(u_income,start = 1, frequency = 12)
ts.plot(ts.income)

#data checking
statrep1<-stat_rep(income[,2]);statrep2<-stat_rep(u_income)
sumdiff <- data.frame(Fields = statrep1[,1], Percent_Diff=(statrep1[,2]-statrep2[,2])/statrep1[,2])
statrep1;statrep2;sumdiff

#Univariate dataset without MOE - Sunspots Dataset
#extract sunspots dataset to a smaller number of data
nrow(sunspots) #with 2820 rows, we will use the last 3 years of observations as simplification
#we cannot use random sampling since this is a time series dataset
#which are time-affected dataset, so it's best to extract like so
n<-36
ex_sunspots <- sunspots[c((nrow(sunspots)-n+1):nrow(sunspots)),]

#convert to time series
t.sunspots <- ts(ex_sunspots[,2],start=1,frequency = 12)
mos <- season(t.sunspots);ts.plot(t.sunspots)

#create a seasonal linear regression (or any regression)
lm.sunspots <- lm(t.sunspots ~ mos-1);sum.lm.sunspots <- summary(lm.sunspots)

#extract the parameter values and standard errors
param.sunspots <- data.frame(Month = c(1:12),
                             ParVal = as.vector(sum.lm.sunspots$coefficients[,1]),
                             ParSE = as.vector(sum.lm.sunspots$coefficients[,2]))

#add moe with 95% of confidence level
param.sunspots <- param.sunspots %>% mutate(Lower=ParVal-qnorm(0.975,0,1)*ParSE,
                                            Upper=ParVal+qnorm(0.975,0,1)*ParSE)

#simulate using the same method as with moe
sim.sunspots <- as.vector(sapply(c(1:3),function(x){
  sapply(c(1:12),function(y){
    param.sunspots$Lower[y]+runif(1)*(param.sunspots$Upper[y]-param.sunspots$Lower[y])
  })
}))

#data checking
statrep1<-stat_rep(ex_sunspots[,2]);statrep2<-stat_rep(sim.sunspots)
sumdiff <- data.frame(Fields = statrep1[,1], Percent_Diff=(statrep1[,2]-statrep2[,2])/statrep1[,2])
statrep1;statrep2;sumdiff

#Multivariate dataset