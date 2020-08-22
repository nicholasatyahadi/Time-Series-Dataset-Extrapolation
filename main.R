# library(tidyverse)
# library(tseries)
# library(TSA)
# library(extraDistr)
# library(fitdistrplus)
# library(actuar)
# library(corrplot)
# library(lattice)
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

stand <- function(data){
  s.data <- NULL
  for(i in 1:length(data)){
    s.data <- c(s.data, (data[i]-min(data))/(max(data)-min(data)))
  }
  return(s.data)
}

rstand <- function(data,ref){
  r.data <- NULL
  for(i in 1:length(data)){
    r.data <- c(r.data,min(ref)+data[i]*(max(ref)-min(ref)))
  }
  return(r.data)
}

#import the dataset
sunspots <- read.table("sunspots.dat",header=T, sep=",")
income <- read.csv("s1902_sum_ann.csv",header = T, sep = ";")
climate <- read.csv('DailyDelhiClimateTrain.csv', header = T, sep = ",")

#making sure both are dataframes, and the the fields in both datasets are appropiate
class(sunspots);class(income);class(climate)
str(sunspots);str(income);str(climate)

#-----------------------------------------------------------------------------------------------------

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
sumdiff <- data.frame(Fields = statrep1[,1], Percent_Diff=(statrep2[,2]-statrep1[,2])/statrep1[,2])
statrep1;statrep2;sumdiff

#-----------------------------------------------------------------------------------------------------

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

#simulate using the same method as with moe and add one year
sim.sunspots <- as.vector(sapply(c(1:4),function(x){
  sapply(c(1:12),function(y){
    param.sunspots$ParVal[y]+runif(1,-1,1)*qnorm(0.975,0,1)*param.sunspots$ParSE[y]
  })
}))

#data checking
statrep1<-stat_rep(ex_sunspots[,2]);statrep2<-stat_rep(sim.sunspots)
sumdiff <- data.frame(Fields = statrep1[,1], Percent_Diff=(statrep2[,2]-statrep1[,2])/statrep1[,2])
statrep1;statrep2;sumdiff

#-----------------------------------------------------------------------------------------------------

#Multivariate dataset - New Delhi Climate Dataset
#Relation check between values
climate1 <- climate %>% mutate(year = format(as.Date(date),"%Y"),month = format(as.Date(date),"%m"))
climate1 <- climate1[-nrow(climate1),]
splom(climate1[,-c(1,(ncol(climate1)-1),ncol(climate1))])

#We do know there's an apparent relation between temperature and humidity
#But looking at the scatterplot, there's a weirdness in the pressure values, maybe a few outliers

#Observe by month, check the average and the standard deviations in handling outliers
sum_climate <- climate1 %>% group_by(year,month) %>% 
  summarise(avg_temp = mean(meantemp),sd_temp = sd(meantemp),
            avg_humid = mean(humidity), sd_humid = sd(humidity),
            avg_wind = mean(wind_speed), sd_wind = sd(wind_speed),
            avg_press = mean(meanpressure), sd_press = sd(meanpressure))

#Turns out, based on the large standard deviation, there must be some outliers in the pressure values.
press_check <- climate1 %>% group_by(year,month) %>% 
  summarise(avg_press = mean(meanpressure), min_press = min(meanpressure), max_press=max(meanpressure))

#Based on https://iridl.ldeo.columbia.edu/dochelp/QA/Basic/atmos_press.html the avg is around 1013.25 millibars
#Hence, if 990<meanpressure<1024 is normal, the rest is not.
#This is derived from the dataset and the website.
climate1 <- climate1 %>% 
  mutate(is_weird = ifelse(meanpressure > 990 & meanpressure < 1024,0,1))

#Then, we make new observation summary specifically for meanpressure using the signed dataset
delrows <- which(climate1$is_weird==1)
press_sum<-climate1[-delrows,] %>% group_by(year,month) %>% 
  summarise(avg_press = mean(meanpressure), sd_press = sd(meanpressure))

#Replace the "weird" values of meanpressure according to the month
for(i in 1:length(delrows)){
  curr_y = climate1$year[delrows[i]];curr_m = climate1$month[delrows[i]]
  sumrow = which(press_sum$year==curr_y & press_sum$month==curr_m)
  climate1$meanpressure[delrows[i]] <- press_sum$avg_press[sumrow]+sample(c(-1,1),1)*press_sum$sd_press[sumrow]
}

#Check the new monthly observation summary
clim_sum<-climate1 %>% group_by(year,month) %>% 
  summarise(avg_temp = mean(meantemp),sd_temp = sd(meantemp),
            avg_humid = mean(humidity), sd_humid = sd(humidity),
            avg_wind = mean(wind_speed), sd_wind = sd(wind_speed),
            avg_press = mean(meanpressure), sd_press = sd(meanpressure))

#Check the new correlation
climate1 <- climate1[,c(1:5)]
splom(climate1[,-1])

#Low: Pressure & Humidity, Pressure & Wind Speed, Temperature & Wind Speed
#Moderate: Wind Speed & Humidity
#High: Temperature & Humidity, Temperature & Pressure

#Since this is a time affected dataset, we might as well understand the correlation between the month and properties
clim_sum$year <- as.numeric(clim_sum$year);clim_sum$month <- as.numeric(clim_sum$month)
splom(clim_sum[,-c(1,2,4,6,8,10)])
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(clim_sum[,-c(1,2,4,6,8,10)]), method = "color",col=col(200),addCoef.col = "black",type ="upper")
#overall, the dataset have quite a high correlated values between each fields.
#most of the relation seem to be quadraticly related, so lets test with both linear and quadratic regression
#then we decide using the R-Squared values.

ts.climatesum <- list(ts(stand(clim_sum$avg_temp),start=2013,frequency=12),
                      ts(stand(clim_sum$avg_humid),start=2013,frequency=12),
                      ts(stand(clim_sum$avg_wind),start=2013,frequency=12),
                      ts(stand(clim_sum$avg_press),start=2013,frequency=12))
mos <- season(ts.climatesum[[1]])[1:48]
reg_df <- data.frame(idx = c(1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5),
                     idy = c(2,3,4,5,3,4,5,2,4,5,2,3,5,2,3,4))
names <- c("time",names(s.clim_sum))
reg_df$xlab <- names[reg_df$idx]
reg_df$ylab <- names[reg_df$idy]
quadr2 <- linearr2 <- adjlinearr2 <- adjquadr2 <- NULL
for(i in 1:nrow(reg_df)){
  y <- ts.climatesum[[(reg_df$idy[i]-1)]]
  if(reg_df$idx[i]==1){
    linearr2 <- c(linearr2,summary(lm(y~mos-1))$r.squared)
    adjlinearr2 <- c(adjlinearr2,summary(lm(y~mos-1))$adj.r.squared)
    quadr2 <- c(quadr2,0);adjquadr2<-c(adjquadr2,0)
  }else{
    x <- ts.climatesum[[(reg_df$idx[i]-1)]];x2 <- x^2
    linearr2 <- c(linearr2,summary(lm(y~x))$r.squared)
    adjlinearr2 <- c(adjlinearr2,summary(lm(y~x))$adj.r.squared)
    quadr2 <- c(quadr2,summary(lm(y~x+x2))$r.squared)
    adjquadr2 <- c(adjquadr2,summary(lm(y~x+x2))$adj.r.squared)
  }
}

reg_df$linearr2 <- linearr2;reg_df$quadr2 <- quadr2;
reg_df$adjlinearr2 <- adjlinearr2;reg_df$adjquadr2<-adjquadr2;reg_df <- reg_df<-reg_df[,-c(1,2)]
#we can see that high rsquared values in relation between temperature and pressure
#which is pressure affects the temperature
#the relation between wind speed with humidity should be noticed
#the regression model might not be a great fit but it aligns with the high correlation coefficient
#whereas the other regression models, doesn't seem to give out the best fit despite a high correlation coefficient
#therefore, we are going to estimate the extrapolated values, using:
#1. Pressure seasonal model
#2. Wind Speed seasonal model
#3. Pressure vs Temperature quadratic model
#4. Wind Speed vs Humidity quadratic model
x1<-ts.climatesum[[4]];x12<-x1^2;mod1 <- lm(x1~(mos-1))
x2<-ts.climatesum[[3]];x22<-x2^2;mod2 <- lm(x2~(mos-1))
mod3 <- lm(ts.climatesum[[1]]~x1+x12)
mod4 <- lm(ts.climatesum[[2]]~x2+x22)

#Estimation setup
param.pressure <- data.frame(Month = c(1:12),
                             ParVal = as.vector(summary(mod1)$coefficients[,1]),
                             ParSE = as.vector(summary(mod1)$coefficients[,2]))
param.wind <- data.frame(Month = c(1:12),
                             ParVal = as.vector(summary(mod2)$coefficients[,1]),
                             ParSE = as.vector(summary(mod2)$coefficients[,2]))
param.temp <- data.frame(ParVal = as.vector(summary(mod3)$coefficients[,1]),
                         ParSE = as.vector(summary(mod3)$coefficients[,2]))
param.humid <- data.frame(ParVal = as.vector(summary(mod4)$coefficients[,1]),
                         ParSE = as.vector(summary(mod4)$coefficients[,2]))

#Extrapolation
sim.pressure <- as.vector(sapply(c(1:4),function(x){
  sapply(c(1:12),function(y){
    param.pressure$ParVal[y]+runif(1,-1,1)*qnorm(0.975,0,1)*param.pressure$ParSE[y]
  })
}))

sim.wind <- as.vector(sapply(c(1:4),function(x){
  sapply(c(1:12),function(y){
    param.wind$ParVal[y]+runif(1,-1,1)*qnorm(0.975,0,1)*param.wind$ParSE[y]
  })
}))

currparams <- param.temp$ParVal+runif(48,-1,1)*qnorm(0.975,0,1)*param.temp$ParSE
sim.temp <- currparams[1]+currparams[2]*sim.pressure+currparams[3]*sim.pressure^2

currparams <- param.humid$ParVal+runif(48,-1,1)*qnorm(0.975,0,1)*param.humid$ParSE
sim.humid <- currparams[1]+currparams[2]*sim.wind+currparams[3]*sim.wind^2

#Converting back to normal values from standardized values
sim.temp <- rstand(sim.temp,clim_sum$avg_temp)
sim.humid <- rstand(sim.humid,clim_sum$avg_humid)
sim.wind <- rstand(sim.wind,clim_sum$avg_wind)
sim.pressure <- rstand(sim.pressure,clim_sum$avg_press)

#Converting back to daily values from monthly agg
days <- data.frame(year=clim_sum$year,
                   month=clim_sum$month)%>%
  mutate(nodays = ifelse(month %in% c(1,3,5,7,8,10,12),31,
                         ifelse(month %in% c(4,6,9,11),30,
                                ifelse(month==2 & as.numeric(year)%%4==0,29,28))))

sim.daily.temp <- sim.daily.humid <- sim.daily.wind <- sim.daily.pressure <- NULL
for(i in 1:4){
  for(j in 1:12){
    row <- i*j
    sim.daily.temp <- c(sim.daily.temp,sim.temp[row]+runif(days$nodays[row],-1,1)*qnorm(0.975,0,1)*clim_sum$sd_temp[row])
    sim.daily.humid <- c(sim.daily.humid,sim.humid[row]+runif(days$nodays[row],-1,1)*qnorm(0.975,0,1)*clim_sum$sd_humid[row])
    sim.daily.wind <- c(sim.daily.wind,sim.wind[row]+runif(days$nodays[row],-1,1)*qnorm(0.975,0,1)*clim_sum$sd_wind[row])
    sim.daily.pressure <- c(sim.daily.pressure,sim.pressure[row]+runif(days$nodays[row],-1,1)*qnorm(0.975,0,1)*clim_sum$sd_press[row])
  }
}

mean(sim.daily.temp)-mean(climate1$meantemp)
mean(sim.daily.humid)-mean(climate1$humidity)
mean(sim.daily.wind)-mean(climate1$wind_speed)
mean(sim.daily.pressure)-mean(climate1$meanpressure)