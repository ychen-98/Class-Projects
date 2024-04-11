
{
  library(MASS)
  library(mfp)
  library(plotrix)
  library(psych)
}
setwd("~/Desktop/IUclasses/stat52501-GLM/Project")

{
  columns <- c("Date","Count","Hour","Temperature","Humidity","Wind Speed", "Visibility", "Dew Point", 
               "Solar Radiation", "Rainfall", "Snowfall", "Season", "Holiday", "Functioning Day")
  data <- read.table(file="SeoulBikeData.csv", header = T, col.names = columns, sep = ",")
  data <- data[data$Functioning.Day == "Yes",]
  data$Functioning.Day <- NULL
  data$Date <- NULL
  data$Hour <- as.factor(data$Hour)
  data$Holiday <- as.factor(data$Holiday)
  data$Holiday <- relevel(data$Holiday, "No Holiday")
  data$Season <- as.factor(data$Season)
  data$Season <- relevel(data$Season, "Winter")
}
{
  data$Rain <- as.factor(cut(data$Rainfall, breaks = c(-Inf, 0.001, 2.5-0.001, 7.6-0.001, Inf), labels=c("None", "Light", "Moderate", "Heavy")))
  data$Snow <- as.factor(cut(data$Snowfall, breaks = c(-Inf, 0.001, 1*13/10-0.001, 2.5*13/10-0.001, Inf), labels=c("None", "Light", "Moderate", "Heavy")))
  data$Wind <- as.factor(cut(data$Wind.Speed, breaks = c(0, 0.5, 1.5, 3.3, 5.5, Inf), include.lowest = T, right = F, labels = seq(0, 4)))
  nb.model1 <- glm.nb(Count~Hour+Temperature+Humidity+Wind+Rain+Snow+Season+Holiday, data)
  
  data$WindGroup <- as.factor(cut(data$Wind.Speed, breaks = c(0, 1.5, 5.5, Inf), include.lowest = T, right = F, labels=c("01", "23", "4")))
  nb.model2 <- update(nb.model1, Count~.-Wind+WindGroup)
}
{
  data$Temperature.FP1 <- I(((data$Temperature+17.9)/10)^3)
  data$Temperature.FP2 <- I(((data$Temperature+17.9)/10)^3*log(((data$Temperature+17.9)/10)))
  data$Humidity.FP1 <- I(((data$Humidity+1)/100)^3)
  data$Humidity.FP2 <- I(((data$Humidity+1)/100)^3*log(((data$Humidity+1)/100)))
  nb.model3 <- glm.nb(Count~Hour+Temperature.FP1+Temperature.FP2+Humidity.FP1+Humidity.FP2+WindGroup+Rain+Snow+Season+Holiday, data)
  
  data$Solar.Radiation.FP1 <- I((data$Solar.Radiation+0.1)^-0.5)
  nb.model4 <- update(nb.model3, Count~.+Solar.Radiation.FP1)
}

nb.int <- update(nb.model4, Count~.+Temperature.FP1*Rain+Temperature.FP1*Snow+Temperature.FP1*(Humidity.FP1+Humidity.FP2)
                 +Temperature.FP1*WindGroup+Temperature.FP1*Solar.Radiation.FP1+(Humidity.FP1+Humidity.FP2)*WindGroup
                 +(Humidity.FP1+Humidity.FP2)*Solar.Radiation.FP1)
summary(nb.int)

# obtain residuals and fitted
eta <- nb.int$linear.predictors 			# x_i^T β ̂ i.e. the linear predictor for everyone, not estimated prob
res.Dev <- residuals(nb.int, type="deviance")
eta.1 <- nb.model4$linear.predictors 			# x_i^T β ̂ i.e. the linear predictor for everyone, not estimated prob
res.Dev.1 <- residuals(nb.model4, type="deviance")
# standardized deviance residuals
res.stdDev <- residuals(nb.int, type="deviance")/sqrt(1 - hatvalues(nb.int))		# regular deviance res/sqrt(1-hi)
lo.stdDev <- loess(res.stdDev~eta)	#loess = Locally Weighted Scatterplot Smoothing, a non-parametric method
# standardized Pearson residuals
res.stdPear <- residuals(nb.int, type="pearson")/sqrt(1 - hatvalues(nb.int))
lo.stdPear <- loess(res.stdPear~eta)
# studentized deleted residuals
res.std <- rstudent(nb.int)
lo.std <- loess(res.std~eta)


#### residual vs fitted
par(mfrow=c(1,3))
plot(y=res.stdDev, x=eta, main="Standardized Deviance Residuals vs linear predictor"
            , xlab=expression(eta), ylab="Standardized Deviance Residuals")
lines(x=eta[order(eta)], y=predict(lo.stdDev)[order(eta)], 
              col='red', lwd=2)	# add loess smoothed curve in red

plot(y=res.stdPear, x=eta, main="Standardized Pearson Residuals vs linear predictor", 
              xlab=expression(eta), ylab="Standardized Pearson Residuals ")
lines(x=eta[order(eta)], y=predict(lo.stdPear)[order(eta)], col='red', lwd=2)

plot(y=res.std, x=eta, main="Studentized deleted Residuals vs linear predictor", 
              xlab=expression(eta), ylab="Studentized deleted Residuals")
lines(x=eta[order(eta)], y=predict(lo.std)[order(eta)], col='red', lwd=2)

res.stdDev[which(abs(res.std)>6)]
data[which(abs(res.std)>6), ]
predict(lo.std)[which(abs(res.std)>6)]

#### index plot
par(mfrow=c(1,3))
data$ObsInd <- seq.int(nrow(data))			# add an obs index var
data[1:5,]
plot(y=res.stdDev, x=data$ObsInd, type="b", 
     main="Index Plot using Standardized D Residual", xlab="Index") 
# type="b": both pts and line
plot(y=res.Dev, x=data$ObsInd, type="b", 
     main="Index Plot using Deviance Residual", xlab="Index")
#### leverage plot
plot(y=hatvalues(nb.int), x=data$ObsInd, type="b", 
     main="Leverage Plot", xlab="Index")
### dfbeta
# gender ref is Female
db1 <- dfbeta(nb.int)		
dim(db1)			# impact of each obs on every parameter estimates beta.hat in the model
#[1] 8465   58

# impact of obs on Age coef
for (i in 2:ncol(db1)) {
  plot(y=db1[,i], x=data$ObsInd, type="b", main=paste0("DfBeta for", colnames(db1)[i], " Coef" )) 	# original: black/white 
}
par(mfrow=c(2,2))
plot(y=db1[,28], x=data$ObsInd, type="b", 
     main=paste0("DfBeta for", colnames(db1)[28], " Coef" ), ylab=" ", xlab="index") 
plot(y=db1[,30], x=data$ObsInd, type="b", 
     main=paste0("DfBeta for", colnames(db1)[30], " Coef" ), ylab=" ", xlab="index") 
plot(y=db1[,54], x=data$ObsInd, type="b", 
     main=paste0("DfBeta for", colnames(db1)[54], " Coef" ), ylab=" ", xlab="index") 
plot(y=db1[,56], x=data$ObsInd, type="b", 
     main=paste0("DfBeta for", colnames(db1)[56], " Coef" ), ylab=" ", xlab="index") 

### dffits
par(mfrow=c(1,2))
df1 <- dffits(nb.int)
length(df1)			# impact of each obs on overall model fitting
plot(y=df1, x=data$ObsInd, type="b", main="DfFits", 
     ylab="DFFITS", xlab="Index")
### Cook's distance
cook1 <- cooks.distance(nb.int)
length(cook1)			# impact of each obs on overall model fitting
plot(y=cook1, x=data$ObsInd, type="b", main="Cook's Distance",
     ylab="Cook's distance", xlab="Index")

par(mfrow=c(2,2))
plot(nb.int)



library(car)
vif(nb.int)
