install.packages()

library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(MASS)
theme_set(theme_pubclean())
library(psych)

setwd("~/Desktop/IUclasses/stat52501-GLM/Project")
mydata <- read.csv("SeoulBikeData.csv", check.names = F)
colnames(mydata)
colnames(mydata) <- c("Date","Count","Hour","Temperature","Humidity","Wind_Speed", "Visibility", "Dew_Point_Temp",
             "Solar_Radiation", "Rainfall", "Snowfall", "Season", "Holiday", "Functioning_Day")
mydata <- as.data.frame(mydata)
str(mydata)
summary(mydata)
table(mydata$Functioning_Day)
mydata <- mydata %>%
  mutate(Count_id = ifelse(Count > "0",1,0))
table(mydata$Functioning_Day, mydata$Count_id) # 295 no and 8465 yes 
dt<-mydata
dt <- mydata %>% 
  filter(dt$Functioning_Day == "Yes") 
dt <- dt[ -c(14:15) ]
sum(is.na(dt)) # no missing values 
#dayofweek
#Sun  Mon  Tue  Wed  Thu  Fri  Sat
#No    24   24   24   72   72   79    0
#Yes 1224 1224 1224 1176 1176 1169 1272
dt$Date<- as.POSIXlt(dt$Date,"%d/%m/%Y")

par(mfrow=c(2,2))

attach(dt)
month<- as.factor(month(dt$Date))
day<-as.factor(day(dt$Date))
dayofweek<- wday(dt$Date, label=TRUE)

data<- dt[-c(1:2)]
data$year <- as.factor(year(dt$Date))
data$month<- as.factor(month(dt$Date))
data$day<-as.factor(day(dt$Date))
data$dayofweek<- wday(dt$Date, label=TRUE)

chrnames<-c("Season", "Holiday", "year", "month", "day", "dayofweek")
data <- data %>% mutate_at(names, as.numeric)
data<-as.matrix(data)
str(data)
cor<-cor(data)
corPlot(cor(data), cex = 1.2, is.corr=FALSE)


ggplot(dt, aes(x = month, y = Count)) +
  geom_line(stat = "identity") +
  geom_boxplot()
dt %>%
  arrange(Count) %>%
  mutate(Season = factor(Season, levels=c("Winter", "Spring", "Summer", "Autumn"))) %>%
  ggplot(aes(x = Season, y = Count)) +
  geom_bar(stat='identity') 
ggplot(dt, aes(x = dayofweek, y = Count)) +
  geom_bar(stat = "identity")
ggplot(data=dt, aes(x=Hour, y=Count, group=1)) +
  geom_bar(stat = "identity")
ggplot(data=dt, aes(x=Hour, y=Count, group=1)) +
  stat_summary(fun=mean, geom="bar", position = "stack") 
ggplot(dt, aes(x = day, y = Count)) +
  geom_bar(stat = "identity")
ggplot(data=dt, aes(x=Solar_Radiation, y=Count, group=1)) +
  geom_bar(stat = "identity")

ggplot(dt)+geom_point(aes(x=as.Date(Date),y=Count,color=as.character(Season)),alpha=0.3)+
  labs(title="Scatter plot of Rent Bike Number vs Date")+
  xlab(label="Date")+
  ylab(label="No. of Rent Bike")+
  scale_x_date(breaks="month",date_labels="%Y-%m")

mod1 <- glm(Count~ dayofweek + Temperature + Humidity + Wind_Speed + Visibility + Dew_Point_Temp + Solar_Radiation +
                             Rainfall + Snowfall + Season + Holiday ,data=dt, family = poisson)
mod2<- glm(Count~ dayofweek + Temperature + Humidity + Wind_Speed + Visibility + Dew_Point_Temp + Solar_Radiation +
             Rainfall + Snowfall + Season + Holiday ,data=dt, family = poisson(link="log"))

mod3 <- glm(Count~ dayofweek + Temperature + Humidity + Wind_Speed + Visibility + Dew_Point_Temp + Solar_Radiation +
              Rainfall + Snowfall + Season + Holiday ,data=dt, family = gaussian(link="log"))

mod4 <- glm.nb(Count~ dayofweek + Temperature + Humidity + Wind_Speed + Visibility + Dew_Point_Temp + Solar_Radiation +
                 Rainfall + Snowfall + Season + Holiday , data = dt, start=c( )) 
coefini=coef(glm(Count~ dayofweek + Temperature + Humidity + Wind_Speed + Visibility + Dew_Point_Temp + Solar_Radiation +
                   Rainfall + Snowfall + Season + Holiday , data=dt))
