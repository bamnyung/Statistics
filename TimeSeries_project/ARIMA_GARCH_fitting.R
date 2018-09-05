##### Time Series Project - Fit bit-coin price data to time-series model #####

library(timeSeries)
library(rugarch)
library(readxl)
library(forecast)
library(ggplot2)
library(cowplot)

##### Plot raw data #####

setwd("C://Users//Minyoung//Desktop//시계열//pj//dataset")
BTC_4h <- read.csv("BTC_4h.csv")
coin_name <- "BTC_4h" #CHANGE THIS!!
coin_raw <- BTC_4h #CHANGE THIS!!!

coin_raw$date <- as.POSIXct(as.character(coin_raw$date), format = "%Y-%m-%d %H:%M:%S") #EXECUTE THIS EXCEPT 1-DAY DATA!!
p_raw <- ggplot(data = coin_raw, aes(x = date, y = close))+ 
  geom_line(color="#00AFBB", size=1)+theme_gray()+
  ggtitle("Close Price of BITCOIN")
p_raw


##### log trasform, 1차 차분 #####

coin_raw$logprice <- log(coin_raw$close)
rownames(coin_raw) <- coin_raw$date
#coin_raw$date <- NULL
p_log <- ggplot(data=coin_raw, aes(x=date, y=logprice))+
  geom_line(color="#E7B800", size=1)+theme_gray()+
  ggtitle("Logprice of BITCOIN")
plot_grid(p_raw, p_log, align = "v",ncol=1)

coin_raw$diff <- c(0, diff(coin_raw$logprice))
ggplot(data=coin_raw, aes(x=date, y=diff))+
  geom_line(color="#00AFBB", size=0.5)+theme_gray()

#Check if diff is stationary
bacf <- acf(coin_raw$logprice[c(6501:7000)], plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
q_log <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
       geom_hline(aes(yintercept = 0)) +
       geom_segment(mapping = aes(xend = lag, yend = 0))+
       ggtitle("ACF of logprice")+
       theme_gray()
bacf <- acf(coin_raw$diff[c(6501:7000)], plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
q_diff <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
       geom_hline(aes(yintercept = 0)) +
       geom_segment(mapping = aes(xend = lag, yend = 0))+
       ggtitle("ACF of diffed logprice")+
       theme_grey()
plot_grid(q_log, q_diff, align = "v",ncol=1)

# ADF-TEST
library(tseries)
adf.test(coin_raw$diff) #0.05수준에서 reject H0=>stationary


##### Fit ARIMA with 500 data (2018-02-07 09:00:00 ~ 2018-05-01 13:00:00) #####

x <- coin_raw$date[6501:7000]
y <- as.xts(subset(coin_raw, select = diff))[6501:7000]
forecast_y <- as.xts(subset(coin_raw, select = diff))[7001:7043]
data <- data.frame(date=x,difflog=y)
arima.order <- arimaorder(auto.arima(y, seasonal=F,allowmean = T)) #auto.arima로 차수결정
arima.fit <- arima(y, order=arima.order, include.mean = T) #arima로 평균이있는 ARIMA fit
arima.order
arima.fit

arima.forecast <- fitted(Arima(forecast_y, model=arima.fit))

ggplot(data=forecast_y, aes(x=index(forecast_y), y=diff))+
  geom_line(color="#00AFBB", size=0.5)+
  geom_line(data=forecast_y, aes(x=index(forecast_y),y=arima.forecast),color='red',size=0.5)+
  theme_gray()
# doesn't explain large fluctuation well


# check heteroscedasticity
library("aTSA")
arima_fit<-arima(y,order=arima.order)
arch.test(arima_fit)

##### Fit GARCH with normal-dist #####
spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder=c(arima.order[1], arima.order[3]), include.mean=T),
  distribution.model="norm"
)
fit = ugarchfit(spec, y, solver = 'hybrid')
# fitted AR(2,0,3)-GARCH(1,1)
coef(fit)
plot(fit)

#normality test - Jarque-Beta test
library(tseries)
jarque.bera.test(fit@fit$residuals/fit@fit$sigma)


#Fit Garch again with t-dist
spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder=c(arima.order[1], arima.order[3]), include.mean=T),
  distribution.model="std"
)
fit = ugarchfit(spec, y, solver = 'hybrid')
# fitted AR(2,0,3)-GARCH(1,1)
coef(fit)

# plot res and standardized_res
res <- residuals(fit)
row.names(res) <- row.names(y)
p_diff <- ggplot(data=res, aes(x=Index,y=res[,1]))+
  geom_line(color="#00AFBB", size=0.5)+
  theme_gray()+
  xlab("date")+ylab("residual")
p_stdiff <- ggplot(data=res, aes(x=Index,y=res[,1]/sigma(fit)))+
  geom_line(color="#E7B800", size=0.5)+
  theme_gray()+
  xlab("date")+ylab("standardized_res")
plot_grid(p_diff, p_stdiff, align = "v", ncol=1)


##### Change-point Detection #####
cov_cusum <- function(dat, h) {
  n <- length(dat)
  xmean <- mean(dat)
  x1 <- dat[1:(n-h)] - xmean
  x2 <- dat[(h+1):n] - xmean
  return(sum(x1 * x2/n))
}

sd_cusum <- function(dat, max_h=length(dat)^(1/3)) {
  n <- length(dat)
  sd2_hat <- cov_cusum(dat,0)
  for(i in 1:max_h) {
    sd2_hat <- sd2_hat + 2*(cov_cusum(dat,i))
  }
  return(sd2_hat)
}

CUSUM_calc <- function(dat) {
  ## return : maximum cusum value and maximum point
  n <- length(dat)
  cusum <- abs((cumsum(dat) - (1:n/n)*sum(dat) ) / ( sqrt(n) * sqrt(sd_cusum(dat))))
  argmax <- which.max(cusum)
  if(max(cusum)>1.345) return(list("CUSUM_statistics" = max(cusum), "change_point"=argmax))
  else return(print("no change"))
}

CUSUM_calc(fit@fit$residuals/fit@fit$sigma)  #No change point

#p_raw +
#  geom_vline(xintercept = as.numeric(coin_raw$date[5500]))+
#  geom_text(aes(x=date[5500], label="\nExpected ChangePoint", y=2800),angle=90)+
#  geom_vline(xintercept = as.numeric(coin_raw$date[5000]), color="red")+
#  geom_vline(xintercept = as.numeric(coin_raw$date[6000]), color="red")+
#  geom_vline(xintercept = as.numeric(coin_raw$date[4500]), color="orange")+
#  geom_vline(xintercept = as.numeric(coin_raw$date[6500]), color="orange")+
#  geom_vline(xintercept = as.numeric(coin_raw$date[3000]), color="yellow")+
#  geom_vline(xintercept = as.numeric(coin_raw$date[7000]), color="yellow")