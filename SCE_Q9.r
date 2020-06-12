library(dslabs)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Lahman)
library(dslabs)
library(AER)
library(tseries)
library(dynlm)
library(stargazer)
library(forecast)
library("readxl")
install.packages("xlsx")
library("xlsx")

# xls files





data <- read_csv("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/SCE_Q9.csv")


data$IQR_Q912 <- as.numeric(as.character(data$IQR_Q912))

data$IQR_Q924 <- as.numeric(as.character(data$IQR_Q924))
#data

#data <- ts(data, start = c(1952,3), frequency = 3)
#head(data)

#data$IQR24minus12 <- data$IQR_Q924 - data$IQR_Q912


data$mean12 <- mean(data$IQR_Q912, na.rm = TRUE)

data$mean24 <- mean(data$IQR_Q924, na.rm = TRUE)

data$mean24

data$mean12

data$IQR24minus12percentage <- ((data$IQR_Q924 - data$IQR_Q912)/(data$IQR_Q924))*100



# replaces Inf by NA Values
is.na(data) <- sapply(data, is.infinite)






data$lambda <- mean(data$IQR24minus12percentage, na.rm = TRUE)


data$Q9_bin1 <-  data$Q9_bin1/100
data$Q9_bin2 <-  data$Q9_bin2/100
data$Q9_bin3 <-  data$Q9_bin3/100
data$Q9_bin4 <-  data$Q9_bin4/100
data$Q9_bin5 <-  data$Q9_bin5/100
data$Q9_bin6 <-  data$Q9_bin6/100
data$Q9_bin7 <-  data$Q9_bin7/100
data$Q9_bin8 <-  data$Q9_bin8/100
data$Q9_bin9 <-  data$Q9_bin9/100
data$Q9_bin10 <-  data$Q9_bin10/100



data$Q9c_bin1 <-  data$Q9c_bin1/100
data$Q9c_bin2 <-  data$Q9c_bin2/100
data$Q9c_bin3 <-  data$Q9c_bin3/100
data$Q9c_bin4 <-  data$Q9c_bin4/100
data$Q9c_bin5 <-  data$Q9c_bin5/100
data$Q9c_bin6 <-  data$Q9c_bin6/100
data$Q9c_bin7 <-  data$Q9c_bin7/100
data$Q9c_bin8 <-  data$Q9c_bin8/100
data$Q9c_bin9 <-  data$Q9c_bin9/100
data$Q9c_bin10 <-  data$Q9c_bin10/100


data$P1 <- data$Q9_bin1

data$P2 <- data$P1 + data$Q9_bin2

data$P3 <- data$P2 + data$Q9_bin3


data$P4 <- data$P3 + data$Q9_bin4


data$P5 <- data$P4 + data$Q9_bin5



data$P6 <- data$P5 + data$Q9_bin6

data$P7 <- data$P6 + data$Q9_bin7




data$P8 <- data$P7 + data$Q9_bin8



data$P9 <- data$P8 + data$Q9_bin9

data$P10 <- data$P9 + data$Q9_bin10



data
data$P1c <- data$Q9c_bin1

data$P2c <- data$P1c + data$Q9c_bin2

data$P3c <- data$P2c + data$Q9c_bin3


data$P4c <- data$P3c + data$Q9c_bin4


data$P5c <- data$P4c + data$Q9c_bin5



data$P6c <- data$P5c + data$Q9c_bin6

data$P7c <- data$P6c + data$Q9c_bin7




data$P8c <- data$P7c + data$Q9c_bin8



data$P9c <- data$P8c + data$Q9c_bin9

data$P10c <- data$P9c + data$Q9c_bin10




data$ERPS12 <- data$P1*(1-data$P1) + data$P2*(1-data$P2) + data$P3*(1-data$P3) + data$P4*(1-data$P4) + data$P5*(1-data$P5) + data$P6*(1-data$P6) + data$P7*(1-data$P7) + data$P8*(1-data$P8) + data$P9*(1-data$P9) + data$P10*(1-data$P10)


data$ERPS24 <- data$P1c*(1-data$P1c) + data$P2c*(1-data$P2c) + data$P3c*(1-data$P3c) + data$P4c*(1-data$P4c) + data$P5c*(1-data$P5c) + data$P6c*(1-data$P6c) + data$P7c*(1-data$P7c) + data$P8c*(1-data$P8c) + data$P9c*(1-data$P9c) + data$P10c*(1-data$P10c)



data$meanERPS12 <- mean(data$ERPS12, na.rm = TRUE)

data$meanERPS24 <- mean(data$ERPS24, na.rm = TRUE)




data$meanERPS12 

data$meanERPS24








