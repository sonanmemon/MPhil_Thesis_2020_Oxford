library(dslabs)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Lahman)
library(dslabs)
library(AER)
library(tseries)
library(dynlm)library(stargazer)
library(forecast)
library("readxl")
library("ggplot2")
install.packages("xlsx")
library("xlsx")
library(ggthemes)
library(tikzDevice)
require(tikzDevice)
install.packages("magrittr") # package installations are only needed the first time you use it
#install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) 

# xls files





data <- read_csv("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/data_masterfile.csv")

data <- ts(data, start = c(1967, 1), frequency = 1)

data



dat <- data[-c(1:24),]





dat <- dat[, c(1,2,4,6)]








lagNewVehicleUnitSales <- lag(dat[,"NewVehicleUnitSales"], k = 1)
dat <- cbind(dat,lagNewVehicleUnitSales)

lagUsedVehicleUnitSales <- lag(dat[,"UsedVehicleUnitSales"], k = 1)

dat <- cbind(dat,lagUsedVehicleUnitSales)


changeUsedVehicleUnitSales <- dat[,"UsedVehicleUnitSales"] - dat[,"lagUsedVehicleUnitSales"]

dat <- cbind(dat,changeUsedVehicleUnitSales)



changeNewVehicleUnitSales <- dat[,"NewVehicleUnitSales"] - dat[,"lagNewVehicleUnitSales"]

dat <- cbind(dat,changeNewVehicleUnitSales)


dat




dat <- data.frame(dat[,"Date"], dat[, "NBER"], dat[,"changeNewVehicleUnitSales"], dat[,"changeUsedVehicleUnitSales"])
dat

dummy <- dat %>% 
  filter(dat....NBER.. == 1)

tikz(file = "newusedvehicles.tex", width = 6, height = 3.7)



r <- ggplot(dat) + geom_rect(data = dummy, 
                                             aes(ymin = -Inf, ymax = Inf, xmin = dat....Date..-0.5, xmax = dat....Date..+0.5), 
                                             fill='gray', alpha = 0.5) + 
  geom_line(aes(x=dat....Date.., y=dat....changeUsedVehicleUnitSales..), color = "red") +
  geom_line(aes(x=dat....Date.., y=dat....changeNewVehicleUnitSales..), color = "blue") +
  scale_y_continuous(breaks = c(-1500, 0, -3000)) +
  labs(x = "Date", y = "Unit Change in Vehicle Sales From Last Year", title = "New and Used Vehicles") + 
  theme(panel.background = element_rect(fill = "azure", colour = "azure", size = 0.5, linetype = "solid"))


print(r)


dev.off()

























