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
ds_theme_set()
library(tikzDevice)
require(tikzDevice)



# Load necessary libraries
library(lpirfs)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(plm)
library(Hmisc)






# Load packages for creating plots
library(gridExtra)
library(ggpubr)
#--- Begin code for panel data
# Load libraries to download and read excel file from the website
library(httr)
library(readxl)



dat <- read_excel("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/R Implementation of Local Projection Methods/master_data_quarterly.xls")





exog_names <- c("INCOME", "EDUC", "INVEST", "YGAP")
endog_names <- c("GGDP", "GOOD")

















dat <- dat %>% 
  group_by(State) %>% 
  mutate(meanGDP_bystate = mean(GGDP)) %>% ungroup()






dat <- subset(dat, State!="NA")


lagpad <- function(x, k) {
  if (k>0) {
    return (c(rep(NA, k), x)[1 : length(x)] );
  }
  else {
    return (c(x[(-k+1) : length(x)], rep(NA, -k)));
  }
}

















firststage <- plm(GOOD ~ Congpres + INCOME + EDUC + INVEST + YGAP, data = dat, model = "within")




coeff1 <- summary(firststage)$coefficients[1, 1]

coeff2 <- summary(firststage)$coefficients[2, 1]


coeff3 <- summary(firststage)$coefficients[3, 1]


coeff4 <- summary(firststage)$coefficients[4, 1]

coeff5 <- summary(firststage)$coefficients[5, 1]


dat$predicted <- coeff1*dat$Congpres + coeff2*dat$INCOME + coeff3*dat$EDUC
+ coeff4*dat$INVEST + coeff5*dat$YGAP


plot(dat$predicted)

dat


m <- dat %>% 
  group_by(time) %>%
  dplyr::summarize(Mean = mean(predicted, na.rm=TRUE)) %>% ungroup()


m



#m <- ts(m, start = c(2004,1), frequency = 4) 

m

recession <- read_csv("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/NBERrecessioncut.csv")

recession <- recession[1:52, 2]
m <- cbind(m, recession)





#m <- ts(m, start = c(2004,1), frequency = 4)





m <- m[, c(2,3)]

m

date = seq(from = as.Date("2004-01-01"), to = as.Date("2016-10-01"), by = 'quarter')
date

m <- cbind(m, date)

m



#dat <- data.frame(data$Goodgr, data$PCEmotorgr)



#tikz(file = "quarterlylpshockseries.tex", width = 6, height = 3.7)

#shockseriesquarterly <- ggplot(m) + geom_line(aes(x = time, y = Mean), color = "red") +ggtitle( paste("Quarterly CS Shock Series Averaged Across States")) +
  #scale_y_continuous(limits = c(0.05, 0.08)) +
  #labs( x = "Time", y = "Shock") +theme_bw()

#This line is only necessary if you want to preview the plot right after compiling
#print(shockseriesquarterly)
#Necessary to close or the tikxDevice .tex file will not be written
#dev.off()


m <- m[-c(49:52),]
m


gdpgrowth <- read_csv("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/gdpgrowth2004q1.csv")

gdpgrowth


m <- cbind(m, gdpgrowth)


m
dummy <- m %>% 
  filter(dummy == 1)

tikz(file = "quarterlyshockseries.tex", width = 6, height = 3.7)


r <- ggplot(m) + geom_rect(data = dummy, aes(ymin = -Inf, ymax = Inf, xmin = date-180, xmax = date+180), 
                                             fill='gray', alpha = 0.2) + 
  geom_line(aes(x=date, y=Mean), color = "red") +
  scale_y_continuous(limits = c(0.05, 0.08)) +
  labs(x = "Date", y = "Consumer Sentiment Shock", title = "Evolution of CS Shock and Recession Window")



print(r)

dev.off()




#ccf(m$Mean, m$gdpgrowth)





