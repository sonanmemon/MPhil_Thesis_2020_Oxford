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



dat <- read_excel("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/R Implementation of Local Projection Methods/Annual_Data_Benhabib_Speigel2018.xls")




exog_names <- c("income", "Educ", "INVEST", "OutputGap")
endog_names <- c("GGDP", "GOOD")


#












# Show irfs

#plot(results_panel_nl)















dat <- dat %>% 
  group_by(state) %>% 
  mutate(meanGDP_bystate = mean(GGDP)) %>% ungroup()



dat <- dat %>% 
  group_by(state) %>% 
  mutate(abovemean = ifelse(GGDP > meanGDP_bystate, 1, 0)) %>% ungroup()


#dat <- dat %>% group_by(state) %>% mutate(Diff_durable = PCE_durables - lag(PCE_durables),lag_durable = lag(PCE_durables),# Difference in durable spending between years durables_gr = (((Diff_durable)/lag(PCE_durables)) * 1))  # growth rate

dat <- subset(dat, state!="NA")


lagpad <- function(x, k) {
  if (k>0) {
    return (c(rep(NA, k), x)[1 : length(x)] );
  }
  else {
    return (c(x[(-k+1) : length(x)], rep(NA, -k)));
  }
}

dat <- dat %>% group_by(state) %>% mutate(ldurables = lagpad(PCE_durables, 1)) %>% ungroup()

dat$ldurables

dat <- dat %>% group_by(state) %>% 
  mutate(Diff_durable = PCE_durables - ldurables,# Difference in durable spending between years 
         durables_gr = (((Diff_durable)/ldurables) * 1)) %>% ungroup()  # growth rate




dat$durables_gr

# using subset function
#datupper <- subset(dat, abovemean == 1, select = c("state", "year", "abovemean", "income", "OutputGap", "INVEST", "Educ", "GGDP", "GOOD", "Congpres"))

#datlower <- subset(dat, abovemean == 0, select = c("state", "year", "abovemean", "income", "OutputGap", "INVEST", "Educ", "GGDP", "GOOD", "Congpres"))




#results_lin_panel_upper <- lp_lin_panel(data_set = datupper, data_sample = "Full",
#                                endog_data = "GGDP", cumul_mult = TRUE,
#                                shock = "GOOD",
#                                diff_shock = FALSE,
#   iv_reg = TRUE, instrum = "Congpres",
#                               panel_model = "within", panel_effect = "individual",
#                               robust_cov = "vcovSCC",
#                              c_exog_data = exog_names,
#                             lags_exog_data = 2,
#                           confint = 1.67, hor = 5)














#datselected <- subset(dat, select = c("state", "year", "income", "OutputGap", "INVEST", "Educ", "GGDP", "GOOD", "Congpres")) # using subset function







#datselected <- subset(datselected, state!="NA") # using subset function








firststage <- plm(GOOD ~ Congpres + income + Educ + INVEST + OutputGap, data = dat, model = "within")

summary(firststage)
predictfirststage <- predict(firststage)


coeff1 <- summary(firststage)$coefficients[1, 1]

coeff2 <- summary(firststage)$coefficients[2, 1]


coeff3 <- summary(firststage)$coefficients[3, 1]


coeff4 <- summary(firststage)$coefficients[4, 1]

coeff5 <- summary(firststage)$coefficients[5, 1]


dat$predicted <- coeff1*dat$Congpres + coeff2*dat$income + coeff3*dat$Educ
+ coeff4*dat$INVEST + coeff5*dat$OutputGap













plot(dat$predicted)

dat


m <- dat %>% 
  group_by(year) %>%
  dplyr::summarize(Mean = mean(predicted, na.rm=TRUE)) %>% ungroup()


m



#m <- ts(m, start = c(2004,1), frequency = 4) 

m

recession <- read_csv("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/NBERrecessioncut.csv")

recession

recession <- recession[-seq(4, NROW(recession), by = 4),]

recession <- recession[-seq(3, NROW(recession), by = 3),]
recession <- recession[-seq(2, NROW(recession), by = 2),]
recession <- recession[1:12,2]
recession

m <- m[1:12,]

m_even = m[seq(2, nrow(m), 2), ]
recession_even = recession[seq(2, nrow(recession), 2), ]

m <- cbind(m, recession)


m


m_even <- cbind(m_even, recession_even)

m_even













dummy <- m %>% 
  filter(dummy == 1)

tikz(file = "annualshockseries.tex", width = 6, height = 3.7)


r <- ggplot(m) + geom_rect(data = dummy, aes(ymin = -Inf, ymax = Inf, xmin = year-18000000, xmax = year+18000000), 
                           fill='gray', alpha = 0.2) + 
  geom_line(aes(x=year, y=Mean), color = "red") +
  scale_y_continuous(limits = c(0.10, 0.15)) +
  labs(x = "Date", y = "Consumer Sentiment Shock", title = "Evolution of CS Shock and Recession Window")



print(r)

dev.off()



dummy <- m_even %>% 
  filter(dummy == 1)

tikz(file = "biennialshockseries.tex", width = 6, height = 3.7)


r_even <- ggplot(m_even) + geom_rect(data = dummy, aes(ymin = -Inf, ymax = Inf, xmin = year-18000000, xmax = year+18000000), 
                           fill='gray', alpha = 0.2) + 
  geom_line(aes(x=year, y=Mean), color = "red") +
  scale_y_continuous(limits = c(0.10, 0.15)) +
  labs(x = "Date", y = "Consumer Sentiment Shock", title = "Evolution of CS Shock and Recession Window")



print(r_even)

dev.off()


#ccf(m$Mean, m$gdpgrowth)





