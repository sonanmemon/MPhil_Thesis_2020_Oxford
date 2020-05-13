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



firststage_g1 <- plm(Good1 ~ Congpres + income + Educ + INVEST + OutputGap, data = dat, model = "within")

summary(firststage_g1)


coeff1_g1 <- summary(firststage_g1)$coefficients[1, 1]

coeff2_g1 <- summary(firststage_g1)$coefficients[2, 1]


coeff3_g1 <- summary(firststage_g1)$coefficients[3, 1]


coeff4_g1 <- summary(firststage_g1)$coefficients[4, 1]

coeff5_g1 <- summary(firststage_g1)$coefficients[5, 1]


dat$predicted_g1 <- coeff1_g1*dat$Congpres + coeff2_g1*dat$income + coeff3_g1*dat$Educ
+ coeff4_g1*dat$INVEST + coeff5_g1*dat$OutputGap




firststage_b1 <- plm(Better1 ~ Congpres + income + Educ + INVEST + OutputGap, data = dat, model = "within")

summary(firststage_b1)


coeff1_b1 <- summary(firststage_b1)$coefficients[1, 1]

coeff2_b1 <- summary(firststage_b1)$coefficients[2, 1]


coeff3_b1 <- summary(firststage_b1)$coefficients[3, 1]


coeff4_b1 <- summary(firststage_b1)$coefficients[4, 1]

coeff5_b1 <- summary(firststage_b1)$coefficients[5, 1]


dat$predicted_b1 <- coeff1_b1*dat$Congpres + coeff2_b1*dat$income + coeff3_b1*dat$Educ
+ coeff4_b1*dat$INVEST + coeff5_b1*dat$OutputGap



firststage_b5 <- plm(Bad5 ~ Congpres + income + Educ + INVEST + OutputGap, data = dat, model = "within")

summary(firststage_b5)


coeff1_b5 <- summary(firststage_b5)$coefficients[1, 1]

coeff2_b5 <- summary(firststage_b5)$coefficients[2, 1]


coeff3_b5 <- summary(firststage_b5)$coefficients[3, 1]


coeff4_b5 <- summary(firststage_b5)$coefficients[4, 1]

coeff5_b5 <- summary(firststage_b5)$coefficients[5, 1]


dat$predicted_b5 <- coeff1_b5*dat$Congpres + coeff5_b1*dat$income + coeff5_b1*dat$Educ
+ coeff4_b5*dat$INVEST + coeff5_b5*dat$OutputGap



results_panel_nl <- lp_nl_panel(data_set = dat,
                                endog_data = "GGDP",
                                cumul_mult = FALSE,
                                shock = "predicted",
                                c_exog_data = exog_names,
                                lags_exog_data = 2,
                                diff_shock = FALSE,
                                panel_model = "within", panel_effect = "individual",
                                robust_cov = "vcovSCC", switching = "OutputGap",
                                lag_switching = TRUE, use_hp = TRUE,lambda = 6.25,
                                gamma = 10,
                                confint = 1.96,
                                hor = 7)



plot(results_panel_nl)





results_panel_nl_g1 <- lp_nl_panel(data_set = dat,
                                   endog_data = "GGDP",
                                   cumul_mult = FALSE,
                                   shock = "predicted_g1",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "OutputGap",
                                   lag_switching = TRUE, use_hp = TRUE,lambda = 6.25,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 7)



plot(results_panel_nl_g1)


results_panel_nl_b1 <- lp_nl_panel(data_set = dat,
                                   endog_data = "GGDP",
                                   cumul_mult = FALSE,
                                   shock = "predicted_b1",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "OutputGap",
                                   lag_switching = TRUE, use_hp = TRUE,lambda = 6.25,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 7)



plot(results_panel_nl_b1)



results_panel_nl_b5 <- lp_nl_panel(data_set = dat,
                                   endog_data = "GGDP",
                                   cumul_mult = FALSE,
                                   shock = "predicted_b5",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "OutputGap",
                                   lag_switching = TRUE, use_hp = TRUE,lambda = 6.25,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 7)



plot(results_panel_nl_b5)






dat$index <- with(dat, ave(state, state, FUN = seq_along))

dat



#dat.odd <- dat.odd %>% group_by(state) %>% mutate(index = seq(1,group_size(state), 1)) %>% ungroup()





dat <- dat[,c(1,117, seq(2,116,1))]






results_panel_nl_c <- lp_nl_panel(data_set = dat,
                                  endog_data = "GPCE",
                                  cumul_mult = FALSE,
                                  shock = "predicted",
                                  c_exog_data = exog_names,
                                  lags_exog_data = 2,
                                  diff_shock = FALSE,
                                  panel_model = "within", panel_effect = "individual",
                                  robust_cov = "vcovSCC", switching = "OutputGap",
                                  lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                  gamma = 10,
                                  confint = 1.67,
                                  hor = 7)





plot(results_panel_nl_c)






results_panel_nl_c_g1 <- lp_nl_panel(data_set = dat,
                                     endog_data = "GPCE",
                                     cumul_mult = FALSE,
                                     shock = "predicted_g1",
                                     c_exog_data = exog_names,
                                     lags_exog_data = 2,
                                     diff_shock = FALSE,
                                     panel_model = "within", panel_effect = "individual",
                                     robust_cov = "vcovSCC", switching = "OutputGap",
                                     lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                     gamma = 10,
                                     confint = 1.67,
                                     hor = 7)





plot(results_panel_nl_c_g1)



results_panel_nl_c_b1 <- lp_nl_panel(data_set = dat,
                                     endog_data = "GPCE",
                                     cumul_mult = FALSE,
                                     shock = "predicted_b1",
                                     c_exog_data = exog_names,
                                     lags_exog_data = 2,
                                     diff_shock = FALSE,
                                     panel_model = "within", panel_effect = "individual",
                                     robust_cov = "vcovSCC", switching = "OutputGap",
                                     lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                     gamma = 10,
                                     confint = 1.67,
                                     hor = 5)





plot(results_panel_nl_c_b1)






results_panel_nl_c_b5 <- lp_nl_panel(data_set = dat,
                                     endog_data = "GPCE",
                                     cumul_mult = FALSE,
                                     shock = "predicted_b5",
                                     c_exog_data = exog_names,
                                     lags_exog_data = 2,
                                     diff_shock = FALSE,
                                     panel_model = "within", panel_effect = "individual",
                                     robust_cov = "vcovSCC", switching = "OutputGap",
                                     lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                     gamma = 10,
                                     confint = 1.67,
                                     hor = 7)





plot(results_panel_nl_c_b5)














results_panel_nl_durable <- lp_nl_panel(data_set = dat,
                                        endog_data = "durables_gr",
                                        cumul_mult = FALSE,
                                        shock = "predicted",
                                        c_exog_data = exog_names,
                                        lags_exog_data = 2,
                                        diff_shock = FALSE,
                                        panel_model = "within", panel_effect = "individual",
                                        robust_cov = "vcovSCC", switching = "OutputGap",
                                        lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                        gamma = 10,
                                        confint = 1.67,
                                        hor = 7)





plot(results_panel_nl_durable)


results_panel_nl_durable_g1 <- lp_nl_panel(data_set = dat,
                                           endog_data = "durables_gr",
                                           cumul_mult = FALSE,
                                           shock = "predicted_g1",
                                           c_exog_data = exog_names,
                                           lags_exog_data = 2,
                                           diff_shock = FALSE,
                                           panel_model = "within", panel_effect = "individual",
                                           robust_cov = "vcovSCC", switching = "OutputGap",
                                           lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                           gamma = 10,
                                           confint = 1.67,
                                           hor = 7)





plot(results_panel_nl_durable_g1)





results_panel_nl_durable_b1 <- lp_nl_panel(data_set = dat,
                                           endog_data = "durables_gr",
                                           cumul_mult = FALSE,
                                           shock = "predicted_b1",
                                           c_exog_data = exog_names,
                                           lags_exog_data = 2,
                                           diff_shock = FALSE,
                                           panel_model = "within", panel_effect = "individual",
                                           robust_cov = "vcovSCC", switching = "OutputGap",
                                           lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                           gamma = 10,
                                           confint = 1.67,
                                           hor = 7)





plot(results_panel_nl_durable_b1)







results_panel_nl_durable_b5 <- lp_nl_panel(data_set = dat,
                                           endog_data = "durables_gr",
                                           cumul_mult = FALSE,
                                           shock = "predicted_b5",
                                           c_exog_data = exog_names,
                                           lags_exog_data = 2,
                                           diff_shock = FALSE,
                                           panel_model = "within", panel_effect = "individual",
                                           robust_cov = "vcovSCC", switching = "OutputGap",
                                           lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                           gamma = 10,
                                           confint = 1.67,
                                           hor = 7)





plot(results_panel_nl_durable_b5)








dat <- subset(dat, state!="NA") # using subset function











dat_even = dat[seq(2, nrow(dat), 2), ]


dat_even$index <- with(dat_even, ave(state, state, FUN = seq_along))






















firststage_even <- plm(GOOD ~ Congpres + income + Educ + INVEST + OutputGap, data = dat_even, model = "within")



coeff1_even <- summary(firststage_even)$coefficients[1, 1]

coeff2_even <- summary(firststage_even)$coefficients[2, 1]


coeff3_even <- summary(firststage_even)$coefficients[3, 1]


coeff4_even <- summary(firststage_even)$coefficients[4, 1]

coeff5_even  <- summary(firststage_even)$coefficients[5, 1]


dat_even$predicted_even <- coeff1_even*dat_even$Congpres + coeff2_even*dat_even$income + coeff3_even*dat_even$Educ
+ coeff4_even*dat_even$INVEST + coeff5_even*dat_even$OutputGap



firststage_g1_even <- plm(Good1 ~ Congpres + income + Educ + INVEST + OutputGap, data = dat_even, model = "within")

summary(firststage_g1)


coeff1_g1_even <- summary(firststage_g1_even)$coefficients[1, 1]

coeff2_g1_even <- summary(firststage_g1_even)$coefficients[2, 1]


coeff3_g1_even <- summary(firststage_g1_even)$coefficients[3, 1]


coeff4_g1_even <- summary(firststage_g1_even)$coefficients[4, 1]

coeff5_g1_even <- summary(firststage_g1_even)$coefficients[5, 1]


dat_even$predicted_g1_even <- coeff1_g1_even*dat_even$Congpres + coeff2_g1_even*dat_even$income + coeff3_g1_even*dat_even$Educ
+ coeff4_g1_even*dat_even$INVEST + coeff5_g1_even*dat_even$OutputGap




firststage_b1_even <- plm(Better1 ~ Congpres + income + Educ + INVEST + OutputGap, data = dat_even, model = "within")




coeff1_b1_even <- summary(firststage_b1_even)$coefficients[1, 1]

coeff2_b1_even <- summary(firststage_b1_even)$coefficients[2, 1]


coeff3_b1_even <- summary(firststage_b1_even)$coefficients[3, 1]


coeff4_b1_even <- summary(firststage_b1_even)$coefficients[4, 1]

coeff5_b1_even <- summary(firststage_b1_even)$coefficients[5, 1]


dat_even$predicted_b1_even <- coeff1_b1_even*dat_even$Congpres + coeff2_b1_even*dat_even$income + coeff3_b1_even*dat_even$Educ
+ coeff4_b1_even*dat_even$INVEST + coeff5_b1_even*dat_even$OutputGap



firststage_b5_even <- plm(Bad5 ~ Congpres + income + Educ + INVEST + OutputGap, data = dat_even, model = "within")

summary(firststage_b5)


coeff1_b5_even <- summary(firststage_b5_even)$coefficients[1, 1]

coeff2_b5_even <- summary(firststage_b5_even)$coefficients[2, 1]


coeff3_b5_even <- summary(firststage_b5_even)$coefficients[3, 1]


coeff4_b5_even <- summary(firststage_b5_even)$coefficients[4, 1]

coeff5_b5_even <- summary(firststage_b5_even)$coefficients[5, 1]


dat_even$predicted_b5_even <- coeff1_b5_even*dat_even$Congpres + coeff2_b5_even*dat_even$income + coeff3_b5_even*dat_even$Educ
+ coeff4_b5_even*dat_even$INVEST + coeff5_b5_even*dat_even$OutputGap




results_panel_nl_gdp_even <- lp_nl_panel(data_set = dat_even,
                                         endog_data = "GGDP",
                                         cumul_mult = FALSE,
                                         shock = "predicted_even",
                                         c_exog_data = exog_names,
                                         lags_exog_data = 2,
                                         diff_shock = FALSE,
                                         panel_model = "within", panel_effect = "individual",
                                         robust_cov = "vcovSCC", switching = "OutputGap",
                                         lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                         gamma = 10,
                                         confint = 1.67,
                                         hor = 3)



plot(results_panel_nl_gdp_even)







results_panel_nl_gdp_g1_even <- lp_nl_panel(data_set = dat_even,
                                            endog_data = "GGDP",
                                            cumul_mult = FALSE,
                                            shock = "predicted_g1_even",
                                            c_exog_data = exog_names,
                                            lags_exog_data = 2,
                                            diff_shock = FALSE,
                                            panel_model = "within", panel_effect = "individual",
                                            robust_cov = "vcovSCC", switching = "OutputGap",
                                            lag_switching = TRUE, use_hp = TRUE,lambda = 6.25,
                                            gamma = 10,
                                            confint = 1.67,
                                            hor = 3)



plot(results_panel_nl_gdp_g1_even)



results_panel_nl_gdp_b1_even <- lp_nl_panel(data_set = dat_even,
                                            endog_data = "GGDP",
                                            cumul_mult = FALSE,
                                            shock = "predicted_b1_even",
                                            c_exog_data = exog_names,
                                            lags_exog_data = 2,
                                            diff_shock = FALSE,
                                            panel_model = "within", panel_effect = "individual",
                                            robust_cov = "vcovSCC", switching = "OutputGap",
                                            lag_switching = TRUE, use_hp = TRUE,lambda = 6.25,
                                            gamma = 10,
                                            confint = 1.67,
                                            hor = 3)



plot(results_panel_nl_gdp_b1_even)



results_panel_nl_gdp_b5_even <- lp_nl_panel(data_set = dat_even,
                                            endog_data = "GGDP",
                                            cumul_mult = FALSE,
                                            shock = "predicted_b5_even",
                                            c_exog_data = exog_names,
                                            lags_exog_data = 2,
                                            diff_shock = FALSE,
                                            panel_model = "within", panel_effect = "individual",
                                            robust_cov = "vcovSCC", switching = "OutputGap",
                                            lag_switching = TRUE, use_hp = TRUE,lambda = 6.25,
                                            gamma = 10,
                                            confint = 1.67,
                                            hor = 3)



plot(results_panel_nl_gdp_b5_even)





results_panel_nl_c_even <- lp_nl_panel(data_set = dat_even,
                                       endog_data = "GPCE",
                                       cumul_mult = FALSE,
                                       shock = "predicted_even",
                                       c_exog_data = exog_names,
                                       lags_exog_data = 2,
                                       diff_shock = FALSE,
                                       panel_model = "within", panel_effect = "individual",
                                       robust_cov = "vcovSCC", switching = "OutputGap",
                                       lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                       gamma = 10,
                                       confint = 1.67,
                                       hor = 3)



plot(results_panel_nl_c_even)


results_panel_nlg1_c_even <- lp_nl_panel(data_set = dat_even,
                                         endog_data = "GPCE",
                                         cumul_mult = FALSE,
                                         shock = "predicted_g1_even",
                                         c_exog_data = exog_names,
                                         lags_exog_data = 2,
                                         diff_shock = FALSE,
                                         panel_model = "within", panel_effect = "individual",
                                         robust_cov = "vcovSCC", switching = "OutputGap",
                                         lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                         gamma = 10,
                                         confint = 1.67,
                                         hor = 3)



plot(results_panel_nlg1_c_even)




results_panel_nlb1_c_even <- lp_nl_panel(data_set = dat_even,
                                         endog_data = "GPCE",
                                         cumul_mult = FALSE,
                                         shock = "predicted_b1_even",
                                         c_exog_data = exog_names,
                                         lags_exog_data = 2,
                                         diff_shock = FALSE,
                                         panel_model = "within", panel_effect = "individual",
                                         robust_cov = "vcovSCC", switching = "OutputGap",
                                         lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                         gamma = 10,
                                         confint = 1.67,
                                         hor = 3)



plot(results_panel_nlb1_c_even)



results_panel_nlb5_c_even <- lp_nl_panel(data_set = dat_even,
                                         endog_data = "GPCE",
                                         cumul_mult = FALSE,
                                         shock = "predicted_b5_even",
                                         c_exog_data = exog_names,
                                         lags_exog_data = 2,
                                         diff_shock = FALSE,
                                         panel_model = "within", panel_effect = "individual",
                                         robust_cov = "vcovSCC", switching = "OutputGap",
                                         lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                         gamma = 10,
                                         confint = 1.67,
                                         hor = 3)



plot(results_panel_nlb5_c_even)




results_panel_nl_durables_even <- lp_nl_panel(data_set = dat_even,
                                              endog_data = "durables_gr",
                                              cumul_mult = FALSE,
                                              shock = "predicted_even",
                                              c_exog_data = exog_names,
                                              lags_exog_data = 2,
                                              diff_shock = FALSE,
                                              panel_model = "within", panel_effect = "individual",
                                              robust_cov = "vcovSCC", switching = "OutputGap",
                                              lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                              gamma = 10,
                                              confint = 1.67,
                                              hor = 3)



plot(results_panel_nl_durables_even)








results_panel_nlg1_durables_even <- lp_nl_panel(data_set = dat_even,
                                                endog_data = "durables_gr",
                                                cumul_mult = FALSE,
                                                shock = "predicted_g1_even",
                                                c_exog_data = exog_names,
                                                lags_exog_data = 2,
                                                diff_shock = FALSE,
                                                panel_model = "within", panel_effect = "individual",
                                                robust_cov = "vcovSCC", switching = "OutputGap",
                                                lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                                gamma = 10,
                                                confint = 1.67,
                                                hor = 3)



plot(results_panel_nlg1_durables_even)




results_panel_nlb1_durables_even <- lp_nl_panel(data_set = dat_even,
                                                endog_data = "durables_gr",
                                                cumul_mult = FALSE,
                                                shock = "predicted_b1_even",
                                                c_exog_data = exog_names,
                                                lags_exog_data = 2,
                                                diff_shock = FALSE,
                                                panel_model = "within", panel_effect = "individual",
                                                robust_cov = "vcovSCC", switching = "OutputGap",
                                                lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                                gamma = 10,
                                                confint = 1.67,
                                                hor = 3)



plot(results_panel_nlb1_durables_even)





results_panel_nlb5_durables_even <- lp_nl_panel(data_set = dat_even,
                                                endog_data = "durables_gr",
                                                cumul_mult = FALSE,
                                                shock = "predicted_b5_even",
                                                c_exog_data = exog_names,
                                                lags_exog_data = 2,
                                                diff_shock = FALSE,
                                                panel_model = "within", panel_effect = "individual",
                                                robust_cov = "vcovSCC", switching = "OutputGap",
                                                lag_switching = TRUE, use_hp = TRUE, lambda = 6.25,
                                                gamma = 10,
                                                confint = 1.67,
                                                hor = 3)



plot(results_panel_nlb5_durables_even)







