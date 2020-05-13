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




results_panel_nl <- lp_nl_panel(data_set = dat,
                                endog_data = "GGDP",
                                cumul_mult = FALSE,
                                shock = "predicted",
                                c_exog_data = exog_names,
                                lags_exog_data = 2,
                                diff_shock = FALSE,
                                panel_model = "within", panel_effect = "individual",
                                robust_cov = "vcovSCC", switching = "YGAP",
                                lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                gamma = 10,
                                confint = 1.67,
                                hor = 20)



plot(results_panel_nl)




firststagegood1 <- plm(Good1 ~ Congpres + INCOME + EDUC + INVEST + YGAP, data = dat, model = "within")




coeff1_g1 <- summary(firststagegood1)$coefficients[1, 1]

coeff2_g1 <- summary(firststagegood1)$coefficients[2, 1]


coeff3_g1 <- summary(firststagegood1)$coefficients[3, 1]


coeff4_g1 <- summary(firststagegood1)$coefficients[4, 1]

coeff5_g1 <- summary(firststagegood1)$coefficients[5, 1]


dat$predicted_g1 <- coeff1_g1*dat$Congpres + coeff2_g1*dat$INCOME + coeff3_g1*dat$EDUC
+ coeff4_g1*dat$INVEST + coeff5_g1*dat$YGAP




results_panel_nl_g1 <- lp_nl_panel(data_set = dat,
                                endog_data = "GGDP",
                                cumul_mult = FALSE,
                                shock = "predicted_g1",
                                c_exog_data = exog_names,
                                lags_exog_data = 2,
                                diff_shock = FALSE,
                                panel_model = "within", panel_effect = "individual",
                                robust_cov = "vcovSCC", switching = "YGAP",
                                lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                gamma = 10,
                                confint = 1.67,
                                hor = 20)



plot(results_panel_nl_g1)






firststagebetter1 <- plm(Better1 ~ Congpres + INCOME + EDUC + INVEST + YGAP, data = dat, model = "within")




coeff1_b1 <- summary(firststagebetter1)$coefficients[1, 1]

coeff2_b1 <- summary(firststagebetter1)$coefficients[2, 1]


coeff3_b1 <- summary(firststagebetter1)$coefficients[3, 1]


coeff4_b1 <- summary(firststagebetter1)$coefficients[4, 1]

coeff5_b1 <- summary(firststagebetter1)$coefficients[5, 1]


dat$predicted_b1 <- coeff1_b1*dat$Congpres + coeff2_b1*dat$INCOME + coeff3_b1*dat$EDUC
+ coeff4_b1*dat$INVEST + coeff5_b1*dat$YGAP




results_panel_nl_b1 <- lp_nl_panel(data_set = dat,
                                   endog_data = "GGDP",
                                   cumul_mult = FALSE,
                                   shock = "predicted_b1",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "YGAP",
                                   lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 20)



plot(results_panel_nl_b1)





firststageb5 <- plm(BAD5 ~ Congpres + INCOME + EDUC + INVEST + YGAP, data = dat, model = "within")




coeff1_b5 <- summary(firststageb5)$coefficients[1, 1]

coeff2_b5 <- summary(firststageb5)$coefficients[2, 1]


coeff3_b5 <- summary(firststageb5)$coefficients[3, 1]


coeff4_b5 <- summary(firststageb5)$coefficients[4, 1]

coeff5_b5 <- summary(firststageb5)$coefficients[5, 1]


dat$predicted_b5 <- coeff1_b5*dat$Congpres + coeff2_b5*dat$INCOME + coeff3_b5*dat$EDUC
+ coeff4_b5*dat$INVEST + coeff5_b5*dat$YGAP




results_panel_nl_b5 <- lp_nl_panel(data_set = dat,
                                   endog_data = "GGDP",
                                   cumul_mult = FALSE,
                                   shock = "predicted_b5",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "YGAP",
                                   lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 20)



plot(results_panel_nl_b5)
























results_panel_nl_c <- lp_nl_panel(data_set = dat,
                                  endog_data = "pce_growth",
                                  cumul_mult = FALSE,
                                  shock = "predicted",
                                  c_exog_data = exog_names,
                                  lags_exog_data = 2,
                                  diff_shock = FALSE,
                                  panel_model = "within", panel_effect = "individual",
                                  robust_cov = "vcovSCC", switching = "YGAP",
                                  lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                  gamma = 10,
                                  confint = 1.67,
                                  hor = 20)





plot(results_panel_nl_c)









results_panel_nl_c_g1 <- lp_nl_panel(data_set = dat,
                                   endog_data = "pce_growth",
                                   cumul_mult = FALSE,
                                   shock = "predicted_g1",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "YGAP",
                                   lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 20)



plot(results_panel_nl_c_g1)




results_panel_nl_c_b1 <- lp_nl_panel(data_set = dat,
                                   endog_data = "pce_growth",
                                   cumul_mult = FALSE,
                                   shock = "predicted_b1",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "YGAP",
                                   lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 20)



plot(results_panel_nl_c_b1)






results_panel_nl_c_b5 <- lp_nl_panel(data_set = dat,
                                   endog_data = "pce_growth",
                                   cumul_mult = FALSE,
                                   shock = "predicted_b5",
                                   c_exog_data = exog_names,
                                   lags_exog_data = 2,
                                   diff_shock = FALSE,
                                   panel_model = "within", panel_effect = "individual",
                                   robust_cov = "vcovSCC", switching = "YGAP",
                                   lag_switching = TRUE, use_hp = TRUE, lambda = 1600,
                                   gamma = 10,
                                   confint = 1.67,
                                   hor = 20)



plot(results_panel_nl_c_b5)













