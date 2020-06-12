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
library(openxlsx)
library(rio)



dat <- read_excel("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/R Implementation of Local Projection Methods/annual_alternative_sentiment.xls")



dat


datnew = dat[seq(1, nrow(dat), 4), ]


datnew

export(datnew, "datnew.xlsx")








