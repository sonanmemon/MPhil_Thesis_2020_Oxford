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





dat <- read_csv("C:/Users/ND.COM/Desktop/michigan_sub.csv")


dat






DATE <- seq(from = as.Date("1978-01-01"), to = as.Date("2019-10-01"), by = 'month')




michigan <- cbind(dat, DATE)

michigan <- michigan %>% 
  rename(
    icshs = ics_ehs,
    icssc = ics_esc,
    icscd = ics_ecd
  )

michigan <- michigan %>% 
  rename(
    icsnc = ics_nc,
    icsne = ics_ne,
    icss = ics_s,
    icsw = ics_w
  )


michigan <- michigan %>% 
  rename(
    icsa1834 = ics_a1834,
    icsa3554 = ics_a3554,
    icsa5597 = ics_a5597)

michigan <- michigan %>% 
  rename(
    icsbottom = ics_y13,
    icsmiddle = ics_y23,
    icstop = ics_y33)



michiganincomegroup <- michigan %>% gather(key = Series, value = Value, icsbottom, icsmiddle, icstop)

michiganagegroup <- michigan %>% gather(key = Series, value = Value, icsa1834, icsa3554, icsa5597)


michiganregion <- michigan %>% gather(key = Series, value = Value, icsnc, icsne, icss, icsw)


michiganeduc <- michigan %>% gather(key = Series, value = Value, icshs, icssc, icscd)



"Sentiment/Growth Rate (Last Year Qr)", title = "Evolution of Consumer Sentiment and GDP (Michigan)")




#michigan <- ggplot(michigan) + geom_line(aes(x = DATE, y = Value, linetype = Series, color = Series, group = Series))




michigan

tikz(file = "ICSbyeduc.tex", width = 6, height = 3.7)
#Simple plot of the dummy data using LaTeX elements







michiganeduc <- ggplot(michiganeduc) + geom_line(aes(x = DATE, y = Value, linetype = Series, color = Series, group = Series))+ggtitle( paste("Evolution of ICS By Education")) +labs( x = "Date", y = "Consumer Sentiment Index") +theme_bw()


#This line is only necessary if you want to preview the plot right after compiling
print(michiganeduc)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()


tikz(file = "ICSbyregion.tex", width = 6, height = 3.7)
#Simple plot of the dummy data using LaTeX elements







michiganregion <- ggplot(michiganregion) + geom_line(aes(x = DATE, y = Value, linetype = Series, color = Series, group = Series))+ggtitle( paste("Evolution of ICS By Region")) +labs( x = "Date", y = "Consumer Sentiment Index") +theme_bw()


#This line is only necessary if you want to preview the plot right after compiling
print(michiganregion)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()



tikz(file = "ICSbyagegroup.tex", width = 6, height = 3.7)
#Simple plot of the dummy data using LaTeX elements







michiganagegroup <- ggplot(michiganagegroup) + geom_line(aes(x = DATE, y = Value, linetype = Series, color = Series, group = Series))+ggtitle( paste("Evolution of ICS By Age Group")) +labs( x = "Date", y = "Consumer Sentiment Index") +theme_bw()


#This line is only necessary if you want to preview the plot right after compiling
print(michiganagegroup)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()






tikz(file = "ICSbyincomegroup.tex", width = 6, height = 3.7)
#Simple plot of the dummy data using LaTeX elements







michiganincomegroup <- ggplot(michiganincomegroup) + geom_line(aes(x = DATE, y = Value, linetype = Series, color = Series, group = Series))+ggtitle( paste("Evolution of ICS By Income Group")) +labs( x = "Date", y = "Consumer Sentiment Index") +theme_bw()


#This line is only necessary if you want to preview the plot right after compiling
print(michiganincomegroup)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()


