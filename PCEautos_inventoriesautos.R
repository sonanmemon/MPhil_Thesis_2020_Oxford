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





data <- read_csv("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/PCE_newautos.csv")


data <- data[-c(1:4),]


data <- ts(data, start = c(1948, 1), frequency = 4)

data

ccf(data[,3], data[,4])









data <- read_csv("C:/Users/ND.COM/Desktop/MPhil Economics/2nd Year/Thesis/Selected Papers/PCE_newautos.csv")

data <- data[-c(1:4),]

data

dat <- data.frame(data$PCE_newautogr, data$autooutput_changeinprivateinventories)






df <- dat

df_x <- eval(substitute(dat$data.PCE_newautogr),df)
df_y <- eval(substitute(dat$data.autooutput_changeinprivateinventories),df)
ccf.object <- ccf(df_x,df_y,plot=FALSE)
output_table <- cbind(lag=ccf.object$lag, x.corr=ccf.object$acf) %>% as_tibble() %>% mutate(cat=ifelse(x.corr>0,"green","red"))
output_table %>% ggplot(aes(x=lag,y=x.corr)) + geom_bar(stat="identity",aes(fill=cat))+scale_fill_manual(values=c("#339933","#cc0000"))+ylab("Cross correlation")+scale_y_continuous(limits=c(-1,1))+theme_bw()+theme(legend.position = "none", plot.title=element_text(size=10))+ggtitle(title) -> p
tikz(file = "PCEautoinventories.tex", width = 6, height = 3.7)
graph <- ggplot(output_table, aes(x=lag,y=x.corr)) + geom_bar(stat="identity",aes(fill=cat))+scale_fill_manual(values=c("#339933","#cc0000"))+ geom_hline(yintercept=0.11) + geom_hline(yintercept = -0.11) + ggtitle("PCE New Autos and Auto Inventories") + ylab("Cross correlations")+scale_y_continuous(limits=c(-1,1))+theme_bw()+theme(legend.position = "none", plot.title=element_text(size=10))

print(graph)
#Necessary to close or the tikxDevice .tex file will not be written





dev.off()





