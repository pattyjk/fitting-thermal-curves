library(readr)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
Corkrey__data <- read_csv("Github/fitting-thermal-curves/Corkrey et al 2016 raw data.csv")
Corkrey__data$rate<-Corkrey__data$rate.min

#select only two strains for testing
#Corkrey__data2<-Corkrey__data[3:16,]

#split dataframe by strain
strains<-split(Corkrey__data, Corkrey__data$strain.code)
#strains<-split(Corkrey__data2, Corkrey__data2$strain.code)


#function to fit each model for all organisms
beta2012<-lapply(strains, function(x) nls_multstart(rate ~ beta_2012(temp = temp, a, b, c, d, e),
                                                           data = x,
                                                           iter = c(6,6,6,6,6),
                                                           start_lower = get_start_vals(x$temp, x$rate, model_name = 'beta_2012') - 10,
                                                           start_upper = get_start_vals(x$temp, x$rate, model_name = 'beta_2012') + 10,
                                                           lower = get_lower_lims(x$temp, x$rate, model_name = 'beta_2012'),
                                                           upper = get_upper_lims(x$temp, x$rate, model_name = 'beta_2012'),
                                                           supp_errors = 'Y',
                                                           convergence_count = FALSE))

#calculate parameters
beta2012_param<-lapply(beta2012, function(x) calc_params(x))

#save data
save(beta2012_param, file = "beta2012_param.RData")
