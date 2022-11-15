`
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

#break dataframe in separate frames for each strain

strains<-split(Corkrey__data, Corkrey__data$strain.code)


sharpeschoolhigh1981<-for (x in strains) {


 nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                     data = x,
                     iter = 500,
                     start_lower = start_vals - 1,
                     start_upper = start_vals + 1,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')
}
`
