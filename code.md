## Analysis of thermal curves

### Install rTPC package
```
download and install devtools executable from CRAN
Run install.packages("devtools") in R
library(devtools)
install_github("padpadpadpad/rTPC")

Should work if R >4.2
```

### Preparing data for analysis
```
library(rTPC)
library(devtools)
library(nls.multstart)
library(broom)
library(tidyverse)
library(readr)
library(ggplot2)
library(dplyr)

#read in taxonomic data
corkrey_taxa_data <- read_delim("C:/Users/patty/OneDrive/Desktop/corkrey_taxa_data.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#read in data
#change path to data
Corkrey__data <- read_csv("C:/Users/patty/OneDrive/Desktop/Corkrey et al 2016 raw data.csv")

#raw data is rate min-1, need in hr-1
Corkrey__data$rate <-Corkrey__data$rate.min*60

# number of unique genus/species
unique(Corkrey__data$binomial.name)
#978

#size of data set
dim(Corkrey__data)
#[1] 10956    10

#models need >7 points, remove any strain that doesn't have 7 points
Corkrey__data<-Corkrey__data[Corkrey__data$strain.code %in% names(which(table(Corkrey__data$strain.code) >7)),]
unique(Corkrey__data$binomial.name)
#269

#size of data set after filtering
dim(Corkrey__data)
#[1] 6167   10

species__counts<-as.data.frame(Corkrey__data %>% group_by(binomial.name) %>% summarize(count=n()))
species__counts$count<-as.numeric(species__counts$count)

#add taxnomic data to species counts
species__counts<-merge(x=species__counts, y=corkrey_taxa_data, by.x="binomial.name", by.y = "bionomial_name")

#count family 
family__counts<-as.data.frame(species__counts %>% group_by(family) %>% summarize(count=n()))
family__counts$count<-as.numeric(family__counts$count)

#num unique fams
unique(family__counts$family)
#104

View(family__counts)
write.table(family__counts, "temp_family_counts.txt", quote=F, row.names = F, sep="\t")
```

### Analyze data/plot
```
# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

#this is optional, comment out if not needed
#test all models on one strain
#select one strain
d <- filter(Corkrey__data, strain.code == 3)

#plot raw data
# show the data
ggplot(d, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

#fit all models to one strain
d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(beta = map(data, ~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                         data = .x,
                                         iter = c(6,6,6,6,6),
                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') - 10,
                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)),
         boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         briere2 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         delong = map(data, ~nls_multstart(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
                                           data = .x,
                                           iter = c(4,4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         flinn = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                          data = .x,
                                          iter = c(5,5,5),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         hinshelwood = map(data, ~nls_multstart(rate~hinshelwood_1947(temp = temp, a, e, b, eh),
                                                data = .x,
                                                iter = c(5,5,5,5),
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') - 1,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') + 1,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
         joehnk = map(data, ~nls_multstart(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4,4, 4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         johnson_lewin = map(data, ~suppressWarnings(nls_multstart(rate~ johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                                                                   data = .x,
                                                                   iter = c(4,4,4,4),
                                                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') - 1,
                                                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') + 1,
                                                                   lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                                   upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                                   supp_errors = 'Y',
                                                                   convergence_count = FALSE))),
         kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                               data = .x,
                                               iter = c(4,4,4,4,4),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') - 10,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') + 10,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin2 = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         lrf = map(data, ~nls_multstart(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                        data = d,
                                        iter = c(3,3,3,3),
                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') - 10,
                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') + 10,
                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
         modifiedgaussian = map(data, ~nls_multstart(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                              data = .x,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') - 10,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') + 10,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         rezende = map(data, ~nls_multstart(rate~rezende_2019(temp = temp, q10, a,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         sharpeschoolfull = map(data, ~nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoollow = map(data, ~nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                                                    data = .x,
                                                    iter = c(4,4,4,4),
                                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
                                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
                                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         spain = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         thomas1 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a,b,c,topt),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') - 1,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') + 2,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         thomas2 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
                                            data = .x,
                                            iter = c(3,3,3,3,3),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = temp, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:weibull)

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# plot
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Fits of every model available in rTPC') +
  geom_hline(aes(yintercept = 0), linetype = 2)
  ```
