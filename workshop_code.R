###code for infectious disease workshop

##clear enviroment
rm(list = ls())

##load packages
library(brms)
library(tidyverse)

##set the seed for the workshop
set.seed(58)

##read in function for data generation
gen_data <- function(n, times) {
  tbl_data = tibble(
    i = factor(1:n),
  ) %>% 
    mutate(
      Ct_peak = rnorm(n(), 17.8, 2.2),
      down_slope = pmax(rnorm(n(), 1.7, 1.7/4), 0.01),
    ) %>% 
    expand_grid(
      t = 0:times,
    ) %>% 
    mutate(
      true_ct = Ct_peak + down_slope * t,
      obs_ct = rnorm(n(), true_ct, 3),
    ) %>% 
    group_by(i) %>% 
    # Keep up until first pos test
    filter(t <= min(c(t[obs_ct > 40], Inf))) %>% 
    ungroup()
  
  tbl_data %>% 
    ggplot(aes(colour = i, fill = i)) +
    geom_line(aes(t, true_ct)) +
    geom_point(aes(t, obs_ct)) +
    theme(legend.position = "none")
  
  tbl_data %>% 
    select(i, t, obs_ct) 
  
  return(tbl_data)}

##generate data with 20 individuals and 15 time points
data <- gen_data(20, 15)

##set up initial model with no priors

##define the family of the data
##here we assume the data is normally distributed by the linear predictor 
##is on the log scale
gauss <- gaussian("identity")

#here we use simple lme4 syntax
#the lack of priors will give everything an improper uniform prior by default
model_no_individual_differences <- brm(obs_ct ~ t, data = data, 
                                       family = gauss,
                                       refresh = 0,
                                       control = list(adapt_delta = 0.95),
                                       cores = 4)
#this should give the same estimate as the MLE

##like a frequentist glm we can get the parameter estimates and 
##(credible) intervals with the summary() function
summary(model_no_individual_differences)
pp_check(model_no_individual_differences)

##we can also view the traces of each individual parameter with the plot()
##function
plot(model_no_individual_differences)

##let's test that this is true
freq_model_no_individual_differences <- glm(obs_ct ~ t, data = data, 
                                            family = gaussian(link = "identity"))

summary(freq_model_no_individual_differences)
summary(model_no_individual_differences)

##now we'll add the random effects were were talking about earlier
##both random slopes and random intercepts
model_individual_differences <- brm(obs_ct ~ t + (1+t|i), data = data,
                                    family = gauss, 
                                    refresh = 0,
                                    control = list(adapt_delta = 0.95),
                                    cores = 4)

##and look at the results
summary(model_individual_differences)
pp_check(model_individual_differences)

##this model may have given you a divergent transition
##what precisely a divergent transition is is beyond the scope of this
##workshop, but what you need to know is, if you have any,
##your estimation is invalid
##you also want to check the Rhat and bulk and tail ESSes
##a rule of thumb is that Rhat should be less than 1.01
##and the ESSes should be at least 600 in both categories
##(though the higher precision you want your numerical estimation
##the high ESSes you require)

##now we'll explore the interaction of priors and sample sizes
##for the sake of time, we'll do two of each
##so let's set up a second dataset with more individuals
larger_data <- gen_data(50, 15)

prior_strong <- c(prior(normal(17, 0.1), class = "Intercept"),
                  prior(normal(1.7, 0.1), class = "b"),
                  prior(exponential(1), class = "sd"),
                  prior(lkj(1), class = "cor"))

##weak prior (i.e. initial), weak data (i.e. initial)
summary(model_individual_differences)

##strong prior, weak data (i.e. initial data)
strong_weak <- brm(obs_ct ~ t + (1+t|i), data = data,
                   family = gauss, 
                   prior = prior_strong,
                   refresh = 0,
                   control = list(adapt_delta = 0.95),
                   cores = 4)
summary(strong_weak)

##weak prior (i.e. initial), stronger data
weak_strong <- brm(obs_ct ~ t + (1+t|i), data = larger_data,
                   family = gauss,
                   refresh = 0,
                   control = list(adapt_delta = 0.95),
                   cores = 4)
summary(weak_strong)

##strong prior, stronger data
strong_strong <- brm(obs_ct ~ t + (1+t|i), data = larger_data,
                   family = gauss,
                   prior = prior_strong,
                   refresh = 0,
                   control = list(adapt_delta = 0.95),
                   cores = 4)
summary(strong_strong)

##what can we conclude from this?
