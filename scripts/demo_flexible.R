# load the model package (or have it built and reloaded as described above)
library(hanojoel)

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4

# creates the odin model
wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model.R",
                                                              package = "hanojoel"),
                                het_brackets = 5,
                                age = init_age,
                                init_EIR = init_EIR,
                                num_int = 2,
                                itn_cov = 0.2,
                                ITN_IRS_on = 2*365,
                                init_ft = prop_treated,
                                country = NULL,#"Uganda",
                                admin2 = NULL)#"Tororo")

# generates model functions with initial state data
mod <- wh$generator(user= wh$state, use_dde = TRUE)

# Runs the model
mod_run <- mod$run(t = 1:time_period)
out <- mod$transform_variables(mod_run)
lines(out$t/365, out$prev, col='red', type='l')
#plot(out$t,out$inc)
