# load the model package (or have it built and reloaded as described above)
library(ICDMM)

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4

# Define time for turning on interventions
ITN_IRS_on <- 5*365

out <-    run_model(            het_brackets = 5,
                                age = init_age,
                                time = time_period,
                                init_EIR = init_EIR,
                                num_int = 2,
                                itn_cov = 0.5,
                                ITN_IRS_on = ITN_IRS_on,
                                init_ft = prop_treated,
                                country = "Kenya",
                                admin2 = "Kisumu")

plot(out$t, out$prev, type='l', ylim=c(0.01, 0.33))
plot(out$t, out$inc, type='l')

