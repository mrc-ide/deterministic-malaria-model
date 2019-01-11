# load the model package (or have it built and reloaded as described above)
library(hanojoel)

# create a vector of age categories
init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide a string for the admin 2 unit that you want to use the seasonality profile for
admin_str <- "Tororo"

# provide the length of time (in days) that you want to run the model for
time_period <- 365*1

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4

# run the model
model_run <- run_model(age=init_age, EIR=init_EIR, ft=prop_treated,
                       admin2=admin_str, time=time_period)

# access model output
output <- model_run$dat

# print the plot of proportion of population in each compartment
print(model_run$plot)
