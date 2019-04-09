library(hanojoel)

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# define the net coverage. Need to define the coverage in the past so can have a coverage using the delay
itn_vector <- c(0, 0.1, 0.3)
t_vector <- c(-25, 2*365, 5*365) # number of days at which the itn coverage changes
ITN_IRS_on <- t_vector[min(which(itn_vector != 0))]

# creates the odin model
wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model_itn.R",
                                                              package = "hanojoel"),
                                num_int = 2,
                                het_brackets = 3,
                                age = init_age,
                                init_EIR = init_EIR,
                                country = NULL,
                                admin2 = NULL,
                                itn_vector = itn_vector,
                                t_vector = t_vector,
                                ITN_IRS_on = ITN_IRS_on)

# Need to define here ITN_IRS_on and pop split

# generates model functions with initial state data
mod <- wh$generator(user= wh$state, use_dde = TRUE)

# Runs the model
mod_run <- mod$run(t = 1:time_period)
out <- mod$transform_variables(mod_run)

plot(out$t/365, out$prev, type='l', col='black')
plot(out$t/365, out$inc, type='l', col='black')
