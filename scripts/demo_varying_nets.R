library(hanojoel)
source("scripts/eqm_soln_varying_nets.R")

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# define the net coverage. Need to define the coverage in the past so can have a coverage using the delay
itn_vector <- c(0, 0.2, 0.4)
t_vector <- c(-22, 2*365, 8*365) # number of days at which the itn coverage changes
ITN_IRS_on <- t_vector[min(which(itn_vector != 0))] # time at which nets switch on

# creates the odin model
wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model_itn.R",
                                                              package = "hanojoel"),
                                num_int = 2,
                                het_brackets = 5,
                                age = init_age,
                                init_EIR = init_EIR,
                                country = NULL,
                                admin2 = NULL,
                                itn_vector = itn_vector,
                                t_vector = t_vector,
                                ITN_IRS_on = ITN_IRS_on)

# edits the equilibrium solution so that it splits the population into equal sized
# groups instead of groups based on the coverage levels.
wh <- edit_equilibrium_varying_nets(wh=wh)

# generates model functions with initial state data
mod <- wh$generator(user= wh$state, use_dde = TRUE)

# runs the model
mod_run <- mod$run(t = 1:time_period)
out <- mod$transform_variables(mod_run)

# plots model outputs
plot(out$t/365, out$prev, type='l', col='black')
plot(out$t/365, out$inc, type='l', col='black')
