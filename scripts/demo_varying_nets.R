library(hanojoel)

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# define the net coverage. Need to define the coverage in the past so can have a coverage using the delay
itn_vector <- c(0, 0.2)
t_vector <- c(-22, 2*365) # number of days at which the itn coverage changes
#t_vector[min(which(itn_vector != 0))] # time at which nets switch on
ITN_IRS_on <- 2*365

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

# Need to add population split and use to edit initial conditions
wh$state$pop_split <- c(0.1, 0.9) # Pop split needs to sum to 1 - expect total population of 1
cov <- wh$mpl$cov

# This edit here only works when IRS is turned off
# Divide the no intervention compartment by cov[1] to get original eqm soln
tmp_init_S <- wh$state$init_S[,,1] / cov[1]
wh$state$init_S[,,1] <- tmp_init_S * wh$state$pop_split[1]
wh$state$init_S[,,2] <- tmp_init_S * wh$state$pop_split[2]
tmp_init_T <- wh$state$init_T[,,1] / cov[1]
wh$state$init_T[,,1] <- tmp_init_T* wh$state$pop_split[1]
wh$state$init_T[,,2] <- tmp_init_T* wh$state$pop_split[2]
tmp_init_D <- wh$state$init_D[,,1] / cov[1]
wh$state$init_D[,,1] <- tmp_init_D* wh$state$pop_split[1]
wh$state$init_D[,,2] <- tmp_init_D* wh$state$pop_split[2]
tmp_init_A <- wh$state$init_A[,,1] / cov[1]
wh$state$init_A[,,1] <- tmp_init_A* wh$state$pop_split[1]
wh$state$init_A[,,2] <- tmp_init_A* wh$state$pop_split[2]
tmp_init_U <- wh$state$init_U[,,1] / cov[1]
wh$state$init_U[,,1] <- tmp_init_U* wh$state$pop_split[1]
wh$state$init_U[,,2] <- tmp_init_U* wh$state$pop_split[2]
tmp_init_P <- wh$state$init_P[,,1] / cov[1]
wh$state$init_P[,,1] <- tmp_init_P* wh$state$pop_split[1]
wh$state$init_P[,,2] <- tmp_init_P* wh$state$pop_split[2]

# generates model functions with initial state data
mod <- wh$generator(user= wh$state, use_dde = TRUE)

# Runs the model
mod_run <- mod$run(t = 1:time_period)
out <- mod$transform_variables(mod_run)

plot(out$t/365, out$prev, type='l', col='black')
plot(out$t/365, out$inc, type='l', col='black')
