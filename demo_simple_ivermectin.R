# loading the package
devtools::load_all()

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10


source("R/mda_ivm_functions.R")
#source("R/ivermectin_params_function.R")
ivm_parms0 = ivm_fun(IVM_start_times = 10000,
                     time_period = time_period,
                     hazard_profile = rep(2, 10),
                     ivm_coverage=0.8,
                     ivm_min_age=5,
                     ivm_max_age = 80)

# system.file("extdata/odin_model_itn_ivermectin.R",package = "ICDMM") 

# creates the odin model
wh <- ICDMM:::create_r_model(odin_model_path = "inst/extdata/odin_model_ivermectin.R",
                             num_int = 2,
                             het_brackets = 5,
                             age = init_age,
                             init_EIR = init_EIR,
                             country = NULL,
                             admin2 = NULL,
                             ttt = ivm_parms0$ttt,
                             eff_len = ivm_parms0$eff_len,
                             haz = ivm_parms0$haz,
                             ivm_cov_par = ivm_parms0$ivm_cov_par,
                             ivm_min_age = ivm_parms0$ivm_min_age,
                             ivm_max_age = ivm_parms0$ivm_max_age,
                             IVRM_start = ivm_parms0$IVRM_start
)


#### USE THIS FUNCTION TO RUN THE MODEL (exactly the same as below, but in a function)
runfun <- function(mod_name){
  #mod_name <- edit_equilibrium_varying_nets(wh=mod_name)
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}


op1 = runfun(wh)



###################################################################################################################
#### run examples of the impact of ivermectin ####


ivm_parms0 = ivm_fun(IVM_start_times = 10000,
                     time_period = time_period,
                     hazard_profile = rep(2, 10),
                     ivm_coverage=0.8,
                     ivm_min_age=5,
                     ivm_max_age = 80)

wh0 <- ICDMM:::create_r_model(odin_model_path = "inst/extdata/odin_model_ivermectin.R",
                              num_int = 2,
                              het_brackets = 5,
                              age = init_age,
                              init_EIR = init_EIR,
                              country = NULL,
                              admin2 = NULL,
                              ttt = ivm_parms0$ttt,
                              eff_len = ivm_parms0$eff_len,
                              haz = ivm_parms0$haz,
                              ivm_cov_par = ivm_parms0$ivm_cov_par,
                              ivm_min_age = ivm_parms0$ivm_min_age,
                              ivm_max_age = ivm_parms0$ivm_max_age,
                              IVRM_start = ivm_parms0$IVRM_start)


ivm_parms1 = ivm_fun(IVM_start_times = c(2190, 2250, 2310),
                     time_period = time_period,
                     hazard_profile = rep(2, 10),
                     ivm_coverage=0.8,
                     ivm_min_age=5,
                     ivm_max_age = 80)

wh1 <- ICDMM:::create_r_model(odin_model_path = "inst/extdata/odin_model_ivermectin.R",
                              num_int = 2,
                              het_brackets = 5,
                              age = init_age,
                              init_EIR = init_EIR,
                              country = NULL,
                              admin2 = NULL,
                              ttt = ivm_parms1$ttt,
                              eff_len = ivm_parms1$eff_len,
                              haz = ivm_parms1$haz,
                              ivm_cov_par = ivm_parms1$ivm_cov_par,
                              ivm_min_age = ivm_parms1$ivm_min_age,
                              ivm_max_age = ivm_parms1$ivm_max_age,
                              IVRM_start = ivm_parms1$IVRM_start)


ivm_parms2 = ivm_fun(IVM_start_times = c(2190, 2250, 2310),
                     time_period = time_period,
                     hazard_profile = rep(2, 28),
                     ivm_coverage=0.8,
                     ivm_min_age=5,
                     ivm_max_age = 80)

wh2 <- ICDMM:::create_r_model(odin_model_path = "inst/extdata/odin_model_ivermectin.R",
                              num_int = 2,
                              het_brackets = 5,
                              age = init_age,
                              init_EIR = init_EIR,
                              country = NULL,
                              admin2 = NULL,
                              ttt = ivm_parms2$ttt,
                              eff_len = ivm_parms2$eff_len,
                              haz = ivm_parms2$haz,
                              ivm_cov_par = ivm_parms2$ivm_cov_par,
                              ivm_min_age = ivm_parms2$ivm_min_age,
                              ivm_max_age = ivm_parms2$ivm_max_age,
                              IVRM_start = ivm_parms2$IVRM_start)



res0 = runfun(wh0)
res1 = runfun(wh1)
res2 = runfun(wh2)

cols = c("grey40", "deeppink2", "deepskyblue3")

par(mfrow=c(1,2), mar=c(5,4,1,1))
plot(res0$t/365, res0$inc*1000*365, type="l", ylab="Annual incidence per 1,000", xlab="Year", lwd=3, col=cols[1], 
     xlim = c(5.7, 8), ylim = c(0, 800), las=1)
lines(res1$t/365, res1$inc*1000*365, lwd=3, col=cols[2])
lines(res2$t/365, res2$inc*1000*365, lwd=3, col=cols[3])
arrows(c(2190, 2250, 2310)/365, -50, c(2190, 2250, 2310)/365, 20, length=0.15, lwd=3, col="goldenrod2")

legend("topright", c("No endectocide", "10 day endectocide with HR=2", "28 day endectocide with HR=2"), 
       col = cols, lwd=3, bty="n", cex=0.8)

plot(res0$t/365, res0$prev*100, type="l", ylab="Slide prevalence (%)", xlab="Year", lwd=3, col=cols[1], 
     xlim = c(5.7, 8), ylim = c(0, 30), las=1)
lines(res1$t/365, res1$prev*100, lwd=3, col=cols[2])
lines(res2$t/365, res2$prev*100, lwd=3, col=cols[3])
arrows(c(2190, 2250, 2310)/365, -50, c(2190, 2250, 2310)/365, 1, length=0.15, lwd=3, col="goldenrod2")

