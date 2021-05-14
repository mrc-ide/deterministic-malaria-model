
## get age categories (currently about 21 when c++ has 36)
init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide a string for the admin 2 unit that you want to use the seasonality profile for
admin_str <- NULL

# provide the length of time (in days) that you want to run the model for
time_period <- 5*365

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4


### run original model
model_run_full <- run_model(age=init_age, init_EIR=init_EIR, init_ft=prop_treated,
                            admin2=admin_str, time=time_period)


### run stripped model no change in EIR
model_run_flat <- run_model_stripped_flat(age=init_age, init_EIR=init_EIR, init_ft=prop_treated,
                                      admin2=admin_str, time=time_period)
### run stripped model with change to EIR=100 at time 0
EIR_from_zero <- 100
model_single_change <- run_model_stripped_flat(age=init_age, init_EIR=init_EIR, init_ft=prop_treated,
                                          admin2=admin_str, time=time_period,EIR_flat=EIR_from_zero)



# generate random walk of EIR
genRandWalk <- function(x,vol,randWalk) {
  if (x == 0)    return (randWalk)
  else return(genRandWalk(x-1,vol,c(randWalk,randWalk[length(randWalk)]*exp(rnorm(1)*vol))))
}
EIR_times<-seq(0,1800,by=30)
EIR_volatility<-0.5
EIR_vals=genRandWalk(length(EIR_times)-1,EIR_volatility,init_EIR)

### run stripped model with random walk on EIR
model_run_strip <- run_model_stripped(age=init_age, init_EIR=init_EIR, init_ft=prop_treated,
                       admin2=admin_str, time=time_period,EIR_times=EIR_times,EIR_vals=EIR_vals)



max_clin=max(model_run_strip$inc,model_single_change$inc,model_run_flat$inc,model_run_full$inc)

plot(model_run_strip$t, model_run_strip$inc,ylim=c(0,max_clin*1.2),col="white")
lines(model_run_full$t, model_run_full$inc,col="black",lwd=4)
lines(model_run_flat$t, model_run_flat$inc,col="grey",lwd=4)
lines(model_single_change$t, model_single_change$inc,col="red",lwd=4)
lines(model_run_strip$t, model_run_strip$inc,col="blue",lwd=4)
















