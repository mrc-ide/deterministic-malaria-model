library(ICDMM)

init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)

#Many vaccine-specific parameters stored here
source("scripts/vaccine_params.R")

#Note: in this version of the model, the vaccine takes the place of IRS.
#This means that 'irs_cov' is the vaccination coverage
#Therefore, at present one cannot model the introduction of the vaccine alongside IRS
#ITNs function as normal.

#Parameters to define
ITN_IRS_on <- 4 * 365 + 100 # When ITNs are introduced
vacc_lag <- 0 #Delay between ITNs being introduced and the start of vaccination campaign
vac_cov <- 0.8 #Vaccine coverage

#Min & Max age for vaccination. Values refer to the position in age vector
age_min_rts <- 3
age_max_rts <- which(init_age > 50)[1] - 1 # index of age vector before age is 50 years
age_min_tbv = 3
age_max_tbv = which(init_age > 50)[1] - 1 # index of age vector before age is 50 years

out <- run_model(model = 'odin_model_TBV',
                 init_EIR = 50,
                 time = 10*365,
                 age = init_age,
                 RTS_switch = 1, #Include RTS,S in vaccination campaign?
                 switch_TBV = 1, #Include TBV in vaccination campaign?
                 switch_TRA_to_TBA = 1, # Set to 1 for the model used in the Nat Comms paper
                 irs_cov = vac_cov,
                 itn_cov = 0,
                 num_int = 4,
                 ITN_IRS_on = ITN_IRS_on,
                 vacc_lag = 0.0,
                 v_interval = v_interval, #frequency of vaccination
                 hill1 = hill1,
                 hill2 = hill2,
                 mu25 = mu25,
                 tau25 = tau25,
                 rho25 = rho25,
                 ds25 = ds25,
                 dl25 = dl25,
                 age_min_tbv = age_min_tbv,
                 age_max_tbv = age_max_tbv,
                 age_min_rts = age_min_rts,
                 age_max_rts = age_max_rts,
                 t_boost_rts = t_boost_rts,
                 CS_peak = CS_peak,
                 CS_boost = CS_boost,
                 p_peak = p_peak,
                 p_boost = p_boost,
                 ds = ds,
                 dl = dl,
                 beta_RTS = beta_RTS,
                 alpha_RTS = alpha_RTS,
                 V_max_RTS = V_max_RTS)
plot(out$t,out$inc05, main= "Incidence in Under 5s", type='l', col = "red")
