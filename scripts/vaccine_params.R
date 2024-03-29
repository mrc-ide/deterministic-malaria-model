#Parameters used in the transmission-blocking vaccine model. Not relevant for other models

#All parameters relating to the vaccine characteristics, that won't change depending on the intervention details
# Equation parameters
CS_peak <-621 #Peak antibody titre (after primary vaccination with RTS,S)
CS_boost <- 277 # Peak antibody titre (after booster dose of RTS,S)
p_peak <- 0.88 #Proportion of antibody response that is short lived (after primary vaccination with RTS,S)
p_boost<-0.7 #Proportion of antibody response that is short lived (after booster dose of RTS,S)
ds <- 45 #decay rate of short-lived anitbody response to RTS,S
dl <- 591 #decay rate of long-lived anitbody response to RTS,S
beta_RTS <- 99.2 #Shape param for dose-response of RTS,S
alpha_RTS  <- 0.74 #Shape param for dose-response of RTS,S
V_max_RTS <- 0.93 # Max RTS,S efficacy

#Scheduling parameters
t_boost_rts <- 365

#Extra params for the TBV vaccine (For Pfs230 estimate)
#v_max_effic <- 0.96
#v_alpha <- 2.13
#v_beta <- 349.58
#v_decay <- 0.022
#v_interval <- 365

#Extra params for the TBV vaccine (For Pfs25 estimate)
hill1 <- 2.50 #Hill parameter for dose-response curve
hill2 <- 0.06 #Hill parameter for dose-response curve
mu25 <- 12.63 #Titre (\mu g/ml) to centre the dose-response curve
tau25 <- 22 #Max. antibody titre
rho25 <- 0.7 #Proportion of antibody response associated with quicker decay
ds25 <- 45 #Short antibody half-life
dl25 <- 590 #Long antibody half-life
ds25worse <- 0.5*ds25#35 #Short antibody half-life
dl25worse <- 0.5*dl25 #450 #Long antibody half-life
ds25better <- 2*ds25#55 #Short antibody half-life
dl25better <- 2*dl25 #730 #Long antibody half-life
v_interval <- 365 #Interval between TBV vaccinations
