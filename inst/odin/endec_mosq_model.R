#chopping out the mosq compartment of Hannah's odin model

##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# defining the effect lenght and hazard of the ivermectin and calculate new daily mortality rates
eff_len = user()

haz[] <- user()
dim(haz) = eff_len

dim(mu_vi) = eff_len
mu_vi[1:eff_len] = haz[i]*mu

spor_len = 10
tau_v = spor_len/10

# initial state values:
init_Sv <- user()
init_Ev <- user()
init_Iv <- user()

initial(Sv) <- init_Sv * mv0
initial(Sv_F1) = 0
initial(Sv_F2) = 0
initial(Sx_F1[]) = 0
initial(Sx_F2[]) = 0

initial(Ev_F1[]) = init_Ev/spor_len * mv0
initial(Ev_F2[]) = 0 #init_Ev/spor_len * mv0/2
initial(Ex_F1[,]) = 0
initial(Ex_F2[,]) = 0

initial(Iv_F1) = init_Iv *mv0
initial(Iv_F2) = 0 #init_Iv *mv0/2
initial(Ix_F1[]) = 0
initial(Ix_F2[]) = 0

dim(Sx_F1) = eff_len
dim(Sx_F2) = eff_len
dim(Ix_F1) = eff_len
dim(Ix_F2) = eff_len

dim(Ev_F1) = spor_len
dim(Ev_F2) = spor_len

dim(Ex_F1) = c(spor_len, eff_len)
dim(Ex_F2) = c(spor_len, eff_len)

## user defined ivermectin age and coverage parameters (coverage = coverage of targeted age group)
ivm_cov_par <- user()
ivm_min_age <- user()
ivm_max_age <- user()
ivm_cov = ivm_cov_par*(exp(-ivm_min_age/21) - exp(-ivm_max_age/21))

Ivtot = Iv_F1 + Iv_F2 + sum(Ix_F1) + sum(Ix_F2)
Evtot = sum(Ev_F1) + sum(Ev_F2) + sum(Ex_F1) + sum(Ex_F2)
Svtot = Sv + Sv_F1 + Sv_F2 + sum(Sx_F1) + sum(Sx_F2)

mv = Svtot + Evtot + Ivtot


# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
#cU <- user() # infectiousness U -> mosq
#cD <- user() # infectiousness D -> mosq
#cT <- user() # T -> mosq
#gamma1 <- user() # fitted value of gamma1 characterises cA function
#dim(cA) <- c(na,nh,num_int)
#cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

# Force of infection from humans to mosquitoes
#dim(FOIvijk) <- c(na,nh,num_int)
#omega <- user() #normalising constant for biting rates
#FOIvijk[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq[k]*foi_age[i]/omega
#lag_FOIv=sum(FOIvijk)

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
#delayGam <- user() # Lag from parasites to infectious gametocytes
#delayMos <- user() # Extrinsic incubation period.
#FOIv <- delay(lag_FOIv, delayGam)
FOIv <- 0.008

# Number of mosquitoes that become infected at each time point
#surv <- exp(-mu*delayMos)
#ince <- FOIv * Sv
#lag_incv <- ince * surv
#incv <- delay(lag_incv, delayMos)
#incv <- lag_incv

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa <- 0.5*PL/dPL
#betaa <- mv0 * mu0 * theta2

# determing the on and off of ivermectin effect
ttt[] <- user()
dim(ttt)<-user()
IVRM_start[]<-user()
dim(IVRM_start) <- length(ttt)

IVRM_sr = interpolate(ttt, IVRM_start, "constant")

deriv(Sv) =               betaa - mu*Sv - avhc*Sv

deriv(Sv_F1) =            if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len))  (avhc-FOIv)*(1-ivm_cov)*Sv - avhc*Sv_F1 - mu*Sv_F1    else (avhc - FOIv)*Sv - avhc*Sv_F1 - mu*Sv_F1

deriv(Sv_F2) =           if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) (avhc - FOIv)*(1-ivm_cov)*Sv_F1 -(FOIv*(1-ivm_cov) + avhc*ivm_cov)*Sv_F2 - mu*Sv_F2     else   (avhc-FOIv)*Sv_F1 - FOIv*Sv_F2 - mu*Sv_F2

deriv(Sx_F1[1:eff_len]) = if (t >= (IVRM_sr + i -1)  &&  t < (IVRM_sr + i) && t < (IVRM_sr + eff_len )) ivm_cov*(avhc-FOIv)*Sv - avhc*Sx_F1[i] - mu_vi[i]*Sx_F1[i]     else - mu_vi[i]*Sx_F1[i] - avhc*Sx_F1[i]

deriv(Sx_F2[1:eff_len]) = if (t >= (IVRM_sr + i -1)  &&  t < (IVRM_sr + i) && t < (IVRM_sr + eff_len )) ivm_cov*(avhc-FOIv)*(Sv_F1 + Sv_F2) + (avhc - FOIv)*Sx_F1[i] - FOIv*Sx_F2[i] - mu_vi[i]*Sx_F2[i]     else  (avhc- FOIv)*Sx_F1[i] - mu_vi[i]*Sx_F2[i] - FOIv*Sx_F2[i]



deriv(Ev_F1[1])  =           if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) FOIv*(1-ivm_cov)*Sv  - mu*Ev_F1[i] - avhc*Ev_F1[i] - 1/tau_v*Ev_F1[i]  else FOIv*Sv - avhc*Ev_F1[i] - 1/tau_v*Ev_F1[i] - mu*Ev_F1[i]

deriv(Ev_F1[2:spor_len]) =   if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len))  1/tau_v*Ev_F1[i-1] - avhc*Ev_F1[i] - 1/tau_v*Ev_F1[i] - mu*Ev_F1[i]  else 1/tau_v*Ev_F1[i-1] - 1/tau_v*Ev_F1[i] - avhc*Ev_F1[i]- mu*Ev_F1[i]


deriv(Ev_F2[1]) =            if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) avhc*(1-ivm_cov)*Ev_F1[i] + (FOIv*(1-ivm_cov))*(Sv_F1+Sv_F2)  - mu*Ev_F2[i] - avhc*ivm_cov*Ev_F2[i] - 1/tau_v*Ev_F2[i]  else avhc*Ev_F1[i] + FOIv*(Sv_F1 + Sv_F2) - 1/tau_v*Ev_F2[i] - mu*Ev_F2[i]

deriv(Ev_F2[2:spor_len]) =   if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len))  1/tau_v*Ev_F2[i-1] + avhc*(1-ivm_cov)*Ev_F1[i] - avhc*ivm_cov*Ev_F2[i] - 1/tau_v*Ev_F2[i] - mu*Ev_F2[i]   else 1/tau_v*Ev_F2[i-1] + avhc*Ev_F1[i] - 1/tau_v*Ev_F2[i] - mu*Ev_F2[i]


deriv(Ex_F1[1, 1:eff_len]) =            if (t >= (IVRM_sr + j -1)  &&  t < (IVRM_sr + j ) && t < (IVRM_sr + eff_len)) FOIv*ivm_cov*Sv  - 1/tau_v*Ex_F1[i,j] - mu_vi[j]*Ex_F1[i,j] - avhc*Ex_F1[i,j] else - 1/tau_v*Ex_F1[i,j] - mu_vi[j]*Ex_F1[i,j]  - avhc*Ex_F1[i,j]

deriv(Ex_F1[2:spor_len, 1:eff_len]) =   if (t >= (IVRM_sr + j -1)  &&  t < (IVRM_sr + j ) && t < (IVRM_sr + eff_len)) 1/tau_v*Ex_F1[i-1,j] - 1/tau_v*Ex_F1[i,j] - mu_vi[j]*Ex_F1[i,j] - avhc*Ex_F1[i,j] else  1/tau_v*Ex_F1[i-1,j]  - (1/tau_v + avhc)*Ex_F1[i,j] - mu_vi[j]*Ex_F1[i,j]


deriv(Ex_F2[1, 1:eff_len]) =            if (t >= (IVRM_sr + j -1)  &&  t < (IVRM_sr + j ) && t < (IVRM_sr + eff_len)) ivm_cov*FOIv*(Sv_F1 + Sv_F2) + FOIv*(Sx_F1[j] + Sx_F2[j]) + avhc*ivm_cov*(Ev_F1[i] + Ev_F2[i]) + avhc*Ex_F1[i,j] - 1/tau_v*Ex_F2[i,j] - mu_vi[j]*Ex_F2[i,j] else  avhc*Ex_F1[i,j] + FOIv*(Sx_F1[j] + Sx_F2[j]) - 1/tau_v*Ex_F2[i,j] - mu_vi[j]*Ex_F2[i,j]

deriv(Ex_F2[2:spor_len, 1:eff_len]) =   if (t >= (IVRM_sr + j -1)  &&  t < (IVRM_sr + j ) && t < (IVRM_sr + eff_len)) 1/tau_v*Ex_F2[i-1,j] + avhc*ivm_cov*(Ev_F1[i] + Ev_F2[i]) + avhc*Ex_F1[i,j] -	1/tau_v*Ex_F2[i,j] - mu_vi[j]*Ex_F2[i,j]  else avhc*Ex_F1[i,j] + 1/tau_v*Ex_F2[i-1,j] - mu_vi[j]*Ex_F2[i,j]  -1/tau_v*Ex_F2[i,j]


deriv(Iv_F1) =  if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) 1/tau_v*Ev_F1[spor_len] - (mu + avhc)*Iv_F1     else  1/tau_v*Ev_F1[spor_len] -   (mu + avhc)*Iv_F1

deriv(Iv_F2) =  if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) 1/tau_v*Ev_F2[spor_len] + avhc*(1-ivm_cov)*Iv_F1 - (mu + ivm_cov*avhc)*Iv_F2     else  1/tau_v*Ev_F2[spor_len] + avhc*Iv_F1 -  mu*Iv_F2

deriv(Ix_F1[1:eff_len]) = if (t >= (IVRM_sr + i -1)  &&  t < (IVRM_sr + i) && t < (IVRM_sr + eff_len)) 1/tau_v*Ex_F1[spor_len,i] - (avhc + mu_vi[i])*Ix_F1[i]   else 1/tau_v*Ex_F1[spor_len,i] - (mu_vi[i] + avhc)*Ix_F1[i]

deriv(Ix_F2[1:eff_len]) = if (t >= (IVRM_sr + i -1)  &&  t < (IVRM_sr + i) && t < (IVRM_sr + eff_len))  avhc*ivm_cov*(Iv_F1 + Iv_F2)  + 1/tau_v*Ex_F2[spor_len,i]  + avhc*Ix_F1[i] - mu_vi[i]*Ix_F2[i]  else 1/tau_v*Ex_F2[spor_len, i] + avhc*Ix_F1[i] - mu_vi[i]*Ix_F2[i]



dim(avhc_i) <- num_int
avhc_i[1:num_int] <- cov[i]*av_mosq[i]
avhc <- sum(avhc_i)   # mean biting rate of mosquitoes on humans in the presence of vector control

