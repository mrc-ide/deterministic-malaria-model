# P. vivax model

# Main differences to P. falciparum model are:
# - human states and transitions
# - immunity
# - stratification into hypnozoite batches (number of batches referred to as k)
# - 2 prevalence outputs: prevLM and prevPCR

## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
nk <- user() # number of hypnozoite batch categories
ft <- user() # proportion of cases treated

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------

# Human states as specified in full transmission model
# Based on supplementary methods and associated code from https://www.nature.com/articles/s41467-018-05860-8

# fitted parameters for human compartments:
eta <- user()          # death rate for exponential population distribution
age_rate[] <- user()   # rate at which humans move through age categories
dim(age_rate) <- na
age[] <- user()        # vector of age groups
dim(age) <- na
het_wt[] <- user()     # weights of the different heterogeneous biting categories
dim(het_wt) <- nh
# Progression rates
rLM <- user()          # rate of movement from I_LM -> I_PCR
rT <- user()           # rate of treatment working: T -> P
rD <- user()           # rate from I_D -> I_LM
rP <- user()           # rate at which prophylaxis wears off P -> S
# dPCR, the duration of subpatent (PCR-detectable) infection (clearance I_PCR -> S)
# depends on immunity

# Hypnozoite parameters:
ff <- user()                   # rate of relapse
gamma_L <- user()              # rate of hypnozoite clearance
kk[] <- user()                 # number of hypnozoite batches for each l index
# in equations below
dim(kk) <- nk
K_max_switch[] <- user()       # switch to set some transitions to 0 when k=K_max
dim(K_max_switch) <- nk
# K_max is the maximum number of hypnozoite batches (model input)
K_max_switch_on[] <- user()    # switch to only apply some transitions when k=K_max
dim(K_max_switch_on) <- nk
K0_switch[] <- user()          # switch to set some transitions to 0 when k=0
dim(K0_switch) <- nk

# S - SUSCEPTIBLE
init_S[,,,] <- user()
dim(init_S) <- c(na,nh,num_int,nk)
initial(S[,,,]) <- init_S[i,j,k,l]
dim(S) <- c(na,nh,num_int,nk)

deriv(S[1, 1:nh, 1:num_int,1]) <- -FOI[i,j,k]*S[i,j,k,l] - ff*kk[l]*S[i,j,k,l] +
  1/dPCR[i,j,k,l]*I_PCR[i,j,k,l] + rP*P[i,j,k,l]  -
  gamma_L*kk[l]*S[i,j,k,l] + gamma_L*(kk[l]+1)*S[i,j,k,l+1] +
  cov[k]*eta*H*het_wt[j] - (eta+age_rate[i])*S[i,j,k,l]
# All babies are born susceptible with 0 hypnozoite batches

deriv(S[1, 1:nh, 1:num_int,2:nk]) <- -FOI[i,j,k]*S[i,j,k,l]- ff*kk[l]*S[i,j,k,l] +
  1/dPCR[i,j,k,l]*I_PCR[i,j,k,l] + rP*P[i,j,k,l]  -
  gamma_L*kk[l]*S[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*S[i,j,k,l+1] -
  (eta+age_rate[i])*S[i,j,k,l]

deriv(S[2:na, 1:nh, 1:num_int, 1:nk]) <- -FOI[i,j,k]*S[i,j,k,l] - ff*kk[l]*S[i,j,k,l] +
  1/dPCR[i,j,k,l]*I_PCR[i,j,k,l] + rP*P[i,j,k,l]  -
  gamma_L*kk[l]*S[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*S[i,j,k,l+1] -
  (eta+age_rate[i])*S[i,j,k,l] + age_rate[i-1]*S[i-1,j,k,l]

# I_PCR - PCR-DETECTABLE (LOW-DENSITY) INFECTION
init_I_PCR[,,,] <- user()
dim(init_I_PCR) <- c(na,nh,num_int,nk)
initial(I_PCR[,,,]) <- init_I_PCR[i,j,k,l]
dim(I_PCR) <- c(na,nh,num_int,nk)

deriv(I_PCR[1, 1:nh, 1:num_int, 1:nk]) <- -FOI[i,j,k]*I_PCR[i,j,k,l] -
  ff*kk[l]*I_PCR[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*(1-phi_LM[i,j,k,l-1])*(S[i,j,k,l-1]+I_PCR[i,j,k,l-1]) +
  ff*kk[l]*(1-phi_LM[i,j,k,l])*(S[i,j,k,l]+I_PCR[i,j,k,l]) - 1/dPCR[i,j,k,l]*I_PCR[i,j,k,l] +
  rLM*I_LM[i,j,k,l] - gamma_L*kk[l]*I_PCR[i,j,k,l] +
  K_max_switch[l]*gamma_L*(kk[l]+1)*I_PCR[i,j,k,l+1] +
  K_max_switch_on[l]*FOI[i,j,k]*(1-phi_LM[i,j,k,l])*(S[i,j,k,l]+I_PCR[i,j,k,l])-
  (eta+age_rate[i])*I_PCR[i,j,k,l]

deriv(I_PCR[2:na, 1:nh, 1:num_int, 1:nk]) <- -FOI[i,j,k]*I_PCR[i,j,k,l] -
  ff*kk[l]*I_PCR[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*(1-phi_LM[i,j,k,l-1])*(S[i,j,k,l-1]+I_PCR[i,j,k,l-1]) +
  ff*kk[l]*(1-phi_LM[i,j,k,l])*(S[i,j,k,l]+I_PCR[i,j,k,l]) - 1/dPCR[i,j,k,l]*I_PCR[i,j,k,l] +
  rLM*I_LM[i,j,k,l] - gamma_L*kk[l]*I_PCR[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*I_PCR[i,j,k,l+1] +
  K_max_switch_on[l]*FOI[i,j,k]*(1-phi_LM[i,j,k,l])*(S[i,j,k,l]+I_PCR[i,j,k,l])-
  (eta+age_rate[i])*I_PCR[i,j,k,l] + age_rate[i-1]*I_PCR[i-1,j,k,l]
# K_max_switch_on expression means that superinfections occurring among
# people with k=K_max do not lead to an increase in hypnozoite batches


# I_LM - LIGHT MICROSCOPY-DETECTABLE (HIGH-DENSITY) INFECTION
init_I_LM[,,,] <- user()
dim(init_I_LM) <- c(na,nh,num_int,nk)
initial(I_LM[,,,]) <- init_I_LM[i,j,k,l]
dim(I_LM) <- c(na,nh,num_int,nk)

deriv(I_LM[1, 1:nh, 1:num_int, 1:nk]) <- -FOI[i,j,k]*I_LM[i,j,k,l] -
  ff*kk[l]*I_LM[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*(1-phi_D[i,j,k,l-1])*
  (phi_LM[i,j,k,l-1]*S[i,j,k,l-1] + phi_LM[i,j,k,l-1]*I_PCR[i,j,k,l-1] + I_LM[i,j,k,l-1]) +
  ff*kk[l]*(1-phi_D[i,j,k,l])*(phi_LM[i,j,k,l]*S[i,j,k,l] + phi_LM[i,j,k,l]*I_PCR[i,j,k,l] + I_LM[i,j,k,l]) -
  rLM*I_LM[i,j,k,l] + rD*I_D[i,j,k,l] -
  gamma_L*kk[l]*I_LM[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*I_LM[i,j,k,l+1]+
  K_max_switch_on[l]*FOI[i,j,k]*(1-phi_D[i,j,k,l])*(phi_LM[i,j,k,l]*S[i,j,k,l]+phi_LM[i,j,k,l]*I_PCR[i,j,k,l]+I_LM[i,j,k,l]) -
  (eta+age_rate[i])*I_LM[i,j,k,l]

deriv(I_LM[2:na, 1:nh, 1:num_int, 1:nk]) <- -FOI[i,j,k]*I_LM[i,j,k,l] -
  ff*kk[l]*I_LM[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*(1-phi_D[i,j,k,l-1])*
  (phi_LM[i,j,k,l-1]*S[i,j,k,l-1] + phi_LM[i,j,k,l-1]*I_PCR[i,j,k,l-1] + I_LM[i,j,k,l-1]) +
  ff*kk[l]*(1-phi_D[i,j,k,l])*(phi_LM[i,j,k,l]*S[i,j,k,l] + phi_LM[i,j,k,l]*I_PCR[i,j,k,l] + I_LM[i,j,k,l]) -
  rLM*I_LM[i,j,k,l] + rD*I_D[i,j,k,l] -
  gamma_L*kk[l]*I_LM[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*I_LM[i,j,k,l+1] +
  K_max_switch_on[l]*FOI[i,j,k]*(1-phi_D[i,j,k,l])*(phi_LM[i,j,k,l]*S[i,j,k,l]+phi_LM[i,j,k,l]*I_PCR[i,j,k,l]+I_LM[i,j,k,l]) -
  (eta+age_rate[i])*I_LM[i,j,k,l] + age_rate[i-1]*I_LM[i-1,j,k,l]

# I_D - CLINICAL DISEASE
init_I_D[,,,] <- user()
dim(init_I_D) <- c(na,nh,num_int,nk)
initial(I_D[,,,]) <- init_I_D[i,j,k,l]
dim(I_D) <- c(na,nh,num_int,nk)

deriv(I_D[1, 1:nh, 1:num_int, 1:nk]) <- -K_max_switch[l]*FOI[i,j,k]*I_D[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*I_D[i,j,k,l-1] +
  (1-ft)*clin_inc[i,j,k,l] - rD*I_D[i,j,k,l] -
  gamma_L*kk[l]*I_D[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*I_D[i,j,k,l+1] -
  (eta+age_rate[i])*I_D[i,j,k,l]

deriv(I_D[2:na, 1:nh, 1:num_int, 1:nk]) <- -K_max_switch[l]*FOI[i,j,k]*I_D[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*I_D[i,j,k,l-1] +
  (1-ft)*clin_inc[i,j,k,l] - rD*I_D[i,j,k,l] -
  gamma_L*kk[l]*I_D[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*I_D[i,j,k,l+1] -
  (eta+age_rate[i])*I_D[i,j,k,l] + age_rate[i-1]*I_D[i-1,j,k,l]

# T- SUCCESSFULLY TREATED
init_T[,,,] <- user()
dim(init_T) <- c(na,nh,num_int,nk)
initial(T[,,,]) <- init_T[i,j,k,l]
dim(T) <- c(na,nh,num_int,nk)

deriv(T[1, 1:nh, 1:num_int, 1:nk]) <-  -K_max_switch[l]*FOI[i,j,k]*T[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*T[i,j,k,l-1] +
  ft*clin_inc[i,j,k,l] - rT*T[i,j,k,l] -
  gamma_L*kk[l]*T[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*T[i,j,k,l+1]-
  (eta+age_rate[i])*T[i,j,k,l]

deriv(T[2:na, 1:nh, 1:num_int, 1:nk]) <- -K_max_switch[l]*FOI[i,j,k]*T[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*T[i,j,k,l-1] +
  ft*clin_inc[i,j,k,l] - rT*T[i,j,k,l] -
  gamma_L*kk[l]*T[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*T[i,j,k,l+1]-
  (eta+age_rate[i])*T[i,j,k,l] + age_rate[i-1]*T[i-1,j,k,l]

# P - PROPHYLAXIS
init_P[,,,] <- user()
dim(init_P) <- c(na,nh,num_int,nk)
initial(P[,,,]) <- init_P[i,j,k,l]
dim(P) <- c(na,nh,num_int,nk)

deriv(P[1, 1:nh, 1:num_int, 1:nk]) <- -K_max_switch[l]*FOI[i,j,k]*P[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*P[i,j,k,l-1] +
  rT*T[i,j,k,l] - rP*P[i,j,k,l] -
  gamma_L*kk[l]*P[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*P[i,j,k,l+1] -
  (eta+age_rate[i])*P[i,j,k,l]

deriv(P[2:na, 1:nh, 1:num_int, 1:nk]) <- -K_max_switch[l]*FOI[i,j,k]*P[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*P[i,j,k,l-1] +
  rT*T[i,j,k,l] - rP*P[i,j,k,l] -
  gamma_L*kk[l]*P[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*P[i,j,k,l+1] -
  (eta+age_rate[i])*P[i,j,k,l] + age_rate[i-1]*P[i-1,j,k,l]

# The number of new clinical episodes at each timestep
# Clinical episodes can occur from new infectious bite (FOI) or relapse (ff)
dim(clin_inc) <- c(na,nh,num_int,nk)
clin_inc[1:na, 1:nh, 1:num_int, 1:nk] <- ff*kk[l]*phi_D[i,j,k,l]*
  (phi_LM[i,j,k,l]*S[i,j,k,l] + phi_LM[i,j,k,l]*I_PCR[i,j,k,l] + I_LM[i,j,k,l]) +
  K0_switch[l]*FOI[i,j,k]*phi_D[i,j,k,l-1]*
  (phi_LM[i,j,k,l-1]*S[i,j,k,l-1] + phi_LM[i,j,k,l-1]*I_PCR[i,j,k,l-1] + I_LM[i,j,k,l-1]) +
  K_max_switch_on[l]*FOI[i,j,k]*phi_D[i,j,k,l]*
  (phi_LM[i,j,k,l]*S[i,j,k,l] + phi_LM[i,j,k,l]*I_PCR[i,j,k,l] + I_LM[i,j,k,l])

# Relapse incidence
dim(clin_inc_relapse) <- c(na,nh,num_int,nk)
clin_inc_relapse[1:na, 1:nh, 1:num_int, 1:nk] <- ff*kk[l]*phi_D[i,j,k,l]*
  (phi_LM[i,j,k,l]*S[i,j,k,l] + phi_LM[i,j,k,l]*I_PCR[i,j,k,l] + I_LM[i,j,k,l])
output(inc_relapse) <- sum(clin_inc_relapse[,,,])

# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,,,])
Th <- sum(T[,,,])
I_Dh <- sum(I_D[,,,])
I_LMh <- sum(I_LM[,,,])
I_PCRh <- sum(I_PCR[,,,])
Ph <- sum(P[,,,])
H <- Sh + Th + I_Dh + I_LMh + I_PCRh + Ph

##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# AP - Anti-parasite immunity due to exposure to previous infection
# AC - Clinical immunity due to exposure to previous infection
# AP_mat - Maternally acquired anti-parasite immunity by babies from mothers (assumed to be proportional to the immunity of a 20 year old woman)
# AC_mat - Maternally acquired clinical immunity by babies from mothers (assumed to be proportional to the immunity of a 20 year old woman)

# Immunity parameters:
# Parameters for acquisition and decay of immunity
d_par <- user()        # duration of anti-parasite immunity
u_par <- user()        # anti-parasite immune boosting refractory period
d_clin <- user()       # duration of clinical immunity
u_clin <- user()       # clinical immune boosting refractory period
p_mat <- user()        # new-born immunity relative to mother’s
d_mat <- user()        # duration of maternal immunity
# Parameters for Hill functions
phi_LM_max <- user()   # probability of LM-detectable infection with no immunity
phi_LM_min <- user()   # probability of LM-detectable infection with full immunity
A_LM_50pc <- user()    # anti-parasite immunity for 50% reduction in LM-detectable infection
K_LM <- user()         # shape parameter for LM-detectable infection
A_PCR_50pc <- user()   # anti-parasite immunity for 50% reduction in duration of PCR-detectable infection
K_PCR <- user()        # shape parameter for duration of PCR-detectable infection
phi_D_max <- user()    # probability of clinical episode with no immunity
phi_D_min <- user()    # probability of clinical episode with full immunity
dPCR_min <- user()     # duration of PCR-detectable infection with full immunity
dPCR_max <- user()     # duration of PCR-detectable infection with no immunity
A_D_50pc <- user()     # anti-parasite immunity for 50% reduction in clinical episode
K_D <- user()          # shape parameter for clinical episode probability
# Other parameters
age20l <- user(integer=TRUE) # lower index of age 20 age compartment
age20u <- user(integer=TRUE) # upper index of age 20 age compartment
age_20_factor <- user()      # factor calculated in equilibrium solution

# TRACKING HYPNOZOITE STATE DISTRIBUTION
init_Hyp[,,,] <- user()
dim(init_Hyp) <- c(na,nh,num_int,nk)
initial(Hyp[,,,]) <- init_Hyp[i,j,k,l]
dim(Hyp) <- c(na,nh,num_int,nk)

deriv(Hyp[1, 1:nh, 1:num_int, 1]) <- -K_max_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1] - gamma_L*kk[l]*Hyp[i,j,k,l] +
  K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1] - (eta+age_rate[i])*Hyp[i,j,k,l] +
  cov[k]*eta*H*het_wt[j]
deriv(Hyp[1, 1:nh, 1:num_int, 2:nk]) <- -K_max_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1] - gamma_L*kk[l]*Hyp[i,j,k,l] +
  K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1] - (eta+age_rate[i])*Hyp[i,j,k,l] -
  reduce_hypnozoites*Hyp[i,j,k,l]

deriv(Hyp[2:na, 1:nh, 1:num_int, 1]) <- -K_max_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1] - gamma_L*kk[l]*Hyp[i,j,k,l] +
  K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1] - (eta+age_rate[i])*Hyp[i,j,k,l] +
  age_rate[i-1]*Hyp[i-1,j,k,l]
deriv(Hyp[2:na, 1:nh, 1:num_int, 2:nk]) <- -K_max_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1] - gamma_L*kk[l]*Hyp[i,j,k,l] +
  K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1] - (eta+age_rate[i])*Hyp[i,j,k,l] +
  age_rate[i-1]*Hyp[i-1,j,k,l] -
  reduce_hypnozoites*Hyp[i,j,k,l]

# Switch to stop rebound in infections due to relapse:
# reduce_hypnozoites is set to 1 and futher depletes hypnozoite reservoir
# if it falls below a certain threshold. This is to avoid rebound in infections
# after decades due to very low levels of hypnozoites remaining
stop_rebound_switch <- user()         # switch to allow option of depleting hypnozoite reservoir
hypnozoite_prev_threshold <- user()   # hypnozoite prevalence in the population at which stop_rebound_switch is triggered
dim(Hyp_prev) <- c(na,nh,num_int)
Hyp_prev[,,] <- sum(Hyp[i,j,k,2:nk])  # calculate hypnozoite prevalence
reduce_hypnozoites <- if (stop_rebound_switch==1 &&
                          sum(Hyp_prev[,,])<hypnozoite_prev_threshold) 1 else 0

# Calculate proportion in each hypnozoite state for each age and heterogeneity group
dim(hyp_wt) <- c(na,nh,num_int,nk)
hyp_wt[,,,] <- Hyp[i,j,k,l]/sum(Hyp[i,j,k,])

# AP - ANTI-PARASITE IMMUNITY
init_AP[,,,] <- user()
dim(init_AP) <- c(na,nh,num_int,nk)
initial(AP[,,,]) <- init_AP[i,j,k,l]
dim(AP) <- c(na,nh,num_int,nk)

deriv(AP[1, 1:nh, 1:num_int, 1:nk]) <- K0_switch[l]*(FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff)/((FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff) * u_par + 1) -
  K_max_switch[l]*FOI[i,j,k]*AP[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]*AP[i,j,k,l-1] - 1/d_par*AP[i,j,k,l] -
  gamma_L*kk[l]*AP[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1]/Hyp[i,j,k,l]*AP[i,j,k,l+1] -
  (eta+age_rate[i])*AP[i,j,k,l]

deriv(AP[2:na, 1:nh, 1:num_int, 1:nk]) <- K0_switch[l]*(FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff)/((FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff) * u_par + 1) -
  K_max_switch[l]*FOI[i,j,k]*AP[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]*AP[i,j,k,l-1] - 1/d_par*AP[i,j,k,l] -
  gamma_L*kk[l]*AP[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1]/Hyp[i,j,k,l]*AP[i,j,k,l+1] -
  (eta+age_rate[i])*AP[i,j,k,l] + age_rate[i-1]*AP[i-1,j,k,l]*Hyp[i-1,j,k,l]/Hyp[i,j,k,l]

# AC - CLINICAL IMMUNITY
init_AC[,,,] <- user()
dim(init_AC) <- c(na,nh,num_int,nk)
initial(AC[,,,]) <- init_AC[i,j,k,l]
dim(AC) <- c(na,nh,num_int,nk)

deriv(AC[1, 1:nh, 1:num_int, 1:nk]) <- K0_switch[l]*(FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff)/((FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff) * u_clin + 1) -
  K_max_switch[l]*FOI[i,j,k]*AC[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]*AC[i,j,k,l-1] - 1/d_clin*AC[i,j,k,l] -
  gamma_L*kk[l]*AC[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1]/Hyp[i,j,k,l]*AC[i,j,k,l+1] -
  (eta+age_rate[i])*AC[i,j,k,l]

deriv(AC[2:na, 1:nh, 1:num_int, 1:nk]) <- K0_switch[l]*(FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff)/((FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]+kk[l]*ff) * u_clin + 1) -
  K_max_switch[l]*FOI[i,j,k]*AC[i,j,k,l] +
  K0_switch[l]*FOI[i,j,k]*Hyp[i,j,k,l-1]/Hyp[i,j,k,l]*AC[i,j,k,l-1] - 1/d_clin*AC[i,j,k,l] -
  gamma_L*kk[l]*AC[i,j,k,l] + K_max_switch[l]*gamma_L*(kk[l]+1)*Hyp[i,j,k,l+1]/Hyp[i,j,k,l]*AC[i,j,k,l+1] -
  (eta+age_rate[i])*AC[i,j,k,l] + age_rate[i-1]*AC[i-1,j,k,l]*Hyp[i-1,j,k,l]/Hyp[i,j,k,l]

# AP_MAT - MATERNALLY ACQUIRED ANTI-PARASITE IMMUNITY
# Maternal immunity is averaged over hypnozoite states of mothers
dim(AP_MAT_mean_pre) <- c(na,nh,num_int,nk)
dim(AP_MAT_mean) <- c(na,nh,num_int)
dim(init_AP_MAT_pre) <- c(nh,num_int)  # proportion of mother's anti-parasite immunity (20 year old woman) averaged over hypnozoite states
dim(AP_MAT) <- c(na,nh,num_int,nk)     # maternal anti-parasite immunity by age
# Note maternal immunity does not vary by hypnozoite state but is stored separately
# for convenience in later calculations

# Calculate AP_MAT averaged over hypnozoite states
AP_MAT_mean_pre[,,,] <- hyp_wt[i,j,k,l]*AP[i,j,k,l]
AP_MAT_mean[,,] <- sum(AP_MAT_mean_pre[i,j,k,])

# Use this for calculation of maternal immunity
init_AP_MAT_pre[1:nh,1:num_int] <- p_mat*(AP_MAT_mean[age20l,i,j] + age_20_factor*(AP_MAT_mean[age20u,i,j]-AP_MAT_mean[age20l,i,j]))
AP_MAT[1:na, 1:nh, 1:num_int, 1:nk] <- init_AP_MAT_pre[j,k]*exp(-age[i]/d_mat)

# AC_MAT - MATERNALLY ACQUIRED CLINICAL IMMUNITY
dim(AC_MAT_mean_pre) <- c(na,nh,num_int,nk)
dim(AC_MAT_mean) <- c(na,nh,num_int)
dim(init_AC_MAT_pre) <- c(nh,num_int)  # proportion of mother's anti-clinical immunity (20 year old woman) averaged over hypnozoite states
dim(AC_MAT) <- c(na,nh,num_int,nk)     # maternal anti-clinical immunity by age

# Calculate AC_MAT averaged over hypnozoite states
AC_MAT_mean_pre[,,,] <- hyp_wt[i,j,k,l]*AC[i,j,k,l]
AC_MAT_mean[,,] <- sum(AC_MAT_mean_pre[i,j,k,])

# Use this for calculation of maternal immunity
init_AC_MAT_pre[1:nh,1:num_int] <- p_mat*(AC_MAT_mean[age20l,i,j] + age_20_factor*(AC_MAT_mean[age20u,i,j]-AC_MAT_mean[age20l,i,j]))
AC_MAT[1:na, 1:nh, 1:num_int, 1:nk] <- init_AC_MAT_pre[j,k]*exp(-age[i]/d_mat)

# HILL FUNCTIONS

# Probability of LM-detectable (vs. PCR-detectable) infection, dependent on anti-parasite immunity
dim(phi_LM) <- c(na,nh,num_int,nk)
phi_LM[1:na,1:nh,1:num_int,1:nk] <- phi_LM_min + (phi_LM_max-phi_LM_min) * 1/(1+((AP[i,j,k,l]+AP_MAT[i,j,k,l])/A_LM_50pc)^K_LM)

# Probability of clinical episode among LM-detectable infections, dependent on clinical immunity
dim(phi_D) <- c(na,nh,num_int,nk)
phi_D[1:na,1:nh,1:num_int,1:nk] <- phi_D_min + (phi_D_max-phi_D_min) * 1/(1+((AC[i,j,k,l]+AC_MAT[i,j,k,l])/A_D_50pc)^K_D)

# Duration of PCR-detectable infection, dependent on anti-parasite immunity
dim(dPCR) <- c(na,nh,num_int,nk)
dPCR[1:na,1:nh,1:num_int,1:nk] <- dPCR_min + (dPCR_max-dPCR_min) * 1/(1+((AP[i,j,k,l]+AP_MAT[i,j,k,l])/A_PCR_50pc)^K_PCR)

##------------------------------------------------------------------------------
##############################################
## FORCE OF INFECTION EXPERIENCED BY HUMANS ##
##############################################
##------------------------------------------------------------------------------

b  <- user()   # mosquito-to-human infection probability
dim(FOI_lag) <- c(na,nh,num_int)
FOI_lag[1:na, 1:nh, 1:num_int] <- EIR[i,j,k] * b
# Unlike in P. falciparum model, b does not depend on level of infection blocking immunity

# Current FOI depends on humans that have been through the latent period
dE <- user() # latent period of human infection
dim(FOI) <- c(na,nh,num_int)
FOI[,,] <- delay(FOI_lag[i,j,k],dE)

# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()
dim(EIR) <- c(na,nh,num_int)
EIR[1:na,1:nh,1:num_int] <- av_human[k] * rel_foi[j] * foi_age[i] * Iv/omega
# av_human is the feeding rate, which depends on the human blood index Q, the
# mean time between feeds gonotrophic cycle tau1+tau2 and vector control interventions.
# rel_foi and foi_age account for biting heterogeneity and the relationship between
# biting and age.

output(Ivout) <- Iv
output(omega) <- omega
output(EIR[,,]) <- EIR

##------------------------------------------------------------------------------
##########################
## SEASONALITY FUNCTION ##
##########################
##------------------------------------------------------------------------------

# Seasonality is added into the model using a Fourier series that was fit to rainfall at every admin 1 level
pi <- user() # weird quirk, need to pass pi

# The parameters for the fourier series
ssa0 <- user()
ssa1 <- user()
ssa2 <- user()
ssa3 <- user()
ssb1 <- user()
ssb2 <- user()
ssb3 <- user()
theta_c <- user()
# Recreation of the rainfall function
theta2 <- if(ssa0 == 0 && ssa1  == 0 && ssa2  == 0 && ssb1  == 0 && ssb2  == 0 && ssb3  == 0 && theta_c  == 0)
  1 else max((ssa0+ssa1*cos(2*pi*t/365)+ssa2*cos(2*2*pi*t/365)+ssa3*cos(3*2*pi*t/365)+ssb1*sin(2*pi*t/365)+ssb2*sin(2*2*pi*t/365)+ ssb3*sin(3*2*pi*t/365) ) /theta_c,0.001)

##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# Largely based on supplementary materials and associated code from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
# Minor modifications as detailed in https://www.nature.com/articles/s41467-018-05860-8
# This is structurally almost identical for the P. vivax and P. falciparum model

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# initial state values:
init_Sv <- user()
init_Ev <- user()
init_Iv <- user()
initial(Sv) <- init_Sv * mv0
initial(Ev) <- init_Ev * mv0
#initial(Ev[1:10]) <- init_Ev/10 * mv0 # Options if not using a delayed delay
#dim(Ev) <- 10
initial(Iv) <- init_Iv * mv0

# Infectiousness of humans to mosquitoes
cPCR <- user()  # infectiousness of PCR-detectable infection
cD <- user()    # infectiousness of untreated clinical infection
cT <- user()    # infectiousness of treated clinical infection
cLM <- user()   # infectiousness of LM-detectable infection

# Force of infection experienced by mosquitoes
dim(FOIvijk) <- c(na,nh,num_int,nk)
omega <- user()  # normalising constant for biting rates
FOIvijk[1:na, 1:nh, 1:num_int, 1:nk] <- (cT*T[i,j,k,l] + cD*I_D[i,j,k,l] + cLM*I_LM[i,j,k,l] +
                                           cPCR*I_PCR[i,j,k,l]) * av_mosq[k] * rel_foi[j] * foi_age[i]/omega
# av_mosq is the feeding rate which depends on the human blood index Q, the
# mean time between feeds/gonotrophic cycle tau1+tau2 and vector control interventions
lag_FOIv <- sum(FOIvijk)

FOIv <- lag_FOIv
# Difference to P. falciparum model: the tI parameter (delayGam) representing the period
# between infection and infectiousness in humans (development of gametocytes) does not
# exist for P. vivax (assumed to be instantaneous).

delayMos <- user() # Extrinsic incubation period

# Number of mosquitoes that become infected at each time point
surv <- exp(-mu*delayMos)
ince <- FOIv * Sv
lag_incv <- ince * surv
incv <- delay(lag_incv, delayMos)
#incv <- lag_incv

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa <- 0.5*PL/dPL
#betaa <- mv0 * mu0 * theta2

deriv(Sv) <- -ince - mu*Sv + betaa
deriv(Ev) <- ince - incv - mu*Ev
deriv(Iv) <- incv - mu*Iv

# Total mosquito population
mv = Sv+Ev+Iv

# model options if don't want to use a delayed delay
#deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
#deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
#mv = Sv+sum(Ev)+Iv


##------------------------------------------------------------------------------
###################
## LARVAL STATES ##
###################
##------------------------------------------------------------------------------

# Model by White et al.
# (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)
# No changes for P. vivax

# EL - early larval instar stage
# LL - late larval instar stage
# PL - pupal stage

# mean carrying capacity from initial mosquito density:
dLL <- user() # development time of larvae
dPL <- user() #development time of pupae
dEL <- user() #development time of early stage
muLL <- user() #daily density dep. mortality rate of larvae
muPL <- user() #daily den. dep. mortality rate of pupae
muEL <- user() #daily den. dep. mortality rate of early stage
gammaL <- user() # eff. of den. dep. on late stage relative to early stage

# fitted entomological parameters:
mv0 <- user()  # initial mosquito density
mu0 <- user()  # baseline mosquito death rate
tau1 <- user() # duration of host-seeking behaviour
tau2 <- user() # duration of resting behaviour
p10 <- user()  # prob of surviving 1 feeding cycle
p2 <- user()   # prob of surviving one resting cycle
betaL <- user() # maximum number of eggs per oviposition per mosq

# Entomological variables:
eov <- betaL/mu*(exp(mu/fv)-1)
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mu0*dLL*(1+dPL*muPL)))
K0 <- 2*mv0*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
KL <- K0*theta2
fv <- 1/( tau1/(1-zbar) + tau2 ) # mosquito feeding rate (zbar from intervention param.)
mu <- -fv*log(p1*p2) # mosquito death rate

# finding equilibrium and initial values for EL, LL & PL
init_PL <- user()
initial(PL) <- init_PL
init_LL <- user()
initial(LL) <- init_LL
init_EL <- user()
initial(EL) <- init_EL

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
deriv(EL) <- beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL
# egg hatching - den. dep. mortality - maturing larvae
deriv(LL) <- EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL
# pupae - mortality - fully developed pupae
deriv(PL) <- LL/dLL - muPL*PL - PL/dPL

##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ##
########################
##------------------------------------------------------------------------------

# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
# No changes for P. vivax

# general parameters
ITN_IRS_on <- user() # days after which interventions begin
num_int <- user() # number of intervention categorys, ITN only, IRS only, neither, both
itn_cov <- user() # proportion of population covered by ITN
irs_cov <- user() # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
dim(cov_) <- 4
cov_[1] <- (1-itn_cov)*(1-irs_cov)  # {No intervention}
cov_[2] <- itn_cov*(1-irs_cov) # 	   {ITN only}
cov_[3] <- (1-itn_cov)*irs_cov	#      {IRS only}
cov_[4] <- itn_cov*irs_cov #	   {Both ITN and IRS}
cov[] <- cov_[i]
dim(cov) <- num_int

IRS_interval <- user() # how long IRS lasts
ITN_interval <- user() # how long ITN lasts
chi <- user() # proportion of vector endophily
Q0 <- user() # proportion of anthropophagy
bites_Bed <- user() # endophagy in bed
bites_Indoors <- user() # endophagy indoors

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN/IRS
# d - probability of dying after hitting ITN/IRS
# s - probability of successful feed after hitting ITN/IRS

# The maximum (and then minimum) r and d values for ITN/IRS on day 0 before they decay
r_ITN0 <- user()
d_ITN0 <- user()
d_IRS0 <- user()
r_IRS0 <- user()
r_ITN1 <- user()
irs_loss <- user()
itn_loss <- user()

# Calculates decay for ITN/IRS
ITN_decay = if(t < ITN_IRS_on) 0 else exp(-((t-ITN_IRS_on)%%ITN_interval) * itn_loss)
IRS_decay = if(t < ITN_IRS_on) 0 else exp(-((t-ITN_IRS_on)%%IRS_interval) * irs_loss)

# The r,d and s values turn on after ITN_IRS_on and decay accordingly
d_ITN <- if(t < ITN_IRS_on) 0 else d_ITN0*ITN_decay
r_ITN <- if(t < ITN_IRS_on) 0 else r_ITN1 + (r_ITN0 - r_ITN1)*ITN_decay
s_ITN <- if(t < ITN_IRS_on) 1 else 1 - d_ITN - r_ITN

r_IRS <- if(t < ITN_IRS_on) 0 else r_IRS0*IRS_decay
d_IRS <- if(t < ITN_IRS_on) 0 else chi*d_IRS0*IRS_decay
s_IRS <- if(t < ITN_IRS_on) 1 else 1 - d_IRS

# probability that mosquito bites and survives for each intervention category
dim(w_) <- 4
w_[1] <- 1
w_[2] <- 1 - bites_Bed + bites_Bed*s_ITN
w_[3] <- 1 - bites_Indoors + bites_Indoors*(1-r_IRS)*s_IRS
w_[4] <- 1 - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN*s_IRS + (bites_Indoors - bites_Bed)*(1-r_IRS)*s_IRS
w[] <- w_[i]
dim(w) <- num_int

# probability that mosq feeds during a single attempt for each int. cat.
dim(yy_) <- 4
yy_[1] <- 1
yy_[2] <- w_[2]
yy_[3] <- 1 - bites_Indoors + bites_Indoors*(1-r_IRS)
yy_[4] <- 1 - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN + (bites_Indoors - bites_Bed)*(1-r_IRS)
yy[] <- yy_[i]
dim(yy) <- num_int

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z_) <- 4
z_[1] <- 0
z_[2] <- bites_Bed*r_ITN
z_[3] <- bites_Indoors*r_IRS
z_[4] <- bites_Bed*(r_IRS+ (1-r_IRS)*r_ITN) + (bites_Indoors - bites_Bed)*r_IRS
z[] <- z_[i]
dim(z) <- num_int

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- num_int
dim(whi) <- num_int
zhi[1:num_int] <- cov[i]*z[i]
whi[1:num_int] <- cov[i]*w[i]
zh <- if(t < ITN_IRS_on) 0 else sum(zhi)
wh <- if(t < ITN_IRS_on) 1 else sum(whi)
# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar <- Q0*zh
# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar <- 1 - Q0 + Q0*wh

# p1 is the updated p10 given that interventions are now in place:
p1 <- wbar*p10/(1-zbar*p10)
Q <- 1-(1-Q0)/wbar # updated anthropophagy given interventions
av <- fv*Q # biting rate on humans
dim(av_mosq) <- num_int
av_mosq[1:num_int] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.
dim(av_human) <- num_int
av_human[1:num_int] <- av*yy[i]/wh # biting rate on humans in each int. cat.

##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------

# Outputs for each compartment across the sum across all ages, biting heterogeneities and intervention categories
output(Sout) <- sum(S[,,,])
output(Tout) <- sum(T[,,,])
output(I_Dout) <- sum(I_D[,,,])
output(I_LMout) <- sum(I_LM[,,,])
output(I_PCRout) <- sum(I_PCR[,,,])
output(Pout) <- sum(P[,,,])

# Outputs for clinical incidence and prevalence on a given day
# population densities for each age category
den[] <- user()
dim(den) <- na
# index of the age vector above 59 months
age59 <- user(integer=TRUE)
# index of the age vector above 5 years
age05 <- user(integer=TRUE)

# Light microscopy prevalence in <59 month year olds
dim(prev_LM_0to59) <- c(age59,nh,num_int,nk)
prev_LM_0to59[1:age59,,,] <- T[i,j,k,l] + I_D[i,j,k,l]  + I_LM[i,j,k,l]
output(prevLM_0to59) <- sum(prev_LM_0to59[,,,])/sum(den[1:age59])

# Light microscopy prevalence across all ages
# This is most commonly measured for P. vivax
dim(prev_LM) <- c(na,nh,num_int,nk)
prev_LM[,,,] <- T[i,j,k,l] + I_D[i,j,k,l]  + I_LM[i,j,k,l]
output(prevLM) <- sum(prev_LM[,,,])/sum(den[])

# PCR prevalence in <59 month year olds
dim(prev_PCR_0to59) <- c(age59,nh,num_int,nk)
prev_PCR_0to59[1:age59,,,] <- T[i,j,k,l] + I_D[i,j,k,l]  + I_LM[i,j,k,l]  + I_PCR[i,j,k,l]
output(prevPCR_0to59) <- sum(prev_PCR_0to59[,,,])/sum(den[1:age59])

# PCR prevalence across all ages
dim(prev_PCR) <- c(na,nh,num_int,nk)
prev_PCR[,,,] <- T[i,j,k,l] + I_D[i,j,k,l]  + I_LM[i,j,k,l] + I_PCR[i,j,k,l]
output(prevPCR) <- sum(prev_PCR[,,,])/sum(den[])

# Clinical incidence in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,num_int,nk)
clin_inc0to5[1:age05,,,] <- clin_inc[i,j,k,l]
output(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])

# Clinical incidence across all ages
output(inc) <- sum(clin_inc[,,,])

# Param checking outputs
output(mu) <- mu
output(beta_larval) <- beta_larval
output(KL) <- KL
output(mv) <- mv
output(Q) <- Q
output(wh) <- wh
output(d_ITN) <- d_ITN
output(r_ITN) <- r_ITN
output(s_ITN) <- s_ITN
output(d_IRS) <- d_IRS
output(r_IRS) <- r_IRS
output(s_IRS) <- s_IRS
output(cov[]) <- TRUE
output(K0) <- K0
