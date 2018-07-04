## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
ft <- user() # proportion of cases treated

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------

# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

# fitted parameters for human compartments:
eta <- user() # death rate for exponential population distribution
age_rate[] <- user() # rate at which humans move through age categories
dim(age_rate) <- na
het_wt[] <- user() # weights of the different heterogeneous biting categories
dim(het_wt) <- nh
rA <- user() # rate of movement from A -> U
rT <- user() # rate of treatment working: T -> P
rD <- user() #  rate from D -> A
rU <- user() # rate of clearance of subpatent infection U -> S
rP <- user() # rate at which prophylaxis wears off P -> S

# S - SUSCEPTIBLE
init_S[,,] <- user()
dim(init_S) <- c(na,nh,num_int)
initial(S[,,]) <- init_S[i,j,k]
dim(S) <- c(na,nh,num_int)

deriv(S[1, 1:nh, 1:num_int]) <- -FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] +
  cov[k]*eta*H*het_wt[j] - (eta+age_rate[i])*S[i,j,k]
deriv(S[2:na, 1:nh, 1:num_int]) <- -FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] -
  (eta+age_rate[i])*S[i,j,k] + age_rate[i-1]*S[i-1,j,k]

# T- SUCCESSFULLY TREATED
init_T[,,] <- user()
dim(init_T) <- c(na,nh,num_int)
initial(T[,,]) <- init_T[i,j,k]
dim(T) <- c(na,nh,num_int)

deriv(T[1, 1:nh, 1:num_int]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k]
deriv(T[2:na, 1:nh, 1:num_int]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k] + age_rate[i-1]*T[i-1,j,k]

# D - CLEAR DISEASE
init_D[,,] <- user()
dim(init_D) <- c(na,nh,num_int)
initial(D[,,]) <- init_D[i,j,k]
dim(D) <- c(na,nh,num_int)

deriv(D[1, 1:nh, 1:num_int]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k]
deriv(D[2:na, 1:nh, 1:num_int]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k] + age_rate[i-1]*D[i-1,j,k]

# A - ASYMPTOMATIC DISEASE
init_A[,,] <- user()
dim(init_A) <- c(na,nh,num_int)
initial(A[,,]) <- init_A[i,j,k]
dim(A) <- c(na,nh,num_int)

deriv(A[1, 1:nh, 1:num_int]) <- (1-phi[i,j,k])*FOI[i,j,k]*Y[i,j,k] - FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k]
deriv(A[2:na, 1:nh, 1:num_int]) <- (1-phi[i,j,k])*FOI[i,j,k]*Y[i,j,k] - FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k] + age_rate[i-1]*A[i-1,j,k]

# U - SUBPATENT DISEASE
init_U[,,] <- user()
dim(init_U) <- c(na,nh,num_int)
initial(U[,,]) <- init_U[i,j,k]
dim(U) <- c(na,nh,num_int)

deriv(U[1, 1:nh, 1:num_int]) <- rA*A[i,j,k] - FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k]
deriv(U[2:na, 1:nh, 1:num_int]) <- rA*A[i,j,k] - FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k] + age_rate[i-1]*U[i-1,j,k]

# P - PROPHYLAXIS
init_P[,,] <- user()
dim(init_P) <- c(na,nh,num_int)
initial(P[,,]) <- init_P[i,j,k]
dim(P) <- c(na,nh,num_int)

deriv(P[1, 1:nh, 1:num_int]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k]
deriv(P[2:na, 1:nh, 1:num_int]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k] +
  age_rate[i-1]*P[i-1,j,k]

# The number of individuals able to acquire clinical malaria
dim(Y) <- c(na,nh,num_int)
Y[1:na, 1:nh, 1:num_int] <- S[i,j,k]+A[i,j,k]+U[i,j,k]

# The number of new cases at this timestep
dim(clin_inc) <- c(na,nh,num_int)
clin_inc[1:na, 1:nh, 1:num_int] <- phi[i,j,k]*FOI[i,j,k]*Y[i,j,k]

# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph

##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

# fitted immunity parameters:
dCM <- user() # decay of maternal immunity
uCA <- user() # scale parameter (see Supplementary mats. 3.1.2)
dCA <- user() # decay for clinical immunity
dB <- user() # decay for infection blocking immunity
uB <- user() # scale param for IB immunity
dID <- user() # decay for detection immunity
uD <- user() # scale param for ID immunity
x_I[] <- user() # intermediate variable for calculating immunity functions
dim(x_I) <- na

# ICM - maternally acquired immunity
init_ICM[,,] <- user()
dim(init_ICM) <- c(na,nh,num_int)
initial(ICM[,,]) <- init_ICM[i,j,k]
dim(ICM) <- c(na,nh,num_int)

deriv(ICM[1, 1:nh, 1:num_int]) <- -1/dCM*ICM[i,j,k] + (init_ICM[i,j,k]-ICM[i,j,k])/x_I[i]
deriv(ICM[2:na, 1:nh, 1:num_int]) <- -1/dCM*ICM[i,j,k] - (ICM[i,j,k]-ICM[i-1,j,k])/x_I[i]

# ICA - exposure driven immunity
init_ICA[,,] <- user()
dim(init_ICA) <- c(na,nh,num_int)
initial(ICA[,,]) <- init_ICA[i,j,k]
dim(ICA) <- c(na,nh,num_int)

deriv(ICA[1, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] -ICA[i,j,k]/x_I[i]
deriv(ICA[2:na, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i]

# clinical immunity is a combination of maternal and exposure-driven immunity
dim(IC) <- c(na,nh,num_int)
IC[,,] <- ICM[i,j,k] + ICA[i,j,k]

# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- user()
phi1 <- user() # these parameters characterise the hill function
IC0 <- user() # for probability of clinical disease
kC <- user() # See supplementary materials 1.1.3
dim(phi) <- c(na,nh,num_int)
phi[1:na,1:nh,1:num_int] <- phi0*((1-phi1)/(1+(IC[i,j,k]/IC0)^kC) + phi1)

# IB - infection blocking immunity
init_IB[,,] <- user()
dim(init_IB) <- c(na,nh,num_int)
initial(IB[,,]) <- init_IB[i,j,k]
dim(IB) <- c(na,nh,num_int)

deriv(IB[1, 1:nh, 1:num_int]) <- EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - IB[i,j,k]/x_I[i]
deriv(IB[2:na, 1:nh, 1:num_int]) <- EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - (IB[i,j,k]-IB[i-1,j,k])/x_I[i]

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- user() # these parameters characterise the hill function for b
b1 <- user() # prob of infection from bite with zero immunity
kB <- user() #
IB0 <- user()
dim(b) <- c(na,nh,num_int)
b[1:na, 1:nh, 1:num_int] <- b0 * ((1-b1)/(1+(IB[i,j,k]/IB0)^kB)+b1)

# detection immunity
init_ID[,,] <- user()
dim(init_ID) <- c(na,nh,num_int)
initial(ID[,,]) <- init_ID[i,j,k]
dim(ID) <- c(na,nh,num_int)

deriv(ID[1, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k]*uD + 1) - ID[i,j,k]/dID - ID[i,j,k]/x_I[i]
deriv(ID[2:na, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k]*uD + 1) - ID[i,j,k]/dID - (ID[i,j,k]-ID[i-1,j,k])/x_I[i]

# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- user() # no idea where these are from
fD0 <- user() # or who fit them
gammaD <- user() #
d1 <- user()
ID0 <- user()
kD <- user()
dim(age) <- na
age[] <- user() # vector of age categories supplied by user

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,num_int)
p_det[,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k]/ID0)^kD)

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int)
FOI_lag[1:na, 1:nh, 1:num_int] <- EIR[i,j,k] * (if(IB[i,j,k]==0) b0 else b[i,j,k])

# Current FOI depends on humans that have been through the latent period and are
# producing gametocytes
dE <- user() # length of time from infection to gametocytogenesis
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
EIR[,,] <- av_human[k] * rel_foi[j] * foi_age[i]/omega*Iv

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
   1 else if(t < (14*365)) max((ssa0+ssa1*cos(2*pi*t/365)+ssa2*cos(2*2*pi*t/365)+ssa3*cos(3*2*pi*t/365)+ssb1*sin(2*pi*t/365)+ssb2*sin(2*2*pi*t/365)+ ssb3*sin(3*2*pi*t/365) ) /theta_c,0.001) else
    theta2_scaled

theta2_raw <- 6.4148 + 3.017*cos(2*pi*t/365) -3.9303*cos(4*pi*t/365) + 0.8912*cos(6*pi*t/365) + 4.209*sin(2*pi*t/365) - 1.2505*sin(4*pi*t/365) - 2.0129*sin(6*pi*t/365)
th_min <- -1.432255
th_max <- 20.67784
theta2_scaled <- max(2*(theta2_raw-th_min)/(th_max-th_min),0.001)
# theta2 <- if(theta2_scaled<0.001) 0.001 else theta2_scaled

##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# initial state values:
init_Sv <- user()
init_Ev <- user()
init_Iv <- user()
initial(Sv) <- init_Sv * mv0
# initial(Ev) <- init_Ev * mv0
initial(Ev[1:10]) <- init_Ev/10 * mv0
dim(Ev) <- 10
initial(Iv) <- init_Iv * mv0

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- user() # infectiousness U -> mosq
cD <- user() # infectiousness D -> mosq
cT <- user() # T -> mosq
gamma1 <- user() # fitted value of gamma1 characterises cA function
dim(cA) <- c(na,nh,num_int)
cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

# Force of infection from humans to mosquitoes
dim(FOIvijk) <- c(na,nh,num_int)
omega <- user() #normalising constant for biting rates
FOIvijk[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq[k]*foi_age[i]/omega
lag_FOIv=sum(FOIvijk)

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- user() # latent period in gametocytogenesis
delayMos <- user() # latent period in humans
FOIv <- delay(lag_FOIv, delayGam)

# Number of mosquitoes that become infected at each time point
surv <- exp(-mu*delayMos)
ince <- FOIv * Sv
lag_incv <- ince * surv
incv <- delay(lag_incv, delayMos)

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa <- 0.5*PL/dPL
#betaa <- mv0 * mu0 * theta2

deriv(Sv) <- -ince - mu*Sv + betaa - feb*Sv + inhib_rate*SvI
#deriv(Ev) <- ince - incv - mu*Ev
deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1] - feb*Ev[1] + inhib_rate*EvI[1]
deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i] - feb*Ev[i] + inhib_rate*EvI[i]
deriv(Iv) <- incv - mu*Iv - feb*Iv + inhib_rate*IvI


# A proportion of what we would consider standard repellency is in fact full feeding inhibition
feb <- if(t < ITN_on) 0 else (cov[2])*bites_Bed*(1-d_ITN-r_ITN)*inhib*inhibition_effect*surv_bioassay


output(f_ITN) <- f_ITN


# PBO
PBO_p <- if(PBO==0) 0 else -0.371
PBO_int <- if(PBO==0) 0 else 1.089
# Washes
Washes <- 0 # 0 washes
# Washes <- 1.31 # 3 washes
# Washes <- -0.379 # 20 washes

p <- 2.363 - 2.57*surv_bioassay + PBO_p + PBO_int*surv_bioassay
inhib <- if(t < ITN_on) 0 else exp(p)/(exp(p)+1)

inhibition_effect <- user()
inhib_length <- user()
inhib_rate <- 1/inhib_length
output(feb) <- feb
output(inhib_rate) <- inhib_rate
dim(EvI) <- 10
initial(SvI) <- 0
initial(EvI[1:10]) <- 0
initial(IvI) <- 0

deriv(SvI) <- feb*Sv - inhib_rate*SvI - mu0*SvI
deriv(EvI[1]) <- feb*Ev[1] - EvI[1] - mu0*EvI[1] - inhib_rate*EvI[1]
deriv(EvI[2:10]) <- feb*Ev[i] + EvI[i-1] - EvI[i] - mu0*EvI[i] - inhib_rate*EvI[i]
deriv(IvI) <- feb*Iv - inhib_rate*IvI - mu0*IvI


# Total mosquito population
#mv = Sv+Ev+Iv
mv = Sv+sum(Ev)+Iv
mvI = Sv+sum(Ev)+Iv + SvI + sum(EvI) + IvI
output(mvI) <- mvI
##------------------------------------------------------------------------------
###################
## LARVAL STATES ##
###################
##------------------------------------------------------------------------------

# Model by White et al.
# (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)

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
mv0 <- user() # initial mosquito density
mu0 <- user() # baseline mosquito death rate
tau1 <- user() # duration of host-seeking behaviour
tau2 <- user() # duration of resting behaviour
p10 <- user() # prob of surviving 1 feeding cycle
p2 <- user() #prob of surviving one resting cycle
betaL <- user() # maximum number of eggs per oviposition per mosq

# Entomological variables:
eov <- betaL/mu*(exp(mu/fv)-1) # eggs of
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mu0*dLL*(1+dPL*muPL)))
K0 <- 2*mv0*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
KL <- K0*theta2
fv <- 1/( tau1/(1-zbar) + tau2 ) # mosquito feeding rate (zbar from intervention param.)
# My change below, p2 is probability of surviving through resting period - now affected by indoor spatial repellents
#mu <- -fv*log(p1*p2) # mosquito death rate
mu <- -fv*log(p1*p2tox)

# finding equilibrium and initial values for EL, LL & PL
init_PL <- user()
initial(PL) <- init_PL
init_LL <- user()
initial(LL) <- init_LL
init_EL <- user()
initial(EL) <- init_EL

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
deriv(EL) <- beta_larval*(1-f_red)*mvI-muEL*(1+(EL+LL)/KL*(1-f_red))*EL - EL/dEL
# egg hatching - den. dep. mortality - maturing larvae
deriv(LL) <- EL/dEL - muLL*(1+gammaL*(EL + LL)/KL *(1-f_red))*LL - LL/dLL
# pupae - mortality - fully developed pupae
deriv(PL) <- LL/dLL - muPL*PL - PL/dPL

##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ##
########################
##------------------------------------------------------------------------------

# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# general parameters
ITN_on <- user() # day when ITNs are first introduced
EM_on <- user() # day when emanators are first introduced
num_int <- user() # number of intervention categorys, ITN only, emanator only, neither, both
itn_cov <- user() # proportion of population covered by ITN
em_cov <- user() # proportion of population covered by emanator
irs_cov <- user() # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
dim(cov) <- num_int
# cov[1] <- (1-itn_cov)*(1-em_cov)  # {No intervention}
# cov[2] <- itn_cov*(1-em_cov) # 	   {ITN only}
# cov[3] <- (1-itn_cov)*em_cov	#      {EM only}
# cov[4] <- itn_cov*em_cov #	   {Both ITN and EM}
# cov[] <- user()
cov[3] <- 0
cov[4] <- 0
cov[2] <- if(t < (365*13)) itn_cov_start else if(t < ((365*13)+60)) itn_cov_start + itn_cov_grad_up*(t-(365*13)) else if(t < (365*14)) itn_cov_max else max(itn_cov_max - (t-(365*14))*itn_cov_grad_down,0.001)
cov[1] <- 1 - cov[2]

itn_cov_start <- if(PBO) 0.32 else 0.39
itn_cov_grad_up <- if(PBO) 0.52/60 else 0.45/60
itn_cov_max <- if(PBO) 0.78 else 0.75
itn_cov_grad_down <- if(PBO) 0.13/365 else 0.24/365

EM_interval <- user() # how long until emanators are refreshed
ITN_interval <- user() # how long until ITNs are repeated
IRS_interval <- user() # how long until IRS is repeated
chi <- user() # proportion of vector endophily
Q0 <- user() # proportion of anthropophagy
bites_Bed <- user() # endophagy in bed
bites_Indoors <- user() # endophagy indoors
bites_Emanator <- user() # endophagy whilst outside, near an emanator
dim(em_human_dist) <- d_len
em_human_dist[] <- user() # distribution of population time spent at each distance from emanator
dim(em_prod_profile) <- d_len
em_prod_profile[] <- user() # proportion of bites averted by emanator at each distance from emanator


## ELLIE'S WORK ##

# ITN parameters with pyrethroid resistance

# Linking cone assay and hut trial work on resistance to d_ITN and r_ITN, done by Ellie
surv_bioassay <- user() # measure of % survival in discriminating dose bioassay
# surv_bioassay <- if(t < 365*19) 0 else 0.5
# surv_bioassay <- if(t < (365*13)) 0 else if(t < ((365*14))) 0.001*(t-(365*13)) else 0.41
# surv_bioassay <- 0.92
output(surv_bioassay) <- surv_bioassay
PBO <- user()
output(PBO) <- PBO
pbo_benefit_a <- 3.407+5.88*((1-surv_bioassay)-0.5)/(1+0.783*((1-surv_bioassay)-0.5))
pbo_benefit <- exp(pbo_benefit_a)/(1+exp(pbo_benefit_a))
mort_assay <- if(PBO==0) 1 - surv_bioassay else pbo_benefit
output(mort_assay) <- mort_assay
output(pbo_benefit) <- pbo_benefit

# Relationship between mortality in bioassay to hut trial, logit scale
mort_hut_a <- 0.6338 + 3.9970 * (mort_assay-0.5)
mort_hut <- exp(mort_hut_a)/(1+exp(mort_hut_a))

# Relationship between hut trial mortality and deterrence
det_hut_a <- 0.07117+1.257*(mort_hut-0.5)-1.517*(mort_hut-0.5)^2
det_hut <- if(det_hut_a < 0) 0 else det_hut_a # censored to stop becoming negative
my_death <- mort_assay*(1-det_hut)
my_kill_det <- 1 - my_death - det_hut
output(my_kill_det) <- my_kill_det
output(my_death) <- my_death
output(det_hut) <- det_hut
# Relationship between hut trial mortality and successful (feed)
suc_hut <- 0.02491*exp(3.317*(1-mort_hut))
rep_hut <- 1-suc_hut-mort_hut

n1n0 <- 1-det_hut
kp1 <- n1n0*suc_hut
jp1 <- n1n0*rep_hut+(n1n0)
lp1 <- n1n0*mort_hut

# New values given resistance
r_ITN0 <- (1-kp1/0.699)*(jp1/(lp1+jp1))
d_ITN0 <- (1-kp1/0.699)*(lp1/(lp1+jp1))

# Insecticide halflife
hut_max_a <- 0.6338 + 3.9970*(1-0.5) # maximum mortality seen in huts
hut_max <- exp(hut_max_a)/(1+exp(hut_max_a))

my_max_washes_a <- -2.360 - 3.048*(hut_max-0.5)
my_max_washes <- log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))

wash_decay_rate_a <- -2.360+-3.048*(mort_hut-0.5)
wash_decay_rate   <- exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a))
itn_half_life     <- (log(2)/wash_decay_rate)/my_max_washes*2.64*365

itn_loss <- log(2)/itn_half_life
r_ITN_min <- 0.24
output(itn_loss) <- itn_loss
output(itn_half_life) <- itn_half_life

## END OF ELLIE'S WORK ##

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN/EM
# d - probability of dying after hitting ITN/EM
# s - probability of successful feed after hitting ITN/EM

# The maximum (and then minimum) r and d values for ITN/EM on day 0 before they decay
# r_ITN0 <- user()
# d_ITN0 <- user()
# r_ITN1 <- user()
em_loss <- user()
# itn_loss <- user()

# Calculates decay for ITN/EM
ITN_decay = if(t < ITN_on) 0 else exp(-((t-ITN_on)%%ITN_interval) * itn_loss)
EM_decay = if(t < EM_on) 0 else exp(-((t-EM_on)%%EM_interval) * em_loss)
output(em_loss) <- em_loss

# The r,d and s values turn on after ITN_EM_on and decay accordingly
d_ITN <- if(t < ITN_on) 0 else d_ITN0*ITN_decay
# r_ITN <- if(t < ITN_on) 0 else r_ITN1 + (r_ITN0 - r_ITN1)*ITN_decay
r_ITN <- if(t < ITN_on) 0 else r_ITN_min + (r_ITN0- r_ITN_min)*ITN_decay
# Feeding inhibition
# f_ITN <- if(t < ITN_on) 0 else (1-d_ITN-r_ITN)*inhib*inhibition_effect*surv_bioassay*ITN_decay
f_ITN <- if(t < ITN_on) 0 else f_ITN_min + (my_kill_det-f_ITN_min)*inhib*inhibition_effect*surv_bioassay*ITN_decay
f_ITN_min <- if(PBO) 0 else 0
output(inhib) <- inhib
# Ellie's + my edit
s_ITN <- if(t < ITN_on) 1 else 1 - d_ITN - r_ITN - f_ITN

# EMANATOR OUTSIDE PARAMETERS
d_len <- user()
dim(rep_EM) <- d_len
# repellency of emanators at different distances
rep_EM[1:d_len] <- if(t < EM_on) 0 else em_prod_profile[i]#*EM_decay
dim(p_EM_vec) <- d_len
# time that humans spend at different distances away from emanator
p_EM_vec[1:d_len] <- em_human_dist[i]*(1-rep_EM[i])

# Values over all distances
r_EM_out0 <- 1-sum(p_EM_vec)
d_EM_out0 <- user()

r_EM_out <- if(t < EM_on) 0 else r_EM_out0*EM_decay

d_EM_out <- if(t < EM_on) 0 else d_EM_out0*EM_decay

s_EM_out <- if(t < EM_on) 1 else 1 - d_EM_out - r_EM_out

# EMANATOR INSIDE parameters

em_in <- user() # Inside effect toggle

r_EM_in0 <- user()
d_EM_in0 <- user()

r_EM_in <- if(t < EM_on) 0 else r_EM_in0*EM_decay

d_EM_in <- if(t < EM_on) 0 else d_EM_in0*EM_decay

s_EM_in <- if(t < EM_on) 1 else 1 - d_EM_in

## experimental fecundity/mortality stuff
f_EM_in0 <- user()
f_EM_in <- if(t < EM_on) 0 else f_EM_in0*EM_decay

# fecundity
# amount to reduce betaa by -> 1 - (% of mosquitoes surviving a feeding attempt in a home with an emanator in * proportion of mosquitoes prevented from laying eggs)
f_red <- if(em_in == 0) 0 else ((cov[3]*w[3]) + (cov[4] * w[4]))*f_EM_in
output(f_red) <- f_red # this is what betaa, the number of new mosquitoes produced is modified by

# toxicity
# affects the probability that the mosquito will survive the subsequent resting period
t_EM_in0 <- user()
t_EM_in <- if(t < EM_on) 0 else t_EM_in0*EM_decay
tox <- cov[1]+cov[2]+(cov[3]*(1-t_EM_in))+(cov[4]*(1-t_EM_in))
p2tox <- p2*tox



# probability that mosquito bites and survives for each intervention category
dim(w) <- num_int
w[1] <- 1
w[2] <- 1 - bites_Bed + bites_Bed*s_ITN
w[3] <- if(em_in == 0) 1 - bites_Emanator + s_EM_out*bites_Emanator else 1 - bites_Emanator + s_EM_out*bites_Emanator - bites_Indoors + s_EM_in*(1-r_EM_in)*bites_Indoors
#w[3] <- 1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + (1-r_EM_in)*bites_Indoors
w[4] <- if(em_in == 0) 1 - bites_Bed - bites_Emanator + bites_Bed*s_ITN + s_EM_out*bites_Emanator else (1 - bites_Indoors + bites_Bed*(1-r_EM_in)*s_EM_in*s_ITN + (bites_Indoors-bites_Bed)*(1-r_EM_in)*s_EM_in - bites_Emanator + s_EM_out*bites_Emanator)
#w[4] <- 1 - bites_Indoors + bites_Bed*(1-r_EM_in)*s_ITN + (bites_Indoors-bites_Bed)*(1-r_EM_in) - bites_Emanator + p_EM*bites_Emanator

# probability that mosq feeds during a single attempt for each int. cat.
dim(yy) <- num_int
yy[1] <- 1
yy[2] <- w[2]
yy[3] <- if(em_in == 0) 1 - bites_Emanator + s_EM_out*bites_Emanator else 1 - bites_Emanator + s_EM_out*bites_Emanator - bites_Indoors + (1-r_EM_in)*bites_Indoors
yy[4] <- if(em_in == 0) 1 - bites_Bed - bites_Emanator + bites_Bed*s_ITN + s_EM_out*bites_Emanator else (1 - bites_Indoors + bites_Bed*(1-r_EM_in)*s_ITN + (bites_Indoors-bites_Bed)*(1-r_EM_in) - bites_Emanator + s_EM_out*bites_Emanator)

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z) <- num_int
z[1] <- 0
z[2] <- bites_Bed*r_ITN
z[3] <- if(em_in == 0) r_EM_out*bites_Emanator else r_EM_out*bites_Emanator + bites_Indoors*r_EM_in
#z[3] <- (1-p_EM)*bites_Emanator + bites_Indoors*r_EM_in
z[4] <- if(em_in == 0) bites_Bed*r_ITN + r_EM_out*bites_Emanator else bites_Bed*(r_EM_in + (1-r_EM_in)*r_ITN) + (bites_Indoors-bites_Bed)*r_EM_in + r_EM_out*bites_Emanator
#z[4] <- bites_Bed*(r_EM_in + (1-r_EM_in)*r_ITN) + (bites_Indoors-bites_Bed)*r_EM_in + (1-p_EM)*bites_Emanator

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- num_int
dim(whi) <- num_int
zhi[1:num_int] <- cov[i]*z[i]
whi[1:num_int] <- cov[i]*w[i]

zh <- if(t < min(ITN_on,EM_on)) 0 else sum(zhi)
wh <- if(t < min(ITN_on,EM_on)) 1 else sum(whi)
# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar <- Q0*zh
# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar <- 1 - Q0 + Q0*wh

# p1 is the updated p10 given that interventions are now in place:
p1 <- wbar*p10/(1-zbar*p10)
Q <- 1-(1-Q0)/wbar # updated anthropophagy given interventions
av <- fv*Q # biting rate on humans
dim(av_mosq) <- num_int
#av_mosq[1:num_int] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.
# ALTERED DUE TO PAGE 6 SUPP MAT 2, the biting rate was previously inflated to account for the fact that some mosquitoes would bite due to IRS and then die
# This essentially meant that as the biting rate on humans covered by emanators dropped, this increased the biting rate on people with no interventions
# av_human[1:num_int] <- av*yy[i]/wh # biting rate on humans in each int. cat.
av_mosq[1:num_int] <- av*w[i]/wh
dim(av_human) <- num_int
av_human[1:num_int] <- av*yy[i]/wh

##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------

# Outputs for each compartment across the sum across all ages, biting heterogeneities and intervention categories
output(Sout) <- sum(S[,,])
output(Tout) <- sum(T[,,])
output(Dout) <- sum(D[,,])
output(Aout) <- sum(A[,,])
output(Uout) <- sum(U[,,])
output(Pout) <- sum(P[,,])

# Outputs for clinical incidence and prevalence on a given day
# population densities for each age category
den[] <- user()
dim(den) <- na
# index of the age vector above 59 months
age59 <- user()
# index of the age vector above 5 years
age05 <- user()

# 2-10 year old outputs
# age2 <- 9
# age10 <- 17
# gap <- age10-age2
age15 <- 18
dim(prev0to15) <- c(age15,nh,num_int)
prev0to15[1:age15,,] <- T[i,j,k] + D[i,j,k] + A[i,j,k] * p_det[i,j,k]
output(prev15) <- sum(prev0to15[,,])/sum(den[1:age15])
# prev2to10[1:gap,,] <- T[age2+i,j,k] + D[age2+i,j,k] + A[age2+i,j,k] + p_det[age2+i,j,k]
# output(prev210) <- sum(prev2to10[,,])/sum(den[age2:age10])

dim(prevITN) <- c(age59,nh)
prevITN[1:age59,] <- T[i,j,2] + D[i,j,2] + A[i,j,2] * p_det[i,j,2]
output(previ) <- sum(prevITN[,])/sum(den[1:age59])

dim(prevnull) <- c(age59,nh)
prevnull[1:age59,] <- T[i,j,1] + D[i,j,1] + A[i,j,1] * p_det[i,j,1]
output(prevn) <- sum(prevnull[,])/sum(den[1:age59])

dim(prev0to59) <- c(age59,nh,num_int)
prev0to59[1:age59,,] <- T[i,j,k] + D[i,j,k] + A[i,j,k] * p_det[i,j,k]
output(prev) <- sum(prev0to59[,,])/sum(den[1:age59])

dim(prevall) <- c(na,nh,num_int)
prevall[,,] <- T[i,j,k] + D[i,j,k] + A[i,j,k] * p_det[i,j,k]
output(allprev) <- sum(prevall[,,])/sum(den[])

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,num_int)
dim(clin_inc_null) <- c(age05,nh)
dim(clin_inc_net) <- c(age05,nh)
clin_inc0to5[1:age05,,] <- clin_inc[i,j,k]
clin_inc_null[1:age05,] <- clin_inc[i,j,1]
clin_inc_net[1:age05,] <- clin_inc[i,j,2]
output(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])
output(inc) <- sum(clin_inc[,,])
output(inc05i) <- sum(clin_inc_net)/sum(den[1:age05])
output(inc05n) <- sum(clin_inc_null)/sum(den[1:age05])
dim(zz) <- num_int

output(Yout) <- sum(Y[,,])
output(phiout) <- sum(phi[,,])
output(FOIout) <- sum(FOI[,,])

output(EIR[,,]) <- EIR

output(clin_inc[,,]) <- clin_inc
output(cov[]) <- cov[i]
# Param checking outputs
output(Y[,,]) <- Y
output(phi[,,]) <- phi
output(KL) <- KL
output(mv) <- mv
output(Q) <- Q
output(wh) <- wh
output(wbar) <- wbar
output(av) <- av
output(av_mosq[]) <- av_mosq[i]
output(av_human[]) <- av_human[i]
output(w[]) <- w[i]
output(yy[]) <- yy[i]
output(zz[]) <- z[i]
output(zbar) <- zbar
output(d_ITN) <- d_ITN
output(r_ITN) <- r_ITN
output(s_ITN) <- s_ITN
output(cov[]) <- cov[i]
output(K0) <- K0
output(theta2) <- theta2
output(EM_decay) <- EM_decay
output(ITN_decay) <- ITN_decay
output(FOI[,,]) <- FOI
output(p1) <- p1
output(p2) <- p2
output(mu) <- mu
output(em_in) <- em_in
output(r_em_in) <- r_EM_in
output(d_em_in) <- d_EM_in
output(s_em_in) <- s_EM_in
output(r_em_out) <- r_EM_out
output(d_em_out) <- d_EM_out
output(s_em_out) <- s_EM_out
output(p2tox) <- p2tox
output(betaa) <- betaa
output(beta_larval) <- beta_larval
