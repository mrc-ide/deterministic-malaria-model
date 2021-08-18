## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
ft <- user() # proportion of cases treated

np <- user() # number of patches

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
init_S[,,,] <- user()
dim(init_S) <- c(na,nh,num_int,np)
initial(S[,,,]) <- init_S[i,j,k,l]
dim(S) <- c(na,nh,num_int,np)

deriv(S[1, 1:nh, 1:num_int, 1:np]) <- -FOI[i,j,k,l]*S[i,j,k,l] + rP*P[i,j,k,l] + rU*U[i,j,k,l] +
  cov[k,l]*eta*H[l]*het_wt[j] - (eta+age_rate[i])*S[i,j,k,l]
deriv(S[2:na, 1:nh, 1:num_int,1:np]) <- -FOI[i,j,k,l]*S[i,j,k,l] + rP*P[i,j,k,l] + rU*U[i,j,k,l] -
  (eta+age_rate[i])*S[i,j,k,l] + age_rate[i-1]*S[i-1,j,k,l]

# T- SUCCESSFULLY TREATED
init_T[,,,] <- user()
dim(init_T) <- c(na,nh,num_int,np)
initial(T[,,,]) <- init_T[i,j,k,l]
dim(T) <- c(na,nh,num_int,np)

deriv(T[1, 1:nh, 1:num_int,1:np]) <- ft*clin_inc[i,j,k,l] - rT*T[i,j,k,l] -
  (eta+age_rate[i])*T[i,j,k,l]
deriv(T[2:na, 1:nh, 1:num_int,1:np]) <- ft*clin_inc[i,j,k,l] - rT*T[i,j,k,l] -
  (eta+age_rate[i])*T[i,j,k,l] + age_rate[i-1]*T[i-1,j,k,l]

# D - CLEAR DISEASE
init_D[,,,] <- user()
dim(init_D) <- c(na,nh,num_int,np)
initial(D[,,,]) <- init_D[i,j,k,l]
dim(D) <- c(na,nh,num_int,np)

deriv(D[1, 1:nh, 1:num_int, 1:np]) <- (1-ft)*clin_inc[i,j,k,l] - rD*D[i,j,k,l] -
  (eta+age_rate[i])*D[i,j,k,l]
deriv(D[2:na, 1:nh, 1:num_int,1:np]) <- (1-ft)*clin_inc[i,j,k,l] - rD*D[i,j,k,l] -
  (eta+age_rate[i])*D[i,j,k,l] + age_rate[i-1]*D[i-1,j,k,l]

# A - ASYMPTOMATIC DISEASE
init_A[,,,] <- user()
dim(init_A) <- c(na,nh,num_int,np)
initial(A[,,,]) <- init_A[i,j,k,l]
dim(A) <- c(na,nh,num_int,np)

deriv(A[1, 1:nh, 1:num_int,1:np]) <- (1-phi[i,j,k,l])*FOI[i,j,k,l]*Y[i,j,k,l] - FOI[i,j,k,l]*A[i,j,k,l] +
  rD*D[i,j,k,l] - rA*A[i,j,k,l] - (eta+age_rate[i])*A[i,j,k,l]
deriv(A[2:na, 1:nh, 1:num_int, 1:np]) <- (1-phi[i,j,k,l])*FOI[i,j,k,l]*Y[i,j,k,l] - FOI[i,j,k,l]*A[i,j,k,l] +
  rD*D[i,j,k,l] - rA*A[i,j,k,l] - (eta+age_rate[i])*A[i,j,k,l] + age_rate[i-1]*A[i-1,j,k,l]

# U - SUBPATENT DISEASE
init_U[,,,] <- user()
dim(init_U) <- c(na,nh,num_int,np)
initial(U[,,,]) <- init_U[i,j,k,l]
dim(U) <- c(na,nh,num_int,np)

deriv(U[1, 1:nh, 1:num_int,1:np]) <- rA*A[i,j,k,l] - FOI[i,j,k,l]*U[i,j,k,l] - rU*U[i,j,k,l] -
  (eta+age_rate[i])*U[i,j,k,l]
deriv(U[2:na, 1:nh, 1:num_int,1:np]) <- rA*A[i,j,k,l] - FOI[i,j,k,l]*U[i,j,k,l] - rU*U[i,j,k,l] -
  (eta+age_rate[i])*U[i,j,k,l] + age_rate[i-1]*U[i-1,j,k,l]

# P - PROPHYLAXIS
init_P[,,,] <- user()
dim(init_P) <- c(na,nh,num_int,np)
initial(P[,,,]) <- init_P[i,j,k,l]
dim(P) <- c(na,nh,num_int,np)

deriv(P[1, 1:nh, 1:num_int,1:np]) <- rT*T[i,j,k,l] - rP*P[i,j,k,l] - (eta+age_rate[i])*P[i,j,k,l]
deriv(P[2:na, 1:nh, 1:num_int,1:np]) <- rT*T[i,j,k,l] - rP*P[i,j,k,l] - (eta+age_rate[i])*P[i,j,k,l] +
  age_rate[i-1]*P[i-1,j,k,l]

# The number of individuals able to acquire clinical malaria
dim(Y) <- c(na,nh,num_int,np)
Y[1:na, 1:nh, 1:num_int,1:np] <- S[i,j,k,l]+A[i,j,k,l]+U[i,j,k,l]

# The number of new cases at this timestep
dim(clin_inc) <- c(na,nh,num_int,np)
clin_inc[1:na, 1:nh, 1:num_int, 1:np] <- phi[i,j,k,l]*FOI[i,j,k,l]*Y[i,j,k,l]

# Sum compartments over all age, heterogeneity and intervention categories
dim(Sh)<- np
Sh[] <- sum(S[,,,i])
dim(Th)<- np
Th[] <- sum(T[,,,i])
dim(Dh)<- np
Dh[] <- sum(D[,,,i])
dim(Ah)<- np
Ah[] <- sum(A[,,,i])
dim(Uh)<- np
Uh[] <- sum(U[,,,i])
dim(Ph)<- np
Ph[] <- sum(P[,,,i])

dim(H)<- np
H[] <- Sh[i] + Th[i]+ Dh[i] + Ah[i] + Uh[i] + Ph[i]

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
age20l <- user(integer=TRUE) # lower index of age 20 age compartment
age20u <- user(integer=TRUE) # upper index of age 20 age compartment
age_20_factor <- user() # factor calculated in equilibrium solution
PM <- user() # immunity constant

# ICM - maternally acquired immunity
init_ICM[,,,] <- user()
dim(init_ICM) <- c(na,nh,num_int,np)
initial(ICM[,,,]) <- init_ICM[i,j,k,l]
dim(ICM) <- c(na,nh,num_int,np)
dim(init_ICM_pre) <- c(nh,num_int,np)
init_ICM_pre[1:nh,1:num_int,1:np] <- PM*(ICA[age20l,i,j,k] + age_20_factor*(ICA[age20u,i,j,k]-ICA[age20l,i,j,k]))

deriv(ICM[1, 1:nh, 1:num_int,1:np]) <- -1/dCM*ICM[i,j,k,l] + (init_ICM_pre[j,k,l]-ICM[i,j,k,l])/x_I[i]
deriv(ICM[2:na, 1:nh, 1:num_int,1:np]) <- -1/dCM*ICM[i,j,k,l] - (ICM[i,j,k,l]-ICM[i-1,j,k,l])/x_I[i]

# ICA - exposure driven immunity
init_ICA[,,,] <- user()
dim(init_ICA) <- c(na,nh,num_int,np)
initial(ICA[,,,]) <- init_ICA[i,j,k,l]
dim(ICA) <- c(na,nh,num_int,np)

deriv(ICA[1, 1:nh, 1:num_int,1:np]) <- FOI[i,j,k,l]/(FOI[i,j,k,l] * uCA + 1) - 1/dCA*ICA[i,j,k,l] -ICA[i,j,k,l]/x_I[i]
deriv(ICA[2:na, 1:nh, 1:num_int,1:np]) <- FOI[i,j,k,l]/(FOI[i,j,k,l] * uCA + 1) - 1/dCA*ICA[i,j,k,l] - (ICA[i,j,k,l]-ICA[i-1,j,k,l])/x_I[i]

# clinical immunity is a combination of maternal and exposure-driven immunity
dim(IC) <- c(na,nh,num_int,np)
IC[,,,] <- ICM[i,j,k,l] + ICA[i,j,k,l]

# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- user()
phi1 <- user() # these parameters characterise the hill function
IC0 <- user() # for probability of clinical disease
kC <- user() # See supplementary materials 1.1.3
dim(phi) <- c(na,nh,num_int,np)
phi[1:na,1:nh,1:num_int,1:np] <- phi0*((1-phi1)/(1+(IC[i,j,k,l]/IC0)^kC) + phi1)

# IB - infection blocking immunity
init_IB[,,,] <- user()
dim(init_IB) <- c(na,nh,num_int,np)
initial(IB[,,,]) <- init_IB[i,j,k,l]
dim(IB) <- c(na,nh,num_int,np)

deriv(IB[1, 1:nh, 1:num_int,1:np]) <- EIR[i,j,k,l]/(EIR[i,j,k,l]* uB + 1) - IB[i,j,k,l]/dB - IB[i,j,k,l]/x_I[i]
deriv(IB[2:na, 1:nh, 1:num_int,1:np]) <- EIR[i,j,k,l]/(EIR[i,j,k,l]* uB + 1) - IB[i,j,k,l]/dB - (IB[i,j,k,l]-IB[i-1,j,k,l])/x_I[i]

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- user() # these parameters characterise the hill function for b
b1 <- user() # prob of infection from bite with zero immunity
kB <- user() #
IB0 <- user()
dim(b) <- c(na,nh,num_int,np)
b[1:na, 1:nh, 1:num_int,1:np] <- b0 * ((1-b1)/(1+(IB[i,j,k,l]/IB0)^kB)+b1)

# detection immunity
init_ID[,,,] <- user()
dim(init_ID) <- c(na,nh,num_int,np)
initial(ID[,,,]) <- init_ID[i,j,k,l]
dim(ID) <- c(na,nh,num_int,np)

deriv(ID[1, 1:nh, 1:num_int,1:np]) <- FOI[i,j,k,l]/(FOI[i,j,k,l]*uD + 1) - ID[i,j,k,l]/dID - ID[i,j,k,l]/x_I[i]
deriv(ID[2:na, 1:nh, 1:num_int,1:np]) <- FOI[i,j,k,l]/(FOI[i,j,k,l]*uD + 1) - ID[i,j,k,l]/dID - (ID[i,j,k,l]-ID[i-1,j,k,l])/x_I[i]

# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- user()
fD0 <- user()
gammaD <- user()
d1 <- user()
ID0 <- user()
kD <- user()
dim(age) <- na
age[] <- user() # vector of age categories supplied by user

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,num_int,np)
p_det[,,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k,l]/ID0)^kD)

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int,np)
FOI_lag[1:na, 1:nh, 1:num_int,1:np] <- EIR[i,j,k,l] * (if(IB[i,j,k,l]==0) b0 else b[i,j,k,l])

# Current FOI depends on humans that have been through the latent period
dE <- user() # latent period of human infection.
dim(FOI) <- c(na,nh,num_int,np)
FOI[,,,] <- delay(FOI_lag[i,j,k,l],dE)


# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
#
dim(mix) <- c(np,np)
mix[,] <- user()                 #  Strength of transmission between populations (patches)
dim(rows) <- c(np,np)
rows[,] <- mix[i,j]*Iv[j]


dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()
dim(EIR) <- c(na,nh,num_int,np)
EIR[,,,] <- av_human[k,l] * rel_foi[j] * foi_age[i] *sum(rows[l,]) /omega

output(Ivout[]) <- Iv[i]
dim(Ivout) <- np

output(omega) <- omega
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

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

dim(Sv) <- np
dim(Ev) <- np
dim(Iv) <- np

dim(init_Sv) <- np
dim(init_Ev) <- np
dim(init_Iv) <- np


# initial state values:
init_Sv[] <- user()
init_Ev[] <- user()
init_Iv[]<- user()
initial(Sv[]) <- init_Sv[i] * mv0[i]
initial(Ev[]) <- init_Ev[i] * mv0[i]
#initial(Ev[1:10]) <- init_Ev/10 * mv0 # Options if not using a delayed delay
#dim(Ev) <- 10
initial(Iv[]) <- init_Iv[i] * mv0[i]

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- user() # infectiousness U -> mosq
cD <- user() # infectiousness D -> mosq
cT <- user() # T -> mosq
gamma1 <- user() # fitted value of gamma1 characterises cA function
dim(cA) <- c(na,nh,num_int,np)
cA[,,,] <- cU + (cD-cU)*p_det[i,j,k,l]^gamma1

# Force of infection from humans to mosquitoes
dim(FOIvijk) <- c(na,nh,num_int,np)
omega <- user() #normalising constant for biting rates
FOIvijk[1:na, 1:nh, 1:num_int,1:np] <- (cT*T[i,j,k,l] + cD*D[i,j,k,l] + cA[i,j,k,l]*A[i,j,k,l] + cU*U[i,j,k,l]) * rel_foi[j] * av_mosq[k,l]*foi_age[i]/omega

dim(lag_FOIv) <- np
lag_FOIv[]=sum(FOIvijk[,,,i])

# Metapopulation force of infection
dim(lag_mix_FOIv) <- np
dim(rows_v) <- c(np,np)
rows_v[,] <- mix[i,j]*lag_FOIv[j]
lag_mix_FOIv[] <- sum(rows_v[i,])

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- user() # Lag from parasites to infectious gametocytes
delayMos <- user() # Extrinsic incubation period.
dim(FOIv) <- np
FOIv[]<- delay(lag_mix_FOIv[i], delayGam)

# Number of mosquitoes that become infected at each time point
dim(surv) <- np
surv[] <- exp(-mu[i]*delayMos)
dim(ince) <- np
ince[] <- FOIv[i] * Sv[i]
dim(lag_incv) <- np
lag_incv[] <- ince[i] * surv[i]
dim(incv) <- np
incv[] <- delay(lag_incv[i], delayMos)
#incv <- lag_incv

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
dim(betaa) <- np
betaa[] <- 0.5*PL[i]/dPL
#betaa <- mv0 * mu0 * theta2

deriv(Sv[1:np]) <- -ince[i] - mu[i]*Sv[i] + betaa[i]
deriv(Ev[1:np]) <- ince[i] - incv[i] - mu[i]*Ev[i]
deriv(Iv[1:np]) <- incv[i] - mu[i]*Iv[i]

# Total mosquito population
dim(mv) <- np
mv[] =Sv[i]+Ev[i] + Iv[i]

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
dim(mv0) <- np
mv0[] <- user() # initial mosquito density
mu0 <- user() # baseline mosquito death rate
tau1 <- user() # duration of host-seeking behaviour
tau2 <- user() # duration of resting behaviour
p10 <- user() # prob of surviving 1 feeding cycle
p2 <- user() #prob of surviving one resting cycle
betaL <- user() # maximum number of eggs per oviposition per mosq

# Entomological variables:
dim(eov) <- np
eov[] <- betaL/mu[i]*(exp(mu[i]/fv[i])-1)
dim(beta_larval) <- np
beta_larval[] <- eov[i]*mu[i]*exp(-mu[i]/fv[i])/(1-exp(-mu[i]/fv[i])) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
dim(lambda) <- np
lambda[] <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval[i]*muLL*dEL/(2*muEL*mu0*dLL*(1+dPL*muPL)))
dim(K0) <- np
K0[] <- 2*mv0[i]*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda[i]+1)/(lambda[i]/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
dim(KL) <- np
KL[] <- K0[i]*theta2
dim(fv) <- np
fv[] <- 1/( tau1/(1-zbar[i]) + tau2 ) # mosquito feeding rate (zbar from intervention param.)
dim(mu) <- np
mu[] <- -fv[i]*log(p1[i]*p2) # mosquito death rate

dim(EL) <- np
dim(LL) <- np
dim(PL) <- np

dim(init_EL) <- np
dim(init_LL) <- np
dim(init_PL) <- np

# finding equilibrium and initial values for EL, LL & PL
init_PL[] <- user()
initial(PL[]) <- init_PL[i]
init_LL[] <- user()
initial(LL[]) <- init_LL[i]
init_EL[] <- user()
initial(EL[]) <- init_EL[i]

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
deriv(EL[1:np]) <- beta_larval[i]*mv[i]-muEL*(1+(EL[i]+LL[i])/KL[i])*EL[i] - EL[i]/dEL
# egg hatching - den. dep. mortality - maturing larvae
deriv(LL[1:np]) <- EL[i]/dEL - muLL*(1+gammaL*(EL[i] + LL[i])/KL[i])*LL[i] - LL[i]/dLL
# pupae - mortality - fully developed pupae
deriv(PL[1:np]) <- LL[i]/dLL - muPL*PL[i] - PL[i]/dPL

##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ##
########################
##------------------------------------------------------------------------------

# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# general parameters
#dim(ITN_IRS_on) <- np
#dim(num_int) <- np
dim(itn_cov) <- np
dim(irs_cov) <- np


ITN_IRS_on <- user() # days after which interventions begin
num_int <- user() # number of intervention categories, ITN only, IRS only, neither, both
itn_cov[] <- user() # proportion of population covered by ITN
irs_cov[] <- user() # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
dim(cov_) <- c(4,np)
cov_[1,] <- (1-itn_cov[j])*(1-irs_cov[j])  # {No intervention}
cov_[2,] <- itn_cov[j]*(1-irs_cov[j]) # 	   {ITN only}
cov_[3,] <- (1-itn_cov[j])*irs_cov[j]	#      {IRS only}
cov_[4,] <- itn_cov[j]*irs_cov[j] #	   {Both ITN and IRS}

dim(cov) <-c(num_int,np)
cov[,] <- cov_[i,j]

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
dim(zhi) <- c(num_int,np)
dim(whi) <- c(num_int,np)
zhi[1:num_int,1:np] <- cov[i,j]*z[i]
whi[1:num_int,1:np] <- cov[i,j]*w[i]
dim(zh) <- np
dim(wh) <- np
zh[1:np] <-if(t < ITN_IRS_on) 0 else sum(zhi[,i])
wh[1:np] <-if(t < ITN_IRS_on) 1 else sum(whi[,i])
# Z (zbar) - average probability of mosquito trying again during single feeding attempt
dim(zbar) <- np
zbar[1:np] <- Q0*zh[i]
# W (wbar) - average probability of mosquito successfully feeding during single attempt
dim(wbar) <- np
wbar[1:np] <- 1 - Q0 + Q0*wh[i]

# p1 is the updated p10 given that interventions are now in place:
dim(p1) <- np
p1[1:np] <- wbar[i]*p10/(1-zbar[i]*p10)
dim(Q) <- np
Q[1:np] <- 1-(1-Q0)/wbar[i] # updated anthropophagy given interventions
dim(av) <- np
av[1:np] <- fv[i]*Q[i] # biting rate on humans

dim(av_mosq) <- c(num_int,np)
av_mosq[1:num_int,1:np] <- av[j]*w[i]/wh[j] # rate at which mosquitoes bite each int. cat.

dim(av_human) <- c(num_int,np)
av_human[1:num_int,1:np] <- av[j]*yy[i]/wh[j] # biting rate on humans in each int. cat.

##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------

# Outputs for each compartment across the sum across all ages, biting heterogeneities and intervention categories
dim(Sout) <- np
dim(Tout) <- np
dim(Dout) <- np
dim(Aout) <- np
dim(Uout) <- np
dim(Pout) <- np


output(Sout[]) <- sum(S[,,,i])
output(Tout[]) <- sum(T[,,,i])
output(Dout[]) <- sum(D[,,,i])
output(Aout[]) <- sum(A[,,,i])
output(Uout[]) <- sum(U[,,,i])
output(Pout[]) <- sum(P[,,,i])

# Outputs for clinical incidence and prevalence on a given day
# population densities for each age category
den[] <- user()
dim(den) <- na
# index of the age vector above 59 months
age59 <- user(integer=TRUE)
# index of the age vector above 5 years
age05 <- user(integer=TRUE)

dim(prev0to59) <- c(age59,nh,num_int,np)
prev0to59[1:age59,,,] <- T[i,j,k,l] + D[i,j,k,l]  + A[i,j,k,l]*p_det[i,j,k,l]
output(prev[]) <- sum(prev0to59[,,,i])/sum(den[1:age59])
dim(prev) <-np

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,num_int,np)
clin_inc0to5[1:age05,,,] <- clin_inc[i,j,k,l]
output(inc05[]) <- sum(clin_inc0to5[,,,i])/sum(den[1:age05])
dim(inc05) <- np
output(inc) <- sum(clin_inc[,,,])

# Param checking outputs
output(mu[]) <- mu
output(beta_larval[]) <- beta_larval
output(KL[]) <- KL
output(mv[]) <- mv
output(Q[]) <- Q
output(wh[]) <- wh
output(d_ITN) <- d_ITN
output(r_ITN) <- r_ITN
output(s_ITN) <- s_ITN
output(d_IRS) <- d_IRS
output(r_IRS) <- r_IRS
output(s_IRS) <- s_IRS
output(cov[,]) <- cov
output(K0[]) <- K0
