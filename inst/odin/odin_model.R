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
age20l <- user(integer=TRUE) # lower index of age 20 age compartment
age20u <- user(integer=TRUE) # upper index of age 20 age compartment
age_20_factor <- user() # factor calculated in equilibrium solution
PM <- user() # immunity constant

# ICM - maternally acquired immunity
init_ICM[,,] <- user()
dim(init_ICM) <- c(na,nh,num_int)
initial(ICM[,,]) <- init_ICM[i,j,k]
dim(ICM) <- c(na,nh,num_int)
dim(init_ICM_pre) <- c(nh,num_int)
init_ICM_pre[1:nh,1:num_int] <- PM*(ICA[age20l,i,j] + age_20_factor*(ICA[age20u,i,j]-ICA[age20l,i,j]))

deriv(ICM[1, 1:nh, 1:num_int]) <- -1/dCM*ICM[i,j,k] + (init_ICM_pre[j,k]-ICM[i,j,k])/x_I[i]
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
aD <- user()
fD0 <- user()
gammaD <- user()
d1 <- user()
ID0 <- user()
kD <- user()
dim(age) <- na
age[] <- user() # vector of age categories supplied by user

scaling_cc_t[] <- user()
dim(scaling_cc_t) <- user()
tt[] <- user()
dim(tt) <- user()
scaling_cc <- interpolate(tt, scaling_cc_t, "linear")
#scaling_cc <-

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,num_int)
p_det[,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k]/ID0)^kD)

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int)
FOI_lag[1:na, 1:nh, 1:num_int] <- EIR[i,j,k] * (if(IB[i,j,k]==0) b0 else b[i,j,k])

# Current FOI depends on humans that have been through the latent period
dE <- user() # latent period of human infection.
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
EIR[,,] <- av_human[k] * rel_foi[j] * foi_age[i] * Iv/omega
output(Ivout) <- Iv

output(omega) <- omega
##------------------------------------------------------------------------------
##########################
## SEASONALITY FUNCTION ## NOW REMOVED AND CREATED WITHIN 'equilibrium-init-create.R'
##########################
##------------------------------------------------------------------------------

# Seasonality is added into the model using a Fourier series that was fit to rainfall at every admin 1 level
#pi <- user() # weird quirk, need to pass pi


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
initial(Ev) <- init_Ev * mv0
#initial(Ev[1:10]) <- init_Ev/10 * mv0 # Options if not using a delayed delay
#dim(Ev) <- 10
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
delayGam <- user() # Lag from parasites to infectious gametocytes
delayMos <- user() # Extrinsic incubation period.
FOIv <- delay(lag_FOIv, delayGam)

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
eov <- betaL/mu*(exp(mu/fv)-1)
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mu0*dLL*(1+dPL*muPL)))
K0 <- 2*mv0*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
KL <- K0*scaling_cc #*scaling_factor for each day to change average carrying capacity and input seasonality.

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
age59 <- user(integer=TRUE)
# index of the age vector above 5 years
age05 <- user(integer=TRUE)

dim(prev0to59) <- c(age59,nh,num_int)
prev0to59[1:age59,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]
output(prev) <- sum(prev0to59[,,])/sum(den[1:age59])

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,num_int)
clin_inc0to5[1:age05,,] <- clin_inc[i,j,k]
output(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])
output(inc) <- sum(clin_inc[,,])

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
