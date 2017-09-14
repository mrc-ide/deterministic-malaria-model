
## MODEL VARIABLES

na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
ft <- user() # proportion of cases treated




##################
## HUMAN STATES ##
##################

# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

# fitted parameters for human compartments:
eta <- user(0.0001305) # death rate for exponential population distribution
age_rate[] <- user() # rate at which humans move through age categories
dim(age_rate) <- na
het_wt[] <- user() # weights of the different heterogeneous biting categories
dim(het_wt) <- nh
rA <- user(0.00512821) # rate of movement from A -> U
rT <- user(0.2) # rate of treatment working: T -> P
rD <- user(0.2) #  rate from D -> A
rU <- user(0.00906627) # rate of clearance of subpatent infection U -> S
rP <- user(0.05) # rate at which prophylaxis wears off P -> S

# S - SUSCEPTIBLE
init_S[,,] <- user()
dim(init_S) <- c(na,nh,ni)
initial(S[,,]) <- init_S[i,j,k]
dim(S) <- c(na,nh,ni)

deriv(S[1, 1:nh, 1:ni]) <- -FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] +
   cov[k]*eta*H*het_wt[j] - (eta+age_rate[i])*S[i,j,k]
deriv(S[2:na, 1:nh, 1:ni]) <- -FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] -
    (eta+age_rate[i])*S[i,j,k] + age_rate[i-1]*S[i-1,j,k]


# T- SUCCESSFULLY TREATED
init_T[,,] <- user()
dim(init_T) <- c(na,nh,ni)
initial(T[,,]) <- init_T[i,j,k]
dim(T) <- c(na,nh,ni)

deriv(T[1, 1:nh, 1:ni]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k]
deriv(T[2:na, 1:nh, 1:ni]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k] + age_rate[i-1]*T[i-1,j,k]

# D - CLEAR DISEASE
init_D[,,] <- user()
dim(init_D) <- c(na,nh,ni)
initial(D[,,]) <- init_D[i,j,k]
dim(D) <- c(na,nh,ni)

deriv(D[1, 1:nh, 1:ni]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k]
deriv(D[2:na, 1:nh, 1:ni]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k] + age_rate[i-1]*D[i-1,j,k]

# A - ASYMPTOMATIC DISEASE
init_A[,,] <- user()
dim(init_A) <- c(na,nh,ni)
initial(A[,,]) <- init_A[i,j,k]
dim(A) <- c(na,nh,ni)

deriv(A[1, 1:nh, 1:ni]) <- (1-phi[i,j,k])*FOI[i,j,k]*Y[i,j,k] - FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k]
deriv(A[2:na, 1:nh, 1:ni]) <- (1-phi[i,j,k])*FOI[i,j,k]*Y[i,j,k] - FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k] + age_rate[i-1]*A[i-1,j,k]

# U - SUBPATENT DISEASE
init_U[,,] <- user()
dim(init_U) <- c(na,nh,ni)
initial(U[,,]) <- init_U[i,j,k]
dim(U) <- c(na,nh,ni)

deriv(U[1, 1:nh, 1:ni]) <- rA*A[i,j,k] - FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k]
deriv(U[2:na, 1:nh, 1:ni]) <- rA*A[i,j,k] - FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k] + age_rate[i-1]*U[i-1,j,k]

# P - PROPHYLAXIS
init_P[,,] <- user()
dim(init_P) <- c(na,nh,ni)
initial(P[,,]) <- init_P[i,j,k]
dim(P) <- c(na,nh,ni)

deriv(P[1, 1:nh, 1:ni]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k]
deriv(P[2:na, 1:nh, 1:ni]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k] +
  age_rate[i-1]*P[i-1,j,k]

# The number of individuals able to acquire clinical malaria
dim(Y) <- c(na,nh,ni)
Y[1:na, 1:nh, 1:ni] <- S[i,j,k]+A[i,j,k]+U[i,j,k]

# The number of new cases at this timestep
dim(clin_inc) <- c(na,nh,ni)
clin_inc[1:na, 1:nh, 1:ni] <- phi[i,j,k]*FOI[i,j,k]*Y[i,j,k]

# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph




#####################
## IMMUNITY STATES ##
#####################

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

# fitted immunity parameters:
dm <- user(67.6952) # decay of maternal immunity
uc <- user(6.06349) # scale parameter (see Supplementary mats. 3.1.2)
dc <- user(10950) # decay for clinical immunity
db <- user(3650) # decay for infection blocking immunity
ub <- user(7.19919) # scale param for IB immunity
dd <- user(3650) # decay for detection immunity
ud <- user(9.44512) # scale param for ID immunity
x_I[] <- user() # intermediate variable for calculating immunity functions
dim(x_I) <- na

# ICM - maternally acquired immunity
init_ICM[,,] <- user()
dim(init_ICM) <- c(na,nh,ni)
initial(ICM[,,]) <- init_ICM[i,j,k]
dim(ICM) <- c(na,nh,ni)

deriv(ICM[1, 1:nh, 1:ni]) <- -1/dm*ICM[i,j,k] + (init_ICM[i,j,k]-ICM[i,j,k])/x_I[i]
deriv(ICM[2:na, 1:nh, 1:ni]) <- -1/dm*ICM[i,j,k] - (ICM[i,j,k]-ICM[i,j,k])/x_I[i]

# ICA - exposure driven immunity
init_ICA[,,] <- user()
dim(init_ICA) <- c(na,nh,ni)
initial(ICA[,,]) <- init_ICA[i,j,k]
dim(ICA) <- c(na,nh,ni)

deriv(ICA[1, 1:nh, 1:ni]) <- FOI[i,j,k]/(FOI[i,j,k] * uc + 1) - 1/dc*ICA[i,j,k] -ICA[i,j,k]/x_I[i]
deriv(ICA[2:na, 1:nh, 1:ni]) <- FOI[i,j,k]/(FOI[i,j,k] * uc + 1) - 1/dc*ICA[i,j,k] - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i]

# clinical immunity is a combination of maternal and exposure-driven immunity
dim(IC) <- c(na,nh,ni)
IC[,,] <- ICM[i,j,k] + ICA[i,j,k]

# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- user(0.791666)
phi1 <- user(0.000737) # these parameters characterise the hill function
IC0 <- user(18.02366) # for probability of clinical disease
kc <- user(2.36949) # See supplementary materials 1.1.3
dim(phi) <- c(na,nh,ni)
phi[1:na,1:nh,1:ni] <- phi0*((1-phi1)/(1+(IC[i,j,k]/IC0)^kc) + phi1)

# IB - infection blocking immunity
init_IB[,,] <- user()
dim(init_IB) <- c(na,nh,ni)
initial(IB[,,]) <- init_IB[i,j,k]
dim(IB) <- c(na,nh,ni)

deriv(IB[1, 1:nh, 1:ni]) <- EIR[i,j,k]/(EIR[i,j,k]* ub + 1) - IB[i,j,k]/db - IB[i,j,k]/x_I[i]
deriv(IB[2:na, 1:nh, 1:ni]) <- EIR[i,j,k]/(EIR[i,j,k]* ub + 1) - IB[i,j,k]/db - (IB[i,j,k]-IB[i-1,j,k])/x_I[i]

# b - probability of disease from infectious bite, depends on infection blocking immunity
bh <- user(0.590076) # these parameters characterise the hill function for b
bmin <- user(0.5) # prob of infection from bite with zero immunity
kb <- user(2.15506) #
IB0 <- user(43.8787)
dim(b) <- c(na,nh,ni)
b[1:na, 1:nh, 1:ni] <- bh * ((1-bmin)/(1+(IB[i,j,k]/IB0)^kb)+bmin)

# detection immunity
init_ID[,,] <- user()
dim(init_ID) <- c(na,nh,ni)
initial(ID[,,]) <- init_ID[i,j,k]
dim(ID) <- c(na,nh,ni)

deriv(ID[1, 1:nh, 1:ni]) <- FOI[i,j,k]/(FOI[i,j,k]*ud + 1) - ID[i,j,k]/dd - ID[i,j,k]/x_I[i]
deriv(ID[2:na, 1:nh, 1:ni]) <- FOI[i,j,k]/(FOI[i,j,k]*ud + 1) - ID[i,j,k]/dd - (ID[i,j,k]-ID[i-1,j,k])/x_I[i]

# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
ad0 <- user(8001.99) # no idea where these are from
fd0 <- user(0.007055) # or who fit them
gammad <- user(4.8183) #
dmin <- user(0.160527)
ID0 <- user(1.577533)
kd <- user(0.476614)
dim(age) <- na
age[] <- user() # vector of age categories supplied by user

dim(fd) <- na
fd[1:na] <- 1-(1-fd0)/(1+(age[i]/ad0)^gammad)
dim(p_det) <- c(na,nh,ni)
p_det[,,] <- dmin + (1-dmin)/(1 + fd[i]*(ID[i,j,k]/ID0)^kd)

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,ni)
FOI_lag[1:na, 1:nh, 1:ni] <- EIR[i,j,k] * (if(IB[i,j,k]==0) bh else b[i,j,k])

# Current FOI depends on humans that have been through the latent period and are
# producing gametocytes
dE <- user(12) # length of time from infection to gametocytogenesis
dim(FOI) <- c(na,nh,ni)
FOI[,,] <- delay(FOI_lag[i,j,k],dE)


# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()
dim(EIR) <- c(na,nh,ni)
EIR[,,] <- av_human[k] * rel_foi[j] * foi_age[i]/omega*Iv





#####################
## MOSQUITO STATES ##
#####################

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# initial state values:
init_Sv <- user()
init_Ev <- user()
init_Iv <- user()
initial(Sv) <- init_Sv * mv0
initial(Ev[1:10]) <- init_Ev/10 * mv0
dim(Ev) <- 10
initial(Iv) <- init_Iv *mv0

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- user(0.006203) # infectiousness U -> mosq
cD <- user(0.0676909) # infectiousness D -> mosq
cT <- user(0.02179647) # T -> mosq
g_inf <- user(1.82425) # fitted value of g_inf characterises cA function
dim(cA) <- c(na,nh,ni)
cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^g_inf


# Force of infection from humans to mosquitoes
dim(FOIvijk) <- c(na,nh,ni)
omega <- user() #normalising constant for biting rates
FOIvijk[1:na, 1:nh, 1:ni] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq[k]*foi_age[i]/omega
lag_FOIv=sum(FOIvijk)
# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
latgam <- user(12.5) # latent period in humans
FOIv <- delay(lag_FOIv, latgam)


# Number of mosquitoes that become infected at each time point
ince <- FOIv * Sv

# Number of mosquitoes born (depends on PL, number of larvae)
betaa <- 0.5*PL/dp
deriv(Sv) <- -ince - mu*Sv + betaa
# Exposed mosquitoes go through 10 compartments to simulate 10 day development of
# sporozoites (and mosquitoes that die off during this wait)
deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
deriv(Iv) <- Ev[10] - mu*Iv
# Total mosquito population
mv = Sv+sum(Ev)+Iv




##########################
## SEASONALITY FUNCTION ##
##########################

# Seasonality is added into the model using a Fourier series that was fit to rainfall at every admin 2 level

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
theta2 <- max((ssa0+ssa1*cos(2*pi*t/365)+ssa2*cos(2*2*pi*t/365)+ssa3*cos(3*2*pi*t/365)+ssb1*sin(2*pi*t/365)+ssb2*sin(2*2*pi*t/365)+ ssb3*sin(3*2*pi*t/365) ) /theta_c,0.001)




###################
## LARVAL STATES ##
###################

# Model by White et al.
# (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)

# EL - early larval instar stage
# LL - late larval instar stage
# PL - pupal stage

# mean carrying capacity from initial mosquito density:
dl <- user(3.72) # development time of larvae
dp <- user(0.643) #development time of pupae
de <- user(6.64) #development time of early stage
mul <- user(0.0348) #daily density dep. mortality rate of larvae
mup <- user(0.249) #daily den. dep. mortality rate of pupae
mue <- user(0.0338) #daily den. dep. mortality rate of early stage
gammal <- user(13.25) # eff. of den. dep. on late stage relative to early stage
K0 <- user() # baseline carrying capacity


# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
KL <- K0*theta2


# fitted entomological parameters:
mu0 <- user(0.132) # daily hazard of death from external causes
mv0 <- user() # initial mosquito density
fv0 <- user(0.33333) # initial mosquito feeding rate
tau1 <- user(0.69) # duration of host-seeking behaviour
tau2 <- user(2.00003) # duration of resting behaviour
p10 <- user(0.9129447) # prob of surviving 1 feeding cycle
p2 <- user(0.7679705) #prob of surviving one resting cycle
beta_larval0 <- user(21.2) # initial number of eggs
eov <- user(78.0345) # maximum number of eggs per oviposition per mosq


# Entomological variables:
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
fv <- 1/( tau1/(1-zbar) + tau2 ) # mosquito feeding rate (zbar from intervention param.)
mu <- -fv*log(p1*p2) # mosquito death rate


# finding equilibrium and initial values for EL, LL & PL
initial(PL) <- 2*dp*mu0*mv0
PL0 <- 2*dp*mu0*mv0
initial(LL) <- dl*(mup+1/dp)*PL0
LL0 <- dl*(mup+1/dp)*PL0
initial(EL) <- (LL0/dl + mul*LL0*(1+gammal*LL0/K0))/(1/de-mul*gammal*LL0/K0)

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
deriv(EL) <- beta_larval*mv-mue*(1+(EL+LL)/KL)*EL - EL/de
# egg hatching - den. dep. mortality - maturing larvae
deriv(LL) <- EL/de - mul*(1+gammal*(EL + LL)/KL)*LL - LL/dl
# pupae - mortality - fully developed pupae
deriv(PL) <- LL/dl - mup*PL - PL/dp



########################
## INTERVENTION MODEL ##
########################

# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# general parameters
ITN_IRS_on <- user(365) # days after which interventions begin
ni <- user(4) # number of intervention categorys, ITN only, IRS only, neither, both
itn_cov <- user(0.8) # proportion of population covered by ITN
irs_cov <- user(0.4) # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
dim(cov) <- ni
cov[1] <- (1-itn_cov)*(1-irs_cov)  # {No intervention}
cov[2] <- itn_cov*(1-irs_cov) # 	   {ITN only}
cov[3] <- (1-itn_cov)*irs_cov	#      {IRS only}
cov[4] <- itn_cov*irs_cov #	   {Both ITN and IRS}


IRS_interval <- user(365) # how long IRS lasts
ITN_interval <- user(1095) # how long ITN lasts
chi <- user(0.86) # proportion of vector endophily
Q0 <- user(0.92) # proportion of anthropophagy
PHI_B <- user(0.89) # endophagy in bed
PHI_I <- user(0.97) # endophagy indoors

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN/IRS
# d - probability of dying after hitting ITN/IRS
# s - probability of successful feed after hitting ITN/IRS

# The maximum (and then minimum) r and d values for ITN/IRS on day 0 before they decay
r_ITN0 <- user(0.56)
d_ITN0 <- user(0.41)
r_IRS0 <- user(0.6)
d_IRS0 <- user(1)
r_ITN_min <- user(0.24)
irs_loss <- user(0.003798067)
itn_loss <- user(0.0007193308)

# Calculates decay for ITN/IRS
ITN_decay = if(t < ITN_IRS_on) 0 else exp(-((t-ITN_IRS_on)%%ITN_interval) * itn_loss)
IRS_decay = if(t < ITN_IRS_on) 0 else exp(-((t-ITN_IRS_on)%%IRS_interval) * irs_loss)

# The r,d and s values turn on after ITN_IRS_on and decay accordingly
d_ITN <- if(t < ITN_IRS_on) 0 else d_ITN0*ITN_decay
r_ITN <- if(t < ITN_IRS_on) 0 else r_ITN_min + (r_ITN0 - r_ITN_min)*ITN_decay
s_ITN <- if(t < ITN_IRS_on) 1 else 1 - d_ITN - r_ITN

r_IRS <- if(t < ITN_IRS_on) 0 else d_IRS0*IRS_decay
d_IRS <- if(t < ITN_IRS_on) 0 else chi*d_IRS0*IRS_decay
s_IRS <- if(t < ITN_IRS_on) 1 else 1 - d_IRS


# probability that mosquito bites and survives for each intervention category
dim(w) <- ni
w[1] <- 1
w[2] <- 1 - PHI_B + PHI_B*s_ITN
w[3] <- 1 - PHI_I + PHI_I*(1-r_IRS)*s_IRS
w[4] <- 1 - PHI_I + PHI_B*(1-r_IRS)*s_ITN*s_IRS + (PHI_I - PHI_B)*(1-r_IRS)*s_IRS
# probability that mosq feeds during a single attempt for each int. cat.
dim(yy) <- ni
yy[1] <- 1
yy[2] <- w[2]
yy[3] <- 1 - PHI_I + PHI_I*(1-r_IRS)
yy[4] <- 1 - PHI_I + PHI_B*(1-r_IRS)*s_ITN + (PHI_I - PHI_B)*(1-r_IRS)
# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z) <- ni
z[1] <- 0
z[2] <- PHI_B*r_ITN
z[3] <- PHI_I*r_IRS
z[4] <- PHI_B*(r_IRS+ (1-r_IRS)*r_ITN) + (PHI_I - PHI_B)*r_IRS

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- ni
dim(whi) <- ni
zhi[1:ni] <- cov[i]*z[i]
whi[1:ni] <- cov[i]*w[i]
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
dim(av_mosq) <- ni
av_mosq[1:ni] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.
dim(av_human) <- ni
av_human[1:ni] <- av*yy[i]/wh # biting rate on humans in each int. cat.


###################
## MODEL OUTPUTS ##
###################

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

dim(prev0to59) <- c(age59,nh,ni)
prev0to59[1:age59,,] <- T[i,j,k] + D[i,j,k] + A[i,j,k] * p_det[i,j,k]
output(prev) <- sum(prev0to59[,,])/sum(den[1:age59])

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,ni)
clin_inc0to5[1:age05,,] <- clin_inc[i,j,k]
output(inc) <- sum(clin_inc0to5)/sum(den[1:age05])
