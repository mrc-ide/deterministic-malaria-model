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
deriv(ICM[2:na, 1:nh, 1:num_int]) <- -1/dCM*ICM[i,j,k] - (ICM[i,j,k]-ICM[i,j,k])/x_I[i]

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
#betaa <- 0.5*PL/dPL
betaa <- mv0 * mu * theta2

deriv(Sv) <- -ince - mu*Sv + betaa
#deriv(Ev) <- ince - incv - mu*Ev
deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
deriv(Iv) <- incv - mu*Iv

# Total mosquito population
#mv = Sv+Ev+Iv
mv = Sv+sum(Ev)+Iv


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
tau1 <- user() # duration of host-seeking behaviour
tau2 <- user() # duration of resting behaviour
p10 <- user() # prob of surviving 1 feeding cycle
p2 <- user() #prob of surviving one resting cycle
betaL <- user() # maximum number of eggs per oviposition per mosq

# Entomological variables:
eov <- betaL/mu*(exp(mu/fv)-1)
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mu*dLL*(1+dPL*muPL)))
K0 <- 2*mv0*dLL*mu*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)

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

# general parameters
ITN_on <- user() # day when ITNs are first introduced
EM_on <- user() # day when emanators are first introduced
IRS_on <- user() # day when IRS begins
num_int <- user() # number of intervention categorys, ITN only, emanator only, neither, both
itn_cov <- user() # proportion of population covered by ITN
em_cov <- user() # proportion of population covered by emanator
irs_cov <- user() # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
dim(cov) <- num_int
cov[1] <- (1-irscov)*(1-itncov)*(1-emcov) # no intervention
cov[2] <- itncov*(1-irscov)*(1-emcov) # itn only
cov[3] <- irscov*(1-itncov)*(1-emcov) # irs only
cov[4] <- emcov*(1-itncov)*(1-irscov) # em only
cov[5] <- itncov*irscov*(1-emcov) # irs and itn
cov[6] <- irscov*(1-itncov)*emcov # irs and em
cov[7] <- itncov*(1-irscov)*emcov # em and itn
cov[8] <- itncov*irscov*emcov # all 3


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
mort_assay <- 1 - surv_bioassay

# Relationship between mortality in bioassay to hut trial, logit scale
mort_hut_a <- 0.6338 + 3.9970 * (mort_assay-0.5)
mort_hut <- exp(mort_hut_a)/(1+exp(mort_hut_a))

# Relationship between hut trial mortality and deterrence
det_hut_a <- 0.07117+1.257*(mort_hut-0.5)-1.517*(mort_hut-0.5)^2
det_hut <- if(det_hut_a < 0) 0 else det_hut_a # censored to stop becoming negative
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

# IRS parameters
sprayChemical <- user()
# 1 = pyrethroid and mort_assay has an impact on efficacy, else
# 2 = Actellic
# 3 = Bendiocarb
# 4 = Sumishield IRS

# Adjusting pyrethroids for resistance
int1 <-  -1.416835
int2 <-  -0.009357982
int3 <- -2.076946
int4 <-  0.002047915
int5 <-  -2.960124
int6 <-  0.0003173997

grad1 <- 2.972465
grad2 <- 0.001040958
grad3 <- 2.352004
grad4 <- 0.009000071
grad5 <- 4.771978
grad6 <- -0.0077501779

alpha1 <-  -1.942755
alpha2 <-  0.03361892
beta1 <-   1.849397
beta2 <-   -0.04921294
thet1 <-  -2.071873
thet2 <-  0.02004906

dim(rhovec2) <- 4
dim(rhovec1) <- 4

rhovec2[1] <- grad1 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * mort_assay))) + int1
rhovec1[1] <- grad2 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * mort_assay))) + int2

rhovec2[2] <-  -0.011
rhovec1[2] <- 2.324 # Actellic

rhovec2[3] <-  -0.036
rhovec1[3] <- 1.491 # Bendiocarb

rhovec2[4] <- -0.009
rhovec1[4] <- 0.951 # Sumishield

rho_IRS2 <- rhovec2[sprayChemical]
rho_IRS1 <- rhovec1[sprayChemical]

dim(tauvec1) <- 4
dim(tauvec2) <- 4

tauvec2[1] <- grad4 * (1 / (1 + exp(-beta1 - beta2 * 100 * mort_assay))) + int4
tauvec1[1] <- grad3 * (1 / (1 + exp(-beta1 - beta2 * 100 * mort_assay))) + int3

tauvec2[2] <- 0.009
tauvec1[2] <- -2.686 # Actellic

tauvec2[3] <- 0.007
tauvec1[3] <- -0.760 # Bendiocarb

tauvec2[4] <- 0.007
tauvec1[4] <- -1.672 # Sumishield

tau_IRS2 <- tauvec2[sprayChemical]
tau_IRS1 <- tauvec1[sprayChemical]

dim(boovec1) <- 4
dim(boovec2) <- 4

boovec2[1] <- grad6 * (1 / (1 + exp(-thet1 - thet2 * 100 * mort_assay))) + int6
boovec1[1] <- grad5 * (1 / (1 + exp(-thet1 - thet2 * 100 *  mort_assay))) + int5

boovec2[2] <- -0.002
boovec1[2] <- -1.266 # Actellic

boovec2[3] <- 0.004
boovec1[3] <- -0.765 # Bendiocarb

boovec2[4] <- -0.008
boovec1[4] <- -0.280 # Sumishield

boo_IRS1 <- boovec1[sprayChemical]
boo_IRS2 <- boovec2[sprayChemical]

mort_hut_IRS <- 1/(1 + exp(-rho_IRS1 - rho_IRS2 * mod(t-IRS_on, IRS_interval)))
suc_hut_IRS <- 1/(1 + exp(-tau_IRS1 - tau_IRS2 * mod(t-IRS_on, IRS_interval)))
det_hut_IRS <- 1/(1 + exp(-boo_IRS1 - boo_IRS2 * mod(t-ITN_IRS_on, IRS_interval)))

rep_hut_IRS   <- 1 - suc_hut_IRS - mort_hut_IRS

kp1_IRS  <- (1 - det_hut_IRS)*suc_hut_IRS
jp1_IRS  <- (1 - det_hut_IRS)*rep_hut_IRS+det_hut_IRS
lp1_IRS  <- (1 - det_hut_IRS)*mort_hut_IRS

r_IRS <- if(t < IRS_on) 0 else (1-kp1_IRS/0.699)*(jp1_IRS/(lp1_IRS+jp1_IRS))	#		{probability of repeating with an encounter with ITN (max)}			; cycle repeating rate
d_IRS = if(t < IRS_on) 0 else (1-kp1_IRS/0.699)*(lp1_IRS/(lp1_IRS+jp1_IRS))	#		{probability of dying behaviour (max)}		; insecticide mortality rate
s_IRS = if(t< IRS_on) 1 else kp1_IRS/0.699 # 		; successful protected human biting

checker_irs = rep_hut_IRS + mort_hut_IRS + suc_hut_IRS

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


# The r,d and s values turn on after ITN_EM_on and decay accordingly
d_ITN <- if(t < ITN_on) 0 else d_ITN0*ITN_decay
# r_ITN <- if(t < ITN_on) 0 else r_ITN1 + (r_ITN0 - r_ITN1)*ITN_decay
r_ITN <- if(t < ITN_on) 0 else r_ITN_min + (r_ITN0- r_ITN_min)*ITN_decay
# Ellie's edit
s_ITN <- if(t < ITN_on) 1 else 1 - d_ITN - r_ITN

dim(r_EM) <- d_len
r_EM[1:d_len] <- if(t < EM_on) 0 else em_prod_profile[i]*EM_decay

# proportion of bites around the emanator that are successful
d_len <- user()
dim(p_EM_vec) <- d_len
p_EM_vec[1:d_len] <- em_human_dist[i]*(1-r_EM[i])
p_EM <- sum(p_EM_vec)
r_EM_in0 <- user()
em_in <- user()
r_EM_in <- if(t < EM_on) 0 else r_EM_in0*EM_decay

# new values since emanators don't kill, just repel
d_EM <- 0
s_EM <- 1 - d_EM

cov[1] <- (1-irscov)*(1-itncov)*(1-emcov) # no intervention
cov[2] <- itncov*(1-irscov)*(1-emcov) # itn only
cov[3] <- irscov*(1-itncov)*(1-emcov) # irs only
cov[4] <- emcov*(1-itncov)*(1-irscov) # em only
cov[5] <- itncov*irscov*(1-emcov) # irs and itn
cov[6] <- irscov*(1-itncov)*emcov # irs and em
cov[7] <- itncov*(1-irscov)*emcov # em and itn
cov[8] <- itncov*irscov*emcov # all 3

# probability that mosquito bites and survives for each intervention category
dim(w) <- num_int
w[1] <- 1 # no intervention
w[2] <- 1 - bites_Bed + bites_Bed*s_ITN # ITN only
w[3] <- 1 - bites_Indoors + bites_Indoors*s_IRS # IRS only
w[4] <- if(em_in == 0) 1 - bites_Emanator + p_EM*bites_Emanator else # EM only outdoors
  1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + (1-r_EM_in)*bites_Indoors*s_EM # EM indoors
w[5] <- 1 - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN*s_IRS + (bites_Indoors-bites_Bed)*s_IRS # ITN and IRS
w[6] <- if(em_in == 0) 1- bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Indoors*s_IRS else # IRS and EM outdoors only
  1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Indoors*(1-r_EM_in)*s_IRS*s_EM # IRS and EM indoors
w[7] <- if(em_in == 0) 1 - bites_Bed - bites_Emanator + bites_Bed*s_ITN + s_EM*bites_Emanator else # EM only outdoors and ITN
  1 - bites_Indoors + bites_Bed*(1-r_EM_in)*s_ITN + (bites_Indoors-bites_Bed)*(1-r_EM_in) - bites_Emanator + p_EM*bites_Emanator # EM indoors and ITN
w[8] <- if(em_in == 0) 1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN*s_IRS + (bites_Indoors-bites_Bed)*s_IRS else # All 3, em outside only
  1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Bed*(1-r_EM_in)*(1-r_IRS)*s_ITN*s_IRS*s_EM + (bites_Indoors-bites_Bed)*(1-r_EM_in)*s_IRS*s_EM # All 3, em inside

# probability that mosq feeds during a single attempt for each int. cat.
dim(yy) <- num_int
yy[1] <- 1 # no intervention
yy[2] <- w[2] # ITN only
yy[3] <- 1 - bites_Indoors + bites_Indoors*(1-r_IRS) # IRS only
yy[4] <- if(em_in == 0) w[4] else # EM only outdoors
  1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + (1-r_EM_in)*bites_Indoors # EM indoors
yy[5] <- 1 - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN + (bites_Indoors-bites_Bed)*(1-r_IRS) # ITN and IRS
yy[6] <- if(em_in == 0) 1- bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Indoors*(1-r_IRS) else # IRS and EM outdoors only
  1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Indoors*(1-r_EM_in)*(1-r_IRS) # IRS and EM indoors
yy[7] <- if(em_in == 0) 1 - bites_Bed - bites_Emanator + bites_Bed*s_ITN + p_EM*bites_Emanator else # EM only outdoors and ITN
  1 - bites_Indoors + bites_Bed*(1-r_EM_in)*s_ITN + (bites_Indoors-bites_Bed)*(1-r_EM_in) - bites_Emanator + p_EM*bites_Emanator # EM indoors and ITN
yy[8] <-  if(em_in == 0) 1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN + (bites_Indoors-bites_Bed)*(1-r_IRS) else # All 3, em outside only
  1 - bites_Emanator + p_EM*bites_Emanator - bites_Indoors + bites_Bed*(1-r_EM_in)*(1-r_IRS)*s_ITN + (bites_Indoors-bites_Bed)*(1-r_EM_in)*(1-r_IRS) # All 3, em inside


# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z) <- num_int
z[1] <- 0 # no intervention
z[2] <- bites_Bed*r_ITN # ITN only
z[3] <- bites_Indoors*r_IRS # IRS only
z[4] <- if(em_in == 0) (1-p_EM)*bites_Emanator else # EM only outdoors
  (1-p_EM)*bites_Emanator + bites_Indoors*r_EM_in # EM only indoors
z[5] <- bites_Bed*(r_IRS+ (1-r_IRS)*r_ITN) + (bites_Indoors - bites_Bed)*r_IRS # ITN and IRS
z[6] <- if(em_in == 0) (1-p_EM)*bites_Emanator + bites_Indoors*r_IRS else # IRS and EM outdoors
  (1-p_EM)*bites_Emanator + bites_Indoors*(r_EM_in+r_IRS*(1-r_EM_in)) # IRS and EM indoors
z[7] <- if(em_in == 0) (1-p_EM)*bites_Emanator + bites_Bed*r_ITN else # ITN and EM outdoors
  (1-p_EM)*bites_Emanator + bites_Bed*(r_EM_in + (1-r_EM_in)*r_ITN) + (bites_Indoors-bites_Bed)*r_EM_in # ITN and EM indoors
z[8] <- if(em_in == 0) (1-p_EM)*bites_Emanator + bites_Bed*(r_IRS+ (1-r_IRS)*r_ITN) + (bites_Indoors - bites_Bed)*r_IRS else # All 3, em outdoors
  (1-p_EM)*bites_Emanator + bites_Bed*(r_EM_in + (1-r_EM_in)*r_IRS + (1-r_EM_in)*(1-r_IRS)*r_ITN) + bites_Indoors*(r_EM_in + (1-r_EM_in)*r_IRS) # All 3, em indoors



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
av_mosq[1:num_int] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.
# note: this is no longer altereted in IRS model
# ALTERED DUE TO PAGE 6 SUPP MAT 2, the biting rate was previously inflated to account for the fact that some mosquitoes would bite due and then die due to IRS
# This essentially meant that as the biting rate on humans covered by emanators dropped, this increased the biting rate on people with no interventions
# av_human[1:num_int] <- av*yy[i]/wh # biting rate on humans in each int. cat.
#av_mosq[1:num_int] <- av*w[i]
dim(av_human) <- num_int
av_human[1:num_int] <- av*yy[i]

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

dim(prev0to59) <- c(age59,nh,num_int)
prev0to59[1:age59,,] <- T[i,j,k] + D[i,j,k] + A[i,j,k] * p_det[i,j,k]
output(prev) <- sum(prev0to59[,,])/sum(den[1:age59])

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,num_int)
clin_inc0to5[1:age05,,] <- clin_inc[i,j,k]
output(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])
output(inc) <- sum(clin_inc[,,])
dim(zz) <- num_int

output(p_EM) <- p_EM
output(Yout) <- sum(Y[,,])
output(phiout) <- sum(phi[,,])
output(FOIout) <- sum(FOI[,,])
output(EIRout) <- sum(EIR[,,])
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
output(r_em_in) <- r_EM_in0
