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
output(avhc) <- avhc

output(IVRM_sr) <- IVRM_sr
