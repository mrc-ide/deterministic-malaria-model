
## MODEL VARIABLES

na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
ft <- user() # proportion of cases treated
ncc <- user() # number of coverage correlation categories (2)

###########################
# MDA set up
###############################

MDA_cov = user()
MDA_grp_prop = user()




ttt[] <- user()
dim(ttt)<-user()



MDA_t1 = user()
MDA_t2 = user()
MDA_t3 = user()
MDA_t4 = user()
MDA_t5 = user()
MDA_t6 = user()
MDA_t7 = user()
MDA_t8 = user()

MDA_sr =       if(t <= MDA_t1 +1) MDA_t1 else if (t > MDA_t1 && t<= MDA_t2+1) MDA_t2 else if (t > MDA_t2 && t<= MDA_t3+1) MDA_t3 else if (t > MDA_t3 && t<= MDA_t4+1) MDA_t4 else if (t > MDA_t4 && t<= MDA_t5+1) MDA_t5 else if (t > MDA_t5 && t<= MDA_t6+1) MDA_t6 else if (t > MDA_t6 && t<= MDA_t7+1) MDA_t7 else MDA_t8

#MDA_sr = 500


dim(MDA_A) <- c(na,nh,num_int,ncc)
dim(MDA_T) <- c(na,nh,num_int,ncc)
dim(MDA_D) <- c(na,nh,num_int,ncc)
dim(MDA_U) <- c(na,nh,num_int,ncc)
dim(MDA_S) <- c(na,nh,num_int,ncc)

# MDA_A[1:na, 1:nh, 1:num_int] = if(t == MDA_sr ) log(1/(1-MDA_cov))*A[i,j,k] else 0
# MDA_T[1:na, 1:nh, 1:num_int] = if(t == MDA_sr ) log(1/(1-MDA_cov))*T[i,j,k] else 0
# MDA_D[1:na, 1:nh, 1:num_int] = if(t == MDA_sr ) log(1/(1-MDA_cov))*D[i,j,k] else 0
# MDA_U[1:na, 1:nh, 1:num_int] = if(t == MDA_sr ) log(1/(1-MDA_cov))*U[i,j,k] else 0
# MDA_S[1:na, 1:nh, 1:num_int] = if(t == MDA_sr ) log(1/(1-MDA_cov))*S[i,j,k] else 0


MDA_succ = 0.95
MDA_eff1 = 1*MDA_succ
MDA_eff2 = MDA_succ * MDA_cov*(1-MDA_grp_prop)/(1-MDA_grp_prop*MDA_cov)

SMCyrs = user()

MDA_A[1:na, 1:nh, 1:num_int, 1:ncc] =  0
MDA_T[1:na, 1:nh, 1:num_int, 1:ncc] =  0
MDA_D[1:na, 1:nh, 1:num_int, 1:ncc] =  0
MDA_U[1:na, 1:nh, 1:num_int, 1:ncc] =  0
MDA_S[1:na, 1:nh, 1:num_int, 1:ncc] =  0

MDA_A[2:SMCyrs, 1:nh, 1:num_int, 1] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff1))*A[i,j,k,l] else 0
MDA_T[2:SMCyrs, 1:nh, 1:num_int, 1] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff1))*T[i,j,k,l] else 0
MDA_D[2:SMCyrs, 1:nh, 1:num_int, 1] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff1))*D[i,j,k,l] else 0
MDA_U[2:SMCyrs, 1:nh, 1:num_int, 1] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff1))*U[i,j,k,l] else 0
MDA_S[2:SMCyrs, 1:nh, 1:num_int, 1] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff1))*S[i,j,k,l] else 0

MDA_A[2:SMCyrs, 1:nh, 1:num_int, 2] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff2))*A[i,j,k,l] else 0
MDA_T[2:SMCyrs, 1:nh, 1:num_int, 2] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff2))*T[i,j,k,l] else 0
MDA_D[2:SMCyrs, 1:nh, 1:num_int, 2] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff2))*D[i,j,k,l] else 0
MDA_U[2:SMCyrs, 1:nh, 1:num_int, 2] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff2))*U[i,j,k,l] else 0
MDA_S[2:SMCyrs, 1:nh, 1:num_int, 2] = if (t >= MDA_sr && t < (MDA_sr + 1)) log(1/(1-MDA_eff2))*S[i,j,k,l] else 0


#
# deriv(Sv) =               if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len  )) betaa - mu*Sv - avB2*Sv*ivm_cov - FOIv*Sv*(1-ivm_cov)                      else betaa - FOIv*Sv - mu*Sv


##################
## HUMAN STATES ##
##################

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
dim(init_S) <- c(na,nh,num_int,ncc)
initial(S[,,,]) <- init_S[i,j,k,l]
dim(S) <- c(na,nh,num_int,ncc)

deriv(S[1, 1:nh, 1:num_int, 1:ncc]) <- -FOI[i,j,k,l]*S[i,j,k,l] + rP*P[i,j,k,l] + rU*U[i,j,k,l] +
   cov[k]*eta*H*het_wt[j]*ccov[l] - (eta+age_rate[i])*S[i,j,k,l]   - MDA_S[i,j,k,l]

# deriv(S[2:na, 1:nh, 1:num_int]) <- -FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] -
#     (eta+age_rate[i])*S[i,j,k] + age_rate[i-1]*S[i-1,j,k]  - MDA_S[i,j,k]

deriv(S[2:na, 1:nh, 1:num_int, 1:ncc]) <-  -FOI[i,j,k,l]*S[i,j,k,l] + rP*P[i,j,k,l] + rU*U[i,j,k,l] -
  (eta+age_rate[i])*S[i,j,k,l] + age_rate[i-1]*S[i-1,j,k,l]  - MDA_S[i,j,k,l]


# T- SUCCESSFULLY TREATED
init_T[,,,] <- user()
dim(init_T) <- c(na,nh,num_int, ncc)
initial(T[,,,]) <- init_T[i,j,k,l]
dim(T) <- c(na,nh,num_int,ncc)

deriv(T[1, 1:nh, 1:num_int, 1:ncc]) <- ft*clin_inc[i,j,k,l] - rT*T[i,j,k,l] -
  (eta+age_rate[i])*T[i,j,k,l]  - MDA_T[i,j,k,l]

deriv(T[2:na, 1:nh, 1:num_int, 1:ncc]) <- ft*clin_inc[i,j,k,l] - rT*T[i,j,k,l] -
  (eta+age_rate[i])*T[i,j,k,l] + age_rate[i-1]*T[i-1,j,k,l] - MDA_T[i,j,k,l]

# D - CLEAR DISEASE
init_D[,,,] <- user()
dim(init_D) <- c(na,nh,num_int, ncc)
initial(D[,,,]) <- init_D[i,j,k,l]
dim(D) <- c(na,nh,num_int,ncc)

deriv(D[1, 1:nh, 1:num_int, 1:ncc]) <- (1-ft)*clin_inc[i,j,k,l] - rD*D[i,j,k,l] -
  (eta+age_rate[i])*D[i,j,k,l] - MDA_D[i,j,k,l]
deriv(D[2:na, 1:nh, 1:num_int, 1:ncc]) <- (1-ft)*clin_inc[i,j,k,l] - rD*D[i,j,k,l] -
  (eta+age_rate[i])*D[i,j,k,l] + age_rate[i-1]*D[i-1,j,k,l] - MDA_D[i,j,k,l]

# A - ASYMPTOMATIC DISEASE
init_A[,,,] <- user()
dim(init_A) <- c(na,nh,num_int,ncc)
initial(A[,,,]) <- init_A[i,j,k,l]
dim(A) <- c(na,nh,num_int,ncc)

deriv(A[1, 1:nh, 1:num_int, 1:ncc]) <- (1-phi[i,j,k,l])*FOI[i,j,k,l]*Y[i,j,k,l] - FOI[i,j,k,l]*A[i,j,k,l] +
  rD*D[i,j,k,l] - rA*A[i,j,k,l] - (eta+age_rate[i])*A[i,j,k,l] - MDA_A[i,j,k,l]
deriv(A[2:na, 1:nh, 1:num_int, 1:ncc]) <- (1-phi[i,j,k,l])*FOI[i,j,k,l]*Y[i,j,k,l] - FOI[i,j,k,l]*A[i,j,k,l] +
  rD*D[i,j,k,l] - rA*A[i,j,k,l] - (eta+age_rate[i])*A[i,j,k,l] + age_rate[i-1]*A[i-1,j,k,l] - MDA_A[i,j,k,l]

# U - SUBPATENT DISEASE
init_U[,,,] <- user()
dim(init_U) <- c(na,nh,num_int,ncc)
initial(U[,,,]) <- init_U[i,j,k,l]
dim(U) <- c(na,nh,num_int,ncc)

deriv(U[1, 1:nh, 1:num_int, 1:ncc]) <- rA*A[i,j,k,l] - FOI[i,j,k,l]*U[i,j,k,l] - rU*U[i,j,k,l] -
  (eta+age_rate[i])*U[i,j,k,l] - MDA_U[i,j,k,l]
deriv(U[2:na, 1:nh, 1:num_int, 1:ncc]) <- rA*A[i,j,k,l] - FOI[i,j,k,l]*U[i,j,k,l] - rU*U[i,j,k,l] -
  (eta+age_rate[i])*U[i,j,k,l] + age_rate[i-1]*U[i-1,j,k,l] - MDA_U[i,j,k,l]

# P - PROPHYLAXIS
init_P[,,,] <- user()
dim(init_P) <- c(na,nh,num_int,ncc)
initial(P[,,,]) <- init_P[i,j,k,l]
dim(P) <- c(na,nh,num_int, ncc)

deriv(P[1, 1:nh, 1:num_int, 1:ncc]) <- rT*T[i,j,k,l] - rP*P[i,j,k,l] - (eta+age_rate[i])*P[i,j,k,l] + MDA_T[i,j,k,l] + MDA_D[i,j,k,l] + MDA_A[i,j,k,l] + MDA_U[i,j,k,l] + MDA_S[i,j,k,l]
deriv(P[2:na, 1:nh, 1:num_int, 1:ncc]) <- rT*T[i,j,k,l] - rP*P[i,j,k,l] - (eta+age_rate[i])*P[i,j,k,l] +
  age_rate[i-1]*P[i-1,j,k,l]  + MDA_T[i,j,k,l] + MDA_D[i,j,k,l] + MDA_A[i,j,k,l] + MDA_U[i,j,k,l] + MDA_S[i,j,k,l]

# The number of individuals able to acquire clinical malaria
dim(Y) <- c(na,nh,num_int, ncc)
Y[1:na, 1:nh, 1:num_int, 1:ncc] <- S[i,j,k,l]+A[i,j,k,l]+U[i,j,k,l]

# The number of new cases at this timestep
dim(clin_inc) <- c(na,nh,num_int, ncc)
clin_inc[1:na, 1:nh, 1:num_int, 1:ncc] <- phi[i,j,k,l]*FOI[i,j,k,l]*Y[i,j,k,l]

dim(clin_inc_ACD) <- c(na,nh,num_int, ncc)
clin_inc_ACD[1:na, 1:nh, 1:num_int, 1:ncc] <- FOI[i,j,k,l]*Y[i,j,k,l]

# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,,,])
Th <- sum(T[,,,])
Dh <- sum(D[,,,])
Ah <- sum(A[,,,])
Uh <- sum(U[,,,])
Ph <- sum(P[,,,])
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
init_ICM[,,,] <- user()
dim(init_ICM) <- c(na,nh,num_int,ncc)
initial(ICM[,,,]) <- init_ICM[i,j,k,l]
dim(ICM) <- c(na,nh,num_int,ncc)

deriv(ICM[1, 1:nh, 1:num_int, 1:ncc]) <- -1/dCM*ICM[i,j,k,l] + (init_ICM[i,j,k,l]-ICM[i,j,k,l])/x_I[i]
deriv(ICM[2:na, 1:nh, 1:num_int, 1:ncc]) <- -1/dCM*ICM[i,j,k,l] - (ICM[i,j,k,l]-ICM[i-1,j,k,l])/x_I[i]

# ICA - exposure driven immunity
init_ICA[,,,] <- user()
dim(init_ICA) <- c(na,nh,num_int,ncc)
initial(ICA[,,,]) <- init_ICA[i,j,k,l]
dim(ICA) <- c(na,nh,num_int,ncc)

deriv(ICA[1, 1:nh, 1:num_int, 1:ncc]) <- FOI[i,j,k,l]/(FOI[i,j,k,l] * uCA + 1) - 1/dCA*ICA[i,j,k,l] -ICA[i,j,k,l]/x_I[i]
deriv(ICA[2:na, 1:nh, 1:num_int, 1:ncc]) <- FOI[i,j,k,l]/(FOI[i,j,k,l] * uCA + 1) - 1/dCA*ICA[i,j,k,l] - (ICA[i,j,k,l]-ICA[i-1,j,k,l])/x_I[i]

# clinical immunity is a combination of maternal and exposure-driven immunity
dim(IC) <- c(na,nh,num_int,ncc)
IC[,,,] <- ICM[i,j,k,l] + ICA[i,j,k,l]

# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- user()
phi1 <- user() # these parameters characterise the hill function
IC0 <- user() # for probability of clinical disease
kC <- user() # See supplementary materials 1.1.3
dim(phi) <- c(na,nh,num_int,ncc)
phi[1:na,1:nh,1:num_int,1:ncc] <- phi0*((1-phi1)/(1+(IC[i,j,k,l]/IC0)^kC) + phi1)

# IB - infection blocking immunity
init_IB[,,,] <- user()
dim(init_IB) <- c(na,nh,num_int, ncc)
initial(IB[,,,]) <- init_IB[i,j,k,l]
dim(IB) <- c(na,nh,num_int,ncc)

deriv(IB[1, 1:nh, 1:num_int, 1:ncc]) <- EIR[i,j,k,l]/(EIR[i,j,k,l]* uB + 1) - IB[i,j,k,l]/dB - IB[i,j,k,l]/x_I[i]
deriv(IB[2:na, 1:nh, 1:num_int, 1:ncc]) <- EIR[i,j,k,l]/(EIR[i,j,k,l]* uB + 1) - IB[i,j,k,l]/dB - (IB[i,j,k,l]-IB[i-1,j,k,l])/x_I[i]

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- user() # these parameters characterise the hill function for b
b1 <- user() # prob of infection from bite with zero immunity
kB <- user() #
IB0 <- user()
dim(b) <- c(na,nh,num_int,ncc)
b[1:na, 1:nh, 1:num_int, 1:ncc] <- b0 * ((1-b1)/(1+(IB[i,j,k,l]/IB0)^kB)+b1)

# detection immunity
init_ID[,,,] <- user()
dim(init_ID) <- c(na,nh,num_int,ncc)
initial(ID[,,,]) <- init_ID[i,j,k,l]
dim(ID) <- c(na,nh,num_int,ncc)

deriv(ID[1, 1:nh, 1:num_int, 1:ncc]) <- FOI[i,j,k,l]/(FOI[i,j,k,l]*uD + 1) - ID[i,j,k,l]/dID - ID[i,j,k,l]/x_I[i]
deriv(ID[2:na, 1:nh, 1:num_int, 1:ncc]) <- FOI[i,j,k,l]/(FOI[i,j,k,l]*uD + 1) - ID[i,j,k,l]/dID - (ID[i,j,k,l]-ID[i-1,j,k,l])/x_I[i]

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
dim(p_det) <- c(na,nh,num_int,ncc)
p_det[,,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k,l]/ID0)^kD)

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int,ncc)
FOI_lag[1:na, 1:nh, 1:num_int, 1:ncc] <- EIR[i,j,k,l]/B2 * (if(IB[i,j,k,l]==0) b0 else b[i,j,k,l])

# Current FOI depends on humans that have been through the latent period and are
# producing gametocytes
dE <- user() # length of time from infection to gametocytogenesis
dim(FOI) <- c(na,nh,num_int,ncc)
FOI[,,,] <- delay(FOI_lag[i,j,k,l],dE)


# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()

dim(ccov) <- ncc
ccov[] <- user() # proportion of population in each of the correlation coverage groups (sums to 1)


dim(EIR) <- c(na,nh,num_int,ncc)
EIR[,,,1] <- av_human[k] * rel_foi[j] * foi_age[i]/omega*Ivtot
EIR[,,,2] <- av_human[k] * rel_foi[j] * foi_age[i]/omega*Ivtot





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
init_Sx[] <- 0
init_Ex[,] <-0
init_Ix[] <- 0

dim(init_Sx) = eff_len
dim(init_Ex) = c(spor_len, eff_len)
dim(init_Ix) = eff_len

# init_ICA[,,] <- user()
# dim(init_ICA) <- c(na,nh,num_int)

# dim(init_Sx) = c(eff_len)
# dim(init_Ex) = c(spor_len, eff_len)
# dim(init_Ix) = c(eff_len)

initial(Sv) <- init_Sv * mv0
initial(Ev[1:spor_len]) <- init_Ev/spor_len * mv0
dim(Ev) <- spor_len
initial(Iv) <- init_Iv *mv0

initial(Sx[1:eff_len]) = 0
initial(Ex[1:spor_len, 1:eff_len]) = 0
initial(Ix[1:eff_len]) = 0



dim(Sx) = eff_len
dim(Ex) = c(spor_len, eff_len)
dim(Ix) = eff_len

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- user() # infectiousness U -> mosq
cD <- user() # infectiousness D -> mosq
cT <- user() # T -> mosq
gamma1<- user() # fitted value of g_inf characterises cA function
dim(cA) <- c(na,nh,num_int,ncc)
cA[,,,] <- cU + (cD-cU)*p_det[i,j,k,l]^gamma1


# Force of infection from humans to mosquitoes
dim(FOIvijk) <- c(na,nh,num_int,ncc)
omega <- user() #normalising constant for biting rates
FOIvijk[1:na, 1:nh, 1:num_int, 1:ncc] <- (cT*T[i,j,k,l] + cD*D[i,j,k,l] + cA[i,j,k,l]*A[i,j,k,l] + cU*U[i,j,k,l]) * rel_foi[j] * av_mosq[k]*foi_age[i]/(omega*B2)
lag_FOIv=sum(FOIvijk)
# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- user() # latent period in gametocytogenesis
FOIv <- delay(lag_FOIv, delayGam)


# # Number of mosquitoes that become infected at each time point
# ince <- FOIv * Sv
#
# # Number of mosquitoes born (depends on PL, number of larvae)
betaa <- 0.5*PL/dPL
# deriv(Sv) <- -ince - mu*Sv + betaa
# # Exposed mosquitoes go through 10 compartments to simulate 10 day development of
# # sporozoites (and mosquitoes that die off during this wait)
# deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
# deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
# deriv(Iv) <- Ev[10] - mu*Iv
# # Total mosquito population
# mv = Sv+sum(Ev)+Iv

eff_len = user()

haz[] <- user()
dim(haz) = eff_len

dim(mu_vi) = eff_len
mu_vi[1:eff_len] = haz[i]*mu

ivm_cov_par <- user()
ivm_min_age <- user()
ivm_cov = ivm_cov_par*(exp(-ivm_min_age/21))
#IVRM = 0
spor_len = 10
tau_v = spor_len/10

Ivtot = Iv + sum(Ix)
Evtot = sum(Ev) + sum(Ex)
Svtot = Sv + sum(Sx)

mv = Svtot + Evtot + Ivtot



#IVRM = if(t >= IVRM_start[1] && t < (IVRM_start[1]+eff_len)) 1 else if(t >= IVRM_start[2] && t < (IVRM_start[2]+eff_len)) 1 else if (t >= IVRM_start[3] && t < (IVRM_start[3]+eff_len)) 1 else if (t >= IVRM_start[4] && t < (IVRM_start[4]+eff_len)) 1 else if (t >= IVRM_start[5] && t < (IVRM_start[5]+eff_len)) 1 else if (t >= IVRM_start[6] && t < (IVRM_start[6]+eff_len)) 1 else 0

#IVRM_sr = if (t > IVRM_start[1] && t <= IVRM_start[2]) IVRM_start[2] else if (t > IVRM_start[2] && t <= IVRM_start[3]) IVRM_start[3] else if (t > IVRM_start[3] && t <= IVRM_start[4]) IVRM_start[4] else if (t > IVRM_start[4] && t <= IVRM_start[5]) IVRM_start[5] else if (t > IVRM_start[5] && t <= IVRM_start[6]) IVRM_start[6] else IVRM_start[1]
#IVRM_sr = 2000

# ttt[] <- user() already defined above
IVRM_start[]<-user()
dim(IVRM_start) <- length(ttt)

IVRM_sr = interpolate(ttt, IVRM_start)


deriv(Sv) =               if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len  )) betaa - mu*Sv - avB2*Sv*ivm_cov - FOIv*Sv*(1-ivm_cov)                      else betaa - FOIv*Sv - mu*Sv

deriv(Ev[1])  =           if (t >= (IVRM_sr )   && t < (IVRM_sr + eff_len )) FOIv*Sv*(1-ivm_cov) -  mu*Ev[1] - (1/tau_v + avB2*ivm_cov)*Ev[1]          else FOIv*Sv - 1/tau_v*Ev[1] - mu*Ev[1]

deriv(Ev[2:spor_len]) =   if (t >= (IVRM_sr )   && t < (IVRM_sr + eff_len )) 1/tau_v*Ev[i-1]  - avB2*ivm_cov*Ev[i] - 1/tau_v*Ev[i] - mu*Ev[i]           else 1/tau_v*Ev[i-1] - 1/tau_v*Ev[i] - mu*Ev[i]

deriv(Iv) =               if (t >= (IVRM_sr )   && t < (IVRM_sr + eff_len )) -avB2*ivm_cov*Iv + 1/tau_v*Ev[spor_len] - mu*Iv                         else 1/tau_v*Ev[spor_len] -  mu*Iv


deriv(Sx[1:eff_len]) =    if (t >= (IVRM_sr + i )  &&  t < (IVRM_sr + i + 1) && t < (IVRM_sr + eff_len)) (ivm_cov*(avB2-FOIv)*Sv)   - FOIv*Sx[i] - mu_vi[i]*Sx[i]     else -mu_vi[i]*Sx[i] - FOIv*Sx[i]

deriv(Ex[1,1:eff_len]) =  if (t >= (IVRM_sr + j )  &&  t < (IVRM_sr + j + 1) && t < (IVRM_sr + eff_len)) (avB2*ivm_cov*Ev[1] + FOIv*Sv*ivm_cov)  - 1/tau_v*Ex[1,j] - mu_vi[j]*Ex[1,j] + FOIv*Sx[j]  else - 1/tau_v*Ex[1,j] - mu_vi[j]*Ex[1,j]  + FOIv*Sx[j]

deriv(Ex[2:spor_len,1:eff_len]) =  if (t >= (IVRM_sr + j )  &&  t < (IVRM_sr + j + 1) && t < (IVRM_sr + eff_len))  avB2*ivm_cov*Ev[i] + 1/tau_v*Ex[i-1,j]  - 1/tau_v*Ex[i,j] - mu_vi[j]*Ex[i,j]  else -1/tau_v*Ex[i,j] + 1/tau_v*Ex[i-1,j] - mu_vi[j]*Ex[i,j]

deriv(Ix[1:eff_len])    =  if (t >= (IVRM_sr + i )  &&  t < (IVRM_sr + i + 1) && t < (IVRM_sr + eff_len)) avB2*ivm_cov*Iv  + 1/tau_v*Ex[spor_len,i] - mu_vi[i]*Ix[i]    else 1/tau_v*Ex[spor_len,i] - mu_vi[i]*Ix[i]





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

B2 <- user()
avB2 = av*B2


# fitted entomological parameters:
mu0 <- user() # baseline mosquito death rate
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



########################
## INTERVENTION MODEL ##
########################

# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# general parameters
ITN_IRS_on <- user() # days after which interventions begin
num_int <- user() # number of intervention categorys, ITN only, IRS only, neither, both
itn_cov <- user() # proportion of population covered by ITN
irs_cov <- user() # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
dim(cov) <- num_int
cov[1] <- (1-itn_cov)*(1-irs_cov)  # {No intervention}
cov[2] <- itn_cov*(1-irs_cov) #      {ITN only}
cov[3] <- (1-itn_cov)*irs_cov #      {IRS only}
cov[4] <- itn_cov*irs_cov #    {Both ITN and IRS}

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
dim(w) <- num_int
w[1] <- 1
w[2] <- 1 - bites_Bed + bites_Bed*s_ITN
w[3] <- 1 - bites_Indoors + bites_Indoors*(1-r_IRS)*s_IRS
w[4] <- 1 - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN*s_IRS + (bites_Indoors - bites_Bed)*(1-r_IRS)*s_IRS

# probability that mosq feeds during a single attempt for each int. cat.
dim(yy) <- num_int
yy[1] <- 1
yy[2] <- w[2]
yy[3] <- 1 - bites_Indoors + bites_Indoors*(1-r_IRS)
yy[4] <- 1 - bites_Indoors + bites_Bed*(1-r_IRS)*s_ITN + (bites_Indoors - bites_Bed)*(1-r_IRS)

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z) <- num_int
z[1] <- 0
z[2] <- bites_Bed*r_ITN
z[3] <- bites_Indoors*r_IRS
z[4] <- bites_Bed*(r_IRS+ (1-r_IRS)*r_ITN) + (bites_Indoors - bites_Bed)*r_IRS

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

###################
## MODEL OUTPUTS ##
###################

# Outputs for each compartment across the sum across all ages, biting heterogeneities and intervention categories
output(Sout) <- sum(S[,,,])
output(Tout) <- sum(T[,,,])
output(Dout) <- sum(D[,,,])
output(Aout) <- sum(A[,,,])
output(Uout) <- sum(U[,,,])
output(Pout) <- sum(P[,,,])
output(Hout) <- sum(S[,,,]) + sum(T[,,,]) + sum(D[,,,]) +  sum(A[,,,]) + sum(U[,,,]) +  sum(P[,,,])

output(H1out) <- sum(S[,,,1]) + sum(T[,,,1]) + sum(D[,,,1]) +  sum(A[,,,1]) + sum(U[,,,1]) +  sum(P[,,,1])
output(H2out) <- sum(S[,,,2]) + sum(T[,,,2]) + sum(D[,,,2]) +  sum(A[,,,2]) + sum(U[,,,2]) +  sum(P[,,,2])

output(EIR1) <- sum(EIR[,,,1])
output(EIR2) <- sum(EIR[,,,2])

output(FOI1) <- sum(FOI[,,,1])
output(FOI2) <- sum(FOI[,,,2])

output(IB1) <- sum(FOI[,,,1])
output(IB2) <- sum(FOI[,,,2])

output(A1out) <- sum(A[,,,1])
output(A2out) <- sum(A[,,,2])

output(S1out) <- sum(S[,,,1])
output(S2out) <- sum(S[,,,2])

output(MDA_sr) <- MDA_sr

output(MDA_Aout1) <- sum(MDA_A[,,,1])
output(MDA_Aout2) <- sum(MDA_A[,,,2])

# Outputs for clinical incidence and prevalence on a given day

# population densities for each age category
den[] <- user()
dim(den) <- na
# index of the age vector above 59 months
age59 <- user()
# index of the age vector above 5 years
age05 <- user()

# index of the age vector above 2 years
age02 <- user()

# index of the age vector above 10 years
age10 <- user()

dim(prev0to59) <- c(age59,nh,num_int,ncc)
prev0to59[1:age59,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l] * p_det[i,j,k,l]
output(prev) <- sum(prev0to59[,,,])/sum(den[1:age59])

age2to10 = age10 - age02 + 1

dim(prev2to10) <- c(age10,nh,num_int,ncc)
prev2to10[1:(age02-1),,,] = 0
prev2to10[age02:age10,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l] * p_det[i,j,k,l]
output(prev2_10) <- sum(prev2to10[,,,])/sum(den[age02:age10])

# dim(prev0to59_pcr) <- c(age59,nh,num_int,ncc)
# prev0to59_pcr[1:age59,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l] + U[i,j,k,l]
# output(prev_pcr) <- sum(prev0to59_pcr[,,,])/sum(den[1:age59])

dim(prev0to5_pcr) <- c(age05,nh,num_int,ncc)
prev0to5_pcr[1:age05,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l] + U[i,j,k,l]
output(prev_pcr_0to5) <- sum(prev0to5_pcr[,,,])/sum(den[1:age05])

# output(age2to10) = age2to10
# output(age10) = age10




dim(prev_all_pcr) <- c(na,nh,num_int,ncc)
prev_all_pcr[,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l] + U[i,j,k,l]
output(prev_pcr) <- sum(prev_all_pcr[,,,])/sum(den[1:na])

dim(prev_all_slide) <- c(na,nh,num_int,ncc)
prev_all_slide[,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l]*p_det[i,j,k,l]
output(prev_slide) <- sum(prev_all_slide[,,,])




dim(prev0to5_slide) <- c(age05,nh,num_int,ncc)
prev0to5_slide[1:age05,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l]*p_det[i,j,k,l]
output(prev_slide_0to5) <- sum(prev0to5_slide[,,,])/sum(den[1:age05])

dim(prev5to10_slide) <- c(age10,nh,num_int,ncc)
prev5to10_slide[1:(age05-1),,,] = 0
prev5to10_slide[age05:age10,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l] * p_det[i,j,k,l]
output(prev_slide_5to10) <- sum(prev5to10_slide[,,,])/sum(den[age05:age10])

dim(prev0to10_slide) <- c(age10,nh,num_int,ncc)
prev0to10_slide[1:age10,,,] <- T[i,j,k,l] + D[i,j,k,l] + A[i,j,k,l]*p_det[i,j,k,l]
output(prev_slide_0to10) <- sum(prev0to10_slide[,,,])/sum(den[1:age10])


#
# dim(pcr_prev_cc1) <- c(na, nh, num_int, 1)
# pcr_prev_cc1[,,,1] <- T[i,j,k,1] + D[i,j,k,1] + A[i,j,k,1] + U[i,j,k,1]
# output(pcrprev_cc1) <- sum(pcr_prev_cc1[,,,])/ccov[1]
#
# dim(pcr_prev_cc2) <- c(na, nh, num_int, 1)
# pcr_prev_cc2[,,,1] <- T[i,j,k,2] + D[i,j,k,2] + A[i,j,k,2] + U[i,j,k,2]
# output(pcrprev_cc2) <- sum(pcr_prev_cc2[,,,])/ccov[2]

# clinical_incidence_treated

output(clin_inc_treated) <- sum(clin_inc)*365*ft

output(clin_inc_all) <- sum(clin_inc)*365


dim(clin_inc0to5) <- c(age05,nh,num_int,ncc)
clin_inc0to5[1:age05,,,] <- clin_inc[i,j,k,l]
output(clin_inc_0to5) <- sum(clin_inc0to5[,,,])/sum(den[1:age05])*365

dim(clin_inc5to10) <- c(age10,nh,num_int,ncc)
clin_inc5to10[1:(age05),,,] = 0
clin_inc5to10[(age05+1):age10,,,] <- clin_inc[i,j,k,l]
output(clin_inc_5to10) <- sum(clin_inc5to10[,,,])/sum(den[(age05+1):age10])*365

dim(clin_inc0to10) <- c(age10,nh,num_int,ncc)
clin_inc0to10[1:age10,,,] <- clin_inc[i,j,k,l]
output(clin_inc_0to10) <- sum(clin_inc0to10[,,,])/sum(den[1:age10])*365


dim(clin_inc0to10_ACD) <- c(age10,nh,num_int,ncc)
clin_inc0to10_ACD[1:age10,,,] <- clin_inc_ACD[i,j,k,l]
output(clin_inc_0to10_ACD) <- sum(clin_inc0to10_ACD[,,,])/sum(den[1:age10])*365





output(dens) = sum(den[])


output(Svtot) <- Svtot
output(Evtot) <- Evtot
output(Ivtot) <- Ivtot
output(Nvtot) <- mv

output(IVRM_sr) <- IVRM_sr
#output(IVRM) <- IVRM
output(mu) <- mu
output(betaa) <- betaa

