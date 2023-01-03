##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

require(odin)


##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes


ento_model_original <- odin::odin({
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
})


# initial state values:
ento_model_simple <- odin::odin({
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
  #cU <- user() # infectiousness U -> mosq
  #cD <- user() # infectiousness D -> mosq
  #cT <- user() # T -> mosq
  #gamma1 <- user() # fitted value of gamma1 characterises cA function
  #dim(cA) <- c(na,nh,num_int)
  #cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

  # Force of infection from humans to mosquitoes
  #dim(FOIvijk) <- c(na,nh,num_int)
  omega <- user() #normalising constant for biting rates
  #FOIvijk[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq[k]*foi_age[i]/omega
  #lag_FOIv=sum(FOIvijk)
  lag_FOIv = 0.6*av0*Q0 #just set to a value and making it a function of Q0 - anthropophagy and av0
  Q0 <- user()
  av0 <- 0.333 #1 bite every 3 days

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
  #betaa <- 0.5*PL/dPL
  #betaa <- mv0 * mu0 * theta2
  betaa <- 0.132 #birth rate: set by NC. for simplicity, have birth rate = death rate
  mu <- 0.132 # death rate: set by NC
  mv0 <- 100 #initial mosquito density

  deriv(Sv) <- -ince - mu*Sv + betaa
  deriv(Ev) <- ince - incv - mu*Ev
  deriv(Iv) <- incv - mu*Iv

  # model options if don't want to use a delayed delay
  #deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
  #deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
  #mv = Sv+sum(Ev)+Iv

  # Total mosquito population
  mv = Sv+Ev+Iv
})






