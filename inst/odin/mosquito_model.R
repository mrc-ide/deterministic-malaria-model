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

ento_model_ivm <- odin::odin({
  init_Sv <- user()
  init_Ev_di_human <- user()
  init_Iv_di_human <- user()
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
  lag_FOIv = 0.6 #just set to a value and making it a function of Q0 - anthropophagy and av0
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

  #susceptible mosquitoes either bite on humans or cows

  deriv(Sv) <- -(av0*Q0*prop_human_ivm) - mu*Sv + betaa - (((1-Q0)*av0)*prop_cow_ivm) #is this line correct?

  #infection for susceptible mosquitoes
  deriv(Sv_di_human) <- (av0*Q0*prop_human_ivm)*Sv - mu_ivm*Sv_di_human - ince
  deriv(Ev_di_human) <- ince - incv - mu_ivm*Ev_di_human
  deriv(Iv_di_human) <- incv - mu*Iv_di_human

  #no infection if the susceptibles mosquitoes feed on cows
  deriv(Sv_di_cow) <- (((1-Q0)*av0)*prop_cow_ivm) - (mu_ivm_cow*Sv_di_cow)

  #cow compartment
  init_calf <- user()
  init_cow <- user()

  initial_calf <- 100
  initial_cow <- 0

  cattle = calf + cow

  deriv(calf) <- (b*cattle) - (d*calf) - (m*calf)
  deriv(cow) <- (m*calf) - (d*cow)

  m <- 0.00136 #maturation into cows which can be treated with IVM (an assumption)
  d <- 0.000152 #say they can live for 18 years
  b <- d #match birth and death rate

  prop_cow_ivm <- ((cow)/(cow+calf))*coverage_cow_ivm #proportion of cows with ivermectin
  prop_human_ivm <- user() #coverage of humans with ivm. age and height restrictions...make up a number

  # model options if don't want to use a delayed delay
  #deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
  #deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
  #mv = Sv+sum(Ev)+Iv

  # Total mosquito population
  mv = Sv+Ev+Iv


})





