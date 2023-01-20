#MOSQUITO MODEL WITH LARVAE AND IVERMECTIN####

#going to keep the same framework wrt the ivermectin compartments
#but going to add in the larvae and the nets

require(odin)

ivm_model_complex <- odin::odin({
  na <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80) #user() # number of age categories
  nh <- 5 #user() # number of biting heterogeneity categories
  ft <- 0.4 #user() # proportion of cases treated

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

  #NO IVM
  init_Sv <- 50000000#user()
  init_Ev <- 0#user()
  init_Iv <- 0#user()
  initial(Sv) <- init_Sv #* mv0
  initial(Ev) <- init_Ev #* mv0
  #initial(Ev[1:10]) <- init_Ev/10 * mv0 # Options if not using a delayed delay
  #dim(Ev) <- 10
  initial(Iv) <- init_Iv #* mv0


  # cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
  # by age/het/int category, infectiousness depends on p_det which depends on detection immunity
  cU <- 0.006203#user() # infectiousness U -> mosq
  cD <- 0.0676909#user() # infectiousness D -> mosq
  cT <- 0.322 * cD #user() # T -> mosq
  gamma1 <- 1.82425 #user() # fitted value of gamma1 characterises cA function. infectiousness state of A
  dim(cA) <- c(na,nh,num_int)
  cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

  # Force of infection from humans to mosquitoes
  dim(FOIvijk) <- c(na,nh,num_int)
  omega <- 1#user() #normalising constant for biting rates
  FOIvijk[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq[k]*foi_age[i] #/omega
  lag_FOIv=sum(FOIvijk)

  # Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
  delayGam <- 12 #user() # Lag from parasites to infectious gametocytes
  delayMos <- 10 #user() # Extrinsic incubation period.
  FOIv <- delay(lag_FOIv, delayGam)

  # Number of mosquitoes that become infected at each time point
  surv <- exp(-mu*delayMos)
  ince <- FOIv * Sv #rate into Ev compartment
  lag_incv <- ince * surv #need to lag the rate into the infectious compartment
  incv <- delay(lag_incv, delayMos)
  #incv <- lag_incv

  # Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
  betaa <- 0.5*PL/dPL #PL is fully developed pupae and dPL is the developmnent time of the pupae. 0.5 because only interested in females
  #betaa <- mv0 * mu0 * theta2

  #NO IVM####
  deriv(Sv) <- -ince - mu*Sv + betaa
  deriv(Ev) <- ince - incv - mu*Ev
  deriv(Iv) <- incv - mu*Iv

  #IVM humans


  #IVM cattle

  # Total mosquito population
  #mv = Sv+Ev+Iv if no ivermectin addition
  #add on the extra compartments e.g. Svih, Evih, Ivih, Svic, Evic, Ivic


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

  # FITTED entomological parameters:
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
})

