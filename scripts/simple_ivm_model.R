#simple ivermectin model

ivm_mod <- odin::odin({
  #NO IVM
  init_Sv <- 5000 #user()
  init_Ev <- 0 #user()
  init_Iv <- 0 #user()
  initial(Sv) <- init_Sv * mv0
  initial(Ev) <- init_Ev * mv0
  #initial(Ev[1:10]) <- init_Ev/10 * mv0 # Options if not using a delayed delay
  #dim(Ev) <- 10
  initial(Iv) <- init_Iv * mv0

  #IVM on humans. MAY NEED TO ADD mv0 terms
  init_Svih <- 0
  init_Evih <- 0
  init_Ivih <- 0
  initial(Svih) <- init_Svih*mv0
  initial(Evih) <- init_Evih*mv0
  initial(Ivih) <- init_Ivih*mv0

  #IVM on cattle. MAY NEED TO ADD mv0 terms
  init_Svic <- 0
  init_Evic <- 0
  init_Ivic <- 0
  initial(Svic) <- init_Svic*mv0
  initial(Evic) <- init_Evic*mv0
  initial(Ivic) <- init_Ivic*mv0

  # cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
  # by age/het/int category, infectiousness depends on p_det which depends on detection immunity
  cU <- user() # infectiousness U -> mosq
  cD <- user() # infectiousness D -> mosq
  cT <- user() # T -> mosq
  gamma1 <- user() # fitted value of gamma1 characterises cA function. infectiousness state of A
  dim(cA) <- c(na,nh,num_int)
  cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

  # Force of infection from humans to mosquitoes
  dim(FOIvijk) <- c(na,nh,num_int)
  omega <- 1#user() #normalising constant for biting rates
  FOIvijk[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq[k]*foi_age[i] /omega
  lag_FOIv=sum(FOIvijk)

  # Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
  delayGam <- 12.5#user() # Lag from parasites to infectious gametocytes
  delayMos <- 10 #user() # Extrinsic incubation period.
  FOIv <- delay(lag_FOIv, delayGam)

  # Number of mosquitoes that become infected at each time point
  surv <- exp(-mu*delayMos) #may have to change this after compare IVM v net death
  ince <- FOIv * Sv #rate into Ev compartment
  lag_incv <- ince * surv #need to lag the rate into the infectious compartment
  incv <- delay(lag_incv, delayMos)
  #incv <- lag_incv

  # Number of mosquitoes that become infected at each time point in human IVM compartments
  ince_ih <- FOIv*Svih
  lag_incv_ih <- ince_ih*surv
  incv_ih <- delay(lag_incv_ih, delayMos)
  #incv_ih <- lag_incv_ih

  #Number of mosquitoes that become infected at each time point in cattle IVm compartments
  ince_ic <- FOIv*Svic
  lag_incv_ic <- ince_ic*surv
  incv_ic <- delay(lag_incv_ic, delayMos)
  #incv_ic <- lag_incv_ic

  # Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
  #betaa <- 0.5*PL/dPL #PL is fully developed pupae and dPL is the developmnent time of the pupae. 0.5 because only interested in females
  betaa <- mu0
  #betaa <- mv0 * mu0 * theta2

  #ivermectin parameters####
  gamma_h <- user() #prop human pop treated with IVM. Will be subject to eligibility criteria so this param is going to become more complicated (will depend on age structure)
  gamma_c < user() #prop cattle treated with IVM (are there eligibility criteria?)

  ivm_human_eff_cov <- fv*Q * gamma_h #effective coverage of IVM on humans given biting rate on humans. don't know if this should be a*Q0 or av_human

  ivm_cow_eff_cov <-fv*(1-Q)*gamma_c

  mu_h <- mu0 + 0.628 #excess mort due to IVM on humans (relative to baseline mort)
  mu_c <- mu0 + 0.628 #excess mort due to IVM on cattle (relative to baseline mort)

  #IVERMECTIN INTEGRATION####

  #no IVM
  deriv(Sv) <- -ince - (ivm_human_eff_cov*Sv) - (ivm_cow_eff_cov*Sv) - mu*Sv + betaa
  deriv(Ev) <- ince - incv -(ivm_human_eff_cov*Ev) - (ivm_cow_eff_cov*Ev) - mu*Ev
  deriv(Iv) <- incv - (ivm_human_eff_cov*Iv) - (ivm_cow_eff_cov*Iv) -mu*Iv

  #IVM humans
  deriv(Svih) <- -ince_ih + (ivm_human_eff_cov*Sv) - (mu_h*Svih)
  deriv(Evih) <- ince_ih - incv_ih + (ivm_human_eff_cov*Ev) - (mu_h*Evih)
  deriv(Ivih) <- incv_ih + (ivm_human_eff_cov*Iv) - (mu_h*Ivih)

  #IVM on cattle
  deriv(Svic) <- - ince_ic + (ivm_cow_eff_cov*Sv) - (mu_c*Svic)
  deriv(Evic) <- ince_ic - incv_ic + (ivm_cow_eff_cov*Ev) - (mu_c*Evic)
  deriv(Ivic) <- incv_ic + (ivm_cow_eff_cov*Iv) - (mu_c*Ivic)

  # Total mosquito population
  mv = Sv+Ev+Iv+Svih+Evih+Ivih+Svic+Evic+Ivic #if no ivermectin addition
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
  #dLL <- user() # development time of larvae
  #dPL <- user() #development time of pupae
  #dEL <- user() #development time of early stage
  #muLL <- user() #daily density dep. mortality rate of larvae
  #muPL <- user() #daily den. dep. mortality rate of pupae
  #muEL <- user() #daily den. dep. mortality rate of early stage
  #gammaL <- user() # eff. of den. dep. on late stage relative to early stage

  # FITTED entomological parameters:
  mv0 <- 100#user() # initial mosquito density (V/H Ratio - I can set this)
  mu0 <- 0.132#user() # baseline mosquito death rate
  tau1 <- 0.69#user() # duration of host-seeking behaviour
  tau2 <- 2.31 #user() # duration of resting behaviour
  p10 <- exp(-mu0 * tau1) # prob of surviving 1 feeding cycle
  p2 <- exp(-mu0 * tau2) #prob of surviving one resting cycle
  betaL <- 21.2#user() # maximum number of eggs per oviposition per mosq

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
  #init_PL <- user()
  #initial(PL) <- init_PL
  #init_LL <- user()
  #initial(LL) <- init_LL
  #init_EL <- user()
  #initial(EL) <- init_EL

  # (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
  #deriv(EL) <- beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL
  # egg hatching - den. dep. mortality - maturing larvae
  #deriv(LL) <- EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL
  # pupae - mortality - fully developed pupae
  #deriv(PL) <- LL/dLL - muPL*PL - PL/dPL

  # model options if don't want to use a delayed delay
  #deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
  #deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
  #mv = Sv+sum(Ev)+Iv

  # Param checking outputs
  output(mu) <- mu
  output(mu_h) <- mu_h
  output(mu_c) <- mu_c
  output(beta_larval) <- beta_larval
  output(KL) <- KL
  output(mv) <- mv
  output(Q) <- Q
  output(Q0) <- Q0
  output(wh) <- wh
  output(lag_incv_ic) <- lag_incv_ic
  output(lag_incv_ih) <- lag_incv_ih
  output(lag_incv) <- lag_incv
  output(d_ITN) <- d_ITN
  output(r_ITN) <- r_ITN
  output(s_ITN) <- s_ITN
  output(d_IRS) <- d_IRS
  output(r_IRS) <- r_IRS
  output(s_IRS) <- s_IRS
  output(cov[]) <- TRUE
  output(K0) <- K0
})
