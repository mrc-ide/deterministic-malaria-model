remove(list=ls())
library(odin)

##### stripped versions ##################################
##runs version that has time-varying EIR
run_model_stripped <- function(het_brackets = 5,
                               age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
                               init_EIR = 10,
                               init_ft = 0.4,
                               country = NULL,
                               admin2 = NULL,
                               time = 100,
                               ...){

  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create_stripped(age_vector=age, EIR=init_EIR,ft=init_ft,
                                            model_param_list = mpl, het_brackets=het_brackets,
                                            country = country,
                                            admin_unit = admin2)
  generator <- odin("./inst/odin/odin_model_stripped.R")
  state_use <- state[names(state) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, 1)

  # run model
  start.time <- Sys.time()
  mod_run <- mod$run(tt)
  print(Sys.time()-start.time)
  # shape output
  out <- mod$transform_variables(mod_run)

  # return mod
  return(out)
}

##runs version that has flat EIR from time zero
run_model_stripped_flat <- function(het_brackets = 5,
                               age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
                               init_EIR = 10,
                               init_ft = 0.4,
                               country = NULL,
                               admin2 = NULL,
                               time = 100,
                               ...){

  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create_stripped(age_vector=age, EIR=init_EIR,ft=init_ft,
                                            model_param_list = mpl, het_brackets=het_brackets,
                                            country = country,
                                            admin_unit = admin2)
  generator <- odin("./inst/odin/odin_model_stripped_flat.R")
  state_use <- state[names(state) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, 1)

  # run model
  start.time <- Sys.time()
  mod_run <- mod$run(tt)
  print(Sys.time()-start.time)
  # shape output
  out <- mod$transform_variables(mod_run)

  # return mod
  return(out)
}


model_param_list_create <- function(
  # age, heterogeneity in exposure,
  eta = 1/(21*365),
  rho = 0.85,
  a0 = 2920,
  sigma2 = 1.67,
  max_age = 100*365,
  #  rate of leaving infection states
  rA = 1/195,
  rT = 0.2,
  rD = 0.2,
  rU = 1/110.299,
  rP = 1/15,
  #  human latent period and time lag from asexual parasites to
  dE  = 12,
  delayGam = 12.5,
  # human infectiousness to mosquitoes
  cD  = 0.0676909,
  cT  =  0.322 * cD,
  cU  = 0.006203,
  gamma1  = 1.82425,
  #  Immunity reducing probability of detection
  d1 = 0.160527,
  dID = 3650,
  ID0 = 1.577533,
  kD = 0.476614,
  uD = 9.44512,
  aD = 8001.99,
  fD0 = 0.007055,
  gammaD = 4.8183,
  alphaA = 0.75735,
  alphaU = 0.185624,
  # Immunity reducing probability of infection
  b0 = 0.590076,
  b1 = 0.5,
  dB = 3650,
  IB0 = 43.8787,
  kB = 2.15506,
  uB = 7.19919,
  # Immunity reducing probability of clinical disease
  phi0 = 0.791666,
  phi1 = 0.000737,
  dCA = 10950,
  IC0 = 18.02366,
  kC = 2.36949,
  uCA = 6.06349,
  PM = 0.774368,
  dCM = 67.6952,
  # entomological parameters
  delayMos = 10,
  tau1 = 0.69,
  tau2 = 2.31,
  mu0 = 0.132,
  Q0 = 0.92,
  chi = 0.86,
  bites_Bed = 0.89,
  bites_Indoors = 0.97,
  # larval parameters daily density dependent mortality rate of egg
  muEL = 0.0338,
  muLL = 0.0348,
  muPL = 0.249,
  dEL = 6.64,
  dLL = 3.72,
  dPL = 0.643,
  gammaL = 13.25,
  km = 11,
  cm = 0.05,
  betaL = 21.2,
  # intervention parameters
  num_int = 1,
  itn_cov = 0,
  irs_cov = 0,
  ITN_IRS_on = -1,
  DY = 365,
  d_ITN0 = 0.41,
  r_ITN0 = 0.56,
  r_ITN1 = 0.24,
  r_IRS0 = 0.6,
  d_IRS0 = 1,
  irs_half_life =   0.5 * DY,
  itn_half_life =   2.64 * DY,
  IRS_interval =   1 * DY,
  ITN_interval =   3 * DY,
  EIR_flat=10,
  EIR_times=c(0),
  EIR_vals=c(10),
  ...

){
  # set up param list
  mp_list <- list()

  # catach extra params and place in list
  extra_param_list <- list(...)
  if(length(extra_param_list)>0){
    if(is.list(extra_param_list[[1]])){
      extra_param_list <- extra_param_list[[1]]
    }
  }
  ## EIR time-series
  mp_list$EIR_times<-EIR_times
 mp_list$EIR_vals<-EIR_vals
  mp_list$EIR_flat<-EIR_flat
  ## DEFAULT PARAMS

  # duration of year
  mp_list$DY <- DY

  # age, heterogeneity in exposure
  mp_list$eta <- eta
  mp_list$rho <- rho
  mp_list$a0 <- a0
  mp_list$sigma2 <- sigma2
  mp_list$max_age <- max_age

  # rate of leaving infection states
  mp_list$rA <- rA
  mp_list$rT <- rT
  mp_list$rD <- rD
  mp_list$rU <- rU
  mp_list$rP <- rP

  # human latent period and time lag from asexual parasites to
  # infectiousness
  mp_list$dE <- dE
  mp_list$delayGam <- delayGam

  # infectiousness to mosquitoes
  mp_list$cD <- cD
  mp_list$cT <- cT
  mp_list$cU <- cU
  mp_list$gamma1 <- gamma1

  # Immunity reducing probability of detection
  mp_list$d1 <- d1
  mp_list$dID <- dID
  mp_list$ID0 <- ID0
  mp_list$kD <- kD
  mp_list$uD <- uD
  mp_list$aD <- aD
  mp_list$fD0 <- fD0
  mp_list$gammaD <- gammaD

  # PCR prevalence parameters
  mp_list$alphaA <- alphaA
  mp_list$alphaU <- alphaU

  # anti-infection immunity
  mp_list$b0 <- b0
  mp_list$b1 <- b1
  mp_list$dB <- dB
  mp_list$IB0 <- IB0
  mp_list$kB <- kB
  mp_list$uB <- uB

  # clinical immunity
  mp_list$phi0 <- phi0
  mp_list$phi1 <- phi1
  mp_list$dCA <- dCA
  mp_list$IC0 <- IC0
  mp_list$kC <- kC
  mp_list$uCA <- uCA
  mp_list$PM <- PM
  mp_list$dCM <- dCM

  # entomological parameters
  mp_list$delayMos <- delayMos
  mp_list$tau1 <- tau1
  mp_list$tau2 <- tau2
  mp_list$mu0 <- mu0
  mp_list$Q0 <- Q0
  mp_list$chi <- chi
  mp_list$bites_Bed <- bites_Bed
  mp_list$bites_Indoors <- bites_Indoors
  mp_list$fv0 <- 1 / (tau1 + tau2)
  mp_list$av0 <- Q0 * mp_list$fv0 # daily feeeding rate on humans
  mp_list$Surv0 <- exp(-mu0 * delayMos) # probability of surviving incubation period
  mp_list$p10 <- exp(-mu0 * tau1)  # probability of surviving one feeding cycle
  mp_list$p2 <- exp(-mu0 * tau2)  # probability of surviving one resting cycle

  # larval parameters
  mp_list$muEL <- muEL
  mp_list$muLL <- muLL
  mp_list$muPL <- muPL
  mp_list$dEL <- dEL
  mp_list$dLL <- dLL
  mp_list$dPL <- dPL
  mp_list$gammaL <- gammaL
  mp_list$km <- km
  mp_list$cm <- cm
  mp_list$betaL <- betaL
  # {White et al. 2011 Parasites and Vectors}
  mp_list$eov <- betaL/mu0 * (exp(mu0/mp_list$fv0) - 1)
  mp_list$b_lambda <- (gammaL * muLL/muEL - dEL/dLL + (gammaL - 1) * muLL * dEL)
  mp_list$lambda <- -0.5 * mp_list$b_lambda +
    sqrt(0.25 * mp_list$b_lambda^2 + gammaL * betaL * muLL * dEL/(2 * muEL * mu0 * dLL * (1 + dPL * muPL)))

  # ITN/IRS parameters
  mp_list$itn_cov <- itn_cov
  mp_list$irs_cov <- irs_cov

  mp_list$num_int <- num_int
  # Catch all: Not defined the correct number of interventions
  if (itn_cov > 0 & num_int == 1){
    stop(message("Incorrect number of interventions for definied ITN coverage. Please ensure you have correctly
                 specified the number of interventions."))
  }
  if (irs_cov > 0 & num_int < 3){
    stop(message("Incorrect number of interventions for definied IRS coverage. Please ensure you have correctly
                 specified the number of interventions."))
  }

  # Sets start time of coverage
  mp_list$ITN_IRS_on <- ITN_IRS_on

  # Sets population split as coverage
  # {No intervention} {ITN only} {IRS only} {Both ITN and IRS}
  cov <- c((1 - itn_cov) * (1 - irs_cov), itn_cov * (1 - irs_cov), (1 - itn_cov) * irs_cov, itn_cov * irs_cov)
  cov <- cov[1:mp_list$num_int]
  mp_list$cov <- cov

  mp_list$d_ITN0 <- d_ITN0
  mp_list$r_ITN0 <- r_ITN0
  mp_list$r_ITN1 <- r_ITN1
  mp_list$r_IRS0 <- r_IRS0
  mp_list$d_IRS0 <- d_IRS0
  mp_list$irs_half_life <- irs_half_life
  mp_list$itn_half_life <- itn_half_life
  mp_list$IRS_interval <- IRS_interval
  mp_list$ITN_interval <- ITN_interval
  mp_list$irs_loss <- log(2)/mp_list$irs_half_life
  mp_list$itn_loss <- log(2)/mp_list$itn_half_life

  # check that none of the spare parameters in the extra
  if(sum(!is.na(match(names(extra_param_list),names(mp_list))))!=0){

    stop (message(cat("Extra params in ... share names with default param names. Please check:\n",
                      names(extra_param_list)[!is.na(match(names(extra_param_list),names(mp_list)))]
    )
    ))
  }

  return(append(mp_list,extra_param_list))
}

equilibrium_init_create_stripped <- function(age_vector, het_brackets,
                                             country = NULL, admin_unit = NULL, ft,
                                             EIR, model_param_list)
{

  # mpl is shorter :)
  mpl <- model_param_list

  ## Check Parameters
  if(!is.numeric(age_vector)) stop("age_vector provided is not numeric")
  if(!is.numeric(het_brackets)) stop("het_brackets provided is not numeric")
  if(!(is.null(country) | is.character(country))) stop("country specified is not character string")
  if(!(is.null(admin_unit) | is.character(admin_unit))) stop("admin_unit specified is not character string")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(EIR)) stop("EIR provided is not numeric")

  ## Handle parameters
  # database for admin units is all in Latin-ASCII for CRAN reasons so must
  # encode parameters accordingly
  if(!is.null(country)) country <- stringi::stri_trans_general(country,"Latin-ASCII")
  if(!is.null(admin_unit)) admin_unit <- stringi::stri_trans_general(admin_unit, "Latin-ASCII")

  ## population demographics
  age <- age_vector * mpl$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")
  age0 <- 2
  age1 <- 10

  age_rate <- age_width <- age_mid_point <- den <- c()
  for (i in 1:(na-1))
  {
    age_width[i] <- age[i+1] - age[i]
    age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group

  }
  age_rate[na] = 0


  den <- 1/(1 + age_rate[1]/mpl$eta)
  for (i in 1:(na-1))
  {
    den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + mpl$eta)  # proportion in each age_vector group
  }

  age59 <- which(age_vector * 12 > 59)[1] - 1  # index of age vector before age is >59 months
  age05 <- which(age_vector > 5)[1] - 1  # index of age vector before age is 5 years

  ## force of infection
  foi_age <- c()
  for (i in 1:na)
  {
    foi_age[i] <- 1 - (mpl$rho * exp(-age[i]/mpl$a0))  #force of infection for each age group
  }
  fden <- foi_age * den
  omega <- sum(fden)  #normalising constant

  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  rel_foi <- exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x)/sum(het_wt * exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x))

  ## EIR
  EIRY_eq <- EIR  # initial annual EIR
  EIRd_eq <- EIRY_eq/mpl$DY
  EIR_eq <- outer(foi_age, rel_foi) * EIRd_eq

  ## Immunity and FOI
  x_I <- den[1]/mpl$eta
  for (i in 2:na)
  {
    x_I[i] <- den[i]/(den[i - 1] * age_rate[i - 1])  #temporary variables
  }
  fd <- 1 - (1 - mpl$fD0)/(1 + (age/mpl$aD)^mpl$gammaD)

  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na)
  {
    age20i[i] <- ifelse(age[i] >= (20 * mpl$DY) & age[i - 1] < (20 * mpl$DY),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * mpl$DY - age[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])

  # finding initial values for all immunity states
  IB_eq <- matrix(0, na, nh)
  FOI_eq <- matrix(0, na, nh)
  ID_eq <- matrix(0, na, nh)
  ICA_eq <- matrix(0, na, nh)
  ICM_init_eq <- vector(length = nh, mode = "numeric")
  ICM_eq <- matrix(0, na, nh)
  cA_eq <- matrix(0, na, nh)
  FOIvij_eq <- matrix(0, na, nh)
  p_det_eq <- matrix(0, na, nh)
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      IB_eq[i, j] <- (ifelse(i == 1, 0, IB_eq[i - 1, j]) +
                        EIR_eq[i,j]/(EIR_eq[i, j] * mpl$uB + 1) * x_I[i])/(1 + x_I[i]/mpl$dB)
      FOI_eq[i, j] <- EIR_eq[i, j] * ifelse(IB_eq[i, j] == 0, mpl$b0,
                                            mpl$b0 * ((1 - mpl$b1)/(1 + (IB_eq[i, j]/mpl$IB0)^mpl$kB) + mpl$b1))
      ID_eq[i, j] <- (ifelse(i == 1, 0, ID_eq[i - 1, j]) +
                        FOI_eq[i, j]/(FOI_eq[i, j] * mpl$uD + 1) * x_I[i])/(1 + x_I[i]/mpl$dID)
      ICA_eq[i, j] <- (ifelse(i == 1, 0, ICA_eq[i - 1, j]) +
                         FOI_eq[i,j]/(FOI_eq[i, j] * mpl$uCA + 1) * x_I[i])/(1 + x_I[i]/mpl$dCA)
      p_det_eq[i, j] <- mpl$d1 + (1 - mpl$d1)/(1 + fd[i] * (ID_eq[i, j]/mpl$ID0)^mpl$kD)
      cA_eq[i, j] <- mpl$cU + (mpl$cD - mpl$cU) * p_det_eq[i, j]^mpl$gamma1
    }
  }
  # needs to be calculated after because it references ICA
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      ICM_init_eq[j] <- mpl$PM * (ICA_eq[age20l, j] + age_20_factor *
                                    (ICA_eq[age20u, j] - ICA_eq[age20l, j]))
      ICM_eq[i, j] <- ifelse(i == 1,
                             ICM_init_eq[j], ICM_eq[i - 1,j])/(1 + x_I[i]/mpl$dCM)
    }
  }

  IC_eq <- ICM_eq + ICA_eq
  phi_eq <- mpl$phi0 * ((1 - mpl$phi1)/(1 + (IC_eq/mpl$IC0)^mpl$kC) + mpl$phi1)


  # human states
  gamma <- mpl$eta + c(age_rate[1:(na - 1)], 0)
  delta <- c(mpl$eta, age_rate[1:(na - 1)])

  betaT <- matrix(rep(mpl$rT + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD <- matrix(rep(mpl$rD + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP <- matrix(rep(mpl$rP + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)

  aT <- FOI_eq * phi_eq * ft/betaT
  aD <- FOI_eq * phi_eq * (1 - ft)/betaD
  aP <- mpl$rT * aT/betaP

  Z_eq <- array(dim = c(na, nh, 4))
  Z_eq[1, , 1] <- den_het[1, ]/(1 + aT[1, ] + aD[1, ] + aP[1, ])
  Z_eq[1, , 2] <- aT[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 3] <- aD[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 4] <- aP[1, ] * Z_eq[1, , 1]

  for (j in 1:nh)
  {
    for (i in 2:na)
    {
      Z_eq[i, j, 1] <- (den_het[i, j] - delta[i] * (Z_eq[i - 1, j, 2]/betaT[i, j] +
                                                      Z_eq[i - 1, j, 3]/betaD[i, j] +
                                                      (mpl$rT *  Z_eq[i - 1, j, 2]/betaT[i, j]
                                                       + Z_eq[i - 1, j, 4])/betaP[i, j]))/(1 + aT[i, j] + aD[i, j] + aP[i, j])
      Z_eq[i, j, 2] <- aT[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 2]/betaT[i, j]
      Z_eq[i, j, 3] <- aD[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 3]/betaD[i, j]
      Z_eq[i, j, 4] <- aP[i, j] * Z_eq[i, j, 1] + delta[i] * (mpl$rT *
                                                                Z_eq[i - 1, j, 2]/betaT[i, j] + Z_eq[i - 1, j, 4])/betaP[i,j]

    }
  }

  Y_eq <- matrix(Z_eq[, , 1], nrow = na, ncol=nh)
  T_eq <- matrix(Z_eq[, , 2], nrow = na, ncol=nh)
  D_eq <- matrix(Z_eq[, , 3], nrow = na, ncol=nh)
  P_eq <- matrix(Z_eq[, , 4], nrow = na, ncol=nh)

  betaS <- apply(FOI_eq, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaA <- apply(FOI_eq * phi_eq + mpl$rA, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaU <- apply(FOI_eq + mpl$rU, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)

  A_eq <- matrix(ncol = nh, nrow = na)
  U_eq <- matrix(ncol = nh, nrow = na)
  S_eq <- matrix(ncol = nh, nrow = na)

  for (i in 1:na)
  {
    for (j in 1:nh)
    {
      A_eq[i, j] <- (delta[i] * ifelse(i == 1, 0, A_eq[i - 1, j]) +
                       FOI_eq[i, j] * (1 - phi_eq[i, j]) * Y_eq[i, j] +
                       mpl$rD * D_eq[i,j])/(betaA[i, j] + FOI_eq[i, j] * (1 - phi_eq[i, j]))
      U_eq[i, j] <- (mpl$rA * A_eq[i, j] + delta[i] * ifelse(i == 1,
                                                             0, U_eq[i - 1, j]))/betaU[i, j]
      S_eq[i, j] <- Y_eq[i, j] - A_eq[i, j] - U_eq[i, j]
      #FOIvij_eq[i, j] <- foi_age[i] * mpl$av0 * (mpl$cT * T_eq[i, j] + mpl$cD *
       #                                            D_eq[i, j] + cA_eq[i, j] * A_eq[i, j] + mpl$cU * U_eq[i, j]) * rel_foi[j]/omega
    }
  }

  # mosquito states
  #FOIv_eq <- sum(FOIvij_eq)
  #Iv_eq <- FOIv_eq * mpl$Surv0/(FOIv_eq + mpl$mu0)
  #Sv_eq <- mpl$mu0 * Iv_eq/(FOIv_eq * mpl$Surv0)
  #Ev_eq <- 1 - Sv_eq - Iv_eq

  # mosquito density needed to give this EIR
  #mv0 <- omega * EIRd_eq/(Iv_eq * mpl$av0)

  # larval states
  #K0 <- 2 * mv0 * mpl$dLL * mpl$mu0 * (1 + mpl$dPL * mpl$muPL) * mpl$gammaL * (mpl$lambda + 1)/(mpl$lambda/(mpl$muLL *
                                                                                                             # mpl$dEL) - 1/(mpl$muLL * mpl$dLL) - 1)
  #PL_eq <- 2 * mpl$dPL * mpl$mu0 * mv0
  #LL_eq <- mpl$dLL * (mpl$muPL + 1/mpl$dPL) * PL_eq
  #EL_eq <- (LL_eq/mpl$dLL + mpl$muLL* LL_eq * (1 + mpl$gammaL * LL_eq/K0))/(1/mpl$dEL - mpl$muLL * mpl$gammaL * LL_eq/K0)

  IB_eq = array(IB_eq, c(na, nh))
  ID_eq = array(ID_eq, c(na, nh))
  ICA_eq = array(ICA_eq, c(na, nh))
  ICM_eq = array(ICM_eq, c(na, nh))
  #ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
  # better het bounds for equilbirum initialisation in individual model
  zetas <- rlnorm(n = 1e5,meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  while(sum(zetas>100)>0){
    zetas[zetas>100] <- rlnorm(n = sum(zetas>100),meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  }

  wt_cuts <- round(cumsum(het_wt)*1e5)
  zeros <- which(wt_cuts==0)
  wt_cuts[zeros] <- 1:length(zeros)
  larges <- which(wt_cuts==1e5)
  wt_cuts[larges] <- (1e5 - (length(larges)-1)):1e5
  wt_cuts <- c(0,wt_cuts)
  het_bounds <- sort(zetas)[wt_cuts]
  het_bounds[length(het_bounds)] <- (mpl$max_age/365)+1

  ## collate init
  res <- list(init_S = S_eq, init_T = T_eq, init_D = D_eq, init_A = A_eq, init_U = U_eq,
              init_P = P_eq, init_Y = Y_eq, init_IB = IB_eq, init_ID = ID_eq, init_ICA = ICA_eq,
              init_ICM = ICM_eq, ICM_init_eq = ICM_init_eq,
              #init_Iv = Iv_eq, init_Sv = Sv_eq,init_Ev = Ev_eq, init_PL = PL_eq, init_LL = LL_eq, init_EL = EL_eq,
              age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
              foi_age = foi_age, rel_foi = rel_foi,
              na = na, nh = nh, x_I = x_I,
              #omega = omega, K0 = K0, mv0 = mv0,theta_c = theta_c, FOIv_eq = FOIv_eq,
              FOI = FOI_eq, EIR_eq = EIR_eq, cA_eq = cA_eq,
              den = den, age59 = age59, age05 = age05,
              # ssa0 = ssa0, ssa1 = ssa1,ssa2 = ssa2, ssa3 = ssa3, ssb1 = ssb1, ssb2 = ssb2, ssb3 = ssb3,pi = pi,
               age = age_vector*mpl$DY, ft = ft,
              betaS = betaS, betaA = betaA, betaU = betaU, FOIvij_eq=FOIvij_eq,
              age_mid_point = age_mid_point, het_bounds = het_bounds,
              age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)

  res <- append(res,mpl)

  return(res)
}

