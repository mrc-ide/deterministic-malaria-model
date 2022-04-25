#------------------------------------------------
#' Model Parameter List Creation: P. falciparum
#'
#' The function \code{model_param_list_create} creates list of model parameters for
#' P. falciparum  to be used within \code{equilibrium_init_create}
#'
#' @param eta Death rate for expoential population distribtuion, i.e. 1/Mean Population Age. Default = 0.0001305
#' @param rho Age-dependent biting parameter. Default = 0.85
#' @param a0 Age-dependent biting parameter. Default = 2920
#' @param sigma2 Variance of the log heterogeneity in biting rates. Default = 1.67
#' @param max_age Maximum age in days. Default = 100*365
#' @param rA Rate of leaving asymptomatic infection. Default = 0.00512821
#' @param rT Rate of leaving treatment. Default = 0.2
#' @param rD Rate of leaving clinical disease. Default = 0.2
#' @param rU Rate of recovering from subpatent infection. Default = 0.00906627
#' @param rP Rate of leaving prophylaxis. Default = 0.06666667
#' @param dE Latent period of human infection. Default = 12
#' @param delayGam Lag from parasites to infectious gametocytes. Default = 12.5
#' @param cD Untreated disease contribution to infectiousness. Default = 0.0676909
#' @param cT Treated disease contribution to infectiousness. Default =   0.322 * cD
#' @param cU Subpatent disease contribution to infectiousness. Default = 0.006203
#' @param gamma1 Parameter for infectiousness of state A. Default = 1.82425
#' @param d1 Minimum probability due to maximum immunity. Default = 0.160527
#' @param dID Inverse of decay rate. Default = 3650
#' @param ID0 Scale parameter. Default = 1.577533
#' @param kD Shape parameter. Default = 0.476614
#' @param uD Duration in which immunity is not boosted. Default = 9.44512
#' @param aD Scale parameter relating age to immunity. Default = 8001.99
#' @param fD0 Time-scale at which immunity changes with age. Default = 0.007055
#' @param gammaD Shape parameter relating age to immunity. Default = 4.8183
#' @param alphaA PCR detection probability parameters state A. Default = 0.757
#' @param alphaU PCR detection probability parameters state U. Default = 0.186
#' @param b0 Maximum probability due to no immunity. Default = 0.590076
#' @param b1 Maximum relative reduction due to immunity. Default = 0.5
#' @param dB Inverse of decay rate. Default = 3650
#' @param IB0 Scale parameter. Default = 43.8787
#' @param kB Shape parameter. Default = 2.15506
#' @param uB Duration in which immunity is not boosted. Default = 7.19919
#' @param phi0 Maximum probability due to no immunity. Default = 0.791666
#' @param phi1 Maximum relative reduction due to immunity. Default = 0.000737
#' @param dCA Inverse of decay rate. Default = 10950
#' @param IC0 Scale parameter. Default = 18.02366
#' @param kC Shape parameter. Default = 2.36949
#' @param uCA Duration in which immunity is not boosted. Default = 6.06349
#' @param PM New-born immunity relative to mother’s. Default = 0.774368
#' @param dCM Inverse of decay rate of maternal immunity. Default = 67.6952
#' @param delayMos Extrinsic incubation period. Default = 10
#' @param tau1 Duration of host seeking, assumed to be constant between species. Default = 0.69
#' @param tau2 Duration of mosquito resting after feed. Default = 2.31
#' @param mu0 Daily mortality of adult mosquitos. Default = 0.132
#' @param Q0 Anthrophagy probability. Default = 0.92
#' @param chi Endophily probability. Default = 0.86
#' @param bites_Bed Percentage of bites indoors and in bed. Default = 0.89
#' @param bites_Indoors Percentage of bites indoors . Default = 0.97
#' @param muEL Per capita daily mortality rate of early stage larvae (low density). Default = 0.0338
#' @param muLL Per capita daily mortality rate of late stage larvae (low density). Default = 0.0348
#' @param muPL Per capita daily mortality rate of pupae. Default = 0.249
#' @param dEL Development time of early stage larvae. Default = 6.64
#' @param dLL Development time of late stage larvae. Default = 3.72
#' @param dPL Development time of pupae. Default = 0.643
#' @param gammaL Relative effect of density dependence on late instars relative to early instars. Default = 13.25
#' @param km Seasonal carrying capacity. Default = 11
#' @param cm Seasonal birth rate. Default = 0.05
#' @param betaL Number of eggs laid per day per mosquito. Default = 21.2
#' @param num_int Number of intervention parameters.  Default = 4
#' @param itn_cov The proportion of people that use an ITN. Default = 0
#' @param irs_cov The proportion of people living in houses that have been sprayed. Default = 0
#' @param ITN_IRS_on Time of ITN and IRS to be activated. Default = -1, i.e. never.
#' @param DY Duration of year (days). Default = 365
#' @param d_ITN0 Probability of dying with an encounter with ITN (max). Default = 0.41
#' @param r_ITN0 Probability of repeating behaviour with ITN (max). Default = 0.56
#' @param r_ITN1 Probability of repeating behaviour with ITN (min). Default = 0.24
#' @param r_IRS0 Probability of repeating behaviour with IRS (min). Default = 0.6
#' @param d_IRS0 Probability of dying with an encounter with IRS (max). Default = 1
#' @param irs_half_life IRS half life. Default =   0.5 * DY
#' @param itn_half_life ITN half life. Default =   2.64 * DY
#' @param IRS_interval How long before IRS is repeated, i.e. when IRS decay = 1. Default =   1 * DY
#' @param ITN_interval How long before ITN is repeated, i.e. when IRS decay = 1.  Default =   3 * DY
#' @param ... Any other parameters needed for non-standard model. If they share the same name
#' as any of the defined parameters \code{model_param_list_create} will stop. You can either write
#' any extra parameters you like individually, e.g. model_param_list_create(extra1 = 1, extra2 = 2)
#' and these parameters will appear appended to the returned list, or you can pass explicitly
#' the ellipsis argument as a list created before, e.g. model_param_list_create(...=list(extra1 = 1, extra2 = 2))
#'
#' @export
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

#' Model Parameter List Creation: P. vivax
#'
#' The function \code{vixax_model_param_list_create} creates list of model parameters for
#' P. vivax to be used within \code{vivax_equilibrium_init_create}
#'
#' @param rho Age-dependent biting parameter. Default = 0.85
#' @param a0 Age-dependent biting parameter. Default = 2920
#' @param sigma2 Variance of the log heterogeneity in biting rates. Default = 1.29
#' @param max_age Maximum age in days. Default = 100*365
#' @param rLM Rate of leaving light microscopy-detectable infection. Default = 0.1
#' @param rT Rate of leaving treatment. Default = 1
#' @param rD Rate of leaving clinical disease. Default = 0.2
#' @param rP Rate of leaving prophylaxis. Default = 0.036
#' @param dPCR_min Minimum duration of PCR-detectable infection. Default = 10
#' @param dPCR_max Maximum duration of PCR-detectable infection. Default = 70
#' @param dE Latent period of human infection. Default = 10
#' @param cD Untreated disease contribution to human-to-mosquito infectiousness. Default = 0.8
#' @param cT Treated disease contribution to human-to-mosquito infectiousness. Default = 0.4
#' @param cPCR PCR-detectable infection contribution to infectiousness. Default = 0.035
#' @param cLM LM-detectable infection contribution to infectiousness. Default = 0.1
#' @param b Mosquito-to-human infection probability. Default = 0.5
#' @param d_par Duration of anti-parasite immunity. Default = 3650
#' @param u_par Anti-parasite immune boosting refractory period. Default = 19.77
#' @param phi_LM_max Probability of LM-detectable infection with no immunity. Default = 0.8918
#' @param phi_LM_min Probability of LM-detectable infection with full immunity. Default = 0.0043
#' @param A_LM_50pc Anti-parasite immunity for 50 percent reduction in LM-detectable infection. Default = 27.52
#' @param K_LM Shape parameter for LM-detectable infection. Default = 2.403
#' @param A_PCR_50pc Anti-parasite immunity for 50 percent reduction in duration of PCR-detectable infection. Default = 9.9
#' @param K_PCR Shape parameter for duration of PCR-detectable infection. Default = 4.602
#' @param d_clin Duration of clinical immunity. Default = 10950
#' @param u_clin Clinical immune boosting refractory period. Default = 7.85
#' @param phi_D_max Probability of clinical episode with no immunity. Default = 0.8605
#' @param phi_D_min Probability of clinical episode with full immunity. Default = 0.018
#' @param A_D_50pc Anti-parasite immunity for 50 percent reduction in clinical episode. Default = 11.538
#' @param K_D Shape parameter for clinical disease. Default = 2.250
#' @param p_mat New-born immunity relative to mother’s. Default = 0.421
#' @param d_mat Duration of maternal immunity. Default = 35.148
#' @param gamma_L Rate of hypnozoite clearance. Default = 0.0026
#' @param ff Rate of relapse. Default = 0.024
#' @param delayMos Extrinsic incubation period/duration of sporogony. Default = 8.4
#' @param tau1 Duration of host seeking, assumed to be constant between species. Default = 0.69
#' @param tau2 Duration of mosquito resting after feed. Default = 2.31
#' @param mu0 Daily mortality of adult mosquitoes. Default = 0.167
#' @param Q0 Anthrophagy probability. Default = 0.5
#' @param chi Endophily probability. Default = 0.9
#' @param bites_Bed Percentage of bites indoors and in bed. Default = 0.9
#' @param bites_Indoors Percentage of bites indoors . Default = 0.9
#' @param muEL Per capita daily mortality rate of early stage larvae (low density). Default = 0.0338
#' @param muLL Per capita daily mortality rate of late stage larvae (low density). Default = 0.0348
#' @param muPL Per capita daily mortality rate of pupae. Default = 0.249
#' @param dEL Development time of early stage larvae. Default = 6.64
#' @param dLL Development time of late stage larvae. Default = 3.72
#' @param dPL Development time of pupae. Default = 0.643
#' @param gammaL Relative effect of density dependence on late instars relative to early instars. Default = 13.25
#' @param betaL Number of eggs laid per day per mosquito. Default = 21.2
#' @param num_int Number of intervention parameters.  Default = 4
#' @param itn_cov The proportion of people that use an ITN. Default = 0
#' @param irs_cov The proportion of people living in houses that have been sprayed. Default = 0
#' @param ITN_IRS_on Time of ITN and IRS to be activated. Default = -1, i.e. never.
#' @param DY Duration of year (days). Default = 365
#' @param d_ITN0 Probability of dying with an encounter with ITN (max). Default = 0.41
#' @param r_ITN0 Probability of repeating behaviour with ITN (max). Default = 0.56
#' @param r_ITN1 Probability of repeating behaviour with ITN (min). Default = 0.24
#' @param r_IRS0 Probability of repeating behaviour with IRS (min). Default = 0.6
#' @param d_IRS0 Probability of dying with an encounter with IRS (max). Default = 1
#' @param irs_half_life IRS half life. Default =   0.5 * DY
#' @param itn_half_life ITN half life. Default =   2.64 * DY
#' @param IRS_interval How long before IRS is repeated, i.e. when IRS decay = 1. Default =   1 * DY
#' @param ITN_interval How long before ITN is repeated, i.e. when IRS decay = 1.  Default =   3 * DY
#' @param stop_rebound_switch Set to 1 to prevent a rebound in cases after decades due to relapse from very small hypnozoite reservoir in the population. Default = 0
#' @param hypnozoite_prev_threshold Hypnozoite prevalence in the population at which stop_rebound_switch is triggered. Default = 1e-05
#' @param ... Any other parameters needed for non-standard model. If they share the same name
#' as any of the defined parameters \code{vivax_model_param_list_create} will stop. You can either write
#' any extra parameters you like individually, e.g. vivax_model_param_list_create(extra1 = 1, extra2 = 2)
#' and these parameters will appear appended to the returned list, or you can pass explicitly
#' the ellipsis argument as a list created before, e.g. vivax_model_param_list_create(...=list(extra1 = 1, extra2 = 2))
#'
#' @export
vivax_model_param_list_create <- function(
  # Parameter values match those in IBM code on Github
  # * indicates the parameter values in Supplementary Table 2 in https://www.nature.com/articles/s41467-018-05860-8

  # Age & heterogeneity in exposure
  eta = 1/(21*365),      # 22.5 years in original IBM
  rho = 0.85,
  a0 = 2920,
  sigma2 = 1.292285,     # *1.53
  max_age = 100*365,     # 80 years in original IBM
  #  Rate of leaving infection states
  rLM = 1/10,
  rT = 1,
  rD = 1/5,
  rP = 1/28,
  dPCR_min = 10,
  dPCR_max = 70,         # *52.6
  # Human latent period  (time lag from asexual parasites does not exist)
  dE  = 10,
  # Human infectiousness to mosquitoes
  cD  = 0.8,
  cT  =  0.4,
  cPCR = 0.035,
  cLM = 0.1,
  # Mosquito-to-human infection probability
  b = 0.5,               # *0.25
  # Anti-parasite immunity reducing probability of detection
  d_par = 10*365,
  u_par = 19.77,         # *42.4
  phi_LM_max = 0.8918,   # *0.93
  phi_LM_min = 0.0043,   # *0.011
  A_LM_50pc = 27.52,     # *18.8
  K_LM = 2.403,          # *3.37
  A_PCR_50pc = 9.9,
  K_PCR = 4.602,         # *3.84
  # Clinical immunity reducing probability of clinical disease
  d_clin = 30*365,
  u_clin = 7.85,         # *4.33
  phi_D_max = 0.8605,    # *0.96
  phi_D_min = 0.018,     # *0.006
  A_D_50pc = 11.538,     # *24.5
  K_D = 2.250,           # *5.63
  # Maternal immunity
  p_mat = 0.421,         # *0.31
  d_mat = 35.148,        # *49.9
  # Hypnozoite parameters
  gamma_L = 1/383,
  ff = 0.02439024,
  # Entomological parameters (Anopheles punctulatus)
  delayMos = 8.4,
  tau1 = 0.69,
  tau2 = 2.31,
  mu0 = 1/6,
  Q0 = 0.5,
  chi = 0.9,
  bites_Bed = 0.9,
  bites_Indoors = 0.9,
  # Larval parameters daily density dependent mortality rate of egg
  muEL = 0.0338,
  muLL = 0.0348,
  muPL = 0.249,
  dEL = 6.64,
  dLL = 3.72,
  dPL = 0.643,
  gammaL = 13.25,
  betaL = 21.19,
  # intervention parameters
  num_int = 1,
  itn_cov = 0,
  irs_cov = 0,
  ITN_IRS_on = -1,
  DY = 365,
  d_ITN0 = 0.3,
  r_ITN0 = 0.6,
  r_ITN1 = 0.2,
  r_IRS0 = 0.4,
  d_IRS0 = 0.5,
  irs_half_life =   0.5 * DY,
  itn_half_life =   2.64 * DY,  # 2.5 years for pyrethroid decay in originalIBM
  IRS_interval =   1 * DY,
  ITN_interval =   3 * DY,
  # Control parameters
  stop_rebound_switch = 0,
  hypnozoite_prev_threshold = 1e-05,
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
  mp_list$rLM <- rLM
  mp_list$rT <- rT
  mp_list$rD <- rD
  mp_list$rP <- rP
  mp_list$dPCR_min <- dPCR_min
  mp_list$dPCR_max <- dPCR_max

  # human latent period
  mp_list$dE <- dE

  # Human infectiousness to mosquitoes
  mp_list$cD <- cD
  mp_list$cT <- cT
  mp_list$cPCR <- cPCR
  mp_list$cLM <- cLM

  # Mosquito-to-human infection probability
  mp_list$b <- b

  # Anti-parasite immunity reducing probability of detection
  mp_list$d_par <- d_par
  mp_list$u_par <- u_par
  mp_list$phi_LM_max <- phi_LM_max
  mp_list$phi_LM_min <- phi_LM_min
  mp_list$A_LM_50pc <- A_LM_50pc
  mp_list$K_LM <- K_LM
  mp_list$A_PCR_50pc <- A_PCR_50pc
  mp_list$K_PCR <-K_PCR

  # Clinical immunity reducing probability of clinical disease
  mp_list$d_clin <- d_clin
  mp_list$u_clin <- u_clin
  mp_list$phi_D_max <- phi_D_max
  mp_list$phi_D_min <- phi_D_min
  mp_list$A_D_50pc <- A_D_50pc
  mp_list$K_D <- K_D

  # Maternal immunity
  mp_list$p_mat <- p_mat
  mp_list$d_mat <- d_mat

  # Hypnozoite parameters
  mp_list$gamma_L <- gamma_L
  mp_list$ff<- ff

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

  # Check that ITN and IRS coverage are not 100%
  if (itn_cov == 1 | irs_cov == 1){
    stop(message("Coverage cannot be 1. Please set to 0.9999 or similar."))
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

  # Control parameters
  mp_list$stop_rebound_switch <- stop_rebound_switch
  mp_list$hypnozoite_prev_threshold <- hypnozoite_prev_threshold

  # check that none of the spare parameters in the extra
  if(sum(!is.na(match(names(extra_param_list),names(mp_list))))!=0){

    stop (message(cat("Extra params in ... share names with default param names. Please check:\n",
                      names(extra_param_list)[!is.na(match(names(extra_param_list),names(mp_list)))]
    )
    ))
  }

  return(append(mp_list,extra_param_list))
}


