#---------------------------------------------------------
#' Equilibrium initialisation list creation: P. falciparum
#'
#' The function \code{equilibrium_init_create} creates an equilibrium initialisation state to be
#' used within later model runs
#'
#' @param age_vector Vector of age brackets.
#' @param het_brackets Integer number of biting heteogenity compartments.
#' @param country String for country of interest. If NULL the seasonal parameters
#' will attempt to be loaded using just the admin unit, however if there is ambiguity
#' in the admin unit an error will be thrown. If both NULL then no seasonality is
#' assumed. Default = NULL.
#' @param admin_unit String for admin unit with country for loading seasonal
#' parameters. If country is NULL, the admin unit will attempt to be located,however
#' if there is ambiguity in the admin unit an error will be thrown. If both country
#' and admin_unit are NULL then no seasonality is assumed. Default = NULL.
#' @param ft Numeric for the frequency of people seeking treatment.
#' @param EIR Numeric for desired annual EIR.
#' @param model_param_list List of epidemiological parameters created by
#'
#' @importFrom stringi stri_trans_general
#' @importFrom statmod gauss.quad.prob
#'
#'
#' @export

equilibrium_init_create <- function(age_vector, het_brackets,
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
  num_int <- mpl$num_int

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

  age02 <- which(age_vector > 2)[1] - 1   # index of age vector at age 2 years
  age10 <- which(age_vector >= 10)[1] - 1  # index of age vector before age is 10 years


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
      FOIvij_eq[i, j] <- foi_age[i] * mpl$av0 * (mpl$cT * T_eq[i, j] + mpl$cD *
                                                   D_eq[i, j] + cA_eq[i, j] * A_eq[i, j] + mpl$cU * U_eq[i, j]) * rel_foi[j]/omega
    }
  }

  # mosquito states
  FOIv_eq <- sum(FOIvij_eq)
  Iv_eq <- FOIv_eq * mpl$Surv0/(FOIv_eq + mpl$mu0)
  Sv_eq <- mpl$mu0 * Iv_eq/(FOIv_eq * mpl$Surv0)
  Ev_eq <- 1 - Sv_eq - Iv_eq

  # mosquito density needed to give this EIR
  mv0 <- omega * EIRd_eq/(Iv_eq * mpl$av0)

  # larval states
  K0 <- 2 * mv0 * mpl$dLL * mpl$mu0 * (1 + mpl$dPL * mpl$muPL) * mpl$gammaL * (mpl$lambda + 1)/(mpl$lambda/(mpl$muLL *
                                                                                                              mpl$dEL) - 1/(mpl$muLL * mpl$dLL) - 1)
  PL_eq <- 2 * mpl$dPL * mpl$mu0 * mv0
  LL_eq <- mpl$dLL * (mpl$muPL + 1/mpl$dPL) * PL_eq
  EL_eq <- (LL_eq/mpl$dLL + mpl$muLL* LL_eq * (1 + mpl$gammaL * LL_eq/K0))/(1/mpl$dEL - mpl$muLL * mpl$gammaL * LL_eq/K0)

  # add in final dimension - interventions
  num_int <- mpl$num_int
  cov <- mpl$cov

  mat <- matrix(0, na, nh)

  S_eq <- vapply(cov, FUN = function(x)
  {
    x * S_eq
  }, mat)
  T_eq <- vapply(cov, FUN = function(x)
  {
    x * T_eq
  }, mat)
  D_eq <- vapply(cov, FUN = function(x)
  {
    x * D_eq
  }, mat)
  A_eq <- vapply(cov, FUN = function(x)
  {
    x * A_eq
  }, mat)
  U_eq <- vapply(cov, FUN = function(x)
  {
    x * U_eq
  }, mat)
  P_eq <- vapply(cov, FUN = function(x)
  {
    x * P_eq
  }, mat)

  IB_eq = array(IB_eq, c(na, nh, num_int))
  ID_eq = array(ID_eq, c(na, nh, num_int))
  ICA_eq = array(ICA_eq, c(na, nh, num_int))
  ICM_eq = array(ICM_eq, c(na, nh, num_int))

  # TODO: Remove this part and put it as an edit to the equilibrium solution
  if(!is.null(mpl$ncc)){
    IB_eq = array(IB_eq, c(na, nh, num_int, mpl$ncc))
    ID_eq = array(ID_eq, c(na, nh, num_int, mpl$ncc))
    ICA_eq = array(ICA_eq, c(na, nh, num_int, mpl$ncc))
    ICM_eq = array(ICM_eq, c(na, nh, num_int, mpl$ncc))

    # add in final dimension - interventions
    all_rounds = mpl$MDA_grp_prop*mpl$MDA_cov
    ccov = c(all_rounds, 1-all_rounds)

    mat2 <- array(0, c(na,nh, num_int))
    S_eq <- vapply(ccov,FUN = function(x){x * S_eq},mat2)
    T_eq <- vapply(ccov,FUN = function(x){x * T_eq},mat2)
    D_eq <- vapply(ccov,FUN = function(x){x * D_eq},mat2)
    A_eq <- vapply(ccov,FUN = function(x){x * A_eq},mat2)
    U_eq <- vapply(ccov,FUN = function(x){x * U_eq},mat2)
    P_eq <- vapply(ccov,FUN = function(x){x * P_eq},mat2)
  }


  admin_units_seasonal <- load_file("admin_units_seasonal.rds")
  admin_matches <- admin_match(admin_unit = admin_unit, country = country,
                               admin_units_seasonal = admin_units_seasonal)

  if(admin_matches == 0){
    ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
  } else {
    ssa0 <- admin_units_seasonal$a0[admin_matches]
    ssa1 <- admin_units_seasonal$a1[admin_matches]
    ssa2 <- admin_units_seasonal$a2[admin_matches]
    ssa3 <- admin_units_seasonal$a3[admin_matches]
    ssb1 <- admin_units_seasonal$b1[admin_matches]
    ssb2 <- admin_units_seasonal$b2[admin_matches]
    ssb3 <- admin_units_seasonal$b3[admin_matches]
    theta_c <- admin_units_seasonal$theta_c[admin_matches]
  }

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
              init_ICM = ICM_eq, ICM_init_eq = ICM_init_eq, init_Iv = Iv_eq, init_Sv = Sv_eq,
              init_Ev = Ev_eq, init_PL = PL_eq, init_LL = LL_eq, init_EL = EL_eq,
              age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
              omega = omega, foi_age = foi_age, rel_foi = rel_foi,
              K0 = K0, mv0 = mv0, na = na, nh = nh, ni = num_int, x_I = x_I,
              FOI = FOI_eq, EIR_eq = EIR_eq, cA_eq = cA_eq,
              den = den, age59 = age59, age05 = age05, age02 = age02, age10 = age10,
              ssa0 = ssa0, ssa1 = ssa1,
              ssa2 = ssa2, ssa3 = ssa3, ssb1 = ssb1, ssb2 = ssb2, ssb3 = ssb3,
              theta_c = theta_c, age = age_vector*mpl$DY, ft = ft, FOIv_eq = FOIv_eq,
              betaS = betaS, betaA = betaA, betaU = betaU, FOIvij_eq=FOIvij_eq,
              age_mid_point = age_mid_point, het_bounds = het_bounds, pi = pi,
              age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)

  res <- append(res,mpl)

  return(res)
}

#' Equilibrium initialisation list creation: P. vivax
#'
#' The function \code{vivax_equilibrium_init_create} creates an equilibrium initialisation state to be
#' used within later model runs
#'
#' @param age_vector Vector of age brackets.
#' @param het_brackets Integer number of biting heteogenity compartments.
#' @param country String for country of interest. If NULL the seasonal parameters
#' will attempt to be loaded using just the admin unit, however if there is ambiguity
#' in the admin unit an error will be thrown. If both NULL then no seasonality is
#' assumed. Default = NULL.
#' @param admin_unit String for admin unit with country for loading seasonal
#' parameters. If country is NULL, the admin unit will attempt to be located,however
#' if there is ambiguity in the admin unit an error will be thrown. If both country
#' and admin_unit are NULL then no seasonality is assumed. Default = NULL.
#' @param ft Numeric for the frequency of people seeking treatment.
#' @param EIR Numeric for desired annual EIR.
#' @param K_max Maximum number of hypnozoite batches. Has to be an integer between
#' 1 and 30. Note a model with stratification into larger number of hypnozoite batches
#' (e.g. K_max > 5) may have to be run on the cluster. For an EIR of 1 or less,
#' K_max = 2 is sufficient.
#' @param model_param_list List of epidemiological parameters created by
#'
#' @importFrom stringi stri_trans_general
#' @importFrom statmod gauss.quad.prob
#'
#'
#' @export

vivax_equilibrium_init_create <- function(age_vector, het_brackets,
                                          country = NULL, admin_unit = NULL, ft,
                                          EIR, K_max, model_param_list)
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
  if(!is.numeric(K_max)) stop("K_max is not provided or not numeric")
  if(!(K_max %in% c(1:30))) stop("K_max is not an integer between 1 and 30")

  ## Handle parameters
  # database for admin units is all in Latin-ASCII for CRAN reasons so must
  # encode parameters accordingly
  if(!is.null(country)) country <- stringi::stri_trans_general(country,"Latin-ASCII")
  if(!is.null(admin_unit)) admin_unit <- stringi::stri_trans_general(admin_unit, "Latin-ASCII")

  ## Population demographics
  age <- age_vector * mpl$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets) # number of heterogeneity groups
  nk <- as.integer(K_max+1)      # number of hypnozoite batch groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")
  age0 <- 2
  age1 <- 10
  num_int <- mpl$num_int

  ## Hypnozoite inputs
  K_max_switch <- c(rep(1,nk-1),0)        # switch to set some transitions to 0 when k=K_max
  K_max_switch_on <- c(rep(0,nk-1),1)     # switch to only apply some transitions when k=K_max
  K0_switch <- c(0, rep(1,nk-1))          # switch to set some transitions to 0 when k=0
  kk_vec <- 0:K_max                       # vector of hypnozoite batches

  age_rate <- age_width <- age_mid_point <- den <- c()
  for (i in 1:(na-1))
  {
    age_width[i] <- age[i+1] - age[i]
    age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group

  }
  age_rate[na] <- 0
  age_mid_point[na] <- 0.5 * (age[na] + mpl$max_age)


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
  ## NOTE: in vivax package foi_age uses age_mid_point instead of age
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


  # Find equilibrium FOI
  FOI_eq <- matrix(NA, nrow=na, ncol=nh)

  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      FOI_eq[i, j] <- EIR_eq[i, j] * mpl$b
    }
  }

  ## Equilibrium state values for humans

  # Isolate demographic component in common to all human equations:
  gamma <- mpl$eta + c(age_rate[1:(na - 1)], 0)   # leaving current age group (death + aging out)
  delta <- c(mpl$eta, age_rate[1:(na - 1)])       # entering current age group (aging in)
  # At age 0, rate of aging in = mortality rate, to replace population

  # All equilibrium state values are calculated using matrix inversion (solve function)
  # Notation: A * X + b = 0
  # where A is a matrix of transitions in differential equations

  ## Equilibrium distribution of hypnozoite batches
  Hyp_eq <- array(NA, dim=c(na,nh,nk))

  b_Hyp <- rep(list(array(0,dim=c(1,nk,nh))), na)

  # Create matrix A for each age and heterogeneity group
  A_Hyp <- rep(list(array(0,dim=c(nk,nk,nh))), na)

  for(j in 1:nh) {
    for(i in 1:na) {
      for (kk in 1:(nk-1)) {
        A_Hyp[[i]][kk,kk,j] <- -FOI_eq[i,j]-kk_vec[kk]*mpl$gamma_L-gamma[i]
        A_Hyp[[i]][kk,kk+1,j] <- kk_vec[kk+1]*mpl$gamma_L
        A_Hyp[[i]][kk+1,kk,j] <- FOI_eq[i,j]
      }
      A_Hyp[[i]][nk,nk,j] <- -K_max*mpl$gamma_L-gamma[i]

    }
  }

  # Calculate equilibrium in each hypnozoite state at age 0:
  for(j in 1:nh) {
    b_Hyp[[1]][,,j] <-  c(-mpl$eta*het_wt[j],rep(0,nk-1))

    # Solve matrices:
    Hyp_eq[1,j,] <- solve(A_Hyp[[1]][,,j], b_Hyp[[1]][,,j])
  }

  # For all other age groups:
  for(j in 1:nh) {
    for(i in 2:na) {
      b_Hyp[[i]][,,j] <-  c(-delta[i]*Hyp_eq[i-1,j,])

      # Solve matrices:
      Hyp_eq[i,j,] <- solve(A_Hyp[[i]][,,j], b_Hyp[[i]][,,j])

    }
  }

  # Calculate proportion in each hypnozoite state for each age and heterogeneity group
  hyp_wt <- array(NA, dim=c(na,nh,nk))
  for (j in 1:nh) {
    for (i in 1:na) {
      for (k in 1:nk) {
        hyp_wt[i,j,k] <- Hyp_eq[i,j,k]/sum(Hyp_eq[i,j,])
      }
    }
  }

  ## Immunity states
  AP_eq <- array(NA, dim=c(na,nh,nk))
  AC_eq <- array(NA, dim=c(na,nh,nk))
  AP_MAT_init_eq <- vector(length = nh, mode = "numeric")
  AP_MAT_eq <- matrix(0, na, nh)
  AC_MAT_init_eq <- vector(length = nh, mode = "numeric")
  AC_MAT_eq <- matrix(0, na, nh)

  b_AP <- list()
  b_AC <- list()

  # Create matrix A for each age and heterogeneity group for anti-parasite (AP) and anti-clinical immunity (AC)
  A_AP <- rep(list(array(0,dim=c(nk,nk,nh))), na)
  A_AC <- rep(list(array(0,dim=c(nk,nk,nh))), na)

  for(j in 1:nh) {
    for(i in 1:na) {
      for (kk in 1:(nk-1)) {
        A_AP[[i]][kk,kk,j] <- -FOI_eq[i,j]-1/mpl$d_par-(kk-1)*mpl$gamma_L-gamma[i]
        A_AP[[i]][kk,kk+1,j] <- mpl$gamma_L*kk * Hyp_eq[i,j,kk+1]/Hyp_eq[i,j,kk]
        A_AP[[i]][kk+1,kk,j] <- FOI_eq[i,j] * Hyp_eq[i,j,kk]/Hyp_eq[i,j,kk+1]

        A_AC[[i]][kk,kk,j] <- -FOI_eq[i,j]-1/mpl$d_clin-(kk-1)*mpl$gamma_L-gamma[i]
        A_AC[[i]][kk,kk+1,j] <- mpl$gamma_L*kk * Hyp_eq[i,j,kk+1]/Hyp_eq[i,j,kk]
        A_AC[[i]][kk+1,kk,j] <- FOI_eq[i,j] * Hyp_eq[i,j,kk]/Hyp_eq[i,j,kk+1]
      }
      A_AP[[i]][nk,nk,j] <- -1/mpl$d_par-K_max*mpl$gamma_L-gamma[i]
      A_AC[[i]][nk,nk,j] <- -1/mpl$d_clin-K_max*mpl$gamma_L-gamma[i]

    }
  }

  # Calculate equilibrium in each AP and AC hypnozoite state at age 0:
  b_AP <- rep(list(array(0,dim=c(1,nk,nh))), na)
  b_AC <- rep(list(array(0,dim=c(1,nk,nh))), na)

  for(j in 1:nh) {
    for(kk in 2:nk) {
      b_AP[[1]][,kk,j] <- -(FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*mpl$ff)/
        ((FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*mpl$ff)*mpl$u_par+1)

      b_AC[[1]][,kk,j] <- -(FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*mpl$ff)/
        ((FOI_eq[1,j]*Hyp_eq[1,j,kk-1]/Hyp_eq[1,j,kk] + kk_vec[kk]*mpl$ff)*mpl$u_clin+1)

    }

    # Solve matrices:
    AP_eq[1,j,] <- solve(A_AP[[1]][,,j], b_AP[[1]][,,j])
    AC_eq[1,j,] <- solve(A_AC[[1]][,,j], b_AC[[1]][,,j])
  }

  # For all other age groups:
  for(j in 1:nh) {
    for(i in 2:na) {

      b_AP[[i]][,1,j] <- -delta[i]*AP_eq[i-1,j,1]*Hyp_eq[i-1,j,1]/Hyp_eq[i,j,1]
      b_AC[[i]][,1,j] <- -delta[i]*AC_eq[i-1,j,1]*Hyp_eq[i-1,j,1]/Hyp_eq[i,j,1]

      for(kk in 2:nk) {
        b_AP[[i]][,kk,j] <- -(FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*mpl$ff)/
          ((FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*mpl$ff)*mpl$u_par+1)-
          delta[i]*AP_eq[i-1,j,kk]*Hyp_eq[i-1,j,kk]/Hyp_eq[i,j,kk]

        b_AC[[i]][,kk,j] <- -(FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*mpl$ff)/
          ((FOI_eq[i,j]*Hyp_eq[i,j,kk-1]/Hyp_eq[i,j,kk] + kk_vec[kk]*mpl$ff)*mpl$u_clin+1)-
          delta[i]*AC_eq[i-1,j,kk]*Hyp_eq[i-1,j,kk]/Hyp_eq[i,j,kk]
      }

      # Solve matrices:
      AP_eq[i,j,] <- solve(A_AP[[i]][,,j], b_AP[[i]][,,j])
      AC_eq[i,j,] <- solve(A_AC[[i]][,,j], b_AC[[i]][,,j])

    }
  }

  # Maternal immunity needs to be calculated after because it references AP and AC

  # Immunity in 20 year old woman is averaged over hypnozoite states
  # (weighted by hypnozoite distribution)

  AP_eq_mean <- apply(AP_eq*hyp_wt, c(1,2), sum)
  AC_eq_mean <- apply(AC_eq*hyp_wt, c(1,2), sum)

  for (j in 1:nh)
  {
    for (i in 1:na)
    {

      AP_MAT_init_eq[j] <- mpl$p_mat * (AP_eq_mean[age20l, j] + age_20_factor *
                                          (AP_eq_mean[age20u, j] - AP_eq_mean[age20l, j]))
      AC_MAT_init_eq[j] <- mpl$p_mat * (AC_eq_mean[age20l, j] + age_20_factor *
                                          (AC_eq_mean[age20u, j] - AC_eq_mean[age20l, j]))

      AP_MAT_eq[i, j] <- AP_MAT_init_eq[j]*exp(-age[i]/mpl$d_mat)
      AC_MAT_eq[i, j] <- AC_MAT_init_eq[j]*exp(-age[i]/mpl$d_mat)
      # Note IBM uses age_mid_point

    }
  }


  # Turn AP_MAT_eq and AC_MAT_eq into an array (with same values for each hypnozoite state)
  # for calculation of phi_LM, phi_D and dPCR
  AP_MAT_eq <- array(AP_MAT_eq, dim=c(na,nh,(K_max+1)))
  AC_MAT_eq <- array(AC_MAT_eq, dim=c(na,nh,(K_max+1)))

  # Calculate probabilities
  phi_LM_eq <- mpl$phi_LM_min + (mpl$phi_LM_max-mpl$phi_LM_min) *
    1/(1+((AP_eq+AP_MAT_eq)/mpl$A_LM_50pc)^mpl$K_LM)
  phi_D_eq <- mpl$phi_D_min + (mpl$phi_D_max-mpl$phi_D_min) *
    1/(1+((AC_eq+AC_MAT_eq)/mpl$A_D_50pc)^mpl$K_D)
  dPCR_eq <- mpl$dPCR_min + (mpl$dPCR_max-mpl$dPCR_min) *
    1/(1+((AP_eq+AP_MAT_eq)/mpl$A_PCR_50pc)^mpl$K_PCR)

  ## Human states (stored in Z)
  n <- 6*nk       # matrix dimensions = epidemiological compartments x hypnozoite groups
  Z <- array(NA, dim=c(na,nh,n))
  b <- list()

  # Define indices:
  iS <- 1
  iI_PCR <- 2
  iI_LM <- 3
  iI_D <- 4
  iT <- 5
  iP <- 6

  # Create matrix A for each age and heterogeneity group
  A <- array(0,dim=c(n,n,nh))
  A <- rep(list(A), na)

  for(j in 1:nh) {
    for(i in 1:na) {

      # Fill in non-zero elements:

      # Transitions for k<K_max
      for (kk in 0:(K_max-1)) {

        # Diagonals
        A[[i]][iS+kk*6,iS+kk*6,j] <- -FOI_eq[i,j]-kk*mpl$ff-kk*mpl$gamma_L-gamma[i]

        A[[i]][iI_PCR+kk*6,iI_PCR+kk*6,j] <- -FOI_eq[i,j]-kk*mpl$ff-1/dPCR_eq[i,j,(kk+1)]-
          kk*mpl$gamma_L-gamma[i]+kk*mpl$ff*(1-phi_LM_eq[i,j,(kk+1)])

        A[[i]][iI_LM+kk*6,iI_LM+kk*6,j] <- -FOI_eq[i,j]-kk*mpl$ff-mpl$rLM-
          kk*mpl$gamma_L-gamma[i]+kk*mpl$ff*(1-phi_D_eq[i,j,(kk+1)])

        A[[i]][iI_D+kk*6,iI_D+kk*6,j] <- -FOI_eq[i,j]-mpl$rD-kk*mpl$gamma_L-gamma[i]

        A[[i]][iT+kk*6,iT+kk*6,j] <- -FOI_eq[i,j]-mpl$rT-kk*mpl$gamma_L-gamma[i]

        A[[i]][iP+kk*6,iP+kk*6,j] <- -FOI_eq[i,j]-mpl$rP-kk*mpl$gamma_L-gamma[i]

        # Additions from gamma_L(k+1):
        A[[i]][iS+kk*6,iS+kk*6+6,j] <- (kk+1)*mpl$gamma_L
        A[[i]][iI_PCR+kk*6,iI_PCR+kk*6+6,j] <- (kk+1)*mpl$gamma_L
        A[[i]][iI_LM+kk*6,iI_LM+kk*6+6,j] <- (kk+1)*mpl$gamma_L
        A[[i]][iI_D+kk*6,iI_D+kk*6+6,j] <- (kk+1)*mpl$gamma_L
        A[[i]][iT+kk*6,iT+kk*6+6,j] <- (kk+1)*mpl$gamma_L
        A[[i]][iP+kk*6,iP+kk*6+6,j] <- (kk+1)*mpl$gamma_L

        # Other transitions:
        A[[i]][iI_PCR+kk*6,iS+kk*6,j] <- kk*mpl$ff*(1-phi_LM_eq[i,j,(kk+1)])
        A[[i]][iI_LM+kk*6,iS+kk*6,j] <- kk*mpl$ff*(1-phi_D_eq[i,j,(kk+1)])*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_LM+kk*6,iI_PCR+kk*6,j] <- kk*mpl$ff*(1-phi_D_eq[i,j,(kk+1)])*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_D+kk*6,iS+kk*6,j] <- kk*mpl$ff*phi_D_eq[i,j,(kk+1)]*(1-ft)*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_D+kk*6,iI_PCR+kk*6,j] <- kk*mpl$ff*phi_D_eq[i,j,(kk+1)]*(1-ft)*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iI_D+kk*6,iI_LM+kk*6,j] <- kk*mpl$ff*phi_D_eq[i,j,(kk+1)]*(1-ft)
        A[[i]][iT+kk*6,iS+kk*6,j] <- kk*mpl$ff*phi_D_eq[i,j,(kk+1)]*ft*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iT+kk*6,iI_PCR+kk*6,j] <- kk*mpl$ff*phi_D_eq[i,j,(kk+1)]*ft*phi_LM_eq[i,j,(kk+1)]
        A[[i]][iT+kk*6,iI_LM+kk*6,j] <- kk*mpl$ff*phi_D_eq[i,j,(kk+1)]*ft

      }

      # Diagonals for K_max:
      A[[i]][iS+K_max*6,iS+K_max*6,j] <- -FOI_eq[i,j]-K_max*mpl$ff-K_max*mpl$gamma_L-gamma[i]
      A[[i]][iI_PCR+K_max*6,iI_PCR+K_max*6,j] <- -FOI_eq[i,j] + FOI_eq[i,j]*(1-phi_LM_eq[i,j,(K_max+1)]) -
        K_max*mpl$ff-1/dPCR_eq[i,j,(K_max+1)]-
        K_max*mpl$gamma_L-gamma[i]+K_max*mpl$ff*(1-phi_LM_eq[i,j,(K_max+1)])
      A[[i]][iI_LM+K_max*6,iI_LM+K_max*6,j] <- -FOI_eq[i,j] + FOI_eq[i,j]*(1-phi_D_eq[i,j,(K_max+1)])-
        K_max*mpl$ff-mpl$rLM-
        K_max*mpl$gamma_L-gamma[i]+K_max*mpl$ff*(1-phi_D_eq[i,j,(K_max+1)])
      A[[i]][iI_D+K_max*6,iI_D+K_max*6,j] <- -mpl$rD-K_max*mpl$gamma_L-gamma[i]
      A[[i]][iT+K_max*6,iT+K_max*6,j] <- -mpl$rT-K_max*mpl$gamma_L-gamma[i]
      A[[i]][iP+K_max*6,iP+K_max*6,j] <- -mpl$rP-K_max*mpl$gamma_L-gamma[i]

      # Other transitions for K_max:
      A[[i]][iI_PCR+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*(1-phi_LM_eq[i,j,(K_max+1)])+
        K_max*mpl$ff*(1-phi_LM_eq[i,j,(K_max+1)])
      A[[i]][iI_LM+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*phi_LM_eq[i,j,(K_max+1)]*(1-phi_D_eq[i,j,(K_max+1)])+
        K_max*mpl$ff*(1-phi_D_eq[i,j,(K_max+1)])*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_LM+K_max*6,iI_PCR+K_max*6,j] <- FOI_eq[i,j]*phi_LM_eq[i,j,(K_max+1)]*(1-phi_D_eq[i,j,(K_max+1)])+
        K_max*mpl$ff*(1-phi_D_eq[i,j,(K_max+1)])*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_D+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)] +
        K_max*mpl$ff*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_D+K_max*6,iI_PCR+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)] +
        K_max*mpl$ff*phi_D_eq[i,j,(K_max+1)]*(1-ft)*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iI_D+K_max*6,iI_LM+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*(1-ft)+
        K_max*mpl$ff*phi_D_eq[i,j,(K_max+1)]*(1-ft)
      A[[i]][iT+K_max*6,iS+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)] +
        K_max*mpl$ff*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iT+K_max*6,iI_PCR+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)] +
        K_max*mpl$ff*phi_D_eq[i,j,(K_max+1)]*ft*phi_LM_eq[i,j,(K_max+1)]
      A[[i]][iT+K_max*6,iI_LM+K_max*6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,(K_max+1)]*ft +
        K_max*mpl$ff*phi_D_eq[i,j,(K_max+1)]*ft

      # Transitions for all k:
      for (kk in c(0:K_max)) {

        A[[i]][iS+kk*6,iI_PCR+kk*6,j] <- 1/dPCR_eq[i,j,(kk+1)]
        A[[i]][iS+kk*6,iP+kk*6,j] <- mpl$rP
        A[[i]][iI_PCR+kk*6,iI_LM+kk*6,j] <- mpl$rLM
        A[[i]][iI_LM+kk*6,iI_D+kk*6,j] <- mpl$rD
        A[[i]][iP+kk*6,iT+kk*6,j] <- mpl$rT

      }


      # Transitions for k>0 (all except S)
      for (kk in c(1:K_max)) {
        A[[i]][iI_PCR+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_LM_eq[i,j,kk])
        A[[i]][iI_PCR+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_LM_eq[i,j,kk])
        A[[i]][iI_LM+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_D_eq[i,j,kk])*phi_LM_eq[i,j,kk]
        A[[i]][iI_LM+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_D_eq[i,j,kk])*phi_LM_eq[i,j,kk]
        A[[i]][iI_LM+kk*6,iI_LM+kk*6-6,j] <- FOI_eq[i,j]*(1-phi_D_eq[i,j,kk])
        A[[i]][iI_D+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*(1-ft)*phi_LM_eq[i,j,kk]
        A[[i]][iI_D+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*(1-ft)*phi_LM_eq[i,j,kk]
        A[[i]][iI_D+kk*6,iI_LM+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*(1-ft)
        A[[i]][iI_D+kk*6,iI_D+kk*6-6,j] <- FOI_eq[i,j]
        A[[i]][iT+kk*6,iS+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*ft*phi_LM_eq[i,j,kk]
        A[[i]][iT+kk*6,iI_PCR+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*ft*phi_LM_eq[i,j,kk]
        A[[i]][iT+kk*6,iI_LM+kk*6-6,j] <- FOI_eq[i,j]*phi_D_eq[i,j,kk]*ft
        A[[i]][iT+kk*6,iT+kk*6-6,j] <- FOI_eq[i,j]
        A[[i]][iP+kk*6,iP+kk*6-6,j] <- FOI_eq[i,j]
      }

    }
  }

  # Calculate equilibrium population in each compartment at age 0:
  for(j in 1:nh) {
    b[[1]]  <-  array(c(-mpl$eta*het_wt[j],rep(0,n-1)), dim=c(1,n,nh))
    # Set dS0/dt at age 0 to births, all other differential equations to 0 (for equilibrium values)

    # Solve matrix:
    Z[1,j,] <- solve(A[[1]][,,j], b[[1]][,,j])
  }

  # For all other age groups:
  for(j in 1:nh) {
    for(i in 2:na) {
      b[[i]] <-  array(c(-delta[i]*Z[i-1,j,1:n]), dim=c(1,n,nh))
      # Set all differential equations to the proportion of the population aging in

      # Solve matrix:
      Z[i,j,] <- solve(A[[i]][,,j], b[[i]][,,j])
    }
  }

  # Fill in equilibrium states
  T_eq <- Z[,,iT+kk_vec*6]
  P_eq <- Z[,,iP+kk_vec*6]
  S_eq <- Z[,,iS+kk_vec*6]
  I_PCR_eq <- Z[,,iI_PCR+kk_vec*6]
  I_LM_eq <- Z[,,iI_LM+kk_vec*6]
  I_D_eq <- Z[,,iI_D+kk_vec*6]

  # Vector model
  FOIvij_eq <- array(dim=c(na,nh,(K_max+1)))
  for (kk in 1:(K_max+1))
  {
    for (j in 1:nh)
    {
      for (i in 1:na)
      {
        FOIvij_eq[i, j, kk] <-  mpl$av0 * foi_age[i] *
          rel_foi[j]/omega * (mpl$cT * T_eq[i, j,kk] +
                                mpl$cD *I_D_eq[i, j,kk] + mpl$cLM * I_LM_eq[i, j,kk] +
                                mpl$cPCR * I_PCR_eq[i, j,kk])

      }
    }
  }

  # Mosquito states
  FOIv_eq <- sum(FOIvij_eq)
  Iv_eq <- FOIv_eq * mpl$Surv0/(FOIv_eq + mpl$mu0)
  Sv_eq <- mpl$mu0 * Iv_eq/(FOIv_eq * mpl$Surv0)
  Ev_eq <- 1 - Sv_eq - Iv_eq

  # Mosquito density needed to give this EIR
  mv0 <- omega * EIRd_eq/(Iv_eq * mpl$av0)

  # Larval states
  K0 <- 2 * mv0 * mpl$dLL * mpl$mu0 * (1 + mpl$dPL * mpl$muPL) * mpl$gammaL * (mpl$lambda + 1)/(mpl$lambda/(mpl$muLL *
                                                                                                              mpl$dEL) - 1/(mpl$muLL * mpl$dLL) - 1)
  PL_eq <- 2 * mpl$dPL * mpl$mu0 * mv0
  LL_eq <- mpl$dLL * (mpl$muPL + 1/mpl$dPL) * PL_eq
  EL_eq <- (LL_eq/mpl$dLL + mpl$muLL* LL_eq * (1 + mpl$gammaL * LL_eq/K0))/(1/mpl$dEL - mpl$muLL * mpl$gammaL * LL_eq/K0)

  # Add in final dimension - interventions
  num_int <- mpl$num_int
  cov <- mpl$cov

  mat <- array(0, dim=c(na,nh,nk))

  S_eq <- vapply(cov, FUN = function(x)
  {
    x * S_eq
  }, mat)
  T_eq <- vapply(cov, FUN = function(x)
  {
    x * T_eq
  }, mat)
  I_D_eq <- vapply(cov, FUN = function(x)
  {
    x * I_D_eq
  }, mat)
  I_LM_eq <- vapply(cov, FUN = function(x)
  {
    x * I_LM_eq
  }, mat)
  I_PCR_eq <- vapply(cov, FUN = function(x)
  {
    x * I_PCR_eq
  }, mat)
  P_eq <- vapply(cov, FUN = function(x)
  {
    x * P_eq
  }, mat)
  Hyp_eq <- vapply(cov, FUN = function(x)
  {
    x * Hyp_eq
  }, mat)

  AP_eq = array(AP_eq, c(na, nh, nk, num_int))
  AC_eq = array(AC_eq, c(na, nh, nk, num_int))
  AP_MAT_eq = array(AP_MAT_eq, c(na, nh, nk, num_int))
  AC_MAT_eq = array(AC_MAT_eq, c(na, nh, nk, num_int))

  # Change array structure to: [age, het groups, intervention groups, hypnozoites]
  S_eq <- aperm(S_eq, c(1,2,4,3))
  I_PCR_eq <- aperm(I_PCR_eq, c(1,2,4,3))
  I_LM_eq <- aperm(I_LM_eq, c(1,2,4,3))
  I_D_eq <- aperm(I_D_eq, c(1,2,4,3))
  T_eq <- aperm(T_eq, c(1,2,4,3))
  P_eq <- aperm(P_eq, c(1,2,4,3))
  Hyp_eq <- aperm(Hyp_eq, c(1,2,4,3))
  AP_eq <- aperm(AP_eq, c(1,2,4,3))
  AC_eq <- aperm(AC_eq, c(1,2,4,3))
  AP_MAT_eq <- aperm(AP_MAT_eq, c(1,2,4,3))
  AC_MAT_eq <- aperm(AC_MAT_eq, c(1,2,4,3))

  # Seasonality
  admin_units_seasonal <- load_file("admin_units_seasonal.rds")
  admin_matches <- admin_match(admin_unit = admin_unit, country = country,
                               admin_units_seasonal = admin_units_seasonal)

  if(admin_matches == 0){
    ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
  } else {
    ssa0 <- admin_units_seasonal$a0[admin_matches]
    ssa1 <- admin_units_seasonal$a1[admin_matches]
    ssa2 <- admin_units_seasonal$a2[admin_matches]
    ssa3 <- admin_units_seasonal$a3[admin_matches]
    ssb1 <- admin_units_seasonal$b1[admin_matches]
    ssb2 <- admin_units_seasonal$b2[admin_matches]
    ssb3 <- admin_units_seasonal$b3[admin_matches]
    theta_c <- admin_units_seasonal$theta_c[admin_matches]
  }

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
  state <- list(init_S = S_eq, init_T = T_eq, init_I_D = I_D_eq, init_I_LM = I_LM_eq,
                init_I_PCR = I_PCR_eq, init_P = P_eq, init_Hyp = Hyp_eq,
                init_AP = AP_eq, init_AC = AC_eq, init_AP_MAT = AP_MAT_eq,
                init_AC_MAT = AC_MAT_eq, init_Iv = Iv_eq, init_Sv = Sv_eq,
                init_Ev = Ev_eq, init_PL = PL_eq, init_LL = LL_eq, init_EL = EL_eq,
                age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
                omega = omega, foi_age = foi_age, rel_foi = rel_foi,
                K0 = K0, mv0 = mv0, na = na, nh = nh,
                nk=nk, kk=kk_vec, K_max = K_max, K_max_switch = K_max_switch,
                K0_switch = K0_switch, K_max_switch_on = K_max_switch_on,
                ni = num_int, hyp_wt = hyp_wt,
                FOI = FOI_eq, EIR_eq = EIR_eq, Hyp_eq = Hyp_eq,
                phi_LM_eq = phi_LM_eq, phi_D_eq = phi_D_eq, dPCR_eq = dPCR_eq,
                den = den, age59 = age59, age05 = age05, ssa0 = ssa0, ssa1 = ssa1,
                ssa2 = ssa2, ssa3 = ssa3, ssb1 = ssb1, ssb2 = ssb2, ssb3 = ssb3,
                theta_c = theta_c, age = age_vector*mpl$DY, ft = ft, FOIv_eq = FOIv_eq,
                FOIvij_eq=FOIvij_eq,
                age_mid_point = age_mid_point, het_bounds = het_bounds, pi = pi,
                age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)

  state <- append(state,mpl)

  return(state)

}
