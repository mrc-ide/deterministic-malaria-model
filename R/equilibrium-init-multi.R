#------------------------------------------------
#' Equilibrium initialisation for multiple populations
#'
#' \code{equilibrium_init_multi} creates an equilibrium initialisation state to be
#' used within later model runs
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
#' @param EIR_vector Vector of numeric for desired annual EIR per patch.
#' @param model_param_list List of epidemiological parameters created by
#'
#' @importFrom stringi stri_trans_general
#' @importFrom statmod gauss.quad.prob
#'
#'
#' @export
#'
#'
#'
equilibrium_init_multi <- function (age_vector, het_brackets,
                                    country = NULL, admin_unit = NULL, ft,
                                    EIR_vector, model_param_list){

  ## population demographics
  age <- age_vector * model_param_list$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  num_int <- model_param_list$num_int
  np <- model_param_list$np

  ## Create 4th dimenstion arrays for humans initial conditions
  S <-  array(0, c(na,nh,num_int,np))
  Tt <- array(0,c(na,nh,num_int,np))
  D <-  array(0, c(na,nh,num_int,np))
  A <-  array(0, c(na,nh,num_int,np))
  U <-  array(0, c(na,nh,num_int,np))
  P <-  array(0, c(na,nh,num_int,np))
  IB <-  array(0, c(na,nh,num_int,np))
  ID <-  array(0, c(na,nh,num_int,np))
  ICA <-  array(0, c(na,nh,num_int,np))
  ICM <-  array(0, c(na,nh,num_int,np))
  # Create 1 dimension array for mosquitos initial condition
  Iv <- vector(mode="numeric", length=np)
  Sv <- vector(mode="numeric", length=np)
  Ev <- vector(mode="numeric", length=np)
  PL <- vector(mode="numeric", length=np)
  EL <- vector(mode="numeric", length=np)
  LL <- vector(mode="numeric", length=np)


  # Parameters that now need an extra dimension
  mv0 <- vector(mode="numeric", length=np)
  EIR_eq <- array(0,c(na,nh,np))
  FOIv_eq <- vector(mode="numeric", length=np)
  FOIvij_eq <- array(0,c(na,nh,np))

  cov <- model_param_list$cov
  temp_param_list <- model_param_list
  # Adding the fourth dimmension
  for (i in 1:length(EIR_vector)){


    temp_cov <- cov[,i]
    temp_param_list$cov <- temp_cov

    init_single <- equilibrium_init_create(age_vector = age_vector, het_brackets=het_brackets,
                                           country = NULL, admin_unit = NULL, ft= ft,
                                           EIR=EIR_vector[i], model_param_list = temp_param_list)



    S[,,,i] <- init_single$init_S
    Tt[,,,i] <- init_single$init_T
    D[,,,i] <- init_single$init_D
    A[,,,i] <- init_single$init_A
    U[,,,i] <- init_single$init_U
    P[,,,i] <- init_single$init_P
    IB[,,,i] <- init_single$init_IB
    ID[,,,i] <- init_single$init_ID
    ICA[,,,i] <- init_single$init_ICA
    ICM[,,,i] <- init_single$init_ICM

    Iv[i] <- init_single$init_Iv
    Sv[i] <- init_single$init_Sv
    Ev[i] <- init_single$init_Ev
    PL[i] <- init_single$init_PL
    LL[i] <- init_single$init_LL
    EL[i] <- init_single$init_EL

    mv0[i] <- init_single$mv0
    EIR_eq[,,i] <- init_single$EIR_eq
    FOIv_eq[i] <- init_single$FOIv_eq
    FOIvij_eq[,,i] <- init_single$FOIvij_eq


  }

  init_single$init_S = S
  init_single$init_T = Tt
  init_single$init_D = D
  init_single$init_A = A
  init_single$init_U = U
  init_single$init_P = P

  init_single$init_IB = IB
  init_single$init_ID = ID
  init_single$init_ICA = ICA
  init_single$init_ICM = ICM

  init_single$init_Iv = Iv
  init_single$init_Sv = Sv
  init_single$init_Ev = Ev

  init_single$init_PL = PL
  init_single$init_LL = LL
  init_single$init_EL = EL


  init_single$mv0 = mv0
  init_single$EIR_eq = EIR_eq
  init_single$FOIv_eq = FOIv_eq
  init_single$FOIvij_eq = FOIvij_eq

  init_single$np = np
  init_single$cov = cov

  return(init_single)

}
