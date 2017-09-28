#' Pipe operator
#'
#' See \code{\link[magrittr]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#------------------------------------------------
#' Grab parameters from odin model
#'
#' \code{grab_user_params} returns vector of all user required variables from an odin
#' model path
#'
#' @param model_path Character path to odin model
#'
#' @return Vector of characters for required variables
#'
#' @examples
#' head(hanojoel:::grab_user_params(system.file("extdata/odin_model.R",package="hanojoel")))
#'

grab_user_params <- function(model_path){

  ## first read lines from model path
  lines <- readLines(model_path)

  ## remove all spaces as makes life easier
  lines <- gsub("[[:space:]]", "", lines)

  ## find all positions that have either <-user() or =user() and dont begin with #
  var_pos <- grep("^[^#].*<-user()|^[^#].*=user()",lines,fixed=FALSE)

  ## remove all characters including and after <-user() or =user() and then remove brackets etc.
  res <- gsub(pattern = "<-user().*$|=user().*$",replacement = "",x = lines)[var_pos] %>%
    sapply(function(x){return(gsub(pattern = "\\[.*\\]|\\(.*\\)",replacement="",x=x))}) %>%
    as.character()


  ## return answer
  return(res)

}


##########
#
# temp <- tempfile()
#
# writeLines(text = "generate_default_model <- function(ft,age,dat,generator,dde = TRUE){
#   mod <- generator(init_S=dat$S,
#                    init_T=dat$T,
#                    init_D=dat$D,
#                    init_A=dat$A,
#                    init_U=dat$U,
#                    init_P=dat$P,
#                    init_ICA = dat$ICA,
#                    init_ICM = dat$ICM,
#                    init_ID = dat$ID,
#                    init_IB = dat$IB,
#                    init_Sv = dat$Sv,
#                    init_Ev = dat$Ev,
#                    init_Iv = dat$Iv,
#                    init_PL = dat$PL,
#                    init_LL = dat$LL,
#                    init_EL = dat$EL,
#                    na = dat$na,
#                    nh = dat$nh,
#                    age_rate = dat$age_rate,
#                    foi_age = dat$foi_age,
#                    het_wt = dat$het_wt,
#                    rel_foi = dat$rel_foi,
#                    omega = dat$omega,
#                    pi = pi,
#                    x_I = dat$x_I,
#                    #eov = dat$eov,
#                    mv0 = dat$mv0,
#                    ssa0 = dat$ssa0,
#                    ssa1 = dat$ssa1,
#                    ssa2 = dat$ssa2,
#                    ssa3 = dat$ssa3,
#                    ssb1 = dat$ssb1,
#                    ssb2 = dat$ssb2,
#                    ssb3 = dat$ssb3,
#                    theta_c = dat$theta_c,
#                    den = dat$den,
#                    age59 = dat$age59,
#                    age05 = dat$age05,
#                    age = age*365,
#                    ft = ft,
#                    use_dde = dde,
#                    eta = dat$eta,
#                    rA = dat$rA,
#                    rT = dat$rT,
#                    rD = dat$rD,
#                    rU = dat$rU,
#                    rP = dat$rP,
#                    dE = dat$dE,
#                    delayGam = dat$delayGam,
#                    delayMos = dat$delayM,
#                    cD = dat$cD,
#                    cT = dat$cT,
#                    cU = dat$cU,
#                    gamma1 = dat$gamma1,
#                    d1 = dat$d1,
#                    dID = dat$dID,
#                    ID0 = dat$ID0,
#                    kD = dat$kD,
#                    uD = dat$uD,
#                    aD = dat$aD,
#                    fD0 = dat$fD0,
#                    gammaD = dat$gammaD,
#                    b0 = dat$b0,
#                    b1 = dat$b1,
#                    dB = dat$dB,
#                    IB0 = dat$IB0,
#                    kB = dat$kB,
#                    uB = dat$uB,
#                    phi0 = dat$phi0,
#                    phi1 = dat$phi1,
#                    dCA = dat$dCA,
#                    IC0 = dat$IC0,
#                    kC = dat$kC,
#                    uCA = dat$uCA,
#                    dCM = dat$dCM,
#                    tau1 = dat$tau1,
#                    tau2 = dat$tau2,
#                    Q0 = dat$Q0,
#                    chi = dat$chi,
#                    bites_Bed = dat$bites_Bed,
#                    bites_Indoors = dat$bites_Indoors,
#                    p10 = dat$p10,
#                    p2 = dat$p2,
#                    muEL = dat$muEL,
#                    muLL = dat$muLL,
#                    muPL = dat$muPL,
#                    dEL = dat$dEL,
#                    dLL = dat$dLL,
#                    dPL = dat$dPL,
#                    gammaL = dat$gammaL,
#                    beta_larval0 = dat$betaL,
#                    num_int = dat$num_int,
#                    itn_cov = dat$itn_cov,
#                    irs_cov = dat$irs_cov,
#                    ITN_IRS_on = dat$ITN_IRS_on,
#                    d_ITN0 = dat$d_ITN0,
#                    r_ITN0 = dat$r_ITN0,
#                    r_ITN1 = dat$r_ITN1,
#                    r_IRS0 = dat$r_IRS0,
#                    d_IRS0 = dat$d_IRS0,
#                    IRS_interval = dat$IRS_interval,
#                    ITN_interval = dat$ITN_interval,
#                    irs_loss = dat$irs_loss,
#                    itn_loss = dat$itn_loss
#   )
# }", temp)
#
# source(temp)
#
