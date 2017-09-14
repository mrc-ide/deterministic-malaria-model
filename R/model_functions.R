#------------------------------------------------
#' Run_Model
#'
#' \code{Run_Model} runs model using declared age, EIR, ft, country, admin
#'
#' @param age Vector of age brackets.
#' @param EIR Numeric of annual EIR
#' @param ft Numeric of proportion of symptomatic cases recieving treatment
#' @param time Numeric of length of time that the model will run for in days
#' @param admin2 Character of admin unit
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous xlab ylab geom_line
#' @importFrom reshape2 melt
#' @export


Run_Model <- function(age,EIR,ft,admin2,time){
  mpl <- Model_Param_List_Create()
  # generate initial state variables from equilibrium solution
  state <- Equilibrium_Init_Create(age=age,EIR=EIR,ft=ft,model.param.list = mpl,het.brackets=5,admin.unit = admin2)
  # create odin generator
  odin_model_path <- system.file("extdata/odin_model2.R",package="hanojoel")
  gen <- odin::odin(odin_model_path,verbose=FALSE,build = TRUE)


  #create model with initial values
  mod <- generate_default_model(ft=ft,age=age,dat=state,generator=gen,dde=TRUE)
  tt <- seq(0,time,1)
  # run model
  mod_run <- mod$run(tt)
  # shape output
  out <- mod$transform_variables(mod_run)
  hum <- list("S"=out$S,"T"=out$T,"D"=out$D,"A"=out$A,"U"=out$U,"P"=out$P)
  humsum <- lapply(hum,function(x){apply(x,sum,MARGIN=1)})
  df <- as.data.frame(matrix(unlist(humsum),nrow=length(tt)))
  colnames(df) <- c("S","T","D","A","U","P")
  df$t <- tt
  brks <- seq(from=0,to=time,by=365)
  ret <- ggplot(melt(df,id.vars="t"),aes(x=t,y=value,col=variable)) + geom_line() +
    scale_colour_manual(values=c("#000000","#CC0000","#339900","#3333FF","#CC9933","#CC3399"),name="Compartment",labels=c("S","T","D","A","U","P")) +
    scale_x_continuous(breaks=brks,labels=0:(length(brks)-1)) +
    xlab("Years") + ylab("Proportion of population")
  return(list("plot"=ret,"dat"=out))
}


## Odin generator function
generate_default_model <- function(ft,age,dat,generator,dde = TRUE){
  mod <- generator(init_S=dat$S,
                   init_T=dat$T,
                   init_D=dat$D,
                   init_A=dat$A,
                   init_U=dat$U,
                   init_P=dat$P,
                   init_ICA = dat$ICA,
                   init_ICM = dat$ICM,
                   init_ID = dat$ID,
                   init_IB = dat$IB,
                   init_Sv = dat$Sv,
                   init_Ev = dat$Ev,
                   init_Iv = dat$Iv,
                   init_PL = dat$PL,
                   init_LL = dat$LL,
                   init_EL = dat$EL,
                   na = dat$na,
                   nh = dat$nh,
                   age_rate = dat$age_rate,
                   foi_age = dat$foi_age,
                   het_wt = dat$het_wt,
                   rel_foi = dat$rel_foi,
                   omega = dat$omega,
                   pi = pi,
                   x_I = dat$x_I,
                   #eov = dat$eov,
                   mv0 = dat$mv0,
                   ssa0 = dat$ssa0,
                   ssa1 = dat$ssa1,
                   ssa2 = dat$ssa2,
                   ssa3 = dat$ssa3,
                   ssb1 = dat$ssb1,
                   ssb2 = dat$ssb2,
                   ssb3 = dat$ssb3,
                   theta_c = dat$theta_c,
                   den = dat$den,
                   age59 = dat$age59,
                   age05 = dat$age05,
                   age = age*365,
                   ft = ft,
                   use_dde = dde,
                   eta = dat$eta,
                   rA = dat$rA,
                   rT = dat$rT,
                   rD = dat$rD,
                   rU = dat$rU,
                   rP = dat$rP,
                   dE = dat$dE,
                   delayGam = dat$delayGam,
                   delayMos = dat$delayM,
                   cD = dat$cD,
                   cT = dat$cT,
                   cU = dat$cU,
                   gamma1 = dat$gamma1,
                   d1 = dat$d1,
                   dID = dat$dID,
                   ID0 = dat$ID0,
                   kD = dat$kD,
                   uD = dat$uD,
                   aD = dat$aD,
                   fD0 = dat$fD0,
                   gammaD = dat$gammaD,
                   b0 = dat$b0,
                   b1 = dat$b1,
                   dB = dat$dB,
                   IB0 = dat$IB0,
                   kB = dat$kB,
                   uB = dat$uB,
                   phi0 = dat$phi0,
                   phi1 = dat$phi1,
                   dCA = dat$dCA,
                   IC0 = dat$IC0,
                   kC = dat$kC,
                   uCA = dat$uCA,
                   dCM = dat$dCM,
                   tau1 = dat$tau1,
                   tau2 = dat$tau2,
                   Q0 = dat$Q0,
                   chi = dat$chi,
                   bites_Bed = dat$bites_Bed,
                   bites_Indoors = dat$bites_Indoors,
                   p10 = dat$p10,
                   p2 = dat$p2,
                   muEL = dat$muEL,
                   muLL = dat$muLL,
                   muPL = dat$muPL,
                   dEL = dat$dEL,
                   dLL = dat$dLL,
                   dPL = dat$dPL,
                   gammaL = dat$gammaL,
                   beta_larval0 = dat$betaL,
                   num_int = dat$num_int,
                   itn_cov = dat$itn_cov,
                   irs_cov = dat$irs_cov,
                   ITN_IRS_on = dat$ITN_IRS_on,
                   d_ITN0 = dat$d_ITN0,
                   r_ITN0 = dat$r_ITN0,
                   r_ITN1 = dat$r_ITN1,
                   r_IRS0 = dat$r_IRS0,
                   d_IRS0 = dat$d_IRS0,
                   IRS_interval = dat$IRS_interval,
                   ITN_interval = dat$ITN_interval,
                   irs_loss = dat$irs_loss,
                   itn_loss = dat$itn_loss
  )
}
