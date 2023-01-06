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
require(tidyverse)

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
  init_Sv_di_human <- 0
  init_Ev_di_human <- 0
  init_Iv_di_human <- 0
  initial(Sv) <- init_Sv * mv0
  initial(Sv_di_human) <- init_Sv_di_human*mv0
  initial(Ev_di_human) <- init_Ev_di_human * mv0
  #initial(Ev[1:10]) <- init_Ev/10 * mv0 # Options if not using a delayed delay
  #dim(Ev) <- 10
  initial(Iv_di_human) <- init_Iv_di_human * mv0

  #ivm cow
  init_Sv_di_cow <- 0

  initial(Sv_di_cow) <- init_Sv_di_cow*mv0

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
  #omega <- user() #normalising constant for biting rates
  #FOIvijk[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq[k]*foi_age[i]/omega
  #lag_FOIv=sum(FOIvijk)
  lag_FOIv = 9 #just set to a value and making it a function of Q0 - anthropophagy and av0
  Q0 <- 0.3
  av0 <- 0.333 #1 bite every 3 days

  # Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
  delayGam <- 12.5 # Lag from parasites to infectious gametocytes
  delayMos <- 10 # Extrinsic incubation period.
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
  init_calf <- 0
  init_cow <- 100

  initial(calf) <- init_calf
  initial(cow) <- init_cow

  cattle = calf + cow

  deriv(calf) <- (b*cattle) - (d*calf) - (m*calf)
  deriv(cow) <- (m*calf) - (d*cow)

  m <- 0.00136 #maturation into cows which can be treated with IVM (an assumption)
  d <- 0.000152 #say they can live for 18 years
  b <- d #match birth and death rate

  prop_cow_ivm <- ((cow)/(cow+calf))*coverage_cow_ivm #proportion of cows with ivermectin
  prop_human_ivm <- coverage_humans_ivm #coverage of humans with ivm. age and height restrictions...make up a number

  coverage_cow_ivm = user() #start off with all cows treated with ivm
  coverage_humans_ivm = user()

  mu_ivm = 0.9 #made up a killing effect of ivm from human treatment
  mu_ivm_cow = 0.5 #made up a killing effect of ivm from the cow treatment

  # model options if don't want to use a delayed delay
  #deriv(Ev[1]) <- ince - Ev[1] - mu*Ev[1]
  #deriv(Ev[2:10]) <- Ev[i-1] - Ev[i] - mu*Ev[i]
  #mv = Sv+sum(Ev)+Iv

  # Total mosquito population
  mv = Sv+Sv_di_human+Ev_di_human+Iv_di_human+Sv_di_cow


})

#define model
params <- list(init_Sv = 1000, coverage_cow_ivm = 0, coverage_humans_ivm = 1)
mod <- ento_model_ivm$new(user = params)

#time points: run for 5 years
t1 <- seq(0, 1825, length.out = 1825)

#run model
yy1 <- mod$run(t1)
df_out <- data.frame(yy1)
df_out

ggplot(df_out) +
  geom_line(aes(x = t, y = Iv_di_human), col = "red")

#then run again with the cow ivermectin
params_int <- list(init_Sv = 1000, coverage_cow_ivm = 0, coverage_humans_ivm = 0.5)
mod_int <- ento_model_ivm$new(user = params_int)

#run int model
yy1_int <- mod_int$run(t1)
df_out_int <- data.frame(yy1_int)

ggplot(df_out)+
  geom_line(aes(x = t, y = Iv_di_human), col = "blue")


#####new mosquito model with ivermectin humans and cattle ####

ivm_model <- odin::odin({
  #no ivermectin

  deriv(Sv) <- R - (mu0*Sv) - (FOIhv*Sv) - (a*gamma_c*(1-Q0))*Sv - (a*gamma_h*Q0)*Sv

  deriv(Ev) <- (FOIhv*Sv) - (mu0*Ev) - (a)*Ev - (g*Ev) - (a*gamma_h*Q0)*Ev

  deriv(Iv) <- (g*Ev) - (mu0*Iv)

  #ivermectin humans

  deriv(Svih) <- -(FOIhv*Svih) + (a*gamma_h*Q0)*Sv - (mu_h*Svih)

  deriv(Evih) <- (FOIhv*Sv) - (g*Evih) + (a*gamma_h*Q0*Evih) - (mu_h*Evih)

  deriv(Ivih) <- (g*Evih) - (mu_h*Ivih)

  #ivermectin cattle

  deriv(Svic) <- (a*gamma_c*(1-Q0))*Sv - (mu_c*Svic) - (FOIhv*Svic)

  deriv(Evic) <- (FOIhv*Svic) - (mu_c*Evic) - (g*Evic) + (a*gamma_h*(1-Q0)*Ev)

  deriv(Ivic) <- (g*Evic) - (mu_c*Ivic)


  #initial conditions

  #no ivermectin####
  initial(Sv) <-init_Sv
  initial(Ev) <- init_Ev
  initial(Iv) <- init_Iv

  init_Sv <- user()
  init_Ev <- user()
  init_Iv <- user()

  #ivermectin humans
  initial(Svih) <- init_Svih
  initial(Evih) <- init_Evih
  initial(Ivih) <- init_Ivih

  init_Svih <- user()
  init_Evih <- user()
  init_Ivih <- user()

  #ivermectin cattle
  initial(Svic) <- init_Svic
  initial(Evic) <- init_Evic
  initial(Ivic) <- init_Ivic

  init_Svic <- user()
  init_Evic <- user()
  init_Ivic <- user()

  #parameters
  FOIhv <- (V/H) * a * Q0 * bv * (Ih/Nh)
  V <- Sv+Ev+Iv+Svih+Evih+Ivih+Svic+Evic+Ivic
  Ih <- user() #I can set prevalence in humans this way
  Nh <- 1000
  H <- Nh
  a <- 0.333 #1 bite every 3 days
  Q0 <- user() #proportion of bites that are on humans
  gamma_c <- user() #proportion of livestock with ivermectin
  gamma_h <- user() #proportion of humans with ivermectin
  mu0 <- 0.132 #baseline mortality rate
  mu_c <- 0.728 #elevated mort rate due to IVM cattle. From Dighe and Elong work. per day
  mu_h <- 0.728 #elevated mort rate due to IVM humans
  bv <- 0.05 #probability of transmission from human to vector
  g <- 10 #latent period. days
  R <- mu0 #setting to baseline death rate for ease
  Nv = Sv+Ev+Iv

  #tracking the EIR
  output(EIR) <- (V/H)*a*Q0*(Iv/Nv)

})

#params no ivermectin on humans or cattle ####
params_noivm <- list(init_Sv = 1000, init_Ev = 0, init_Iv = 0,
               init_Svih = 0, init_Evih = 0, init_Ivih = 0,
               init_Svic = 0, init_Evic = 0, init_Ivic = 0,
               gamma_c = 0, gamma_h = 0, Q0 = 0.7, Ih = 600)

mod_noivm <- ivm_model$new(user = params_noivm)

#time points: run for 5 years
t1_noivm <- seq(0, 90, length.out = 90)

#run model
yy1_noivm <- mod_noivm$run(t1_noivm)
df_out_noivm <- data.frame(yy1_noivm)
df_out_noivm$EIR
max(df_out_noivm$EIR)
plot1 <- ggplot(df_out_noivm) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("No ivermectin treatment, Q0 = 0.7")+
  ylim(0, 0.003)+
  geom_hline(aes(yintercept = max(df_out_noivm$EIR)), linetype = "dashed")


#params 100% coverage humans and cattle ####
params_ivm_ch <- list(init_Sv = 1000, init_Ev = 0, init_Iv = 0,
               init_Svih = 0, init_Evih = 0, init_Ivih = 0,
               init_Svic = 0, init_Evic = 0, init_Ivic = 0,
               gamma_c = 1, gamma_h = 1, Q0 = 0.7, Ih = 600)

mod_ivm_ch <- ivm_model$new(user = params_ivm_ch)

#time points: run for 5 years
t1_ivm_ch <- seq(0, 90, length.out = 90)

#run model
yy1_ivm_ch <- mod_ivm_ch$run(t1_ivm_ch)
df_out_ivm_ch <- data.frame(yy1_ivm_ch)
df_out_ivm_ch$EIR

plot2 <- ggplot(df_out_ivm_ch) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage human and cattle IVM, Q0 = 0.7")+
  ylim(0, 0.003)+
  geom_hline(aes(yintercept = max(df_out_noivm$EIR)), linetype = "dashed")


#params no ivermectin and low Q0 ####

params_lowQ0 <- list(init_Sv = 1000, init_Ev = 0, init_Iv = 0,
                      init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                      init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                      gamma_c = 0, gamma_h = 0, Q0 = 0.4, Ih = 600)

mod_lowQ0 <- ivm_model$new(user = params_lowQ0)

#time points: run for 5 years
t1_lowQ0 <- seq(0, 90, length.out = 90)

#run model
yy1_lowQ0 <- mod_lowQ0$run(t1_lowQ0)
df_out_lowQ0 <- data.frame(yy1_lowQ0)
df_out_lowQ0$EIR

plot3 <- ggplot(df_out_lowQ0) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage human and cattle IVM, Q0 = 0.4")+
  ylim(0, 0.003)+
  geom_hline(aes(yintercept = max(df_out_lowQ0$EIR)), linetype = "dashed")


#treat humans and low Q0####
params_lowQ0_ivh <- list(init_Sv = 1000, init_Ev = 0, init_Iv = 0,
                     init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                     init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                     gamma_c = 0, gamma_h = 1, Q0 = 0.4, Ih = 600)

mod_lowQ0_ivh <- ivm_model$new(user = params_lowQ0_ivh)

#time points: run for 5 years
t1_lowQ0_ivh <- seq(0, 90, length.out = 90)

#run model
yy1_lowQ0_ivh <- mod_lowQ0_ivh$run(t1_lowQ0_ivh)
df_out_lowQ0_ivh <- data.frame(yy1_lowQ0_ivh)
df_out_lowQ0_ivh$EIR

plot4 <- ggplot(df_out_lowQ0_ivh) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage humans, no cattle IVM, Q0 = 0.4")+
  ylim(0, 0.003)+
  geom_hline(aes(yintercept = max(df_out_lowQ0$EIR)), linetype = "dashed")+
  geom_hline(aes(yintercept = max(df_out_lowQ0_ivh$EIR)), linetype = "dashed")

#then the additional benefit of spraying cattle
params_5 <- list(init_Sv = 1000, init_Ev = 0, init_Iv = 0,
                         init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                         init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                         gamma_c = 1, gamma_h = 1, Q0 = 0.4, Ih = 600)

mod_5 <- ivm_model$new(user = params_5)

#time points: run for 5 years
t1_5 <- seq(0, 90, length.out = 90)

#run model
yy1_5 <- mod_5$run(t1_5)
df_out_5 <- data.frame(yy1_5)
df_out_5$EIR

plot5 <- ggplot(df_out_5) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage human and cattle IVM, Q0 = 0.4")+
  ylim(0, 0.003)+
  geom_hline(aes(yintercept = max(df_out_lowQ0$EIR)), linetype = "dashed")+
  geom_hline(aes(yintercept = max(df_out_lowQ0_ivh$EIR)), linetype = "dashed")


require(cowplot)

summary_plots <- plot_grid(plot1, plot2, plot3, plot4, plot5)



#with larval....needs more work####
ivm_model_larval <- odin::odin({
  #no ivermectin

  deriv(Sv) <- R - (mu0*Sv) - (FOIhv*Sv) - (a*gamma_c*(1-Q0))*Sv - (a*gamma_h*Q0)*Sv

  deriv(Ev) <- (FOIhv*Sv) - (mu0*Ev) - (a)*Ev - (g*Ev) - (a*gamma_h*Q0)*Ev

  deriv(Iv) <- (g*Ev) - (mu0*Iv)

  #ivermectin humans

  deriv(Svih) <- -(FOIhv*Svih) + (a*gamma_h*Q0)*Sv - (mu_h*Svih)

  deriv(Evih) <- (FOIhv*Sv) - (g*Evih) + (a*gamma_h*Q0*Evih) - (mu_h*Evih)

  deriv(Ivih) <- (g*Evih) - (mu_h*Ivih)

  #ivermectin cattle

  deriv(Svic) <- (a*gamma_c*(1-Q0))*Sv - (mu_c*Svic) - (FOIhv*Svic)

  deriv(Evic) <- (FOIhv*Svic) - (mu_c*Evic) - (g*Evic) + (a*gamma_h*(1-Q0)*Ev)

  deriv(Ivic) <- (g*Evic) - (mu_c*Ivic)


  #initial conditions

  #no ivermectin####
  initial(Sv) <-init_Sv
  initial(Ev) <- init_Ev
  initial(Iv) <- init_Iv

  init_Sv <- user()
  init_Ev <- user()
  init_Iv <- user()

  #ivermectin humans
  initial(Svih) <- init_Svih
  initial(Evih) <- init_Evih
  initial(Ivih) <- init_Ivih

  init_Svih <- user()
  init_Evih <- user()
  init_Ivih <- user()

  #ivermectin cattle
  initial(Svic) <- init_Svic
  initial(Evic) <- init_Evic
  initial(Ivic) <- init_Ivic

  init_Svic <- user()
  init_Evic <- user()
  init_Ivic <- user()

  #parameters
  FOIhv <- (V/H) * a * Q0 * bv * (Ih/Nh)
  V <- Sv+Ev+Iv+Svih+Evih+Ivih+Svic+Evic+Ivic
  Ih <- user() #I can set prevalence in humans this way
  Nh <- 1000
  H <- Nh
  a <- 0.333 #1 bite every 3 days
  Q0 <- user() #proportion of bites that are on humans
  gamma_c <- user() #proportion of livestock with ivermectin
  gamma_h <- user() #proportion of humans with ivermectin
  mu0 <- 0.132 #baseline mortality rate
  mu_c <- 0.728 #elevated mort rate due to IVM cattle. From Dighe and Elong work. per day
  mu_h <- 0.728 #elevated mort rate due to IVM humans
  bv <- 0.05 #probability of transmission from human to vector
  g <- 10 #latent period. days
  R <- (1/2)*PL/dPL
  Nv = Sv+Ev+Iv


  #tracking the EIR
  output(EIR) <- (V/H)*a*Q0*(Iv/Nv)
})


