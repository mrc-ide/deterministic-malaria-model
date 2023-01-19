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



#####new mosquito model with ivermectin humans and cattle ####

ivm_model <- odin::odin({
  #no ivermectin

  deriv(Sv) <- R - (mu0*Sv) - (FOIhv*Sv) - (a*gamma_c*(1-Q0))*Sv - (a*gamma_h*Q0)*Sv

  deriv(Ev) <- (FOIhv*Sv) - (mu0*Ev) - (g*Ev) - (a*gamma_h*Q0)*Ev

  deriv(Iv) <- (g*Ev) - (mu0*Iv) - (a*gamma_h*Q0*Iv) - (a*gamma_c*(1-Q0)*Iv)

  #ivermectin humans

  deriv(Svih) <- -(FOIhv*Svih) + (a*gamma_h*Q0)*Sv - (mu_h*Svih)

  deriv(Evih) <- (FOIhv*Sv) - (g*Evih) + (a*gamma_h*Q0*Evih) - (mu_h*Evih)

  deriv(Ivih) <- (g*Evih) - (mu_h*Ivih) + (a*gamma_h*Q0*Iv)

  #ivermectin cattle

  deriv(Svic) <- (a*gamma_c*(1-Q0))*Sv - (mu_c*Svic) - (FOIhv*Svic)

  deriv(Evic) <- (FOIhv*Svic) - (mu_c*Evic) - (g*Evic) + (a*gamma_h*(1-Q0)*Ev)

  deriv(Ivic) <- (g*Evic) - (mu_c*Ivic) - (a*gamma_c*(1-Q0)*Iv)


  #initial conditions

  #no ivermectin####
  initial(Sv) <- init_Sv
  initial(Ev) <- init_Ev
  initial(Iv) <- init_Iv

  init_Sv <- 50000000
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
  FOIhv <- a * Q0 * bv * 0.2
  V <- Sv+Ev+Iv+Svih+Evih+Ivih+Svic+Evic+Ivic
  #Ih <- user() #I can set prevalence in humans this way
  Nh <- 1000
  H <- Nh
  a <- 0.333 #1 bite every 3 days
  Q0 <- user() #proportion of bites that are on humans
  gamma_c <- user() #proportion of livestock with ivermectin
  gamma_h <- user() #proportion of humans with ivermectin
  mu0 <- 0.1 #baseline mortality rate
  mu_c_elev <- 0.628
  mu_h_elev <- 0.628
  mu_c <- mu0 + mu_c_elev #elevated mort rate due to IVM cattle. From Dighe and Elong work. per day
  mu_h <- mu0 + mu_h_elev #elevated mort rate due to IVM humans
  bv <- 0.05 #probability of transmission from human to vector
  g <- 0.1 #1/latent period. 10 days
  R <- mu0*Nv #setting to baseline death rate for ease
  Nv = Sv+Ev+Iv

  #tracking the EIR
  output(EIR) <- (V/H)*a*Q0*((Iv+Ivih+Ivic)/Nv)
  output(total_mosq) <- Nv

})

#params no ivermectin on humans or cattle ####
params_1 <- list(init_Ev = 0, init_Iv = 0,
               init_Svih = 0, init_Evih = 0, init_Ivih = 0,
               init_Svic = 0, init_Evic = 0, init_Ivic = 0,
               gamma_c = 0, gamma_h = 0, Q0 = 0.7)

mod_1 <- ivm_model$new(user = params_1)

#time points: run for 90 days
t1_1 <- seq(0, 90, length.out = 90)

#run model
yy1_1 <- mod_1$run(t1_1)
df_1 <- data.frame(yy1_1)
df_1$total_mosq #getting a stable mosquito population this way

ylim_upper <- max(df_1$EIR)
plot1 <- ggplot(df_1) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("No ivermectin treatment, Q0 = 0.7")+
  ylim(0, ylim_upper)

#params 100% coverage humans and cattle ####
params_2 <- list(init_Ev = 0, init_Iv = 0,
               init_Svih = 0, init_Evih = 0, init_Ivih = 0,
               init_Svic = 0, init_Evic = 0, init_Ivic = 0,
               gamma_c = 1, gamma_h = 1, Q0 = 0.7, Ih = 600)

mod_2 <- ivm_model$new(user = params_2)

#time points: run for 90 days
t1_2 <- seq(0, 90, length.out = 90)

#run model
yy1_2 <- mod_2$run(t1_2)
df_2 <- data.frame(yy1_2)

plot2 <- ggplot(df_2) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage human and cattle IVM, Q0 = 0.7")+
  ylim(0, ylim_upper)

#params 100% IVM humans, no cattle
params_3 <- list(init_Ev = 0, init_Iv = 0,
                      init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                      init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                      gamma_c = 0, gamma_h = 1, Q0 = 0.7, Ih = 600)

mod_3 <- ivm_model$new(user = params_3)

#time points: run for 90 days
t1_3 <- seq(0, 90, length.out = 90)

#run model
yy1_3 <- mod_3$run(t1_3)
df_3 <- data.frame(yy1_3)

plot_3 <- ggplot(df_3) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage human, no cattle IVM, Q0 = 0.7")+
  ylim(0, ylim_upper)

#params 100% cov cattle, no humans
params_4 <- list(init_Ev = 0, init_Iv = 0,
                 init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                 init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                 gamma_c = 1, gamma_h = 0, Q0 = 0.7, Ih = 600)

mod_4 <- ivm_model$new(user = params_4)

#time points: run for 90 days
t1_4 <- seq(0, 90, length.out = 90)

#run model
yy1_4 <- mod_4$run(t1_4)
df_4 <- data.frame(yy1_4)

plot_4 <- ggplot(df_4) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage cattle, no human IVM, Q0 = 0.7")+
  ylim(0, ylim_upper)


#params no ivermectin and low Q0 ####

params_5 <- list(init_Ev = 0, init_Iv = 0,
                      init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                      init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                      gamma_c = 0, gamma_h = 0, Q0 = 0.4, Ih = 600)

mod_5 <- ivm_model$new(user = params_5)

#time points: run for 90 days
t1_5 <- seq(0, 90, length.out = 90)

#run model
yy1_5 <- mod_5$run(t1_5)
df_5 <- data.frame(yy1_5)

plot5 <- ggplot(df_5) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("No IVM, Q0 = 0.4")+
  ylim(0, ylim_upper)

#treat humans and cattle and low Q0####
params_6 <- list(init_Ev = 0, init_Iv = 0,
                     init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                     init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                     gamma_c = 1, gamma_h = 1, Q0 = 0.4, Ih = 600)

mod_6 <- ivm_model$new(user = params_6)

#time points: run for 90 days
t1_6 <- seq(0, 90, length.out = 90)

#run model
yy1_6 <- mod_6$run(t1_6)
df_6 <- data.frame(yy1_6)

plot6 <- ggplot(df_6) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage humans and cattle IVM, Q0 = 0.4")+
  ylim(0, ylim_upper)

#then the additional benefit of spraying cattle
params_7 <- list(init_Ev = 0, init_Iv = 0,
                         init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                         init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                         gamma_c = 0, gamma_h = 1, Q0 = 0.4, Ih = 600)

mod_7 <- ivm_model$new(user = params_7)

#time points: run for 5 years
t1_7 <- seq(0, 90, length.out = 90)

#run model
yy1_7 <- mod_7$run(t1_7)
df_7 <- data.frame(yy1_7)


plot7 <- ggplot(df_7) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% coverage human, no cattle IVM, Q0 = 0.4")+
  ylim(0, ylim_upper)

#100% cov cattle, no human IVM. Q0 = 0.4

params_8 <- list( init_Ev = 0, init_Iv = 0,
                 init_Svih = 0, init_Evih = 0, init_Ivih = 0,
                 init_Svic = 0, init_Evic = 0, init_Ivic = 0,
                 gamma_c = 1, gamma_h = 0, Q0 = 0.4, Ih = 600)

mod_8 <- ivm_model$new(user = params_8)

#time points: run for 90 days
t1_8 <- seq(0, 90, length.out = 90)

#run model
yy1_8 <- mod_5$run(t1_8)
df_8 <- data.frame(yy1_8)

plot8 <- ggplot(df_8) +
  geom_line(aes(x = t, y = EIR), col = "red")+
  ggtitle("100% cattle IVM and no human, Q0 = 0.4")+
  ylim(0, ylim_upper)

require(cowplot)

summary_plots <- plot_grid(plot1, plot2, plot_3, plot_4, plot5, plot6, plot7, plot8,
                           ncol = 4, nrow = 4)

#viz of the different scenarios####

#no ivermectin treatment in both
df_1 <- df_1 %>%
  mutate(scenario = "No IVM, Q0 = 0.7")

df_5 <- df_5 %>%
  mutate(scenario = "No IVM, Q0 = 0.4")

no_ivm_dat <- rbind(df_1, df_5)

no_ivm_plot <- ggplot(no_ivm_dat) +
  geom_line(aes(x  = t, y = EIR, col = as.factor(scenario)))+
  theme_minimal()+
  labs(col = "Scenario")+
  ylim(0, ylim_upper)


#100% IVM in both
df_2 <- df_2 %>%
  mutate(scenario = "100% IVM human and cattle, Q0 = 0.7")

df_6 <- df_6 %>%
  mutate(scenario = "100% IVM human and cattle, Q0 = 0.4")

ivm_both_dat <- rbind(df_2, df_6)

ivm_both_plot <- ggplot(ivm_both_dat)+
  geom_line(aes(x  = t, y = EIR, col = as.factor(scenario)))+
  theme_minimal()+
  labs(col = "Scenario")+
  ylim(0, ylim_upper)

#100% cov humans, no cattle
df_3 <- df_3 %>%
  mutate(scenario = "100% IVM humans, 0% cattle, Q0 = 0.7")

df_7 <- df_7 %>%
  mutate(scenario = "100% IVM humans, 0% cattle, Q0 = 0.4")

ivm_human_only_dat <- rbind(df_3, df_7)

ivm_human_only_plot <- ggplot(ivm_human_only_dat)+
  geom_line(aes(x = t, y = EIR, col = as.factor(scenario)))+
  theme_minimal()+
  labs(col = "Scenario")+
  ylim(0, ylim_upper)

#100% cov cattle, no human
df_4 <- df_4 %>%
  mutate(scenario = "100% cattle, 0% human, Q0 = 0.7")

df_8 <- df_8 %>%
  mutate(scenario = "100% cattle, 0% human, Q0 = 0.4")

ivm_cattle_only_dat <- rbind(df_4, df_8)

ivm_cattle_only_plot <- ggplot(ivm_cattle_only_dat)+
  geom_line(aes(x = t, y = EIR, col = as.factor(scenario)))+
  theme_minimal()+
  labs(col = "Scenario")+
  ylim(0, ylim_upper)

scenario_plots <- plot_grid(no_ivm_plot, ivm_both_plot, ivm_human_only_plot, ivm_cattle_only_plot,
                            labels = c("A", "B", "C", "D"))


#do the same but anthropophagy plots so two panels

#high Q0 df
high_Q0_df <- do.call("rbind", list(df_1, df_2, df_3, df_4))

high_Q0_plot <- ggplot(high_Q0_df)+
  geom_line(aes(x = t, y = EIR, col = as.factor(scenario)))+
  theme_minimal()+
  labs(col = "Scenario")+
  ggtitle("Q0 = 0.7")+
  ylim(0, ylim_upper)

#low Q0 df
low_Q0_df <- do.call("rbind", list(df_5, df_6, df_7, df_8))

low_Q0_df %>%
  filter(scenario == "100% cattle, 0% human, Q0 = 0.4") %>%
  ggplot()+
  geom_line(aes(x = t, y = EIR, col = as.factor(scenario)))+
  theme_minimal()+
  #labs(col = "Scenario")+
  ggtitle("Q0 = 0.4")+
  ylim(0, ylim_upper)

low_Q0_plot <- ggplot(low_Q0_df)+
  geom_line(aes(x = t, y = EIR, col = as.factor(scenario)))+
  theme_minimal()+
  labs(col = "Scenario")+
  ggtitle("Q0 = 0.4")+
  ylim(0, ylim_upper)

unique(low_Q0_df$scenario)
plot_grid(high_Q0_plot, low_Q0_plot)


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


