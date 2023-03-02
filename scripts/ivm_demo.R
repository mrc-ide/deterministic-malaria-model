library(ICDMM)
require(tidyverse)

out <- run_model(model = "mosquito_ivermectin_model",
                 init_EIR = 100,
                  #increasing this EIR gets rid of the dde error
                 gamma_c = 0,
                 gamma_h = 0, time = 730) #need to run a bit longer to get to equilibrium
out$ince
out$delayGam
out$delayMos
out$incv


out$Sv
out$Ev
out$Iv

out$Svih
out$Evih
out$Ivih

out$Svic
out$Evic
out$Ivic

plot(out$t,out$prev, main= "Prevalance", ylim = c(0, 1), type='l')
plot(out$t,out$mv, main= "Total mosq", type='l')
plot(out$t,out$Sv, main= "Total mosq", type='l')


output_mod <- run_model(model = "odin_model", init_EIR = 100, time = 730)
output_mod$delayGam
output_mod$delayMos
plot(output_mod$t,output_mod$prev, main= "Prevalance", ylim = c(0, 1), type='l')
plot(output_mod$t,output_mod$mv, main= "Total mosq", type='l')
plot(output_mod$t,output_mod$Sv, main= "Total mosq", type='l')


output_mod$mv
output_mod$Sv
output_mod$Ev
output_mod$Iv




#MODEL WITHOUT DELAY#####
out <- run_model(model = "mosquito_ivermectin_model_no_delay",
                 init_EIR = 1,
                 gamma_c = 0,
                 gamma_h = 0, time = 500)
plot(out$t,out$prev, main= "Prevalance", ylim = c(0, 1), type='l') #why is the prevalence dropping (ever so slightly)?
out$mv #this is now v stable
plot(out$t,out$mv, main= "Total Mosquitoes", ylim = c(0, 50), type='l')

out_no_ivm <- data.frame(t = out$t, Sv = out$Sv, Ev = out$Ev, Iv = out$Iv)

ggplot(out_no_ivm)+
  geom_line(aes(x = t, y = Sv), col = "blue")+
  geom_line(aes(x = t, y = Ev.47), col = "red")+
  geom_line(aes(x = t, y = Iv), col = "green")+
  ylab("Total Mosq")

out$Sv
out$Ev
out$Iv

out$Svih
out$Evih
out$Ivih

out$Svic
out$Evic
out$Ivic


#now run this model with nets: so far with IVM only, no mort from nets because IVM mort is expressed in terms of mu0
#to check this is correct, try outputting mu and mu_c and mu_h and compare to mu
#if mu_h and mu_c are much bigger than mu, can express mu_h in terms of mu0
out_nets <- run_model(model = "mosquito_ivermectin_model_no_delay",
                      init_EIR = 10,
                      gamma_c = 0,
                      gamma_h = 0,
                      time = 500,
                      ITN_IRS_on = 365,
                      itn_cov = 0.75,
                      num_int = 2) #getting dde error

out_nets$mu
out_nets$mu_c
out_nets$mu_h

#here can see that deaths from IVM would happen faster than death from nets, so probably okay to express IVM mort in terms of mu0

out_no_delay <- run_model(model = "odin_model_no_delay", init_EIR = 100, time = 500)

out_no_del <- data.frame(Sv = out_no_delay$Sv, Ev = out_no_delay$Ev, Iv = out_no_delay$Iv)

out_no_delay$mv
plot(out_no_delay$t,out_no_delay$prev, main= "Prevalance", ylim = c(0, 1), type='l')
plot(out_no_delay$t,out_no_delay$mv, main= "Total Mosquitoes", ylim = c(0, 6), type='l')

out_delay <- run_model(model = "odin_model", time = 500)
out_delay$mv
