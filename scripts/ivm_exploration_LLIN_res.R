#exploring the mosquito_ivermectin_model odin model with insecticide resistance and LLINs####

#thinking: LLINs and resistance will both (separately and together?)

out_1 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 0,
                   gamma_h = 0, time = 1000)

out_2 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 0,
                   gamma_h = 0, time = 1000,
                   ITN_IRS_on = 365, itn_cov = 0.75, num_int = 2) #nets only at 75% coverage


out_3 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 0,
                   gamma_h = 1, time = 1000,
                   ITN_IRS_on = 365, itn_cov = 0.75, num_int = 2) #nets only at 75% coverage and human ivm


out_4 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 0,
                   gamma_h = 1, time = 1000) #no nets

out_5 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 1,
                   gamma_h = 1, time = 1000, ITN_IRS_on = 365, itn_cov = 0.75, num_int = 2) #nets (75% cov), human and cattle ivm


plot(out_1$t,out_1$prev, main= "Prevalance", type='l', ylim = c(0, 1))
lines(out_2$t, out_2$prev, col = "blue")
lines(out_3$t, out_3$prev, col = "red")
lines(out_4$t, out_4$prev, col = "green")
lines(out_5$t, out_5$prev, col = "pink")
abline(v = 365, lty = 2)

#reach a lower prevalence when have both nets and IVM, but ignoring the correlation in ivm and llin coverage

#first, need to actually quantify the proportion of bites that are taken on cattle due to nets. Let's consider a situation where net cov is 60%, 80% and 100%
#explore Q0 of 40% and 60%

mod_1 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   gamma_c = 0, gamma_h = 0,
                   time = 2000,
                   ITN_IRS_on = 365, itn_cov = 0.6, num_int = 2,
                   Q0 = 0.4)

mod_2 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   gamma_c = 0, gamma_h = 0,
                   time = 2000,
                   ITN_IRS_on = 365, itn_cov = 0.8, num_int = 2,
                   Q0 = 0.4)

mod_3 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   gamma_c = 0, gamma_h = 0,
                   time = 2000,
                   ITN_IRS_on = 365, itn_cov = 1, num_int = 2,
                   Q0 = 0.4)

plot(mod_1$t,mod_1$prev, main= "Prevalence", type='l', ylim = c(0, 1))
lines(mod_2$t, mod_2$prev, col = "blue")
lines(mod_3$t, mod_3$prev, col = "green")

#track Q
plot(mod_1$t,mod_1$Q, main= "Q", type='l', ylim = c(0, 1))
lines(mod_2$t, mod_2$Q, col = "blue")
lines(mod_3$t, mod_3$Q, col = "green")

#look at Q at day 1, 366 and 1000 (60% net coverage and Q0 = 0.4)
mod_1$Q[1] #0.4
mod_1$Q[366] #0.23
mod_1$Q[1460] #0.31
((0.31-0.23)/0.23)*100 #34.78%

#80% net coverage, Q0 = 0.4
mod_2$Q[1] #0.4
mod_2$Q[366] #0.1709
mod_2$Q[1460] #0.283
(0.283-0.1709)/0.1709 #65.59% increase in Q

#100% net coverage, Q0 = 0.4
mod_3$Q[1] #0.4
mod_3$Q[366] #0.0835
mod_3$Q[1460] #0.2464719
(0.2465-0.0835)/0.0835 #195.21% increase

mod_4 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   gamma_c = 0, gamma_h = 0,
                   time = 2000,
                   ITN_IRS_on = 365, itn_cov = 0.6, num_int = 2,
                   Q0 = 0.6)

mod_4$Q[1] #0.6
mod_4$Q[366] #0.4196271
mod_4$Q[1460] #0.5101814
(0.51018-0.41962)/0.41962 #21.58% increase in Q

mod_5 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   gamma_c = 0, gamma_h = 0,
                   time = 2000,
                   ITN_IRS_on = 365, itn_cov = 0.8, num_int = 2,
                   Q0 = 0.6)
mod_5$Q[1] #0.6
mod_5$Q[366] #0.31695
mod_5$Q[1460] #0.4705529
(0.4705529-0.31695)/0.31695 #48.46% increase in Q


mod_6 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   gamma_c = 0, gamma_h = 0,
                   time = 2000,
                   ITN_IRS_on = 365, itn_cov = 1, num_int = 2,
                   Q0 = 0.6)
mod_6$Q[1] #0.6
mod_6$Q[366] #0.1701589
mod_6$Q[1460] #0.4239477
(0.4239477 - 0.1701589)/0.1701589 #149.15% increase
