# load the model package
library(ICDMM)
library (ggplot2)
library(reshape2)


init_age <- c(0,1,5,10,20,30,50)
init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)

init_age <- c(0,0.5,1,5,10,20,30,60)

# 1. Normal model 3 different EIR  ####

# run the model
model_run1 <- run_model(init_EIR =1,age =init_age)#,admin2=admin_str)
plot(rowSums(model_run1$S), t = "l", ylim = c(0, 1))

model_run2 <- run_model(init_EIR =10,age =init_age)#,admin2=admin_str)
lines(rowSums(model_run2$S), col = "red")

model_run3 <- run_model(init_EIR =100,age =init_age)#,admin2=admin_str)
lines(rowSums(model_run3$S), col = "blue")


# 1.1 Meta pop model with 1 patch and 3 different EIR
mixing <- matrix (1, nrow=1)
model_met1 <- run_model_metapop(init_EIR = 1,
                                age =init_age,
                                mix=mixing)
plot(rowSums(model_met1$S), t = "l", ylim = c(0, 1))

model_met2 <- run_model_metapop(init_EIR = 10,
                                age =init_age,
                                mix=mixing)
lines(rowSums(model_met2$S), col = "red")

model_met3 <- run_model_metapop(init_EIR = 100,
                                age =init_age,
                                mix=mixing)
lines(rowSums(model_met3$S), col = "blue")


#1.2 Three parches isolated with different EIR
mixing <- matrix (c(1,0,0,
                    0,1,0,
                    0,0,1),
                  nrow=3)

model_met <- run_model_metapop(init_EIR = c(1,10,100),
                               age =init_age,
                               mix=mixing,
                               time = 1000,
                               irs_cov = c(0,0,0),
                               itn_cov= c(0,0,0))

Ssum <-   data.frame(model_met$Sout)
#data.frame(model_met$Sout+ model_met$Tout+model_met$Dout+ model_met$Aout+model_met$Uout+ model_met$Pout) #
#data.frame(lapply(list(model_met$S),function(x){apply(x,sum,MARGIN=c(1,5))}))
Ssum$t <- model_met$t


plot(Ssum$X1, t = "l", ylim = c(0, 1))
lines(Ssum$X2, col = "red")
lines(Ssum$X3, col = "blue")



#2 Interventions ####


# 2.1 Simple model:

# ITN only
out <- run_model(init_EIR = 50,
                 age =init_age,
                 time = 10000,
                 ITN_IRS_on = 365,
                 itn_cov=0.75,
                 num_int = 2)

plot(out$t,rowSums(out$S), main= "Susceptible", type='l', ylim = c(0, 1))

# ITN + IRS
out2 <- run_model(init_EIR = 50,
                  age =init_age,
                  time = 1000,
                  ITN_IRS_on = 365,
                  itn_cov = 1,
                  irs_cov = 1,
                  num_int = 4)


lines(out2$t,rowSums(out2$S), col = "blue")
abline(v = 365, lty = 2)




# 2.2 Metapop model with 1 patch

mixing <- matrix(1,nrow=1)

out_metapop <- run_model_metapop(init_EIR = 50,
                                 age =init_age,
                                 time = 10000,
                                 ITN_IRS_on = 365,
                                 itn_cov = 0.75,
                                 num_int = 2,
                                 mix = mixing)


plot(out_metapop$t,rowSums(out_metapop$S), main= "Susceptible-Meta", type='l', ylim = c(0, 1))

out_metapop2 <- run_model_metapop(init_EIR = 50,
                                  age =init_age,
                                  time = 1000,
                                  ITN_IRS_on = 365,
                                  itn_cov = 0.75,
                                  irs_cov = 0.4,
                                  num_int = 4,
                                  mix = mixing)

lines(out_metapop2$t,rowSums(out_metapop2$S), col = "blue")
abline(v = 365, lty = 2)



# 2.3 Metapop with two independent patches


d2<-matrix(c(0.5,0.5,0.5,0.5),nrow=2)

d3<-matrix(c(1,0,0,
             0,1,0,
             0,0,1),nrow = 3)

mixing <- d2




#

out_metapop <- run_model_metapop(init_EIR =c(50,50),
                                 age =init_age,
                                 time = 1000,
                                 ITN_IRS_on = 365,
                                 itn_cov = c(0.75,0.75),
                                 irs_cov = c (0,0),
                                 num_int = 2,
                                 mix = mixing)

Ssum <- data.frame(lapply(list(out_metapop$S),function(x){apply(x,sum,MARGIN=c(1,5))}))
Ssum$t <-out_metapop$t

x<-data.frame(out_metapop$Sout+ out_metapop$Tout+out_metapop$Dout+ out_metapop$Aout+out_metapop$Uout+ out_metapop$Pout) #



plot(Ssum$t,Ssum$X1, type='l', ylim = c(0, 1))
lines(Ssum$X2, col = "red")
lines(Ssum$X3, col = "green")


# ITN + IRS
out_meta_dob2 <- run_model_metapop(init_EIR = c(10,50),
                                   age =init_age,
                                   time =1000,
                                   ITN_IRS_on = 365,
                                   itn_cov = c(0.75,0.75),
                                   irs_cov = c(0.4,0.4),
                                   num_int = 4,
                                   mix =mixing)

Ssum2 <- data.frame(lapply(list(out_meta_dob2$S),function(x){apply(x,sum,MARGIN=c(1,5))}))
Ssum2$t <-out_meta_dob2$t

plot(Ssum2$t,Ssum2$X1, type='l', ylim = c(0, 1))

lines(Ssum2$X1, col = "blue")
lines(Ssum2$X2, col = "yellow")
lines(Ssum2$X3, col = "black")



