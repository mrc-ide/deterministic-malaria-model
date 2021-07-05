## ------- Debug checking points   ------   ####


# S <- data.frame(lapply(list(state_use$init_S),function(x){apply(x,sum,MARGIN=c(4))}))
# T <- data.frame(lapply(list(state_use$init_T),function(x){apply(x,sum,MARGIN=c(4))}))
# A <- data.frame(lapply(list(state_use$init_A),function(x){apply(x,sum,MARGIN=c(4))}))
# U <- data.frame(lapply(list(state_use$init_U),function(x){apply(x,sum,MARGIN=c(4))}))
# D <- data.frame(lapply(list(state_use$init_D),function(x){apply(x,sum,MARGIN=c(4))}))
# P <- data.frame(lapply(list(state_use$init_P),function(x){apply(x,sum,MARGIN=c(4))}))
#
#
# Total_hum<- S+T+A+U+D+P
#
# Total_mosq <- state_use$init_Iv + state_use$init_Ev + state_use$init_Sv
# Total_larv <- state_use$init_PL + state_use$init_LL + state_use$init_EL
#
#
# # Immunity status
#
# IB <- data.frame(lapply(list(state_use$init_IB),function(x){apply(x,sum,MARGIN=c(4))}))
# ID <- data.frame(lapply(list(state_use$init_ID),function(x){apply(x,sum,MARGIN=c(4))}))
# ICA <- data.frame(lapply(list(state_use$init_ICA),function(x){apply(x,sum,MARGIN=c(4))}))
# ICM <- data.frame(lapply(list(state_use$init_ICM),function(x){apply(x,sum,MARGIN=c(4))}))
#
#
## ---- Comparing 2 vs 4 interventions ----- ####

# x<- np2$init_P - np4$init_P[,,1:2,]
# any (x !=0)
#
#
# S2 <- data.frame(lapply(list(np2$init_S),function(x){apply(x,sum,MARGIN=c(1,3))}))
# S4 <- data.frame(lapply(list(np4$init_S),function(x){apply(x,sum,MARGIN=c(1,3))}))
#
#
# T2 <- data.frame(lapply(list(np2$init_T),function(x){apply(x,sum,MARGIN=c(1,3))}))
# T4 <- data.frame(lapply(list(np4$init_T),function(x){apply(x,sum,MARGIN=c(1,3))}))
#
#
# A2 <- data.frame(lapply(list(np2$init_A),function(x){apply(x,sum,MARGIN=c(3,4))}))
# A4 <- data.frame(lapply(list(np4$init_A),function(x){apply(x,sum,MARGIN=c(3,4))}))
#
# U2 <- data.frame(lapply(list(np2$init_U),function(x){apply(x,sum,MARGIN=c(3,4))}))
# U4 <- data.frame(lapply(list(np4$init_U),function(x){apply(x,sum,MARGIN=c(3,4))}))
#
# D2 <- data.frame(lapply(list(np2$init_D),function(x){apply(x,sum,MARGIN=c(3,4))}))
# D4 <- data.frame(lapply(list(np4$init_D),function(x){apply(x,sum,MARGIN=c(3,4))}))
#
# P2 <- data.frame(lapply(list(np2$init_P),function(x){apply(x,sum,MARGIN=c(3,4))}))
# P4 <- data.frame(lapply(list(np4$init_P),function(x){apply(x,sum,MARGIN=c(3,4))}))
#
# ICM2 <- data.frame(lapply(list(np2$init_ICM),function(x){apply(x,sum,MARGIN=c(3,4))}))
# ICM4<- data.frame(lapply(list(np4$init_ICM),function(x){apply(x,sum,MARGIN=c(3,4))}))

#
