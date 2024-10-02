#compute cases averted when we take hourly vs segmented data collection methods to estimate phi

#compare this to the cases averted for comparing different net types (use a low resistance)

#eg 10% for pyrethroid, 0% for PBO and 0% for IG2

#deterministic malaria model
#looking to see the impact of different phi estimates on LLIN impact
require(ICDMM)
require(tidyverse)
# load the model package (or have it built and reloaded as described above)

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 100

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0

# Define time for turning on interventions
ITN_IRS_on <- 365

#from SeyoumB: study furthest away from the straight line
seyoumB <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/SeyoumB.csv")

#phiBs
phi_B_seg <-  round(seyoumB$new_PHI_B_seg,2)
phi_B_seg_L <- round(seyoumB$new_PHI_B_seg_L, 2)
phi_B_seg_U <- round(seyoumB$new_PHI_B_seg_U, 2)

phi_B_hourly <- round(seyoumB$new_PHI_B,2)
phi_B_hourly_L <- round(seyoumB$new_PHI_B_L, 2)
phi_B_hourly_U <- round(seyoumB$new_PHI_B_U, 2)

#phiIs
phi_I_seg <-  round(seyoumB$new_PHI_I_seg,2)
phi_I_seg_L <- round(seyoumB$new_PHI_I_seg_L,2)
phi_I_seg_U <- round(seyoumB$new_PHI_I_seg_U,2)

phi_I_hourly <- round(seyoumB$new_PHI_I,2)
phi_I_hourly_L <- round(seyoumB$new_PHI_I_L,2)
phi_I_hourly_U <- round(seyoumB$new_PHI_I_U,2)


bites_in_vec <- c(phi_I_seg, phi_I_seg_L, phi_I_seg_U, phi_I_hourly, phi_I_hourly_L, phi_I_hourly_U)
bites_bed_vec <- c(phi_B_seg, phi_B_seg_L, phi_B_seg_U, phi_B_hourly, phi_B_hourly_L, phi_B_hourly_U)
itn_cov_vec <- c(0.2, 0.8)

#get net efficacies
pyr_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_only_nets.csv")
head(pyr_df)

pyr_df <- pyr_df %>%
  filter(resistance == 0.50)
pyr_df$g

dn0_vec_pyr <- c(pyr_df$dn0_lo10, pyr_df$dn0_med, pyr_df$dn0_up90)
rn0_vec_pyr <- c(pyr_df$rn0_lo10, pyr_df$rn0_med, pyr_df$rn0_up90)

#ELLIE SAID USE THE GAMMAN PYR-ONLY HALF LIFE FOR ALL OTHER NET TYPES
gamman_vec_pyr <- c(pyr_df$gamman_lo10, pyr_df$gamman_med, pyr_df$gamman_up90)
ento_df_pyr <- expand.grid(bites_Bed = bites_bed_vec,
                       itn_cov = itn_cov_vec,
                       d_ITN0 = dn0_vec_pyr)

ento_df_pyr <- ento_df_pyr %>%
  mutate(r_ITN0 = case_when(d_ITN0 == pyr_df$dn0_lo10 ~ rn0_vec_pyr[1],
                            d_ITN0 == pyr_df$dn0_med ~ rn0_vec_pyr[2],
                            d_ITN0 == pyr_df$dn0_up90 ~ rn0_vec_pyr[3],
                            TRUE ~ NA_real_),
         itn_half_life = case_when(d_ITN0 == pyr_df$dn0_lo10 ~ gamman_vec_pyr[1]*365,
                                   d_ITN0 == pyr_df$dn0_med ~ gamman_vec_pyr[2]*365,
                                   d_ITN0 == pyr_df$dn0_up90 ~ gamman_vec_pyr[3]*365,
                                   TRUE ~ NA_real_))

ento_df_pyr <- ento_df_pyr %>%
  mutate(bites_Indoors = case_when(bites_Bed == bites_bed_vec[1] ~ bites_in_vec[1],
                                   bites_Bed == bites_bed_vec[2] ~ bites_in_vec[2],
                                   bites_Bed == bites_bed_vec[3] ~ bites_in_vec[3],
                                   bites_Bed == bites_bed_vec[4] ~ bites_in_vec[4],
                                   bites_Bed == bites_bed_vec[5] ~ bites_in_vec[5],
                                   bites_Bed == bites_bed_vec[6] ~ bites_in_vec[6]))


ento_param_list_pyr <- list()
for (i in seq_len(nrow(ento_df_pyr))){
  ento_param_list_pyr[[i]] <- as.numeric(ento_df_pyr[i,])
}

out_nets_pyr <- function(itn_input){
  bites_Bed_in <- itn_input[1]
  itn_cov_in <- itn_input[2]
  d_ITN0_in <- itn_input[3]
  r_ITN0_in <- itn_input[4]
  itn_half_life_in <- itn_input[5]
  bites_In_in<- itn_input[6]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 2,
                      itn_cov = itn_cov_in,
                      ITN_IRS_on = ITN_IRS_on,
                      init_ft = prop_treated,
                      bites_Bed = bites_Bed_in,
                      bites_Indoors = bites_In_in,
                      itn_half_life = itn_half_life_in,
                      d_ITN0 = d_ITN0_in,
                      r_ITN0 = r_ITN0_in)
  return(output)
}


#my_sim_antag_ITN <- function(){
#  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
#  pyr_out_list_antag_ITN <- lapply(pyr_param_list, antag_ITN_cov_loop)
#  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
#  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
#    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
#    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
#                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
#    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
#  return(pyr_out_df_antag_ITN)
#
#}



my_sim_nets_pyr <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  ento_out <- lapply(ento_param_list_pyr, out_nets_pyr)
  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_df_pyr)), function(x){
    df <- as.data.frame(ento_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, itn_cov, EIR_tot, prev,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, bites_Indoors, Q0,inc05,s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
  return(ento_out_mod)

}

ento_nets_pyr <- my_sim_nets_pyr()



ento_nets_pyr <- ento_nets_pyr %>%
  mutate(segments = case_when(bites_Indoors == bites_in_vec[1] & bites_Bed == bites_bed_vec[1] ~ "segment_M",
                              bites_Indoors == bites_in_vec[2] & bites_Bed == bites_bed_vec[2] ~ "segment_L",
                              bites_Indoors == bites_in_vec[3] & bites_Bed == bites_bed_vec[3] ~ "segment_U",
                              bites_Bed == bites_in_vec[4] & bites_Bed == bites_bed_vec[4] ~ "hourly_M",
                              bites_Bed == bites_in_vec[5] & bites_Bed == bites_bed_vec[5] ~ "hourly_L",
                              bites_Bed == bites_in_vec[6] & bites_Bed == bites_bed_vec[6] ~ "hourly_L",
                              TRUE ~ NA_character_))

ento_nets_pyr <- ento_nets_pyr %>%
  mutate(net_eff = case_when(d_ITN0 == dn0_vec_pyr[1] & r_ITN0 == rn0_vec_pyr[1] ~ "estim_lower",
                             d_ITN0 == dn0_vec_pyr[2] & r_ITN0 == rn0_vec_pyr[2] ~ "estim_med",
                             d_ITN0 == dn0_vec_pyr[3] & r_ITN0 == rn0_vec_pyr[3] ~ "estim_upper",
                             TRUE ~ NA_character_))
#repeat for pbo
pbo_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_pbo_nets.csv")
head(pbo_df)

pbo_df <- pbo_df %>%
  filter(resistance == 0.00)

dn0_vec_pbo <- c(pbo_df$dn0_lo10, pbo_df$dn0_med, pbo_df$dn0_up90)
rn0_vec_pbo <- c(pbo_df$rn0_lo10, pbo_df$rn0_med, pbo_df$rn0_up90)
ento_df_pbo <- expand.grid(bites_Bed = bites_bed_vec,
                           itn_cov = itn_cov_vec,
                           d_ITN0 = dn0_vec_pbo)

ento_df_pbo <- ento_df_pbo %>%
  mutate(r_ITN0 = case_when(d_ITN0 == pbo_df$dn0_lo10 ~ rn0_vec_pbo[1],
                            d_ITN0 == pbo_df$dn0_med ~ rn0_vec_pbo[2],
                            d_ITN0 == pbo_df$dn0_up90 ~ rn0_vec_pbo[3],
                            TRUE ~ NA_real_),
         itn_half_life = case_when(d_ITN0 == pbo_df$dn0_lo10 ~ gamman_vec_pyr[1]*365,
                                   d_ITN0 == pbo_df$dn0_med ~ gamman_vec_pyr[2]*365,
                                   d_ITN0 == pbo_df$dn0_up90 ~ gamman_vec_pyr[3]*365,
                                   TRUE ~ NA_real_))

ento_df_pbo <- ento_df_pbo %>%
  mutate(bites_Indoors = case_when(bites_Bed == bites_bed_vec[1] ~ bites_in_vec[1],
                                   bites_Bed == bites_bed_vec[2] ~ bites_in_vec[2]))



ento_param_list_pbo <- list()
for (i in seq_len(nrow(ento_df_pbo))){
  ento_param_list_pbo[[i]] <- as.numeric(ento_df_pbo[i,])
}

out_nets_pbo <- function(itn_input){
  bites_Bed_in <- itn_input[1]
  itn_cov_in <- itn_input[2]
  d_ITN0_in <- itn_input[3]
  r_ITN0_in <- itn_input[4]
  itn_half_life_in <- itn_input[5]
  bites_In_in<- itn_input[6]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 2,
                      itn_cov = itn_cov_in,
                      ITN_IRS_on = ITN_IRS_on,
                      init_ft = prop_treated,
                      bites_Bed = bites_Bed_in,
                      bites_Indoors = bites_In_in,
                      itn_half_life = itn_half_life_in,
                      d_ITN0 = d_ITN0_in,
                      r_ITN0 = r_ITN0_in)
  return(output)
}


#my_sim_antag_ITN <- function(){
#  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
#  pyr_out_list_antag_ITN <- lapply(pyr_param_list, antag_ITN_cov_loop)
#  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
#  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
#    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
#    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
#                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
#    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
#  return(pyr_out_df_antag_ITN)
#
#}
my_sim_nets_pbo <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  ento_out <- lapply(ento_param_list_pbo, out_nets_pbo)
  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_df_pbo)), function(x){
    df <- as.data.frame(ento_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, itn_cov, EIR_tot, prev,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, bites_Indoors, Q0,inc05,s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pbo"))}, simplify = F))
  return(ento_out_mod)

}

ento_nets_pbo <- my_sim_nets_pbo() #this has been run



ento_nets_pbo <- ento_nets_pbo %>%
  mutate(segments = case_when(bites_Indoors == bites_in_vec[1] & bites_Bed == bites_bed_vec[1] ~ "segment_M",
                              bites_Indoors == bites_in_vec[2] & bites_Bed == bites_bed_vec[2] ~ "segment_L",
                              bites_Indoors == bites_in_vec[3] & bites_Bed == bites_bed_vec[3] ~ "segment_U",
                              bites_Bed == bites_in_vec[4] & bites_Bed == bites_bed_vec[4] ~ "hourly_M",
                              bites_Bed == bites_in_vec[5] & bites_Bed == bites_bed_vec[5] ~ "hourly_L",
                              bites_Bed == bites_in_vec[6] & bites_Bed == bites_bed_vec[6] ~ "hourly_L",
                              TRUE ~ NA_character_))

ento_nets_pbo <- ento_nets_pbo %>%
  mutate(net_eff = case_when(d_ITN0 == dn0_vec_pbo[1] & r_ITN0 == rn0_vec_pbo[1] ~ "estim_lower",
                             d_ITN0 == dn0_vec_pbo[2] & r_ITN0 == rn0_vec_pbo[2] ~ "estim_med",
                             d_ITN0 == dn0_vec_pbo[3] & r_ITN0 == rn0_vec_pbo[3] ~ "estim_upper",
                             TRUE ~ NA_character_))

#IG2
IG2_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_pyrrole_nets.csv")
head(IG2_df)

IG2_df <- IG2_df %>%
  filter(resistance == 0.00)

dn0_vec_IG2 <- c(IG2_df$dn0_lo10, IG2_df$dn0_med, IG2_df$dn0_up90)
rn0_vec_IG2 <- c(IG2_df$rn0_lo10, IG2_df$rn0_med, IG2_df$rn0_up90)
ento_df_IG2 <- expand.grid(bites_Bed = bites_bed_vec,
                           itn_cov = itn_cov_vec,
                           d_ITN0 = dn0_vec_IG2)

ento_df_IG2 <- ento_df_IG2 %>%
  mutate(r_ITN0 = case_when(d_ITN0 == IG2_df$dn0_lo10 ~ rn0_vec_IG2[1],
                            d_ITN0 == IG2_df$dn0_med ~ rn0_vec_IG2[2],
                            d_ITN0 == IG2_df$dn0_up90 ~ rn0_vec_IG2[3],
                            TRUE ~ NA_real_),
         itn_half_life = case_when(d_ITN0 == IG2_df$dn0_lo10 ~ gamman_vec_pyr[1]*365,
                                   d_ITN0 == IG2_df$dn0_med ~ gamman_vec_pyr[2]*365,
                                   d_ITN0 == IG2_df$dn0_up90 ~ gamman_vec_pyr[3]*365,
                                   TRUE ~ NA_real_))

ento_df_IG2 <- ento_df_IG2 %>%
  mutate(bites_Indoors = case_when(bites_Bed == bites_bed_vec[1] ~ bites_in_vec[1],
                                   bites_Bed == bites_bed_vec[2] ~ bites_in_vec[2]))



ento_param_list_IG2 <- list()
for (i in seq_len(nrow(ento_df_IG2))){
  ento_param_list_IG2[[i]] <- as.numeric(ento_df_IG2[i,])
}

out_nets_IG2 <- function(itn_input){
  bites_Bed_in <- itn_input[1]
  itn_cov_in <- itn_input[2]
  d_ITN0_in <- itn_input[3]
  r_ITN0_in <- itn_input[4]
  itn_half_life_in <- itn_input[5]
  bites_In_in<- itn_input[6]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 2,
                      itn_cov = itn_cov_in,
                      ITN_IRS_on = ITN_IRS_on,
                      init_ft = prop_treated,
                      bites_Bed = bites_Bed_in,
                      bites_Indoors = bites_In_in,
                      itn_half_life = itn_half_life_in,
                      d_ITN0 = d_ITN0_in,
                      r_ITN0 = r_ITN0_in)
  return(output)
}


#my_sim_antag_ITN <- function(){
#  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
#  pyr_out_list_antag_ITN <- lapply(pyr_param_list, antag_ITN_cov_loop)
#  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
#  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
#    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
#    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
#                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
#    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
#  return(pyr_out_df_antag_ITN)
#
#}
my_sim_nets_IG2 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  ento_out <- lapply(ento_param_list_IG2, out_nets_IG2)
  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_df_IG2)), function(x){
    df <- as.data.frame(ento_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, itn_cov, EIR_tot, prev,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, bites_Indoors, Q0,inc05,s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "IG2"))}, simplify = F))
  return(ento_out_mod)

}

ento_nets_IG2 <- my_sim_nets_IG2() #this has been run



ento_nets_IG2 <- ento_nets_IG2 %>%
  mutate(segments = case_when(bites_Indoors == bites_in_vec[1] & bites_Bed == bites_bed_vec[1] ~ "segment_M",
                              bites_Indoors == bites_in_vec[2] & bites_Bed == bites_bed_vec[2] ~ "segment_L",
                              bites_Indoors == bites_in_vec[3] & bites_Bed == bites_bed_vec[3] ~ "segment_U",
                              bites_Bed == bites_in_vec[4] & bites_Bed == bites_bed_vec[4] ~ "hourly_M",
                              bites_Bed == bites_in_vec[5] & bites_Bed == bites_bed_vec[5] ~ "hourly_L",
                              bites_Bed == bites_in_vec[6] & bites_Bed == bites_bed_vec[6] ~ "hourly_L",
                              TRUE ~ NA_character_))

ento_nets_IG2 <- ento_nets_IG2 %>%
  mutate(net_eff = case_when(d_ITN0 == dn0_vec_IG2[1] & r_ITN0 == rn0_vec_IG2[1] ~ "estim_lower",
                             d_ITN0 == dn0_vec_IG2[2] & r_ITN0 == rn0_vec_IG2[2] ~ "estim_med",
                             d_ITN0 == dn0_vec_IG2[3] & r_ITN0 == rn0_vec_IG2[3] ~ "estim_upper",
                             TRUE ~ NA_character_))

#baseline scenario
ento_param_list_baseline <- list()
for (i in seq_len(nrow(ento_df_IG2))){
  ento_param_list_baseline[[i]] <- as.numeric(ento_df_IG2[i,]) #some of this gets overwritten
}

out_baseline <- function(itn_input){
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 1,
                      #itn_cov = itn_cov_in,
                      #ITN_IRS_on = ITN_IRS_on,
                      init_ft = prop_treated)
                      #bites_Bed = bites_Bed_in,
                      #bites_Indoors = bites_In_in,
                      #d_ITN0 = d_ITN0_in,
                      #r_ITN0 = r_ITN0_in)
  return(output)
}

my_sim_baseline <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  ento_out <- lapply(ento_param_list_baseline, out_baseline)
  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_df_IG2)), function(x){
    df <- as.data.frame(ento_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, EIR_tot, prev,
                                     Q0, inc05))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "baseline"))}, simplify = F))
  return(ento_out_mod)

}

ento_baseline <- my_sim_baseline()

ento_baseline <- ento_baseline %>%
  mutate(net_eff = "estim_med",
         itn_cov = 0.8) #just for ease to help with facet plot later

ento_nets_data <- list(ento_nets_pyr, ento_nets_pbo, ento_nets_IG2)
ento_nets_data <- do.call("rbind", ento_nets_data)

saveRDS(ento_nets_data, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/ento_nets_data_sampling_methods.rds")

ento_nets_data <- readRDS("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/ento_nets_data_sampling_methods.rds")

ento_nets_data2 <- dplyr::bind_rows(ento_nets_data, ento_baseline)

unique(ento_nets_data2$net_type)
require(tidyverse)
ento_nets_data2 %>%
  filter(net_eff == "estim_med" & itn_cov == 0.8) %>%
  ggplot()+
  aes(x = t/365, y = prev, col = as.factor(segments))+
  geom_line()+
  facet_grid(cols = vars(net_type))


ento_nets_data2 %>%
  filter(net_eff == "estim_med" & itn_cov == 0.8 & segments == "hourly") %>%
  ggplot()+
  aes(x = t/365, y = prev)+
  geom_line()+
  facet_grid(cols = vars(net_type))

ento_nets_data2 %>%
  filter(net_eff == "estim_med" & itn_cov == 0.8 & segments == "hourly") %>%
  ggplot()+
  aes(x = t/365, y = inc05)+
  geom_line(aes(col = as.factor(net_type)))

#for showing the CIs too
#ento_nets_CI_prev_df <- ento_nets_data %>%
#  select(t, itn_cov, bites_Bed, bites_Indoors, segments, prev, net_eff)
#
#ento_nets_CI_prev_df_wide <- spread(ento_nets_CI_prev_df, key = net_eff, value = prev)

#show that the cases averted by the different types of bednet are more different than the estimate from each sampling collection method

cases_baseline <- ento_baseline %>%
  filter(between(t, 365, 365*10)) %>%
  reframe(tot_cases = sum(inc05))

cases_IG2 <- ento_nets_IG2 %>%
  group_by(segments, net_eff) %>%
  filter(between(t, 365, 365*10),
         net_eff == "estim_med") %>%
  summarise(tot_cases = sum(inc05)) %>%
  mutate(total_case_baseline = cases_baseline$tot_cases,
         cases_averted = total_case_baseline - tot_cases,
         cases_averted_percent = cases_averted/total_case_baseline,
         net_type = "IG2")

cases_pbo <- ento_nets_pbo %>%
  group_by(segments, net_eff) %>%
  filter(between(t, 365, 365*10),
         net_eff == "estim_med") %>%
  summarise(tot_cases = sum(inc05))%>%
  mutate(total_case_baseline = cases_baseline$tot_cases,
         cases_averted = total_case_baseline - tot_cases,
         cases_averted_percent = cases_averted/total_case_baseline,
         net_type = "PBO")

cases_pyr <- ento_nets_pyr %>%
  group_by(segments, net_eff) %>%
  filter(between(t, 365, 365*10),
         net_eff == "estim_med") %>% #just taking the median point for efficacy
  summarise(tot_cases = sum(inc05))%>%
  mutate(total_case_baseline = cases_baseline$tot_cases,
         cases_averted = total_case_baseline - tot_cases,
         cases_averted_percent = cases_averted/total_case_baseline,
         net_type = "pyrethroid only")

case_summary_list <- list(cases_IG2, cases_pbo, cases_pyr)
case_summary <- do.call("rbind", case_summary_list)

ggplot(case_summary, aes(x = net_type, y = cases_averted_percent, fill = as.factor(segments)))+
  geom_bar(stat = "identity", position = position_dodge())+
  ylim(0, 1)+
  ylab("Relative difference cases averted due to LLINs \n compared to baseline scenario with no intervention (%)")+
  theme_bw()+
  xlab("Net type")+
  labs(fill = "Mosquito sampling method")+
  theme(legend.position = c(0.8, 0.92))

write.csv(case_summary, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/cases_averted_summary.csv")

ggplot(case_summary, aes(x = net_type, y = cases_averted, fill = as.factor(segments)))+
  geom_bar(stat = "identity", position = position_dodge())+
  ylab("Absolute cases averted due to LLINs compared to baseline scenario with no intervention")

#differences in the cases averted by the different net types

#what is the difference in cases averted (compared to a baseline scenario of no net) of IG2 and PBO relative to pyr, for each data collection method

#1) Compute within difference of cases averted for each net type, acc each sampling method

sampling_summary <- case_summary %>%
  select(segments, cases_averted, net_type) %>%
  pivot_wider(names_from = segments, values_from = cases_averted) %>%
  mutate(abs_diff_cases_av_sampling_method = hourly-segments,
         rel_diff_cases_av_sampling_method = ((hourly-segments)/hourly)*100)

net_summary <- case_summary %>%
  select(segments, cases_averted, net_type) %>%
  filter(net_type != "PBO") %>%
  pivot_wider(names_from = net_type, values_from = cases_averted) %>%
  mutate(abs_diff_net_type = IG2-`pyrethroid only`,
         rel_diff_net_type = ((IG2-`pyrethroid only`)/IG2)*100)

