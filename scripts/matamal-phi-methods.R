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
init_EIR <- 20 # low transmission setting for Guinea-Bissau

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0

# Define time for turning on interventions
ITN_IRS_on <- 365

#from geissbulherB: study furthest away from the straight line y = x
study_dev <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/GeissbulherB.csv")

#phiBs
phi_B_seg <-  round(study_dev$new_PHI_B_seg,3)
phi_B_seg_L <- round(study_dev$new_PHI_B_seg_L, 3)
phi_B_seg_U <- round(study_dev$new_PHI_B_seg_U, 3)

phi_B_hourly <- round(study_dev$new_PHI_B,3)
phi_B_hourly_L <- round(study_dev$new_PHI_B_L, 3)
phi_B_hourly_U <- round(study_dev$new_PHI_B_U, 3)

#phiIs
phi_I_seg <-  round(study_dev$new_PHI_I_seg,3)
phi_I_seg_L <- round(study_dev$new_PHI_I_seg_L,3)
phi_I_seg_U <- round(study_dev$new_PHI_I_seg_U,3)

phi_I_hourly <- round(study_dev$new_PHI_I,3)
phi_I_hourly_L <- round(study_dev$new_PHI_I_L,3)
phi_I_hourly_U <- round(study_dev$new_PHI_I_U,3)


bites_in_vec <- c(phi_I_seg, phi_I_seg_L, phi_I_seg_U, phi_I_hourly, phi_I_hourly_L, phi_I_hourly_U)
bites_bed_vec <- c(phi_B_seg, phi_B_seg_L, phi_B_seg_U, phi_B_hourly, phi_B_hourly_L, phi_B_hourly_U)

#get net efficacies
pyr_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_only_nets.csv")
head(pyr_df)

pyr_df <- pyr_df %>%
  filter(resistance == 0.55)
pyr_df$g

dn0_vec_pyr <- c(pyr_df$dn0_lo10, pyr_df$dn0_med, pyr_df$dn0_up90)
rn0_vec_pyr <- c(pyr_df$rn0_lo10, pyr_df$rn0_med, pyr_df$rn0_up90)

#ELLIE SAID USE THE GAMMAN PYR-ONLY HALF LIFE FOR ALL OTHER NET TYPES
gamman_vec_pyr <- c(pyr_df$gamman_lo10, pyr_df$gamman_med, pyr_df$gamman_up90)
ento_df_pyr <- expand.grid(bites_Bed = bites_bed_vec,
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
  d_ITN0_in <- itn_input[2]
  r_ITN0_in <- itn_input[3]
  itn_half_life_in <- itn_input[4]
  bites_In_in<- itn_input[5]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 2,
                      itn_cov = 0.8,
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
  mutate(segments = case_when(bites_Bed == bites_bed_vec[1] ~ "segment_M",
                              bites_Bed == bites_bed_vec[2] ~ "segment_L",
                              bites_Bed == bites_bed_vec[3] ~ "segment_U",
                              bites_Bed == bites_bed_vec[4] ~ "hourly_M",
                              bites_Bed == bites_bed_vec[5] ~ "hourly_L",
                              bites_Bed == bites_bed_vec[6] ~ "hourly_U"))


ento_nets_pyr <- ento_nets_pyr %>%
  mutate(net_eff = case_when(d_ITN0 == dn0_vec_pyr[1] & r_ITN0 == rn0_vec_pyr[1] ~ "estim_lower",
                             d_ITN0 == dn0_vec_pyr[2] & r_ITN0 == rn0_vec_pyr[2] ~ "estim_med",
                             d_ITN0 == dn0_vec_pyr[3] & r_ITN0 == rn0_vec_pyr[3] ~ "estim_upper",
                             TRUE ~ NA_character_))
#repeat for pbo
pbo_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_pbo_nets.csv")
head(pbo_df)

pbo_df <- pbo_df %>%
  filter(resistance == 0.55)

dn0_vec_pbo <- c(pbo_df$dn0_lo10, pbo_df$dn0_med, pbo_df$dn0_up90)
rn0_vec_pbo <- c(pbo_df$rn0_lo10, pbo_df$rn0_med, pbo_df$rn0_up90)
ento_df_pbo <- expand.grid(bites_Bed = bites_bed_vec,

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
                                   bites_Bed == bites_bed_vec[2] ~ bites_in_vec[2],
                                   bites_Bed == bites_bed_vec[3] ~ bites_in_vec[3],
                                   bites_Bed == bites_bed_vec[4] ~ bites_in_vec[4],
                                   bites_Bed == bites_bed_vec[5] ~ bites_in_vec[5],
                                   bites_Bed == bites_bed_vec[6] ~ bites_in_vec[6]))



ento_param_list_pbo <- list()
for (i in seq_len(nrow(ento_df_pbo))){
  ento_param_list_pbo[[i]] <- as.numeric(ento_df_pbo[i,])
}

out_nets_pbo <- function(itn_input){
  bites_Bed_in <- itn_input[1]
  d_ITN0_in <- itn_input[2]
  r_ITN0_in <- itn_input[3]
  itn_half_life_in <- itn_input[4]
  bites_In_in<- itn_input[5]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 2,
                      itn_cov = 0.8,
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
  mutate(segments = case_when(bites_Bed == bites_bed_vec[1] ~ "segment_M",
                              bites_Bed == bites_bed_vec[2] ~ "segment_L",
                              bites_Bed == bites_bed_vec[3] ~ "segment_U",
                              bites_Bed == bites_bed_vec[4] ~ "hourly_M",
                              bites_Bed == bites_bed_vec[5] ~ "hourly_L",
                              bites_Bed == bites_bed_vec[6] ~ "hourly_U"))

ento_nets_pbo <- ento_nets_pbo %>%
  mutate(net_eff = case_when(d_ITN0 == dn0_vec_pbo[1] & r_ITN0 == rn0_vec_pbo[1] ~ "estim_lower",
                             d_ITN0 == dn0_vec_pbo[2] & r_ITN0 == rn0_vec_pbo[2] ~ "estim_med",
                             d_ITN0 == dn0_vec_pbo[3] & r_ITN0 == rn0_vec_pbo[3] ~ "estim_upper",
                             TRUE ~ NA_character_))

#IG2
IG2_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_pyrrole_nets.csv")
head(IG2_df)

IG2_df <- IG2_df %>%
  filter(resistance == 0.55)

dn0_vec_IG2 <- c(IG2_df$dn0_lo10, IG2_df$dn0_med, IG2_df$dn0_up90)
rn0_vec_IG2 <- c(IG2_df$rn0_lo10, IG2_df$rn0_med, IG2_df$rn0_up90)
ento_df_IG2 <- expand.grid(bites_Bed = bites_bed_vec,

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
                                   bites_Bed == bites_bed_vec[2] ~ bites_in_vec[2],
                                   bites_Bed == bites_bed_vec[3] ~ bites_in_vec[3],
                                   bites_Bed == bites_bed_vec[4] ~ bites_in_vec[4],
                                   bites_Bed == bites_bed_vec[5] ~ bites_in_vec[5],
                                   bites_Bed == bites_bed_vec[6] ~ bites_in_vec[6]))



ento_param_list_IG2 <- list()
for (i in seq_len(nrow(ento_df_IG2))){
  ento_param_list_IG2[[i]] <- as.numeric(ento_df_IG2[i,])
}

out_nets_IG2 <- function(itn_input){
  bites_Bed_in <- itn_input[1]
  d_ITN0_in <- itn_input[2]
  r_ITN0_in <- itn_input[3]
  itn_half_life_in <- itn_input[4]
  bites_In_in<- itn_input[5]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 2,
                      itn_cov = 0.8,
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
  mutate(segments = case_when(bites_Bed == bites_bed_vec[1] ~ "segment_M",
                              bites_Bed == bites_bed_vec[2] ~ "segment_L",
                              bites_Bed == bites_bed_vec[3] ~ "segment_U",
                              bites_Bed == bites_bed_vec[4] ~ "hourly_M",
                              bites_Bed == bites_bed_vec[5] ~ "hourly_L",
                              bites_Bed == bites_bed_vec[6] ~ "hourly_U"))

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
         itn_cov = 0.8,
         segment_type = "baseline") #just for ease to help with facet plot later

ento_nets_data <- list(ento_nets_pyr, ento_nets_pbo, ento_nets_IG2)
ento_nets_data <- do.call("rbind", ento_nets_data)

saveRDS(ento_nets_data, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/ento_nets_data_sampling_methods.rds")

ento_nets_data <- readRDS("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/ento_nets_data_sampling_methods.rds")

#ento_nets_data2 <- dplyr::bind_rows(ento_nets_data, ento_baseline)

unique(ento_nets_data2$segments)
require(tidyverse)

ento_nets_data2 <- ento_nets_data %>%
  filter(net_eff == "estim_med") %>%
  mutate(segment_type = case_when(grepl("segment", segments) ~ "segment",
                                  grepl("hourly", segments) ~ "hourly"),
         segment_uncertainty = case_when(grepl("_M", segments) ~ "median",
                                         grepl("_L", segments) ~ "lower",
                                          grepl("_U", segments) ~ "upper"),
         net_type = case_when(net_type == "IG2" ~ "IG2",
                              net_type == "pbo" ~ "PBO",
                              net_type == "pyr_only" ~ "Pyrethroid-only"))

#ento_nets_data3 <- ento_nets_data2 %>%
#  select(t, inc05, segment_type, segment_uncertainty, net_type) %>%
#  pivot_wider(names_from = "segment_uncertainty", values_from  = "inc05"
             # )
#write.csv(ento_nets_data3, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/incidence_dynamics_study_dev.csv")

ento_nets_data3 <- ento_nets_data2 %>%
  filter(segment_uncertainty == "median", segment_type == "segment") %>% #take epi estimates derived from the median phi-B estimate
  select(t, inc05,net_type)

ento_baseline_rbind <- ento_baseline %>%
  mutate(net_type = "baseline") %>%
  select(t, inc05, net_type)

nets_data <- rbind(ento_nets_data3, ento_baseline_rbind)
write.csv(nets_data,"C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/net_efficacy_phi_sampling.csv")

#incidence for each net type, as derived from the segmented collection method
inc_dynamics_plot <- ggplot(nets_data, aes(x = t/365, y = inc05, col = net_type))+
  geom_line()+
  theme_bw()+
  labs(col = "Intervention scenario")+
  ylab("Clinical incidence in children under 5 years-old")+
  xlab("Time (years)")+
  theme(legend.position = c(0.3, 0.8))

#inc_dynamics_plot <- ggplot(ento_nets_data3, aes(x = t/365, y = median, col = as.factor(segment_type),
#                                    ymin = lower, ymax = upper))+
#  geom_line(size = 1)+
#  facet_grid(~ net_type)+
#  geom_smooth(stat = "identity", aes(fill = as.factor(segment_type)))+
#  #guides(fill= "none")+
#  theme_classic()+
#  #theme(legend.position = c(0.3, 0.3))+
#  #scale_x_continuous(breaks = seq(1, 10, 3))+
#  labs(col = "Mosquito collections", fill = "Mosquito collections")+
#  #guides(fill = "none", colour = "none")+
#  ylab("Clinical incidence in \n under 5 year-olds")+
#  xlab("Time (years)")+
#  guides(col = "none", fill = "none")


#ggplot(ento_nets_data3, aes(x = t/365, y = median, col = as.factor(segment_type),
#                            ymin = lower, ymax = upper))+
#  geom_line(size = 1)+
#  facet_wrap(vars(net_type, segment_type))+
#  geom_smooth(stat = "identity", aes(fill = as.factor(segment_type)))+
#  #guides(fill= "none")+
#  theme_classic()+
#  #theme(legend.position = c(0.3, 0.3))+
#  #scale_x_continuous(breaks = seq(1, 10, 3))+
#  labs(col = "Mosquito collections", fill = "Mosquito collections")+
#  #guides(fill = "none", colour = "none")+
#  ylab("incidence in \n under 5 year-olds")+
#  xlab("Time (years)")

#ggplot(ento_nets_data4, aes(x = t, y = median))+
#  geom_point()
#
#ento_nets_data2 %>%
#  filter(net_eff == "estim_med" & itn_cov == 0.8) %>%
#  ggplot()+
#  aes(x = t/365, y = prev, col = as.factor(segments))+
#  geom_line()+
#  facet_grid(cols = vars(net_type))
#
#
#ento_nets_data2 %>%
#  filter(net_eff == "estim_med" & itn_cov == 0.8 & segments == "hourly") %>%
#  ggplot()+
#  aes(x = t/365, y = prev)+
#  geom_line()+
#  facet_grid(cols = vars(net_type))
#
#ento_nets_data2 %>%
#  filter(net_eff == "estim_med" & itn_cov == 0.8 & segments == "hourly") %>%
#  ggplot()+
#  aes(x = t/365, y = inc05)+
#  geom_line(aes(col = as.factor(net_type)))

#for showing the CIs too
#ento_nets_CI_prev_df <- ento_nets_data %>%
#  select(t, itn_cov, bites_Bed, bites_Indoors, segments, prev, net_eff)
#
#ento_nets_CI_prev_df_wide <- spread(ento_nets_CI_prev_df, key = net_eff, value = prev)

#show that the cases averted by the different types of bednet are more different than the estimate from each sampling collection method

cases_baseline <- ento_baseline %>%
  filter(between(t, 365, 365*10)) %>%
  reframe(tot_cases = sum(inc05)) #384 cases at baseline

cases_IG2 <- ento_nets_IG2 %>%
  group_by(segments, net_eff) %>%
  filter(between(t, 365, 365*10),
         net_eff == "estim_med",
         itn_cov == 0.8) %>%
  summarise(tot_cases = sum(inc05)) %>%
  mutate(total_case_baseline = cases_baseline$tot_cases,
         cases_averted = total_case_baseline - tot_cases,
         cases_averted_percent = (cases_averted/total_case_baseline)*100,
         net_type = "IG2")

cases_pbo <- ento_nets_pbo %>%
  group_by(segments, net_eff) %>%
  filter(between(t, 365, 365*10),
         net_eff == "estim_med",
         itn_cov == 0.8) %>%
  summarise(tot_cases = sum(inc05))%>%
  mutate(total_case_baseline = cases_baseline$tot_cases,
         cases_averted = total_case_baseline - tot_cases,
         cases_averted_percent = (cases_averted/total_case_baseline)*100,
         net_type = "PBO")

cases_pyr <- ento_nets_pyr %>%
  group_by(segments, net_eff) %>%
  filter(between(t, 365, 365*10),
         net_eff == "estim_med",
         itn_cov == 0.8) %>% #just taking the median point for efficacy
  summarise(tot_cases = sum(inc05))%>%
  mutate(total_case_baseline = cases_baseline$tot_cases,
         cases_averted = total_case_baseline - tot_cases,
         cases_averted_percent = (cases_averted/total_case_baseline)*100,
         net_type = "pyrethroid only")

case_summary_list <- list(cases_IG2, cases_pbo, cases_pyr)
case_summary <- do.call("rbind", case_summary_list)

case_summary_wide <- case_summary %>%
  mutate(segment_type = case_when(grepl("segment", segments) ~ "segment",
                                  grepl("hourly", segments) ~ "hourly"),
         segment_uncertainty = case_when(grepl("_M", segments) ~ "median",
                                         grepl("_L", segments) ~ "lower",
                                         grepl("_U", segments) ~ "upper"),
         net_type = case_when(net_type == "IG2" ~ "IG2",
                              net_type == "PBO" ~ "PBO",
                              net_type == "pyrethroid only" ~ "Pyrethroid-only"))
case_summary_rel<- case_summary_wide %>%
  ungroup() %>%
  select(segment_type, segment_uncertainty,cases_averted_percent, net_type) %>%
  pivot_wider(names_from = segment_uncertainty, values_from = cases_averted_percent)

write.csv(case_summary_rel, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/case_summary_rel.csv")
cases_averted_percent <- ggplot(case_summary_rel,aes(x = net_type, y = median,col = as.factor(segment_type)))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5),
                width = 0.2)+
  ylab("Cases averted (%) due to LLINs \n compared to baseline scenario of no intervention")+
  theme_bw()+
  xlab("Net type")+
  labs(col = "Mosquito sampling method")+
  theme(legend.position = c(0.8,0.8))

case_summary_abs<- case_summary_wide %>%
  ungroup() %>%
  select(segment_type, segment_uncertainty,cases_averted, net_type) %>%
  pivot_wider(names_from = segment_uncertainty, values_from = cases_averted)

ggplot(case_summary_abs,aes(x = net_type, y = median,col = as.factor(segment_type)))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5),
                width = 0.2)+
  ylab("Abs difference cases averted due to LLINs \n compared to baseline scenario with no intervention")+
  theme_bw()+
  xlab("Net type")+
  labs(fill = "Mosquito sampling method")





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

