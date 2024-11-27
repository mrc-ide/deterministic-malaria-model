#modelling impact of LLINs in different phi-B settings in the Bijagos islands
#differences in phi-B are written by different assumptions on human sleeping patterns

# A) make a dynamics with the median estimate from each scenario (incidence in under 5s)
# B) make a plot with cases averted from each scenario

#show the extent to which assumptions about sleeping patterns (and how this can change throughout the year) may influence the results

require(ICDMM)
require(tidyverse)
phi_vals <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/phi_calc_stan.csv")


# load the model package (or have it built and reloaded as described above)

# create a vector of age categories
init_age <- c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)

# provide a value of the annual EIR for this model run
init_EIR <- 30 # low transmission setting for Guinea-Bissau

# provide the length of time (in days) that you want to run the model for
time_period <- 365*15

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0

# Define time for turning on interventions
ITN_IRS_on <- (365*5) + (6*30) #5.5 years


#get net efficacies
pyr_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_only_nets.csv")
head(pyr_df)

pyr_df <- pyr_df %>%
  filter(resistance == 0.55)
pyr_df$g

dn0_vec_pyr <- c(pyr_df$dn0_lo10, pyr_df$dn0_med, pyr_df$dn0_up90)
rn0_vec_pyr <- c(pyr_df$rn0_lo10, pyr_df$rn0_med, pyr_df$rn0_up90)

bites_bed_vec <- c(phi_vals$phi_mean_round, phi_vals$phi_lower_round, phi_vals$phi_upper_round)


#first get model in the right space

init_EIR_in <- c(15, 20, 30)
bites_bed_in <- phi_vals$phi_mean_round[1]
gamman_in <- pyr_df$gamman_med*365
rn0_in <- pyr_df$rn0_med
dn0_in <- pyr_df$dn0_med

space_mod <- data.frame(init_EIR_in, bites_bed_in, gamman_in, rn0_in, dn0_in)
names(space_mod) <- c("init_EIR", "bites_Bed", "itn_half_life", "r_ITN0", "d_ITN0")

space_mod_list <- list()

for (i in seq_len(nrow(space_mod))){
  space_mod_list[[i]] <- as.numeric(space_mod[i,])
}

out_space_mod <- function(itn_input){
  init_EIR_in <- itn_input[1]
  bites_Bed_in <- itn_input[2]
  itn_half_life_in <- itn_input[3]
  r_ITN0_in <- itn_input[4]
  d_ITN0_in <- itn_input[5]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR_in,
                      num_int = 2,
                      itn_cov = 0.86,
                      ITN_IRS_on = ITN_IRS_on,
                      init_ft = prop_treated,
                      bites_Bed = bites_Bed_in,
                      itn_half_life = itn_half_life_in,
                      d_ITN0 = d_ITN0_in,
                      r_ITN0 = r_ITN0_in,
                      country = "Guinea_Bissau",
                      admin2 = "Boloma")
  return(output)
}

my_sim_space_mod <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  ento_out <- lapply(space_mod_list, out_space_mod)
  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
  ento_out_mod <- do.call(rbind, sapply(1:(nrow(space_mod)), function(x){
    df <- as.data.frame(ento_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, itn_cov, EIR_tot, prev,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, bites_Indoors, Q0,inc05,s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "space_mod"))}, simplify = F))
  return(ento_out_mod)

}

ento_space_mod <- my_sim_space_mod()
#
#
ggplot(ento_space_mod, aes(x = t, y  = prev))+
  geom_point()+
  facet_wrap(vars(ref))+
  geom_vline(aes(xintercept = (12*365)+60), col = "red")+ #time of survey
  geom_vline(aes(xintercept = ITN_IRS_on), col = "blue")+ #time nets go on
  ylim(0, 0.8)+
  geom_hline(aes(yintercept = 0.06))

#ggplot(ento_space_mod, aes(x = t, y  = mv))+
#  geom_point()+
#  facet_wrap(vars(ref))+
#  geom_vline(aes(xintercept = 365*8.5), col = "red")+
#  geom_vline(aes(xintercept = 365), col = "blue")#+
#  #ylim(0, 1)+
#  #geom_hline(aes(yintercept = 0.06))
#

###based on this, use EIR of 30
init_EIR <- 30


#ELLIE SAID USE THE GAMMAN PYR-ONLY HALF LIFE FOR ALL OTHER NET TYPES
gamman_vec_pyr <- c(pyr_df$gamman_lo10, pyr_df$gamman_med, pyr_df$gamman_up90)
ento_df_pyr <- expand.grid(bites_Bed = bites_bed_vec,
                           d_ITN0 = dn0_vec_pyr[2])

ento_df_pyr <- ento_df_pyr %>%
  mutate(r_ITN0 = case_when(d_ITN0 == pyr_df$dn0_med ~ rn0_vec_pyr[2],
                            TRUE ~ NA_real_),
         itn_half_life = case_when(d_ITN0 == pyr_df$dn0_med ~ gamman_vec_pyr[2]*365,
                                   TRUE ~ NA_real_))
ento_param_list_pyr <- list()
for (i in seq_len(nrow(ento_df_pyr))){
  ento_param_list_pyr[[i]] <- as.numeric(ento_df_pyr[i,])
}

out_nets_pyr <- function(itn_input){
  bites_Bed_in <- itn_input[1]
  d_ITN0_in <- itn_input[2]
  r_ITN0_in <- itn_input[3]
  itn_half_life_in <- itn_input[4]
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 2,
                      itn_cov = 0.86,
                      ITN_IRS_on = ITN_IRS_on,
                      init_ft = prop_treated,
                      bites_Bed = bites_Bed_in,
                      itn_half_life = itn_half_life_in,
                      d_ITN0 = d_ITN0_in,
                      r_ITN0 = r_ITN0_in,
                      country = "Guinea-Bissau",
                      admin2 = "Bolama")
  return(output)
}

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
  mutate(scenario = case_when(bites_Bed == bites_bed_vec[1] ~ "data_sleep_M",
                              bites_Bed == bites_bed_vec[2] ~ "80_sleep_M",
                              bites_Bed == bites_bed_vec[3] ~ "60_sleep_M",
                              bites_Bed == bites_bed_vec[4] ~ "data_sleep_L",
                              bites_Bed == bites_bed_vec[5] ~ "80_sleep_L",
                              bites_Bed == bites_bed_vec[6] ~ "60_sleep_L",
                              bites_Bed == bites_bed_vec[7] ~ "data_sleep_U",
                              bites_Bed == bites_bed_vec[8] ~ "80_sleep_U",
                              bites_Bed == bites_bed_vec[9] ~ "60_sleep_U"))

##pbo
#pbo_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_pbo_nets.csv")
#
#pbo_df <- pbo_df %>%
#  filter(resistance == 0.55)
#
#dn0_vec_pbo <- c(pbo_df$dn0_lo10, pbo_df$dn0_med, pbo_df$dn0_up90)
#rn0_vec_pbo <- c(pbo_df$rn0_lo10, pbo_df$rn0_med, pbo_df$rn0_up90)
#
#bites_bed_vec <- c(phi_vals$median, phi_vals$lower, phi_vals$upper)
#
##ELLIE SAID USE THE GAMMAN PYR-ONLY HALF LIFE FOR ALL OTHER NET TYPES
#ento_df_pbo <- expand.grid(bites_Bed = bites_bed_vec,
#                           d_ITN0 = dn0_vec_pbo[2])
#
#ento_df_pbo <- ento_df_pbo %>%
#  mutate(r_ITN0 = case_when(d_ITN0 == pbo_df$dn0_med ~ rn0_vec_pbo[2],
#                            TRUE ~ NA_real_),
#         itn_half_life = case_when(d_ITN0 == pbo_df$dn0_med ~ gamman_vec_pyr[2]*365,
#                                   TRUE ~ NA_real_))
#ento_param_list_pbo <- list()
#for (i in seq_len(nrow(ento_df_pbo))){
#  ento_param_list_pbo[[i]] <- as.numeric(ento_df_pbo[i,])
#}
#
#out_nets_pbo <- function(itn_input){
#  bites_Bed_in <- itn_input[1]
#  d_ITN0_in <- itn_input[2]
#  r_ITN0_in <- itn_input[3]
#  itn_half_life_in <- itn_input[4]
#  output <- run_model(het_brackets = 5,
#                      age = init_age,
#                      time = time_period,
#                      init_EIR = init_EIR,
#                      num_int = 2,
#                      itn_cov = 0.8,
#                      ITN_IRS_on = ITN_IRS_on,
#                      init_ft = prop_treated,
#                      bites_Bed = bites_Bed_in,
#                      itn_half_life = itn_half_life_in,
#                      d_ITN0 = d_ITN0_in,
#                      r_ITN0 = r_ITN0_in,
#                      country = "Guinea-Bissau",
#                      admin2 = "Bolama")
#  return(output)
#}
#
#my_sim_nets_pbo <- function(){
#  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
#  ento_out <- lapply(ento_param_list_pbo, out_nets_pbo)
#  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
#  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_df_pbo)), function(x){
#    df <- as.data.frame(ento_out[[x]])
#    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, itn_cov, EIR_tot, prev,
#                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, bites_Indoors, Q0,inc05,s_ITN, d_ITN, r_ITN))
#    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "PBO"))}, simplify = F))
#  return(ento_out_mod)
#
#}
#
#ento_nets_pbo <- my_sim_nets_pbo()
#
#ento_nets_pbo <- ento_nets_pbo %>%
#  mutate(scenario = case_when(bites_Bed == bites_bed_vec[1] ~ "data_sleep_M",
#                              bites_Bed == bites_bed_vec[2] ~ "80_sleep_M",
#                              bites_Bed == bites_bed_vec[3] ~ "60_sleep_M",
#                              bites_Bed == bites_bed_vec[4] ~ "data_sleep_L",
#                              bites_Bed == bites_bed_vec[5] ~ "80_sleep_L",
#                              bites_Bed == bites_bed_vec[6] ~ "60_sleep_L",
#                              bites_Bed == bites_bed_vec[7] ~ "data_sleep_U",
#                              bites_Bed == bites_bed_vec[8] ~ "80_sleep_U",
#                              bites_Bed == bites_bed_vec[9] ~ "60_sleep_U"))
#
##IG2
#IG2_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_pyrrole_nets.csv")
#
#IG2_df <- IG2_df %>%
#  filter(resistance == 0.55)
#
#dn0_vec_IG2 <- c(IG2_df$dn0_lo10, IG2_df$dn0_med, IG2_df$dn0_up90)
#rn0_vec_IG2 <- c(IG2_df$rn0_lo10, IG2_df$rn0_med, IG2_df$rn0_up90)


##ELLIE SAID USE THE GAMMAN PYR-ONLY HALF LIFE FOR ALL OTHER NET TYPES
#ento_df_IG2 <- expand.grid(bites_Bed = bites_bed_vec,
#                           d_ITN0 = dn0_vec_IG2[2])
#
#ento_df_IG2 <- ento_df_IG2 %>%
#  mutate(r_ITN0 = case_when(d_ITN0 == IG2_df$dn0_med ~ rn0_vec_IG2[2],
#                            TRUE ~ NA_real_),
#         itn_half_life = case_when(d_ITN0 == IG2_df$dn0_med ~ gamman_vec_pyr[2]*365,
#                                   TRUE ~ NA_real_))
#ento_param_list_IG2 <- list()
#for (i in seq_len(nrow(ento_df_IG2))){
#  ento_param_list_IG2[[i]] <- as.numeric(ento_df_IG2[i,])
#}
#
#out_nets_IG2 <- function(itn_input){
#  bites_Bed_in <- itn_input[1]
#  d_ITN0_in <- itn_input[2]
#  r_ITN0_in <- itn_input[3]
#  itn_half_life_in <- itn_input[4]
#  output <- run_model(het_brackets = 5,
#                      age = init_age,
#                      time = time_period,
#                      init_EIR = init_EIR,
#                      num_int = 2,
#                      itn_cov = 0.8,
#                      ITN_IRS_on = ITN_IRS_on,
#                      init_ft = prop_treated,
#                      bites_Bed = bites_Bed_in,
#                      itn_half_life = itn_half_life_in,
#                      d_ITN0 = d_ITN0_in,
#                      r_ITN0 = r_ITN0_in,
#                      country = "Guinea-Bissau",
#                      admin2 = "Bolama")
#  return(output)
#}
#
#my_sim_nets_IG2 <- function(){
#  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
#  ento_out <- lapply(ento_param_list_IG2, out_nets_IG2)
#  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
#  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_df_IG2)), function(x){
#    df <- as.data.frame(ento_out[[x]])
#    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, itn_cov, EIR_tot, prev,
#                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, bites_Indoors, Q0,inc05,s_ITN, d_ITN, r_ITN))
#    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "IG2"))}, simplify = F))
#  return(ento_out_mod)
#
#}
#
#ento_nets_IG2 <- my_sim_nets_IG2()
#
#ento_nets_IG2 <- ento_nets_IG2 %>%
#  mutate(scenario = case_when(bites_Bed == bites_bed_vec[1] ~ "data_sleep_M",
#                              bites_Bed == bites_bed_vec[2] ~ "80_sleep_M",
#                              bites_Bed == bites_bed_vec[3] ~ "60_sleep_M",
#                              bites_Bed == bites_bed_vec[4] ~ "data_sleep_L",
#                              bites_Bed == bites_bed_vec[5] ~ "80_sleep_L",
#                              bites_Bed == bites_bed_vec[6] ~ "60_sleep_L",
#                              bites_Bed == bites_bed_vec[7] ~ "data_sleep_U",
#                              bites_Bed == bites_bed_vec[8] ~ "80_sleep_U",
#                              bites_Bed == bites_bed_vec[9] ~ "60_sleep_U"))
#
#baseline
#baseline scenario
init_EIR_in <- 30
ento_baseline_df <- data.frame(init_EIR_in)
ento_param_list_baseline <- list()
for (i in seq_len(nrow(ento_baseline_df))){
  ento_param_list_baseline[[i]] <- as.numeric(ento_baseline_df[i,]) #some of this gets overwritten
}

out_baseline <- function(itn_input){
  output <- run_model(het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = init_EIR,
                      num_int = 1,
                      #itn_cov = itn_cov_in,
                      #ITN_IRS_on = ITN_IRS_on,
                      init_ft = prop_treated,
                      country = "Guinea-Bissau",
                      admin2 = "Bolama")
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
  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_baseline_df)), function(x){
    df <- as.data.frame(ento_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, EIR_tot, prev,
                                       Q0, inc05))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "baseline"))}, simplify = F))
  return(ento_out_mod)

}

ento_baseline <- my_sim_baseline()

cases_baseline <- ento_baseline %>%
  filter(between(t, 365*11.5, 365*14.5)) %>%
  reframe(tot_cases = sum(inc05)) #4.210176 per person

ggplot(ento_baseline, aes(x = t, y = prev))+
  geom_line()+
  ylim(0, 1)

ggplot(ento_baseline, aes(x = t, y = inc05))+
  geom_line()


ento_baseline <- ento_baseline %>%
  mutate(
         itn_cov = 0.86,
         sleeping_scenario = "No interventions",
         net_type = "Baseline (no interventions)") #just for ease to help with facet plot later

ento_nets_data <- ento_nets_pyr

saveRDS(ento_nets_data, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/bijagos_scenarios.rds")
require(tidyverse)
ento_nets_data <- readRDS("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/bijagos_scenarios.rds")

ento_nets_data2 <- ento_nets_data %>%
  filter(net_type == "pyr_only") %>%
  mutate(sleeping_scenario = case_when(grepl("data", scenario) ~ "Observed data (100%) on population in bed at night",
                                  grepl("80", scenario) ~ "80% population in bed at night",
                                  grepl("60", scenario) ~ "60% population in bed at night"),
         phi_uncertainty = case_when(grepl("_M", scenario) ~ "mean",
                                         grepl("_L", scenario) ~ "lower",
                                         grepl("_U", scenario) ~ "upper"),
         net_type = case_when(
                              net_type == "pyr_only" ~ "Pyrethroid-only"))

ento_nets_data3 <- ento_nets_data2 %>%
  filter(phi_uncertainty == "mean") %>% #take epi estimates derived from the median phi-B estimate
  select(t, inc05,net_type, sleeping_scenario)



ento_baseline_rbind <- ento_baseline %>%
  mutate(
    itn_cov = 0.86,
    sleeping_scenario = "No interventions",
    net_type = "Baseline (no interventions)") %>% #just for ease to help with facet plot later
  select(t, inc05, net_type, sleeping_scenario)

nets_data_phi <- rbind(ento_nets_data3, ento_baseline_rbind) %>%
  rename(LLIN = net_type)

sleeping_pals <- c('#bebada','#fb8072','#80b1d3','#fdb462')
time_labs <- c("0.5", "1.5", "2.5")
breaks <- c(12, 13, 14)

write.csv(nets_data_phi, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/nets_data_phi.csv", row.names = FALSE)

nets_data_phi <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/nets_data_phi.csv")
phi_scenario_plots <- nets_data_phi %>%
  filter(between(t, 365*11.5, 365*14.5)) %>%
  ggplot()+
  aes(x = t/365, y = inc05*1000, col = sleeping_scenario)+
  geom_line(size = 2)+
  theme_bw()+
  scale_colour_manual(values = sleeping_pals, name = "Scenario")+
  ylab("Clinical incidence in children under 5 years-old per 1000 persons")+
  scale_x_continuous(labels = time_labs, name = "Time (years) since previous ITN distribution",
                     breaks = breaks)+
  xlab("Time (years)")+
  theme(legend.position = c(0.5, 0.8),
        legend.direction = "vertical",
        text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=14), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=14), #change font size of legend text
        legend.title=element_text(size=14),
        axis.ticks.length = unit(5, "pt"))+
  ylim(0, 25)

ggsave(phi_scenario_plots, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/plots/phi_scenario_bijagos.png")


cases_baseline <- ento_baseline %>%
  filter(between(t, 365*11.5, 365*14.5)) %>%
  reframe(tot_cases = sum(inc05)) #4.210176

#cases_IG2 <- ento_nets_IG2 %>%
#  group_by(scenario) %>%
#  filter(between(t, 365*11.5, 365*14.5),
#         itn_cov == 0.8) %>%
#  summarise(tot_cases = sum(inc05)) %>%
#  mutate(total_case_baseline = cases_baseline$tot_cases,
#         cases_averted = total_case_baseline - tot_cases,
#         cases_averted_percent = (cases_averted/total_case_baseline)*100,
#         net_type = "IG2")
#
#cases_pbo <- ento_nets_pbo %>%
#  group_by(scenario) %>%
#  filter(between(t, 365*11.5, 365*14.5),
#         itn_cov == 0.8) %>%
#  summarise(tot_cases = sum(inc05))%>%
#  mutate(total_case_baseline = cases_baseline$tot_cases,
#         cases_averted = total_case_baseline - tot_cases,
#         cases_averted_percent = (cases_averted/total_case_baseline)*100,
#         net_type = "PBO")

cases_pyr <- ento_nets_data2 %>%
  group_by(scenario) %>%
  filter(between(t, 365*11.5, 365*14.5), #between y8 and y10
         itn_cov == 0.86) %>%
  summarise(tot_cases = sum(inc05))%>%
  mutate(total_case_baseline = cases_baseline$tot_cases,
         cases_averted = total_case_baseline - tot_cases,
         cases_averted_percent = (cases_averted/total_case_baseline)*100,
         net_type = "pyrethroid only")

#case_summary_list <- list(cases_IG2, cases_pbo, cases_pyr)
case_summary <-cases_pyr

case_summary_wide <- case_summary %>%
  mutate(sleeping_scenario = case_when(grepl("data", scenario) ~ "Observed data on population in bed at night",
                                       grepl("80", scenario) ~ "80% population in bed at night",
                                       grepl("60", scenario) ~ "60% population in bed at night"),
         phi_uncertainty = case_when(grepl("_M", scenario) ~ "mean",
                                     grepl("_L", scenario) ~ "lower",
                                     grepl("_U", scenario) ~ "upper"))

case_summary_rel<- case_summary_wide %>%
  group_by(sleeping_scenario, net_type) %>%
  select(sleeping_scenario, phi_uncertainty,cases_averted_percent, net_type) %>%
  pivot_wider(names_from = phi_uncertainty, values_from = cases_averted_percent)

write.csv(case_summary_rel, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/case_summary_rel_bijagos_scenario.csv")

sleeping_pals2 <- c(sleeping_pals[1], sleeping_pals[2], sleeping_pals[4])

case_summary_rel <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/case_summary_rel_bijagos_scenario.csv")

cases_averted_percent <- ggplot(case_summary_rel,aes(x = sleeping_scenario, y = mean,fill = as.factor(sleeping_scenario)))+
  #geom_point(position = position_dodge(width = 0.5),
  #           shape = 4,
  #           size = 5)+
  geom_bar(stat = "identity", position = "dodge", width = 0.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5),
                width = 0.2, size = 1)+
  scale_fill_manual(values = sleeping_pals2, name = "Sleeping pattern scenario")+
  ylab("Cases averted (%) in children under 5-years-old due to ITNs \n compared to baseline scenario of no intervention")+
  theme_bw()+
  xlab("Net type")+
  scale_x_discrete(labels = c("60% population \nin bed at night", "80% population \nin bed at night", "100% population \nin bed at night"), name = "Scenario")+
  guides(fill = "none")+
  ylim(0, 100)+
  theme(legend.position = c(0.3, 0.8),
        legend.direction = "vertical",
        text=element_text(size=12), #change font size of all text
        axis.text=element_text(size=14,), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        plot.title=element_text(size=14), #change font size of plot title
        legend.text=element_text(size=14), #change font size of legend text
        legend.title=element_text(size=14),
        axis.ticks.length = unit(5, "pt"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

phi_scenario_plots_panel <- cowplot::plot_grid(phi_scenario_plots, cases_averted_percent,
                                         labels = c("A", "B"))

ggsave(phi_scenario_plots_panel, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/plots/phi_scenario_bijagos_panel.pdf")
