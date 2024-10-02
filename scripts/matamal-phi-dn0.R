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

#from GeissbulherB: study furthest away from the straight line
bites_in_vec <- c(0.641, 0.512)
bites_bed_vec <- c(0.554, 0.417)
itn_cov_vec <- c(0.2, 0.8)

#get net efficacies
net_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_only_nets.csv")
head(net_df)

net_df <- net_df %>%
  filter(resistance == 0.00)

dn0_vec <- c(net_df$dn0_lo10, net_df$dn0_med, net_df$dn0_up90)
rn0_vec <- c(net_df$rn0_lo10, net_df$rn0_med, net_df$rn0_up90)
ento_df <- expand.grid(bites_Bed = bites_bed_vec,
                       itn_cov = itn_cov_vec,
                       d_ITN0 = dn0_vec)

ento_df <- ento_df %>%
  mutate(r_ITN0 = case_when(d_ITN0 == net_df$dn0_lo10 ~ rn0_vec[1],
                            d_ITN0 == net_df$dn0_med ~ rn0_vec[2],
                            d_ITN0 == net_df$dn0_up90 ~ rn0_vec[3],
                            TRUE ~ NA_real_))

ento_df <- ento_df %>%
  mutate(bites_Indoors = case_when(bites_Bed == bites_bed_vec[1] ~ bites_in_vec[1],
                                   bites_Bed == bites_bed_vec[2] ~ bites_in_vec[2]))






ento_param_list <- list()
for (i in seq_len(nrow(ento_df))){
  ento_param_list[[i]] <- as.numeric(ento_df[i,])
}

out_nets <- function(itn_input){
  bites_Bed_in <- itn_input[1]
  itn_cov_in <- itn_input[2]
  d_ITN0_in <- itn_input[3]
  r_ITN0_in <- itn_input[4]
  bites_In_in<- itn_input[5]
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



my_sim_nets <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  ento_out <- lapply(ento_param_list, out_nets)
  #res_ento_out <- lapply(ento_out, runfun) #put these values into the model
  ento_out_mod <- do.call(rbind, sapply(1:(nrow(ento_df)), function(x){
    df <- as.data.frame(ento_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, itn_cov, EIR_tot, prev,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, bites_Indoors, Q0,s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "nets"))}, simplify = F))
  return(ento_out_mod)

}

ento_nets <- my_sim_nets() #this has been run

head(ento_nets)

ento_nets <- ento_nets %>%
  mutate(segments = case_when(bites_Indoors == bites_in_vec[1] & bites_Bed == bites_bed_vec[1] ~ "segments",
                              TRUE ~ "hourly"))

ento_nets <- ento_nets %>%
  mutate(net_eff = case_when(d_ITN0 == dn0_vec[1] & r_ITN0 == rn0_vec[1] ~ "estim_lower",
                             d_ITN0 == dn0_vec[2] & r_ITN0 == rn0_vec[2] ~ "estim_med",
                             d_ITN0 == dn0_vec[3] & r_ITN0 == rn0_vec[3] ~ "estim_upper",
                             TRUE ~ NA_character_))

#go from wide to long, then can use the lower and upper as bound for uncertainty estimates
ento_nets_med <- ento_nets %>%
  filter(net_eff == "estim_med")

ento_nets_CIs_EIR <- ento_nets %>%
  select(t, itn_cov, bites_Bed, bites_Indoors, segments, EIR_tot, net_eff)

ento_nets_CIs_EIR_wide <- spread(ento_nets_CIs_EIR, key = net_eff, value = EIR_tot)

ento_nets_CIs_prev <- ento_nets %>%
  select(t, itn_cov, bites_Bed, segments, prev, net_eff)

ento_nets_CIs_prev %>%
  ggplot(aes(x = t, y = prev, col = as.factor(segments)))+
  geom_line()+
  facet_wrap(vars(net_eff, itn_cov))

ento_nets_CIs_prev_wide <- spread(ento_nets_CIs_prev, key = net_eff, value = prev)


head(ento_nets)
ento_nets %>%
  filter(between(t, 1800, 2000)) %>%
  select(t, EIR_tot, prev, itn_cov) %>%
  ggplot()+
  aes(x = t, y = EIR_tot)+
  geom_line()+
  facet_grid(~itn_cov, labeller = label_both)

net.labs <- c("ITN coverage = 20%", "ITN coverage = 80%")
names(net.labs) <- c("0.2", "0.8")

#eir spike at 35.505368, t = 2919 for hourly segments and itn_cov = 0.8

ch <- ento_nets_CIs_EIR_wide %>%
  filter(segments == "segments" & itn_cov == 0.8) %>%
  filter(between(t, 2900, 3000))

ggplot(ch, aes(x = t, y = estim_med))+
  geom_line()

check <- ento_nets_CIs_EIR_wide %>%
  filter(itn_cov == 0.8 & segments == "segments") %>%
  filter(between(t, 7.5*365, 9*365)) %>%
  filter(estim_med <= 41.19240)


ggplot(check, aes(x = t, y = estim_med))+
  geom_line()

ento_nets_CIs_EIR_wide_spikerm <- ento_nets_CIs_EIR_wide %>%
  #filter(itn_cov == 0.8 & segments == "hourly") %>%
  filter(estim_med <= 100) %>% # spike here
  mutate(estim_med = case_when(itn_cov == 0.8 & segments == "hourly" & t >2000 & t <3000 & estim_med >35.505368 ~ 35.505368,
                               itn_cov == 0.8 & segments == "segments" & t > 7.5*365 & t < 9*365 & estim_med >= 41.19240 ~ 41.19240, #needs more work
                               TRUE ~ estim_med))

ento_nets_CIs_EIR_wide %>%
  filter(itn_cov == 0.8, segments == "hourly") %>%
  ggplot()+
  aes(x = t/365, y = estim_med)+
  geom_line()+
  ylab("Average annual EIR")

#write.csv(ento_nets_CIs_EIR_wide_spikerm, "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/model_eir_review.csv", row.names = FALSE)

eir_plots <- ggplot(ento_nets_CIs_EIR_wide, aes(x = t/365, y = estim_med, col = as.factor(segments),
                                           ymin = estim_lower, ymax = estim_upper))+
  geom_line(size = 1)+
  facet_grid(~ itn_cov, labeller = labeller(itn_cov = net.labs))+
  geom_smooth(stat = "identity", aes(fill = as.factor(segments)))+
  #guides(fill= "none")+
  theme_classic()+
  theme(legend.position = c(0.25, 0.15))+
  scale_x_continuous(breaks = seq(1, 10, 3))+
  labs(col = "Mosquito collections", fill = "Mosquito collections")+
  ylab("Average annual EIR in adults and children")+
  xlab("Time (years)")


#eir_plots <- ggplot(ento_nets_CIs_EIR_wide_spikerm, aes(x = t/365, y = estim_med, col = as.factor(segments),
#                                                        ymin = estim_lower, ymax = estim_upper))+
#  geom_line(size = 1)+
#  facet_grid(~ itn_cov, labeller = labeller(itn_cov = net.labs))+
#  geom_smooth(stat = "identity", aes(fill = as.factor(segments)))+
#  #guides(fill= "none")+
#  theme_classic()+
#  theme(legend.position = c(0.25, 0.15))+
#  scale_x_continuous(breaks = seq(1, 10, 3))+
#  labs(col = "Mosquito collections", fill = "Mosquito collections")+
#  ylab("Average annual EIR in adults and children")+
#  xlab("Time (years)")


#ento_nets_prev %>%
#  filter(t == 1) %>%
#  select(itn_cov, estim_lower, estim_upper, estim_med) %>%
#  print(n = 100)
#

write.csv(ento_nets_CIs_prev_wide, "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/model_prev_review.csv", row.names = FALSE)


prev_plots <- ggplot(ento_nets_CIs_prev_wide, aes(x = t/365, y = estim_med, col = as.factor(segments),
                                                  ymin = estim_lower, ymax = estim_upper))+
  geom_line(size = 1)+
  facet_grid(~ itn_cov, labeller = labeller(itn_cov = net.labs))+
  geom_smooth(stat = "identity", aes(fill = as.factor(segments)))+
  #guides(fill= "none")+
  theme_classic()+
  theme(legend.position = c(0.3, 0.3))+
  scale_x_continuous(breaks = seq(1, 10, 3))+
  labs(col = "Mosquito collections", fill = "Mosquito collections")+
  guides(fill = "none", colour = "none")+
  ylab("Slide prevalence in \n under 5 year-olds")+
  xlab("Time (years)")+
  ylim(0, 1)

epi_plots <- cowplot::plot_grid(eir_plots, prev_plots, labels = c("A)", "B)"))
ggsave(epi_plots, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/plots/model_plots.png")




#at the start of a new distribution, what is the absolute and relative difference in prevalence and EIR
#no interventions, prev = 0.693 and EIR_tot = 99.99663

no_int <- ento_nets %>%
  filter(t == 1) %>%
  select(EIR_tot, prev) %>%
  rename(eir_no_int = EIR_tot,
         prev_no_int = prev)

eir_no_int <- unique(no_int$eir_no_int)
prev_no_int <- unique(no_int$prev_no_int)

int <- ento_nets %>%
  filter(t == (5*365)+30) %>% #the day after the ITN distrib
  group_by(ITN_cov, segments) %>%
  summarise(EIR_tot = EIR_tot,
            prev = prev) %>%
  mutate(EIR_no_int = eir_no_int,
         PREV_no_int = prev_no_int)

seg_diff <- int %>%
  filter(segments == "segments") %>%
  summarise(abs_diff_eir = EIR_no_int - EIR_tot,
            rel_diff_eir = ((EIR_no_int - EIR_tot)/EIR_no_int)*100,
            abs_diff_prev = PREV_no_int - prev,
            rel_diff_prev = ((PREV_no_int - prev)/PREV_no_int)*100)

hourly_diff <- int %>%
  filter(segments == "hourly") %>%
  summarise(abs_diff_eir = EIR_no_int - EIR_tot,
            rel_diff_eir = ((EIR_no_int - EIR_tot)/EIR_no_int)*100,
            abs_diff_prev = PREV_no_int - prev,
            rel_diff_prev = ((PREV_no_int - prev)/PREV_no_int)*100)

#segmented --> underestimate bites_Indoors and bites_Bed --> underestimate the impact of LLINs. with hourly, nets do more so less room for endectocide.
