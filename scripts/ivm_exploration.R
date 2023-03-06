#Exploring the ivermectin model####

#have created an IVM model with x3 SEI compartments, depending on whether the mosquito gets IVM from human/cattle/none at all
require(ICDMM)
require(tidyverse)
require(cowplot)

#base model
out_0 <- run_model(model = "odin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                    time = 730)

out_0_df <- as.data.frame(out_0) %>%
  select(t, Sv, Ev, Iv) %>%
  rowwise() %>%
  mutate(total_mosq = sum(Sv, Ev, Iv))

out_0_df_prop <- out_0_df %>%
  rowwise() %>%
  mutate(Sv = Sv/total_mosq,
         Ev = Ev/total_mosq,
         Iv = Iv/total_mosq)

out_0_df_long <- gather(out_0_df_prop, state_var, prop_mosq, Sv:Iv, factor_key=TRUE)

s_comp <- c("Sv", "Svic", "Svih")
e_comp <- c("Ev", "Evic", "Evih")
i_comp <- c("Iv", "Ivic", "Ivih")

out_0_df_long <- out_0_df_long %>%
  mutate(state_categ = case_when(state_var %in% s_comp ~ "Susceptible",
                                 state_var %in% e_comp ~ "Exposed",
                                 state_var %in% i_comp ~ "Infectious"))


ggplot(out_0_df_long, aes(x = t, y = prop_mosq, col = as.factor(state_categ)))+
  geom_line()+
  ylim(0, 1)

#no ivermectin####
out_1 <- run_model(model = "mosquito_ivermectin_model",
                 init_EIR = 100,
                 #increasing this EIR gets rid of the dde error
                 gamma_c = 0,
                 gamma_h = 0, time = 730) #need to run a bit longer to get to equilibrium

out_1_df <- as.data.frame(out_1) %>%
  select(t, Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic) %>%
  rowwise() %>%
  mutate(total_mosq = sum(Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic))

out_1_df_prop <- out_1_df %>%
  rowwise() %>%
  mutate(Sv = Sv/total_mosq,
         Ev = Ev/total_mosq,
         Iv = Iv/total_mosq,
         Svih = Svih/total_mosq,
         Evih = Evih/total_mosq,
         Ivih = Ivih/total_mosq,
         Svic = Svic/total_mosq,
         Evic = Evic/total_mosq,
         Ivic = Ivic/total_mosq)

out_1_df_long <- gather(out_1_df_prop, state_var, prop_mosq, Sv:Ivic, factor_key=TRUE)


s_comp <- c("Sv", "Svic", "Svih")
e_comp <- c("Ev", "Evic", "Evih")
i_comp <- c("Iv", "Ivic", "Ivih")

out_1_df_long <- out_1_df_long %>%
  mutate(state_categ = case_when(state_var %in% s_comp ~ "Susceptible",
                                 state_var %in% e_comp ~ "Exposed",
                                 state_var %in% i_comp ~ "Infectious"),
         ivm_categ = case_when(grepl("ih", state_var) ~ "IVM from humans",
                               grepl("ic", state_var) ~ "IVM from cattle",
                               TRUE ~ "No IVM"))

out_1_plot <- ggplot(out_1_df_long, aes(x = t, y = prop_mosq, col = as.factor(state_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(ivm_categ,"No IVM", "IVM from humans","IVM from cattle"))+
  ylim(0, 1)+
  ggtitle("No ivermectin treatment")+
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes \n in compartment")+
  theme_minimal()

#colour-code this so coloured if no ivm, cattle or human

#100% human coverage####
out_2 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 0,
                   gamma_h = 1, time = 730) #need to run a bit longer to get to equilibrium

out_2_df <- as.data.frame(out_2) %>%
  select(t, Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic) %>%
  rowwise() %>%
  mutate(total_mosq = sum(Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic))

out_2_df_prop <- out_2_df %>%
  rowwise() %>%
  mutate(Sv = Sv/total_mosq,
         Ev = Ev/total_mosq,
         Iv = Iv/total_mosq,
         Svih = Svih/total_mosq,
         Evih = Evih/total_mosq,
         Ivih = Ivih/total_mosq,
         Svic = Svic/total_mosq,
         Evic = Evic/total_mosq,
         Ivic = Ivic/total_mosq)

out_2_df_long <- gather(out_2_df_prop, state_var, prop_mosq, Sv:Ivic, factor_key=TRUE)


s_comp <- c("Sv", "Svic", "Svih")
e_comp <- c("Ev", "Evic", "Evih")
i_comp <- c("Iv", "Ivic", "Ivih")

out_2_df_long <- out_2_df_long %>%
  mutate(state_categ = case_when(state_var %in% s_comp ~ "Susceptible",
                                 state_var %in% e_comp ~ "Exposed",
                                 state_var %in% i_comp ~ "Infectious"),
         ivm_categ = case_when(grepl("ih", state_var) ~ "IVM from humans",
                               grepl("ic", state_var) ~ "IVM from cattle",
                               TRUE ~ "No IVM"))


out_2_plot <- ggplot(out_2_df_long, aes(x = t, y = prop_mosq, col = as.factor(state_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(ivm_categ,"No IVM", "IVM from humans","IVM from cattle"))+
  ylim(0, 1)+
  ggtitle("100% human IVM coverage")+
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes \n in compartment")+
  theme_minimal()


#do work from here
#100% cow coverage####
out_3 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 1,
                   gamma_h = 0, time = 730) #need to run a bit longer to get to equilibrium

out_3_df <- as.data.frame(out_3) %>%
  select(t, Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic) %>%
  rowwise() %>%
  mutate(total_mosq = sum(Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic))

out_3_df_prop <- out_3_df %>%
  rowwise() %>%
  mutate(Sv= Sv/total_mosq,
         Ev = Ev/total_mosq,
         Iv = Iv/total_mosq,
         Svih = Svih/total_mosq,
         Evih = Evih/total_mosq,
         Ivih = Ivih/total_mosq,
         Svic = Svic/total_mosq,
         Evic = Evic/total_mosq,
         Ivic = Ivic/total_mosq)

out_3_df_long <- gather(out_3_df_prop, state_var, prop_mosq, Sv:Ivic, factor_key=TRUE)

s_comp <- c("Sv", "Svic", "Svih")
e_comp <- c("Ev", "Evic", "Evih")
i_comp <- c("Iv", "Ivic", "Ivih")

out_3_df_long <- out_3_df_long %>%
  mutate(state_categ = case_when(state_var %in% s_comp ~ "Susceptible",
                                 state_var %in% e_comp ~ "Exposed",
                                 state_var %in% i_comp ~ "Infectious"),
         ivm_categ = case_when(grepl("ih", state_var) ~ "IVM from humans",
                               grepl("ic", state_var) ~ "IVM from cattle",
                               TRUE ~ "No IVM"))
out_3_plot <- ggplot(out_3_df_long, aes(x = t, y = prop_mosq, col = as.factor(state_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(ivm_categ,"No IVM", "IVM from humans","IVM from cattle"))+
  ylim(0, 1)+
  ggtitle("100% cattle IVM coverage")+
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes \n in compartment")+
  theme_minimal()

#50% coverage in each
out_4 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 0.5,
                   gamma_h = 0.5, time = 730) #need to run a bit longer to get to equilibrium

out_4_df <- as.data.frame(out_4) %>%
  select(t, Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic) %>%
  rowwise() %>%
  mutate(total_mosq = sum(Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic))

out_4_df_prop <- out_4_df %>%
  rowwise() %>%
  mutate(Sv= Sv/total_mosq,
         Ev = Ev/total_mosq,
         Iv = Iv/total_mosq,
         Svih = Svih/total_mosq,
         Evih = Evih/total_mosq,
         Ivih = Ivih/total_mosq,
         Svic = Svic/total_mosq,
         Evic = Evic/total_mosq,
         Ivic = Ivic/total_mosq)

out_4_df_long <- gather(out_4_df_prop, state_var, prop_mosq, Sv:Ivic, factor_key=TRUE)


s_comp <- c("Sv", "Svic", "Svih")
e_comp <- c("Ev", "Evic", "Evih")
i_comp <- c("Iv", "Ivic", "Ivih")

out_4_df_long <- out_4_df_long %>%
  mutate(state_categ = case_when(state_var %in% s_comp ~ "Susceptible",
                                 state_var %in% e_comp ~ "Exposed",
                                 state_var %in% i_comp ~ "Infectious"),
         ivm_categ = case_when(grepl("ih", state_var) ~ "IVM from humans",
                               grepl("ic", state_var) ~ "IVM from cattle",
                               TRUE ~ "No IVM"))

out_4_plot <- ggplot(out_4_df_long, aes(x = t, y = prop_mosq, col = as.factor(state_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(ivm_categ,"No IVM", "IVM from humans","IVM from cattle"))+
  ylim(0, 1)+
  ggtitle("50% cattle, 50% IVM coverage")+
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes \n in compartment")+
  theme_minimal()

#plot showing prop of mosquitoes, for each ivm intervention, in each compartment
ivermectin_coverage_plots <- plot_grid(out_1_plot, out_2_plot, out_3_plot, out_4_plot)
ggsave(ivermectin_coverage_plots, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/plots_ivm_model/odin_model/ivm_cov_plots.svg")
ggsave(ivermectin_coverage_plots, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/plots_ivm_model/odin_model/ivm_cov_plots.pdf",
       device = cairo_pdf,
       width = 297,
       height = 210,
       units = "mm")

#show proportions the other way
out_1_plot_2 <- ggplot(out_1_df_long, aes(x = t, y = prop_mosq, col = as.factor(ivm_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(state_categ,"Susceptible", "Exposed","Infectious"))+
  ggtitle("No IVM treatment")+
  labs(col = "IVM treatment", y = "Proportion of mosquitoes in compartment \n (for given treatment)")+
  theme_minimal()

out_2_plot_2 <- ggplot(out_2_df_long, aes(x = t, y = prop_mosq, col = as.factor(ivm_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(state_categ,"Susceptible", "Exposed","Infectious"))+
  ggtitle("100% human IVM coverage")+
  labs(col = "IVM treatment", y = "Proportion of mosquitoes in compartment \n (for given treatment)")+
  theme_minimal()

out_3_plot_2 <- ggplot(out_3_df_long, aes(x = t, y = prop_mosq, col = as.factor(ivm_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(state_categ,"Susceptible", "Exposed","Infectious"))+
  ggtitle("100% cattle IVM coverage")+
  labs(col = "IVM treatment", y = "Proportion of mosquitoes in compartment \n (for given treatment)")+
  theme_minimal()

out_4_plot_2 <- ggplot(out_4_df_long, aes(x = t, y = prop_mosq, col = as.factor(ivm_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(state_categ,"Susceptible", "Exposed","Infectious"))+
  ggtitle("100% human IVM coverage")+
  labs(col = "IVM treatment", y = "Proportion of mosquitoes in compartment \n (for given treatment)")+
  theme_minimal()

ivm_cov_plots_2 <- plot_grid(out_1_plot_2, out_2_plot_2, out_3_plot_2, out_4_plot_2)
ggsave(ivm_cov_plots_2, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/plots_ivm_model/odin_model/ivm_cov_plots2.pdf",
       device = cairo_pdf,
       width = 297,
       height = 210,
       units = "mm")

#look at sporozoite prevalence in each mosquito compartment

#no ivm
out_1_sporo <- as.data.frame(out_1)
out_1_sporo_dat <- out_1_sporo %>%
  select(t, sporo_rate_no_ivm, sporo_rate_ivm_human, sporo_rate_ivm_cattle, sporo_rate_total)

out_1_sporo_dat_long <- gather(out_1_sporo_dat, ivm_categ, sporo_rate, sporo_rate_no_ivm:sporo_rate_total, factor_key = TRUE)

no_ivm_SR <- ggplot(out_1_sporo_dat_long, aes(x = t, y = sporo_rate))+
  geom_line()+
  facet_wrap(~ivm_categ)+
  ggtitle("No Ivermectin")


#ivm in humans
out_2_sporo <- as.data.frame(out_2)
out_2_sporo_dat <- out_2_sporo %>%
  select(t, sporo_rate_no_ivm, sporo_rate_ivm_human, sporo_rate_ivm_cattle, sporo_rate_total)

out_2_sporo_dat_long <- gather(out_2_sporo_dat, ivm_categ, sporo_rate, sporo_rate_no_ivm:sporo_rate_total, factor_key = TRUE)

ivm_h_SR <- ggplot(out_2_sporo_dat_long, aes(x = t, y = sporo_rate))+
  geom_line()+
  facet_wrap(~ivm_categ)+
  ggtitle("Human IVM")

#ivm in cattle
out_3_sporo <- as.data.frame(out_3)
out_3_sporo_dat <- out_3_sporo %>%
  select(t, sporo_rate_no_ivm, sporo_rate_ivm_human, sporo_rate_ivm_cattle, sporo_rate_total)

out_3_sporo_dat_long <- gather(out_3_sporo_dat, ivm_categ, sporo_rate, sporo_rate_no_ivm:sporo_rate_total, factor_key = TRUE)

ivm_c_SR <- ggplot(out_3_sporo_dat_long, aes(x = t, y = sporo_rate))+
  geom_line()+
  facet_wrap(~ivm_categ)+
  ggtitle("Cattle IVM")

#visualising this differently: just need the sporozoite rate (across all mosq pops)

out_1_sporo_total <- as.data.frame(out_1)
out_1_sporo_total <- out_1_sporo_total %>%
  select(t, sporo_rate_total, prev) %>%
  mutate(ivm_categ = "No IVM")

out_2_sporo_total <- as.data.frame(out_2)
out_2_sporo_total <- out_2_sporo_total %>%
  select(t, sporo_rate_total, prev) %>%
  mutate(ivm_categ = "Human IVM")

out_3_sporo_total <- as.data.frame(out_3)
out_3_sporo_total <- out_3_sporo_total %>%
  select(t, sporo_rate_total, prev) %>%
  mutate(ivm_categ = "Cattle IVM")

out_4_sporo_total <- as.data.frame(out_4)
out_4_sporo_total <- out_4_sporo_total %>%
  select(t, sporo_rate_total, prev) %>%
  mutate(ivm_categ = "50% Cattle, 50% Human IVM")

df_list <- list(out_1_sporo_total, out_2_sporo_total, out_3_sporo_total, out_4_sporo_total)

sporo_total_output <- do.call(rbind, df_list)

sporo_rate_plot <- ggplot(sporo_total_output, aes(x = t, y = sporo_rate_total, col = as.factor(ivm_categ)))+
  geom_line()+
  ylab("Sporozoite Rate")+
  xlab("Time")+
  labs(col = "IVM category")+
  theme_minimal()
ggsave(sporo_rate_plot, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/plots_ivm_model/odin_model/sporo_rate_plot.svg")

#malaria prevalence
malaria_prev <- ggplot(sporo_total_output, aes(x = t, y = prev, col = as.factor(ivm_categ)))+
  geom_line()+
  ylab("Prevalence")+
  xlab("Time")+
  labs(col = "IVM category")+
  theme_minimal()
ggsave(malaria_prev, file = "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/plots_ivm_model/odin_model/prev_plot.svg")


#under what scenarios for Q does cattle ivm have a bigger impact####. Need Q0 < 0.5
out_Q_no_ivm <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   time = 730,
                   Q0 = 0.5, gamma_h = 0, gamma_c = 0)
out_Q_ivm_h <- run_model(model = "mosquito_ivermectin_model",
                          init_EIR = 100,
                          #increasing this EIR gets rid of the dde error
                          time = 730,
                          Q0 = 0.5, gamma_h = 1, gamma_c = 0)

out_Q_ivm_c <- run_model(model = "mosquito_ivermectin_model",
                         init_EIR = 100,
                         #increasing this EIR gets rid of the dde error
                         time = 730,
                         Q0 = 0.5, gamma_h = 0, gamma_c = 1)

out_Q_ivm_c_h <- run_model(model = "mosquito_ivermectin_model",
                         init_EIR = 100,
                         #increasing this EIR gets rid of the dde error
                         time = 730,
                         Q0 = 0.5, gamma_h = 0.5, gamma_c = 0.5)

out_Q_no_ivm_df <- as.data.frame(out_Q_no_ivm)
out_Q_ivm_h_df <- as.data.frame(out_Q_ivm_h)
out_Q_ivm_c_df <- as.data.frame(out_Q_ivm_c)
out_Q_ivm_c_h_df <- as.data.frame(out_Q_ivm_c_h)

out_Q_no_ivm_df <- out_Q_no_ivm_df %>%
  mutate(ivm_categ = "No IVM")

out_Q_ivm_h_df <- out_Q_ivm_h_df %>%
  mutate(ivm_categ = "IVM human")

out_Q_ivm_c_df <- out_Q_ivm_c_df %>%
  mutate(ivm_categ = "IVM cattle")

out_Q_ivm_c_h_df <- out_Q_ivm_c_h_df %>%
  mutate(ivm_categ = "50% Cattle, 50% Human IVM")


df_list2 <- list(out_Q_no_ivm_df, out_Q_ivm_h_df, out_Q_ivm_c_df, out_Q_ivm_c_h_df)

sporo_total_output2 <- do.call(rbind, df_list2)

ggplot(sporo_total_output2, aes(x = t, y = prev, col = as.factor(ivm_categ)))+
  geom_line()+
  ylab("Prevalence")+
  xlab("Time")+
  labs(col = "IVM category")

#outputting Q under different net coverages, then repeating analysis to see how cattle ivermectin will work####
