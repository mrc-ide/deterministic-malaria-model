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
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes in compartment")+
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
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes in compartment")+
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
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes in compartment")+
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
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes in compartment")+
  theme_minimal()

#plot showing prop of mosquitoes, for each ivm intervention, in each compartment
ivermectin_coverage_plots <- plot_grid(out_1_plot, out_2_plot, out_3_plot, out_4_plot)
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



#outputting Q under different net coverages####
