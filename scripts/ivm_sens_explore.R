#performing sensitivity analyses on the ivm model###

#MORTALITY

#1: excess mortality modelled purely with a switch
require(ICDMM)
require(tidyverse)

out_1 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c_0 = 0,
                   gamma_h_0 = 1, time = 730,
                   gamma_c_min_age = 5,
                   gamma_h_min_age = 5,
                   ivm_c_on = 731,
                   ivm_h_on = 20, #time turn on
                   mu_c_0 = 0,
                   mu_h_0 = 0.628)

out_1_df <- as.data.frame(out_1) %>%
  select(t, Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic, mv) %>%
  rowwise() %>%
  mutate(total_mosq = mv)

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
  #xlim(0,100)+
  ggtitle("No ivermectin treatment")+
  labs(col = "Mosquito infection state", y = "Proportion of mosquitoes \n in compartment")+
  theme_minimal()+
  geom_hline(yintercept = 0.940, lty = "dashed")

##looking at the number of mosquitoes
out_1_df_long_num <- gather(out_1_df, state_var, num_mosq, Sv:mv, factor_key=TRUE)
out_1_df_long_num <- out_1_df_long_num %>%
  mutate(state_categ = case_when(state_var %in% s_comp ~ "Susceptible",
                                 state_var %in% e_comp ~ "Exposed",
                                 state_var %in% i_comp ~ "Infectious",
                                 TRUE ~ "Total Mosq"),
         ivm_categ = case_when(grepl("ih", state_var) ~ "IVM from humans",
                               grepl("ic", state_var) ~ "IVM from cattle",
                               TRUE ~ "No IVM"))

ggplot(out_1_df_long_num, aes(x = t, y = num_mosq, col = as.factor(state_categ)))+
  geom_line()+
  facet_wrap(~fct_relevel(ivm_categ,"No IVM", "IVM from humans","IVM from cattle"))+
  #ylim(0, 1)+
  #xlim(0,100)+
  ggtitle("No ivermectin treatment")+
  labs(col = "Mosquito infection state", y = "Number of mosquitoes \n in compartment")+
  theme_minimal()+
  geom_hline(yintercept = 39.1, lty = "dashed")

#2: allow excess mortality to change over time (using a daily hazard)






#compare 1 and 2 in terms of prevalence. Is there a way to simplify incorporation of excess mort?


#DELIVERY

#3: all-in-one vs staggered distribution. Difference between the treatment on d1, 2, 3 every month of entire pop and doing this more gradually
#i.e. subpopA d1,d2,d3 then subpopB d1,2,3 then subpopC d1,d2,d3.
