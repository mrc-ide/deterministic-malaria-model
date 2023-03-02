#Exploring the ivermectin model####

#have created an IVM model with x3 SEI compartments, depending on whether the mosquito gets IVM from human/cattle/none at all
require(ICDMM)
require(tidyverse)

#no ivermectin####
out_1 <- run_model(model = "mosquito_ivermectin_model",
                 init_EIR = 100,
                 #increasing this EIR gets rid of the dde error
                 gamma_c = 0,
                 gamma_h = 0, time = 730) #need to run a bit longer to get to equilibrium

out_1_df <- as.data.frame(out_1) %>%
  select(t, Sv, Ev, Iv, Svih, Evih, Ivih, Svic, Evic, Ivic)

out_1_df_long <- gather(out_1_df, state_var, count, Sv:Ivic, factor_key=TRUE)


s_comp <- c("Sv", "Svic", "Svih")
e_comp <- c("Ev", "Evic", "Evih")
i_comp <- c("Iv", "Ivic", "Ivih")

out_1_df_long <- out_1_df_long %>%
  mutate(state_categ = case_when(state_var %in% s_comp ~ "Susceptible",
                                 state_var %in% e_comp ~ "Exposed",
                                 state_var %in% i_comp ~ "Infectious"),
         ivm_categ = case_when(grepl("ih", state_var) ~ "Ivermectin-treated human",
                               grepl("ic", state_var) ~ "Ivermectin-treated cattle",
                               TRUE ~ "No ivermectin treatment"))


ggplot(out_1_df_long, aes(x = t, y = count, col = as.factor(ivm_categ)))+
  geom_line()+
  facet_wrap(~state_categ)

#colour-code this so coloured if no ivm, cattle or human

#100% human coverage####
out_2 <- run_model(model = "mosquito_ivermectin_model",
                   init_EIR = 100,
                   #increasing this EIR gets rid of the dde error
                   gamma_c = 0,
                   gamma_h = 1, time = 730) #need to run a bit longer to get to equilibrium

