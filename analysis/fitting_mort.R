#fitting a curve to the ivermectin hazards data --> ivermectin-induced mortality rates

require(tidyverse)
hazards <- read.table("C:/Users/nc1115/Documents/github/ivRmectin/IVM_derivation/ivermectin_hazards.txt", header = TRUE)
hazards$d300
hazards$d400
hazards_long <- gather(hazards, dose, hazard, d400:d300, factor_key = TRUE)
ggplot(hazards_long, aes(x = day, y = hazard, col = as.factor(dose)))+
  geom_point()+
  theme_minimal()

#start plotting the functions on


f_haz <- function(x) 1*x^2*exp(-0.2*x)

hazards_long %>%
  filter(dose == "d300") %>%
  ggplot(aes(x = day, y = hazard))+
  geom_point()+
  xlim(0, 30)+
  geom_function(fun = f_haz, colour = "red")

