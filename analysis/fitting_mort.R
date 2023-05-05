#fitting a curve to the ivermectin hazards data --> ivermectin-induced mortality rates

require(tidyverse)
require(optim)

hazards <- read.table("C:/Users/nc1115/Documents/github/ivRmectin/IVM_derivation/ivermectin_hazards.txt", header = TRUE)
hazards$d300
hazards$d400
hazards_long <- gather(hazards, dose, hazard, d400:d300, factor_key = TRUE)

hazards_data <- hazards %>%
  select(day, d300) %>%
  drop_na()
#fitting with optim
#minimise the sum sqs for each parameter in optim
#https://www.r-bloggers.com/2013/03/how-to-use-optim-in-r/


#writing for two parameters
f_haz <- function(x, d, c, a)d*x^2*exp(-c*x) + a

hazards_long %>%
  filter(dose == "d300") %>%
  ggplot(aes(x = day, y = hazard))+
  geom_point()+
  ylim(0, 15)+
  geom_function(fun = f_haz, args = list(a = 1.5, c = 1 , d= 13 ), colour = "red")

hazards_long %>%
  filter(dose == "d300") %>%
  ggplot(aes(x = day, y = hazard))+
  geom_point()+
  xlim(0, 30)+
  geom_function(fun = f_haz, colour = "red")

fitMe <- function(params) {
  a <- params[1]
  c <- params[2]
  d <- params[3]
  t <- 1:28
  hazard_out <- d*t^2*exp(-c*t)+a
  error <- sum((hazards_data$d300 - hazard_out)^2)
  #return(hazard_out)
}

#fitMe(params = c(1.5, 1, 13))

o <- optim(c(a = 1.5 , c = 1, d = 13), fitMe)
o$par

f_haz2 <- function(x, d, c, a)d*x^2*exp(-c*x) + a
hazards_long %>%
  filter(dose == "d300") %>%
  ggplot(aes(x = day, y = hazard))+
  xlim(1, 28)+
  geom_point()+
  geom_function(fun = f_haz2, args = list(a = o$par[1], c = o$par[2], d= o$par[3]), colour = "red")

#try recreating like Isaac has said.
fake_df <- o$par
