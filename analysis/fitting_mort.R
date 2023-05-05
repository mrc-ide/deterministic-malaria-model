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
  #n <- params[4]
  n <- 0.5 #square root test
  t <- 1:28
  hazard_out <- d*(t^n)*exp(-c*t)+a
  error <- sum((hazards_data$d300 - hazard_out)^2)
  return(error)
}

#fitMe(params = c(1.5, 1, 13))

#o <- optim(c(a = 1.5 , c = 1, d = 13, n = 1), fitMe)
o <- optim(c(a = 1.5 , c = 1, d = 13), fitMe)
o$par

f_haz2 <- function(x, d, c, a, n)d*(x^n)*exp(-c*x) + a
ori_plot <- hazards_long %>%
  filter(dose == "d300") %>%
  ggplot(aes(x = day, y = hazard))+
  xlim(1, 28)+
  geom_point()+
  #geom_function(fun = f_haz2, args = list(a = o$par[1], c = o$par[2], d= o$par[3], n = o$par[4]), colour = "red")
  geom_function(fun = f_haz2, args = list(a = o$par[1], c = o$par[2], d= o$par[3], n = 0.5), colour = "red")

#try recreating like Isaac has said.
t <- seq(1:28)
a <- rep(o$par[1], 28) #changing it slightly with +0.5
c <- rep(o$par[2], 28) #http://127.0.0.1:29903/graphics/a02c4310-2703-4f21-8a82-6254dc923956.http://127.0.0.1:29903/graphics/a02c4310-2703-4f21-8a82-6254dc923956.pngpng
d <- rep(o$par[3], 28)
#n <- rep(o$par[4], 28)
n <- rep(0.5, 28)
fake_df <- data.frame(t, a, c, d, n)
fake_df <- fake_df %>%
  mutate(daily_haz = d*(t^n)*exp(-c*t) + a)

#then
fitMe_again <- function(params) {
  a <- params[1]
  c <- params[2]
  d <- params[3]
  #n <- params[4]
  n <- 0.5
  t <- 1:28
  hazard_out <- d*(t^n)*exp(-c*t)+a
  error <- sum((fake_df$daily_haz - hazard_out)^2)
  return(error)
}

p <- optim(c(a = 1.5 , c = 1, d = 13), fitMe_again)
p$par
o$par

#can get the curve back when use the parameters, as expected
new_plot <- ggplot(fake_df, aes(x = t, y = daily_haz))+
  geom_point()+
  xlim(1, 28)+
  ylim(0, 10)+
  #geom_function(fun = f_haz2, args = list(a = p$par[1], c = p$par[2], d= p$par[3], n = p$par[4]), colour = "red")
  geom_function(fun = f_haz2, args = list(a = p$par[1], c = p$par[2], d= p$par[3], n = 0.5), colour = "red")


#great, all works
