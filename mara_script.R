library(ICDMM)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# create a vector of age categories
init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3,4,5,6,7.5,10,15,20,30,40,50,60,70,80)

# provide the length of time (in days) that you want to run the model for
time_period <- 365*6

# provide a value for the proportion of cases that are treated
prop_treated <- 0.4

# provide a value of the annual EIR for this model run
init_EIR <- 110
ITN_IRS_on <- 3 * 365

country_str <- "Burkina Faso"
admin_str <- "Sahel"

out0 <- run_model(model = 'odin_model',
                  het_brackets = 5,
                  age = init_age,
                  num_int = 4,
                  init_EIR = init_EIR,
                  time = time_period,
                  init_ft = prop_treated,
                  country = country_str,
                  admin2 = admin_str,
                  irs_cov = 0,
                  itn_cov = 0,
                  ITN_IRS_on = ITN_IRS_on)

# 1. Model with intervention

out <- run_model(model = 'odin_model',
                 het_brackets = 5,
                 age = init_age,
                 num_int = 4,
                 time = time_period,
                 init_EIR = init_EIR,
                 init_ft = prop_treated,
                 country = country_str,
                 admin2 = admin_str,
                 irs_cov = 0,
                 itn_cov = 0.7,
                 ITN_IRS_on = ITN_IRS_on)

df <- data.frame('t' = out0$t, 'pr0' = out0$prev, 'pr' = out$prev)
dfm <- melt(df, id.vars = 't')
ggplot(dfm, aes(x=t, y = value, color = variable)) + geom_line() + theme_bw() + xlab('Time (days)') +
  ylab('prevalence') + labs(color = '')

plot(out0$mv[1*365:time_period])
min(out0$mv[1*365:time_period])
max(out0$mv[1*365:time_period])
mean(out0$mv[1*365:time_period])

dm <- data.frame('t' = out0$t[1*365:time_period], 'mv' = out0$mv[1*365:time_period])
ggplot(dm, aes(x=t,y=mv)) + geom_line() + theme_bw() + ylim(c(0,1))

#For what proportion of the year are there more humans than mosquitoes?
table(out0$mv[1*365:time_period] < 1)
345/(345+1481) #nearly 19%
