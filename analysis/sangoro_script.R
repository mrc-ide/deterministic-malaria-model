library(ICDMM)
library(tidyverse)
library(cowplot) #this package allows you to plot figures side by side
model_0 <- run_model(model = "odin_model",
                     init_EIR = 100, # set initial EIR
                     time = 730 #run model for two years
                     )
plot(model_0$t,model_0$prev, main= "Prevalance (init_EIR = 100)", ylim = c(0, 1), type='l')

model_1 <- run_model(model = "odin_model",
                     init_EIR = 50,
                     time = 730)
plot(model_1$t, model_1$prev, main = "Prevalence (init_EIR = 50)", ylim = c(0, 1), type = "l") #note that the prevalence on this plot is lower, because had a lower EIR


#now use ggplot to get the same plots
model_0_df <- as.data.frame(model_0) #need to convert model output to dataframe for ggplot
prev_plot_0 <- ggplot(model_0_df, aes(x = t, y = prev))+
  geom_line()+ #this tells ggplot we are going to use lines
  xlim(0, 800) + #sets limit of x axis
  ylim(0, 1)+ #sets limit of y axis
  ggtitle("Prevalence of malaria (init_EIR = 100)") #set a title for the plot
prev_plot_0 #run this to show the plot


model_1_df <- as.data.frame(model_1)
prev_plot_1 <- ggplot(model_1_df, aes(x = t, y = prev))+
  geom_line()+
  xlim(0, 800)+
  ylim(0, 1)+
  ggtitle("Prevalence of malaria (init_EIR = 50)")
prev_plot_1 #run this to show the plot

#now plot these side by side
plot_grid(prev_plot_0, prev_plot_1, labels = c("A", "B"))


#can also combine the data into one plot and colour code by EIR
#this %>% is a pipe - it connects pieces of code which use tidyverse functions
model_0_df <- model_0_df %>%
  mutate(EIR = "init_EIR = 100") # I am creating a column which tells us the EIR in this df

model_1_df <- model_1_df %>%
  mutate(EIR = "init_EIR = 50")

model_output <- rbind(model_0_df, model_1_df) #combine the dataframes

ggplot(model_output, aes(x = t, y = prev, col = as.factor(EIR)))+
  geom_line()+
  scale_color_manual(values = c("green", "blue"), name = "EIR", labels = c("100", "50"))+
  theme_minimal() + # this takes away the grey background
  labs(x = "Time", y = "Prevalence")+
  ylim(0, 1)+
  xlim(0, 800)
