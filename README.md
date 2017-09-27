# deterministic-R-model
This repository contains the R package version of the Imperial College malaria full transmission model.

## Learning the model
Assuming that you are new to the model, the first step is to read the papers that describe the model to get a good understanding of the theory behind it and how it works. To start with you should read:
* This [Griffin 2010 paper](http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324) is the original paper explaining the model and shows the kind of work produced by it. Make sure to read the supplementary material, this explains most of the model in detail.
* The [Griffin 2014 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923296/) explains further refinements to the model and again the supplementary material is a good place to start to get your head around how the model works.
* This [White 2011 paper](https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153) details the larval section of the model, which you will also need to be familiar with.

## Using Odin
To run the model forward in time we need to solve a system of differential equations, we will do this using the R package that Rich Fitzjohn made called [Odin](https://github.com/richfitz/odin). The Odin package lets us write down the differential equations in a domain-specific language that is similar to R, it then converts the model that we have written down into C code and uses a numerical method to get an approximation of the solutions to the system of differential equations.

* [Install the Odin package](https://github.com/richfitz/odin#installation) and [read the vignettes](https://richfitz.github.io/odin/vignettes/odin.html) to make sure you understand how to use Odin. 
* Install a couple more R packages that you will need to get the model running:
`install.packages("statmod","ggplot2","reshape2","dplyr","magrittr")`

## Installing the model
Easiest way to install this model (until we can make it publicly available) is probably to clone this repository to your local machine, open the `hanojoel.Rproj` file in RStudio and then Ctrl+Shift+B to build and reload the package

## Running the model
Below is some code that shows how you might go about performing simple runs of the model
```
  # load the model package (or have it built and reloaded as described above)
  library(hanojoel)

  # create a vector of age categories
  init_age <- c(0,0.25,0.5,0.75,1,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)
  
  # provide a value of the annual EIR for this model run
  init_EIR <- 10
  
  # provide a string for the admin 2 unit that you want to use the seasonality profile for
  admin_str <- "Katanga"
  
  # provide the length of time (in days) that you want to run the model for
  time_period <- 365*2
  
  # provide a value for the proportion of cases that are treated (referred to as ft in the paper)
  prop_treated <- 0.6
  
  # run the model
  model_run <- Run_Model(age=init_age, EIR=init_EIR, ft = prop_treated, admin2 = admin_str, time = time_period)
  
  # After this has finished, `model_run` will have 2 components. model_run$plot is a ggplot of the main compartments of the model
  # over time (summed over all age, biting & intervention categories). model_run$dat is a list containing 
  # the estimated values for all of the compartments in the model. For example, model_run$dat$S will 
  # be a big 3-dimensional array with the indices relating to age, biting heterogeneity and intervention categories.
  ```
