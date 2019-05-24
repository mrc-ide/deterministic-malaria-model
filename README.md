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
`install.packages("statmod","ggplot2","reshape2","dplyr","magrittr", "RecordLinkage")`

## Installing the model
Easiest way to install this model (until we can make it publicly available) is probably to clone this repository to your local machine, open the `hanojoel.Rproj` file in RStudio and then Ctrl+Shift+B to build and reload the package

## Running the model
Below is some code that shows how you might go about performing simple runs of the model
```
  # load the model package (or have it built and reloaded as described above)
  library(hanojoel)

  # create a vector of age categories
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)
  
  # provide a value of the annual EIR for this model run
  init_EIR <- 10
  
  # provide a string for the admin 2 unit that you want to use the seasonality profile for
  admin_str <- "Katanga"
  
  # provide the length of time (in days) that you want to run the model for
  time_period <- 365*2
  
  # provide a value for the proportion of cases that are treated (referred to as ft in the paper)
  prop_treated <- 0.6
  
  # run the model
  model_run <- run_model(age=init_age, EIR=init_EIR, ft = prop_treated, admin2 = admin_str, time = time_period)
  
  # After this has finished, `model_run` will have 2 components. model_run$plot is a ggplot of the main compartments of the model
  # over time (summed over all age, biting & intervention categories). model_run$dat is a list containing 
  # the estimated values for all of the compartments in the model. For example, model_run$dat$S will 
  # be a big 3-dimensional array with the indices relating to age, biting heterogeneity and intervention categories.
  
```
## Running alternative models

The basic model used above has been adapted by various people, to add SMC, ivermectin 
courses, MDA etc. As such these model changes will have a number of diffrent parameters
required within the actual odin model code (e.g. see "inst/extdata/odin_model_IVM_SMChet.R").
As such it becomes annoying having to do housekeeping about what parameters are required 
from the user, and creating the model becomes quite drawn out in terms of specifying what 
all the user paramaters are. As such there is an altenrative function `create_r_model` 
that will take some of the burden away, as demonstrated with this model that assumes a 
proportion of hrp2 deleted parasites exist.

```
# first we use the function, by specifying very similar parameters such as 
# the odin model path, and the initial EIR, ft, age categories, region etc.
wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model_hrp2.R",package = "hanojoel"),
                                het_brackets = 5,
                                age = init_age,
                                init_EIR = init_EIR,
                                init_ft = 0.4,
                                country = "Uganda",
                                admin2 = "Tororo")

```
The outpput of `create_r_model` returns a list with 4 elements. 
1. The first is an automatically created function called `generate_default_model`, which has scanned over the odin model and
identified all the user required parameters before writing to a temp file the required function
that pulls all the user required variables from your initial state solution. 
2. The initial state solution is the second element in the list. This is a list that contains
the equilibirum solution for the default parameters and also contains all the parameters used
to create the initial solution. These parameters will also then be used when the model is run 
forwards, and thus if you want to run the model for the next year with a different value for `ft` 
you can change that here. This will no affect the initial solution, but simply the parameters 
that are used when running the deterministic model forward from the initial solution.
3. The actual odin model generated using the odin package.
4. The model parameter list that was used to create the initial solution.

We can now use the following to create our model:

```
mod <- wh$generate_model_function(dat = wh$state,generator = wh$generator,dde = TRUE)
```

We'll see the model errors intially as we didn't declare a value for `hrp2_prop`. This
variable is not needed in the basic model, but was found to be needed in this extended 
model. Thus we should add it to out model state as an extra variable.

```
wh$state$hrp2_prop <- 0.1

mod <- wh$generate_model_function(dat = wh$state,generator = wh$generator,dde = TRUE)
mod_run <- mod$run(t = 1:90)
out <- mod$transform_variables(mod_run)
plot(out$t,out$prev)

# If we wanted to increase the treatment for example:
wh$state$ft <- 0.7
mod <- wh$generate_model_function(dat = wh$state,generator = wh$generator,dde = TRUE)
mod_run <- mod$run(t = 1:90)
out2 <- mod$transform_variables(mod_run)

# And then compare the change in the prevalence
par(mfrow=c(2,1))
plot(out$t,out$prev)
plot(out2$t,out2$prev)


```

This will now work and should show how once you have the model set up you can then directly
change model variables within the state list object to alter how the deterministic model is 
run. 

## Further thoughts and rambles

Although the second way of running the model takes some of the burden from writing out the 
generate_model_function, it does start to become a bit messy still. In the above example it 
was only one variable that would affect how the model was run. We could have also specified 
this in the intial call to `create_r_model` like this:

```
wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model_hrp2.R",package = "hanojoel"),
                                het_brackets = 5,
                                age = init_age,
                                init_EIR = init_EIR,
                                init_ft = 0.4,
                                country = "Uganda",
                                admin2 = "Tororo",
                                hrp2_prop = 0.1)

```

This will append the parameter to the model parameter list and it will appear in the intial state list
object. If we had lots of variables we need to add then we could also specify these as a list:

```
wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model_hrp2.R",package = "hanojoel"),
                                het_brackets = 5,
                                age = init_age,
                                init_EIR = init_EIR,
                                init_ft = 0.4,
                                country = "Uganda",
                                admin2 = "Tororo",
                                ... = list("hrp2_prop"=0.1, "IVM_start"=1)
                                )

```

However, this does not change how the initial solution was created as ultimately the intial 
solution doesn't change at the moment, regardless of what model you want to run. This will be an issue if,
for example your model has a 4th dimension, such as the IVM_SMChet model. As such there is now a catch in the 
`equilibirum_init_create` function that looks for a variable called `ncc`, so that if you declared `ncc`
as an extra variable within `create_r_model` the `equilibirum_init_create` function will catch this. 

This starts to become quite ugly and not reliably expandable, as you will need to change the `equilibirum_init_create`
function every time a new model is needed that sets up the initial solution differently. It's not perhaps the worst 
thing in the world, but it might be nice here to then specify a path name as a new variable to where your initial solution
function is if it becomes very different. 

The above is hoped to possibly explain how we might start thinking about how to make the malaria models more 
modular and operational in the future, but the above might not be the correct way yet. 
