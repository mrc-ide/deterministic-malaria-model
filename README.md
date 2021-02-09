# Imperial College Deterministic Malaria Model
This repository contains the R package version of the Imperial College malaria full transmission model. This is an age-stratified model, which incorporates the acquisition and loss of immunity experienced by human populations living in malaria-endemic settings. The parameterisation of the model is informed by data from a range of settings across sub-Saharan Africa. The model facilitates the introduction of public health interventions against malaria, such as insecticide-treated nets or indoor residual spraying. The seasonality of malaria transmission, driven by rainfall patterns, can also be specified in the model. Characterising malaria transmission observed in a particular setting can help inform decision-making on which interventions against malaria could be most impactful. 

## The background to the model
This model has been developed over a number of years, and has been described in multiple publications. We recommend studying these papers to gain an understanding of the structure of the model and the data used to parameterise it. The key publications are: 

* [Griffin et al. PloS Medicine, 2010](https://doi.org/10.1371/journal.pmed.1000324) is the original paper that describes the model and illustrates the type of output it produces. The supplementary materials contain extensive detail of the model structure.
* [Griffin et al. Nature Communications, 2014](https://doi.org/10.1038/ncomms4136) explains further refinements to the model, which are described in the supplementary materials
* [White et al. Parasites & Vectors, 2011](https://doi.org/10.1186/1756-3305-4-153) details the larval stages of the mosquito model, which you will also need to be familiar with.
* [Griffin et al. Lancet, 2016](https://doi.org/10.1016/S1473-3099(15)00423-5) explores the impact of multiple intervention scenarios upon case incidence and malaria mortality rates.

## Using Odin
To run the model forward in time we need to solve a system of differential equations, we will do this using the R package that Rich Fitzjohn made called [Odin](https://github.com/richfitz/odin). The Odin package lets us write down the differential equations in a domain-specific language that is similar to R, it then converts the model that we have written down into C code and uses a numerical method to get an approximation of the solutions to the system of differential equations.

* [Install the Odin package](https://github.com/richfitz/odin#installation) and [read the vignettes](https://richfitz.github.io/odin/vignettes/odin.html) to make sure you understand how to use Odin. Note: to be able to run this model, you will need to install the developmental version of the odin package directly from Github, using the devtools package.
* Install a couple more R packages that you will need to get the model running:
`install.packages("statmod","ggplot2","reshape2","dplyr","magrittr", "RecordLinkage")`

## Installing the model
Easiest way to install this model is probably to clone this repository to your local machine, open the `ICDMM.Rproj` file in RStudio and then Ctrl+Shift+B to build and reload the package. Alternatively, run the command `devtools::install_github("mrc-ide/ICDMM")` in RStudio.

## Running the model

There are vignettes that explain how to run the [simple model](https://mrc-ide.github.io/deterministic-malaria-model/articles/run_model.html) and the more [flexible versions](https://mrc-ide.github.io/deterministic-malaria-model/articles/create_r_model.html) of the model. 
