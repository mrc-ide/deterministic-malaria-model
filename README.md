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
There are vignettes that explain how to run the [simple model](https://mrc-ide.github.io/deterministic-malaria-model/articles/run_model.html) and the more [flexible versions](https://mrc-ide.github.io/deterministic-malaria-model/articles/create_r_model.html) of the model. 
