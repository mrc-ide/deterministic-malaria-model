# Imperial College Deterministic Malaria Model - Ivermectin integration is NOT complete. Please refer to main for the full malaria model
This repository contains the R package version of the Imperial College malaria full transmission model. This is an age-stratified model, which incorporates the acquisition and loss of immunity experienced by human populations living in malaria-endemic settings. The parameterisation of the model is informed by data from a range of settings across sub-Saharan Africa. The model facilitates the introduction of public health interventions against malaria, such as insecticide-treated nets or indoor residual spraying. The seasonality of malaria transmission, driven by rainfall patterns, can also be specified in the model. Characterising malaria transmission observed in a particular setting can help inform decision-making on which interventions against malaria could be most impactful. 

## The background to the model
This model has been developed over a number of years, and has been described in multiple publications. We recommend studying these papers to gain an understanding of the structure of the model and the data used to parameterise it. The key publications are: 

* [Griffin et al. PloS Medicine, 2010](https://doi.org/10.1371/journal.pmed.1000324) is the original paper that describes the model and illustrates the type of output it produces. The supplementary materials contain extensive detail of the model structure.
* [Griffin et al. Nature Communications, 2014](https://doi.org/10.1038/ncomms4136) explains further refinements to the model, which are described in the supplementary materials
* [White et al. Parasites & Vectors, 2011](https://doi.org/10.1186/1756-3305-4-153) details the larval stages of the mosquito model, which you will also need to be familiar with.
* [Griffin et al. Lancet, 2016](https://doi.org/10.1016/S1473-3099(15)00423-5) explores the impact of multiple intervention scenarios upon case incidence and malaria mortality rates.

## Installing the model
Install the [Odin](https://github.com/mrc-ide/odin) package from Github (we need a newer version than the one available on CRAN):
```
devtools::install_github("mrc-ide/odin")
```
Then install the ICDMM package from github:
```
devtools::install_github("mrc-ide/deterministic-malaria-model")
```
Some users have experienced problems running the package on R version 3.X. Please update to R version 4.X if you have any issues.

## Running the model

There are vignettes that explain how to run the [simple model](https://mrc-ide.github.io/deterministic-malaria-model//articles/run_model_example.html) and the more [flexible versions](https://mrc-ide.github.io/deterministic-malaria-model//articles/run_model.html) of the model. Alternative versions of the model, exploring the impact of other public health interventions against malaria, have also been developed and are discussed [here](https://mrc-ide.github.io/deterministic-malaria-model//articles/run_model_alternative.html).
