% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_model.R
\name{run_model}
\alias{run_model}
\title{Run named odin model}
\usage{
run_model(
  model = "odin_model",
  het_brackets = 5,
  age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40,
    50, 60, 70, 80),
  init_EIR = 10,
  init_ft = 0.4,
  country = NULL,
  admin2 = NULL,
  time = 100,
  ...
)
}
\arguments{
\item{model}{Name of odin model used. Default = \code{"odin_model"}}

\item{het_brackets}{Numeric for heterogeneity brackets}

\item{age}{Age vector}

\item{init_EIR}{Initial EIR for initial solution}

\item{init_ft}{Initial ft for initial solution}

\item{country}{Character of country in which admin unit exists}

\item{admin2}{Character for admin level 2}

\item{time}{Integer for length of simulation in days. Default = 100}

\item{...}{Any other parameters needed for non-standard model. If they share the same name
as any of the defined parameters \code{model_param_list_create} will stop. You can either write
any extra parameters you like individually, e.g. create_r_model(extra1 = 1, extra2 = 2)
and these parameters will appear appended to the returned list, or you can pass explicitly
the ellipsis argument as a list created before, e.g. create_r_model(...=list(extra1 = 1, extra2 = 2))}
}
\value{
list of generator function, initial state, model parameters and generator
}
\description{
\code{run_model} returns list with generator function automatically created given the odin model
specified.
}
