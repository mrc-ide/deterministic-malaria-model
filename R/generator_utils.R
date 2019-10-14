#' Create generator for model creation
#'
#' \code{create_r_model} returns list with generator function automatically created given the odin model
#' specified.
#'
#' @param odin_model_path Character path to odin model
#' @param het_brackets Numeric for heterogeneity brackets
#' @param age Age vector
#' @param init_EIR Initial EIR for initial solution
#' @param init_ft Initial ft for initial solution
#' @param country Character of country in which admin unit exists
#' @param admin2 Character for admin level 2
#' @param ... Any other parameters needed for non-standard model. If they share the same name
#' as any of the defined parameters \code{model_param_list_create} will stop. You can either write
#' any extra parameters you like individually, e.g. create_r_model(extra1 = 1, extra2 = 2)
#' and these parameteres will appear appended to the returned list, or you can pass explicitly
#' the ellipsis argument as a list created before, e.g. create_r_model(...=list(extra1 = 1, extra2 = 2))
#'
#' @return list of generator function, initial state, model parameters and generator
#'
#' @export

create_r_model <- function(odin_model_path = system.file("extdata/odin_model.R",package="hanojoel"),
                           het_brackets = 5,
                           age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
                           init_EIR = 10,
                           init_ft = 0.4,
                           country = NULL,
                           admin2 = NULL,
                           ...){

  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector=age, EIR=init_EIR,ft=init_ft,
                                   model_param_list = mpl, het_brackets=het_brackets,
                                   country = country,
                                   admin_unit = admin2)

  # create odin generator
  gen <- odin::odin(odin_model_path, verbose=FALSE)
  state <- state[names(state) %in% names(formals(gen))]

  # return mod
  return(list("state"=state,"generator"=gen,"mpl"=mpl))
}

