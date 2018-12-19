#' Pipe operator
#'
#' See \code{\link[magrittr]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#------------------------------------------------
#' Grab parameters from odin model
#'
#' \code{grab_user_params} returns vector of all user required variables from an odin
#' model path
#'
#' @param model_path Character path to odin model
#'
#' @return Vector of characters for required variables
#'
#'

grab_user_params <- function(model_path){

  ## first read lines from model path
  lines <- readLines(model_path)

  ## remove all spaces as makes life easier
  lines <- gsub("[[:space:]]", "", lines)

  ## find all positions that have either <-user() or =user() and dont begin with #
  var_pos <- grep("^[^#].*<-user()|^[^#].*=user()",lines,fixed=FALSE)

  ## remove all characters including and after <-user() or =user() and then remove brackets etc.
  res <- gsub(pattern = "<-user().*$|=user().*$",replacement = "",x = lines)[var_pos] %>%
    sapply(function(x){return(gsub(pattern = "\\[.*\\]|\\(.*\\)",replacement="",x=x))}) %>%
    as.character()


  ## return answer
  return(res)

}


#------------------------------------------------
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
#'

create_r_model <- function(odin_model_path = system.file("extdata/odin_model.R",package="hanojoel"),
                           het_brackets = 5,
                           age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
                           init_EIR = 10,
                           init_ft = 0.4,
                           country = NULL,
                           admin2 = NULL,
                           ...){

#  browser()
  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector=age,EIR=init_EIR,ft=init_ft,
                                   model_param_list = mpl,het_brackets=5,
                                   country = country,
                                   admin_unit = admin2)
  # create odin generator
  gen <- odin::odin(odin_model_path,verbose=FALSE,build = TRUE)

  ## now autogenerate the model including all user allocation
  # first grab all user required params
  model_formals <- formals(gen)

  user_vars <- names(model_formals$user)[-match(c("","age"),names(model_formals$user))]

  # create temp tp store generator
  temp <- tempfile(fileext = ".txt")


  # create the generator
  writeLines(text = paste0("generate_default_model <- function(dat,generator,dde = TRUE){\n",
                           "mod <- generator(",
                           paste0(sapply(user_vars, function(x){(paste0(x,"=dat$",x,",\n") )}),collapse = ""),
                           "age = dat$age*365,
                           use_dde=dde)\n}"),con = temp)

  # source generator
  #source(temp)


  #create model with initial values
  #mod <- generate_default_model(dat=state,generator=gen,dde=TRUE)

  # return mod
  return(list("generate_model_function"=source(temp)$value,"state"=state,"generator"=gen,"mpl"=mpl))
}

