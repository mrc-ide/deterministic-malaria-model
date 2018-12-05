#------------------------------------------------
#' Run_Model
#'
#' \code{Run_Model} runs model using declared age, EIR, ft, country, admin
#'
#' @param age Vector of age brackets.
#' @param EIR Numeric of annual EIR
#' @param ft Numeric of proportion of symptomatic cases recieving treatment
#' @param time Numeric of length of time that the model will run for in days
#' @param admin2 Character of admin unit
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous xlab ylab geom_line ggplot
#' @importFrom reshape2 melt
#' @importFrom odin odin
#'
#' @export


Run_Model <- function(age,EIR,ft,admin2,time){
  mpl <- model_param_list_create()
  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector=age,EIR=EIR,ft=ft,
                                   model_param_list = mpl,het_brackets=5,
                                   admin_unit = admin2)
  # create odin generator
  odin_model_path <- system.file("extdata/odin_model.R",package="hanojoel")
  gen <- odin::odin(odin_model_path,verbose=FALSE,build = TRUE)


  #create model with initial values
  mod <- generate_default_model(ft=ft,age=age,dat=state,generator=gen,dde=TRUE)
  tt <- seq(0,time,1)
  # run model
  mod_run <- mod$run(tt)
  # shape output
  out <- mod$transform_variables(mod_run)
  hum <- list("S"=out$S,"T"=out$T,"D"=out$D,"A"=out$A,"U"=out$U,"P"=out$P)
  humsum <- lapply(hum,function(x){apply(x,sum,MARGIN=1)})
  df <- as.data.frame(matrix(unlist(humsum),nrow=length(tt)))
  colnames(df) <- c("S","T","D","A","U","P")
  df$t <- tt
  brks <- seq(from=0,to=time,by=365)
  ret <- ggplot(melt(df,id.vars="t"),aes(x=t,y=value,col=variable)) + geom_line() +
    scale_colour_manual(values=c("#000000","#CC0000","#339900","#3333FF","#CC9933","#CC3399"),name="Compartment",labels=c("S","T","D","A","U","P")) +
    scale_x_continuous(breaks=brks,labels=0:(length(brks)-1)) +
    xlab("Years") + ylab("Proportion of population")
  return(list("plot"=ret,"dat"=out))
}


## Odin generator function
generate_default_model <- function(ft,age,dat,generator,dde = TRUE){
  dat$pi <- pi
  dat$ft <- ft
  dat$age <- age
  mod <- generator(user=dat, use_dde=dde)
}
