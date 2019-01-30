#------------------------------------------------
#' run_model
#'
#' \code{run_model} runs model using declared age, EIR, ft, country, admin
#'
#' @param age Vector of age brackets.
#'   Default=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60)
#' @param EIR Numeric of annual EIR. Default=10
#' @param ft Numeric of proportion of symptomatic cases recieving treatment.
#'   Default=0.4
#' @param time Numeric of length of time that the model will run for in days.
#'   Default=365
#' @param admin2 Character of admin unit. Default = "Tororo"
#'
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous
#'   xlab ylab geom_line ggplot
#' @importFrom reshape2 melt
#' @importFrom odin odin
#'
#' @export


run_model <- function(age=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60),
                      EIR=10,
                      ft=0.4,
                      admin2="Tororo",
                      time=365){

  mpl <- model_param_list_create()

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector=age,EIR=EIR,ft=ft,
                                   model_param_list = mpl,het_brackets=5,
                                   admin_unit = admin2)

  # create odin generator
  odin_model_path <- system.file("extdata/odin_model.R",package="hanojoel")
  gen <- odin::odin(odin_model_path,verbose=FALSE, build = TRUE)

  # create model with initial values
  mod <- gen(user=state, use_dde=TRUE)
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
  ret <- ggplot(melt(df,id.vars="t"),
                aes(x=t,y=.data$value,col=.data$variable)) +
    geom_line() +
    scale_colour_manual(
      values=c("#000000","#CC0000","#339900","#3333FF","#CC9933","#CC3399"),
      name="Compartment",
      labels=c("S","T","D","A","U","P")
    ) +
    scale_x_continuous(breaks=brks,labels=0:(length(brks)-1)) +
    xlab("Years") + ylab("Proportion of population")
  return(list("plot"=ret,"dat"=out))

}
