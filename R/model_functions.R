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

  mpl <- model_param_list_create(num_int = 1)

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector=age,EIR=EIR,ft=ft,
                                   model_param_list = mpl,het_brackets=5,
                                   admin_unit = admin2)

  # create odin generator
  odin_model_path <- system.file("extdata/odin_model.R",package="hanojoel")
  gen <- odin::odin(odin_model_path,verbose=FALSE)

  # There are many parameters used that should not be passed through
  # to the model.
  state_use <- state[names(state) %in% names(formals(gen))]

  # create model with initial values
  mod <- gen(user=state_use, use_dde=TRUE)
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

#------------------------------------------------
#' load_file
#'
#' \code{load_file} loads package file
#'
#' @description Load a file from within the inst/extdata folder of the
#'   hanojoel package. File extension must be one of .csv, .txt, or .rds.
#'
#' @param name the name of a file within the inst/extdata folder.
#'
#' @importFrom utils read.csv
#'
#' @export

load_file <- function(name) {

  # check that valid file extension
  ext <- strsplit(name, "\\.")[[1]]
  ext <- ext[length(ext)]
  if(is.element(ext, c("csv", "rds")) == FALSE){
    stop("file extension not valid")
  }

  # get full file path
  name_full <- system.file("extdata/", name, package="hanojoel", mustWork = TRUE)

  # read in file
  if (ext == "rds") {
    ret <- readRDS(name_full)
  } else {
    ret <-  read.csv(file=name_full, header=TRUE, sep=",")
  }

  return(ret)
}

#------------------------------------------------
#' match admin region
#'
#' @param country Character for country within which admin2 is in.
#'   Default = NULL
#' @param admin Character for admin region. Some fuzzy logic will be used to
#'   match. If not provided then no seasonality is introduced. Default = NULL

admin_match <- function(admin_unit = NULL, country = NULL,
                        admin_units_seasonal = admin_units_seasonal){

  # intiialise admin match as no match
  admin_matches <- 0

  if(!is.null(admin_unit)){

    # if there is no country given then search for the admin unit
    if(is.null(country)){

      # find exact match
      admin_matches <- grep(paste("^",admin_unit,"\\b",sep=""),admin_units_seasonal$admin1)
      # if exact does not match try fuzzy match up to dist of 4 which should catch having nop spaces or separators etc
      if(length(admin_matches)==0){
        admin_matches <- which(adist(admin_units_seasonal$admin1,admin_unit)<=4)
      }
      if(length(admin_matches)>1) stop("Admin unit string specified is ambiguous without country")

      # if we do have a country though find that match first and then find admin
    } else {

      # first find an exact match
      country_matches <- grep(paste("^",country,"\\b",sep=""), admin_units_seasonal$country)
      if(length(unique(admin_units_seasonal$country[country_matches]))==1){
        chosen_country <- unique(admin_units_seasonal$country[country_matches])
      } else if(length(unique(admin_units_seasonal$country[country_matches]))==0){
        # if exact does not match try fuzzy match up to dist of 2 which should catch having no spaces or separators etc
        country_matches <- which(adist(admin_units_seasonal$country,y = country)<=2)
        if(length(unique(admin_units_seasonal$country[country_matches]))==1){
          chosen_country <- unique(admin_units_seasonal$country[country_matches])
        } else if(length(unique(admin_units_seasonal$country[country_matches]))==0) stop ("Country string specified not close enough to those in database")
      }

      # find exact match
      admin_sub_matches <- grep(paste("^",admin_unit,"\\b",sep=""),admin_units_seasonal$admin1[country_matches])
      # if exact does not match try fuzzy match up to dist of 4 which should catch having nop spaces or separators etc
      if(length(admin_sub_matches)==0){
        admin_sub_matches <- which(adist(admin_units_seasonal$admin1[country_matches],admin_unit)<=1)
      }
      if(length(admin_sub_matches)>1) stop("Admin unit string specified is not close enougth to those in the database")

      admin_matches <- which(admin_units_seasonal$admin1 == admin_units_seasonal$admin1[country_matches][admin_sub_matches])
    }

  }

  if(admin_matches == 0){
    admin_matches = NULL
  }

  return (admin_matches)
}
