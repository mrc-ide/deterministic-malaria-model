#------------------------------------------------
#' run_model_example
#'
#' \code{run_model_example} runs model using declared age, EIR, ft, country, admin
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
#' @importFrom stats coef
#' @useDynLib ICDMM
#'
#' @export


run_model_example <- function(age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60),
                      EIR = 10,
                      ft = 0.4,
                      admin2 = "Tororo",
                      time = 365){

  mpl <- model_param_list_create(num_int = 1)

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector = age,
                                   EIR = EIR,
                                   ft = ft,
                                   model_param_list = mpl,
                                   het_brackets = 5,
                                   admin_unit = admin2)

  # create odin generator
  generator <- odin_model

  # There are many parameters used that should not be passed through
  # to the model.
  state_use <- state[names(state) %in% coef(generator)$name]

  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, 1)

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
#'   ICMDMM package. File extension must be one of .csv, .txt, or .rds.
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
  name_full <- system.file("extdata/", name, package="ICDMM", mustWork = TRUE)

  # read in file
  if (ext == "rds") {
    ret <- readRDS(name_full)
  } else {
    ret <-  read.csv(file=name_full, header=TRUE, sep=",")
  }

  return(ret)
}

#------------------------------------------------
#' match clean
#'
#' @param a First string to compare. Default = NULL
#' @param b Second string to compare. Default = NULL
#'
#' @importFrom RecordLinkage levenshteinSim
#'
#' @export

match_clean <- function(a, b){

  a <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(b, "latin-ascii")))

  ret <- which(b %in% a)

  if(length(ret) == 0){
    distance <- levenshteinSim(a, b)
    ret <- which.max(distance)
  }
  return(ret)
}

#------------------------------------------------
#' match admin region
#'
#' \code{admin_match} Matches the user input admin unit and country with data
#'
#' @param country Character for country within which admin unit is in.
#'   Default = NULL
#' @param admin_unit Character for admin region. Some fuzzy logic will be used to
#'   match. If not provided then no seasonality is introduced. Default = NULL
#' @param admin_units_seasonal Dataframe of seasonality data for country and admin unit
#'
#' @export

admin_match <- function(admin_unit = NULL, country = NULL,
                        admin_units_seasonal){

  # intialise admin match as no match
  admin_matches <- 0

  if (!is.null(admin_unit)) {

    # if there is no country given then search for the admin unit
    if (is.null(country)) {

      # find exact match
      admin_matches <- which(tolower(admin_units_seasonal$admin1) %in% tolower(admin_unit))

      # if exact does not match try closest match
      if (length(admin_matches) < 1) {
        admin_matches <- match_clean(admin_unit, admin_units_seasonal$admin1)
      } else if (length(admin_matches) > 1){
        stop('Please specify the country of admin unit.  There are multiple with same name.')
      }

      # if we do have a country though find that match first and then find admin
    } else {

      # first find an exact match
      country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country))

      if (length(country_matches) < 1) {
        country_name <- admin_units_seasonal$country[match_clean(country, admin_units_seasonal$country)]
        country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country_name))
      }

      sub_admin_units_seasonal <- admin_units_seasonal[country_matches,]

      # find exact match
      admin_sub_matches <- which(tolower(sub_admin_units_seasonal$admin1) %in% tolower(admin_unit))

      # if exact does not match try closest match
      if (length(admin_sub_matches) != 1) {
        admin_sub_matches <- match_clean(admin_unit,
                                         sub_admin_units_seasonal$admin1)
      } else if (length(admin_sub_matches) > 1){
        stop('There are multiple admins with same name - check the data file!')
      }

      admin_matches <- country_matches[admin_sub_matches]
    }

    message("Requested: ", admin_unit,
            "\nReturned: ", admin_units_seasonal$admin1[admin_matches], ", ",
            admin_units_seasonal$country[admin_matches])
  }

  return(admin_matches)
}

#------------------------------------------------
#' Generate Fourier
#'
#' \code{fourier} computes fourier series to be used later for seasonality in carrying capacity
#' @param x Vector of times
#' @param ss Vector of fourier parameters in this order a0, a1, b1, a2, b2, a3, b3
#'
#' @export
#'
#'
fourier <- function(x, ss) {
  two_pi <- 6.2831853071796
  raw <- (ss[1] + ss[2]*cos(two_pi*x/365) + ss[4]*cos(2*two_pi*x/365) + ss[6]*cos(3*two_pi*x/365)
          + ss[3]*sin(two_pi*x/365) + ss[5]*sin(2*two_pi*x/365) + ss[7]*sin(3*two_pi*x/365))

  return (raw)
}

#------------------------------------------------
#' Seasonality
#'
#' \code{seasonality} Generate seasonality scaling vector for one year, to be applied to carrying capacity.
#' @param ss Vector of fourier parameters in this order a0, a1, b1, a2, b2, a3, b3
#'
#' @export
#'
#'
seasonality <- function(ss){

  # define vector of times spanning one year
  tvec = 1:365

  # calculate Fourier series curve
  seasonality <- sapply(tvec, fourier, ss=ss)
  theta_c <- sum(seasonality)/365
  seasonality <- seasonality/theta_c

  # ensure that scaling factor never goes below zero (this can happen in practice
  # because we are only using the first few terms in an infinite series)
  seasonality[seasonality<0.001] <- 0.001

  seasonality_list <- list(seasonality, theta_c)

  return (seasonality_list)

}
