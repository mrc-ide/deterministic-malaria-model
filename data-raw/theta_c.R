## This code computes the theta_c parameter needed for the site seasonality parameterisation
## HJT Unwin

#' \code{fourier} computes fourier series to be used later for seasonality in carrying capacity
#' @param x Vector of times
#' @param ss Vector of fourier parameters in this order a0, a1, b1, a2, b2, a3, b3
#'
#' @export
#'
## Fourier series
fourier <- function(x, ss) {
  two_pi <- 6.2831853071796
  raw <- (ss[1] + ss[2]*cos(two_pi*x/365) + ss[4]*cos(2*two_pi*x/365) + ss[6]*cos(3*two_pi*x/365)
           + ss[3]*sin(two_pi*x/365) + ss[5]*sin(2*two_pi*x/365) + ss[7]*sin(3*two_pi*x/365))

  return (raw)
  }

#---------------------------------------------------------------------------------------------------
## function to create annual season patterns
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

#---------------------------------------------------------------------------------------------------
library(readr)
library(dplyr)

data <- read_csv("data-raw/admin_units_seasonal.csv", col_types = "dcccddddddd")

# Initialises empty vector for theta c
theta_c <- vector(length=nrow(data))

# Computes theta c for each admin1 region
for (i in 1:nrow(data)){

  # Make list of fourier coefficients
  ss <- c(data$seasonal_a0[i], data$seasonal_a1[i], data$seasonal_b1[i],
          data$seasonal_a2[i], data$seasonal_b2[i], data$seasonal_a3[i],
          data$seasonal_b3[i])
  seasonality_list <- seasonality(ss)

  # Add theta_c to a vector
  theta_c[i] <- seasonality_list[[2]]

  # Plots seasonality curve
  #plot(1:365, seasonality_list[[1]])
}

# Extracts data from CSV to form data for deterministic malaria model
country <- data$NAME_0
admin1 <- data$NAME_1
ID_1 <- data$DIDE_CODE
a0 <- data$seasonal_a0
a1 <- data$seasonal_a1
a2 <- data$seasonal_a2
a3 <- data$seasonal_a3
b1 <- data$seasonal_b1
b2 <- data$seasonal_b2
b3 <- data$seasonal_b3

# Saves rda file of data.
admin_units_seasonal <- data.frame(country, admin1, ID_1, a0, a1, b1, a2, b2, a3, b3, theta_c)
saveRDS(admin_units_seasonal, file = "inst/extdata/admin_units_seasonal.rds")
