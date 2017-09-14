##'  Admin level 1 africa seasonal parameters
##'
##'  These datasets represent the data fitted within the Imperial College Malaria model for
##'  relating seasonal profiles to malaria transmission intensity at level 1 admin regions
##'  across Africa
##'
##' @docType data
##'
##' @format A dataframe of 576 observations of 13 variables:
##'
##'   \code{$admin_units_seasonal}: A dataframe of admin units and their seasonal parameters
##'     \itemize{
##'       \item country: Country string
##'       \item admin1: Admin 1 string
##'       \item map_prev_2010: 2010 Atlas map microscopy prevalence in 2-10 year olds
##'       \item id: Numeric vector of 1:576
##'       \item ID_1: Numeric vector referencing shape file geoshape ids
##'       \item a0: Average value of fourier series
##'       \item a1: First of partial series cos terms
##'       \item b1: First of partial series sine terms
##'       \item a2: Second of partial series cos terms
##'       \item b2: Second of partial series sine terms
##'       \item a3: Third of partial series cos terms
##'       \item b3: Third of partial series sine terms
##'       \item theta_c: Rainfall normalising constant
##'     }
##'
##' @rdname admin_units_seasonal
##'
##' @aliases admin_units_seasonal
##'
##'
"admin_units_seasonal"
