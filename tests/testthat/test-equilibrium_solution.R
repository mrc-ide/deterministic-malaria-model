context("test-equilibrium_solution.R")
library(hanojoel)

age_vector = c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)
mpl <- model_param_list_create()
eqm_soln <- equilibrium_init_create(age_vector=age_vector, het_brackets=5, country = NULL, admin_unit = NULL,
                                    ft=0.4, EIR=10, model_param_list=mpl)

bm_age <- unlist(read.delim("C:/Users/hunwin/Desktop/age.txt")[1, 2:14], use.names=FALSE)
bm_age_width <- unlist(read.delim("C:/Users/hunwin/Desktop/age_width.txt")[1, 2:13], use.names=FALSE)
bm_age_rate <- unlist(read.delim("C:/Users/hunwin/Desktop/age_rate.txt")[1, 2:14], use.names=FALSE)

#TODO: change this test to reflect how data is spat out at end of equilibrium soltion
test_that("equilibrium age data", {
  expect_identical(eqm_soln$age, bm_age/365)
  expect_identical(eqm_soln$age_width, bm_age_width)
  expect_equal(eqm_soln$age_rate, bm_age_rate, tolerance=1e-8)
})
