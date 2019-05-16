context("test-compare_with_BM.R")
library(hanojoel)

test_that("compare_no_ints_bm", {
  bm_no_ints <- unlist(read.delim("bm_data/no_ints_prev.txt"), use.names=FALSE)
  bm_no_ints <- matrix(bm_no_ints, nrow=11, ncol=2)

  init_age <- c(0, 10, 20, 30)
  init_EIR <- 12
  time_period <- 11
  prop_treated <- 0.4
  # TODO: Fix IRN_IRS_on problem with numbers of interventions or num_ints
  wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model.R",
                                                                package = "hanojoel"),
                                  num_int = 1,
                                  het_brackets = 1,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  init_ft = prop_treated,
                                  country = NULL,
                                  admin2 = NULL)
  mod <- wh$generator(user= wh$state, use_dde = TRUE)
  mod_run <- mod$run(t = 1:time_period)
  out <- mod$transform_variables(mod_run)
  expect_equal(out$prev, bm_no_ints[,2], tolerance=1e-6)
})

