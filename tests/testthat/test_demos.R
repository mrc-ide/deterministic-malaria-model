context("test-demos_run.R")
source("eqm_soln_varying_nets.R")

test_that("model_run demo runs", {
  set.seed(1234)
  # define input parameters
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,
                7.5,10,15,20,30,40,50,60,70,80)
  init_EIR <- 10
  admin_str <- NULL
  time_period <- 30*1
  prop_treated <- 0.4
  # run the model
  model_run <- run_model_example(age=init_age, EIR=init_EIR, ft = prop_treated,
                         admin2 = admin_str, time = time_period)
  # objects out correct
  expect_equal(c("plot", "dat"), names(model_run))
  expect_true(is.list(model_run$dat))
  expect_true(all(class(model_run$plot) == c("gg", "ggplot")))
  # time handled right
  expect_equal(length(model_run$dat$t), 31)
  # equilibrium init check indirectly (though this could be risky)
  expect_equal(model_run$dat$prev[10]-model_run$dat$inc[15], 0.2287737)
})

test_that("create_r_model demo runs", {
  # define input parameters
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,
                7.5,10,15,20,30,40,50,60,70,80)
  init_EIR <- 10
  time_period <- 30*1
  prop_treated <- 0.4
  # creates the odin model
  out <- run_model(odin_model_path = system.file("extdata/odin_model.R",
                                                                package = "ICDMM"),
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  init_ft = prop_treated,
                                  country = NULL,
                                  admin2 = NULL)

  expect_equal(out$prev[10]-out$inc[15], 0.2287737)
})

test_that("compare model outputs", {
  # define input parameters
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,
                7.5,10,15,20,30,40,50,60,70,80)
  init_EIR <- 10
  time_period <- 30*1
  prop_treated <- 0.4
  # creates the odin model
  out <- run_model(odin_model_path = system.file("extdata/odin_model.R",
                                                                package = "ICDMM"),
                                  num_int = 1,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  init_ft = prop_treated,
                                  country = NULL,
                                  admin2 = NULL)
  expect_equal(out$prev[10]-out$inc[15], 0.2287737)
})

