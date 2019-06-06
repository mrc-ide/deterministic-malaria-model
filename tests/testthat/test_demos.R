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
  model_run <- run_model(age=init_age, EIR=init_EIR, ft = prop_treated,
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
  wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model.R",
                                                                package = "hanojoel"),
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  init_ft = prop_treated,
                                  country = NULL,
                                  admin2 = NULL)

  # generates model functions with initial state data
  mod <- wh$generator(user= wh$state, use_dde = TRUE)
  # Runs the model
  mod_run <- mod$run(t = 1:time_period)
  out <- mod$transform_variables(mod_run)
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
  wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model.R",
                                                                package = "hanojoel"),
                                  num_int = 1,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  init_ft = prop_treated,
                                  country = NULL,
                                  admin2 = NULL)
  # generates model functions with initial state data
  mod <- wh$generator(user= wh$state, use_dde = TRUE)
  # Runs the model
  mod_run <- mod$run(t = 1:(time_period+1))
  out <- mod$transform_variables(mod_run)
  model_run <- run_model(age=init_age, EIR=init_EIR, ft = prop_treated,
                         admin2 = NULL, time = time_period)
  expect_equal(model_run$dat$prev, out$prev, tolerance= 1e-8)
  expect_equal(out$prev[10]-out$inc[15], 0.2287737)
})

test_that("compare varying itns and not", {
  # define input parameters
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,
                7.5,10,15,20,30,40,50,60,70,80)
  init_EIR <- 10
  time_period <- 30

  # Specify coverage as a coverage after a time
  wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model.R",
                                                                package = "hanojoel"),
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  itn_cov = 0.3,
                                  ITN_IRS_on = 20,
                                  num_int = 2,
                                  country = NULL,
                                  admin2 = NULL)
  mod <- wh$generator(user= wh$state, use_dde = TRUE)
  mod_run <- mod$run(t = 1:(time_period))
  out <- mod$transform_variables(mod_run)

  # Specify coverage as a vector
  wh2 <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model_itn.R",package = "hanojoel"),
                                   het_brackets = 5,
                                   age = init_age,
                                   init_EIR = init_EIR,
                                   num_int = 2,
                                   t_vector = c(-25, 20),
                                   itn_vector = c(0, 0.3),
                                   ITN_IRS_on = 20,
                                   pop_split = c(0.5, 0.5),
                                   country = NULL,
                                   admin2 = NULL)
  wh2 <- edit_equilibrium_varying_nets(wh=wh2)
  mod2 <- wh2$generator(user= wh2$state, use_dde = TRUE)
  mod_run2 <- mod2$run(t = 1:(time_period))
  out2 <- mod2$transform_variables(mod_run2)

  # Specify coverage as old coverage vector
  wh3 <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model_itn.R",package = "hanojoel"),
                                   het_brackets = 5,
                                   age = init_age,
                                   init_EIR = init_EIR,
                                   num_int = 2,
                                   t_vector = c(-25, 20),
                                   itn_vector = c(0, 0.3),
                                   ITN_IRS_on = 20,
                                   pop_split = c(0.7, 0.3),
                                   country = NULL,
                                   admin2 = NULL)
  wh3 <- edit_equilibrium_varying_nets(wh=wh3)
  mod3 <- wh3$generator(user= wh3$state, use_dde = TRUE)
  mod_run3 <- mod3$run(t = 1:(time_period))
  out3 <- mod3$transform_variables(mod_run3)

  expect_equal(out$prev, out2$prev, tolerance=1e-5)
  expect_equal(out$prev, out3$prev, tolerance=1e-5)
  expect_equal(out$inc, out2$inc, tolerance=1e-5)
  expect_equal(out$inc, out3$inc, tolerance=1e-5)

})
