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
  expect_equal(length(model_run$dat$t), 31L)
  # equilibrium init check indirectly (though this could be risky)
  expect_equal(1.196338e-05, model_run$dat$inc05[30] - model_run$dat$inc05[1])
})

test_that("create_r_model demo runs", {
  # define input parameters
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,
                7.5,10,15,20,30,40,50,60,70,80)
  init_EIR <- 10
  time_period <- 30*1
  prop_treated <- 0.4
  # creates the odin model
  wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model.R",package = "hanojoel"),
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
  expect_equal(1.196338e-05, out$inc05[30] - out$inc05[1])
})

test_that("compare model outputs", {
  # define input parameters
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,
                7.5,10,15,20,30,40,50,60,70,80)
  init_EIR <- 10
  time_period <- 30*1
  prop_treated <- 0.4
  # creates the odin model
  wh <- hanojoel:::create_r_model(odin_model_path = system.file("extdata/odin_model.R",package = "hanojoel"),
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
  expect_equal(model_run$dat$prev, out$prev, tolerance= 1e-12)
})
