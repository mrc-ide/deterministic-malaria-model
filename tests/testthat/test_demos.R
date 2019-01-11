test_that("model_run demo runs", {

  set.seed(1234)

  # create a vector of age categories
  init_age <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,
                7.5,10,15,20,30,40,50,60,70,80)

  # provide a value of the annual EIR for this model run
  init_EIR <- 10

  # provide a string for the admin 2 unit
  admin_str <- "Tororo"

  # provide the length of time (in days) that you want to run the model for
  time_period <- 30*1

  # provide a value for the proportion of cases that are treatedr)
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
  expect_equal(-0.0001534032, model_run$dat$inc05[30] - model_run$dat$inc05[1])
})
