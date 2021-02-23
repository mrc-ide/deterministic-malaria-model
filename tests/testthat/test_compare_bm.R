context("test-compare_with_BM.R")

test_that("compare_no_ints_bm", {
  bm_no_ints <- unlist(read.delim(system.file("testdata/bm_data/no_ints_prev.txt", package="ICDMM")), use.names=FALSE)
  bm_no_ints <- matrix(bm_no_ints, nrow=11, ncol=2)

  init_age <- c(0, 10, 20, 30)
  init_EIR <- 12
  time_period <- 11
  prop_treated <- 0.4
  # TODO: Fix IRN_IRS_on problem with numbers of interventions or num_ints
  wh <- ICDMM::run_model(num_int = 1,time = 10,
                                  het_brackets = 1,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  init_ft = prop_treated,
                                  country = NULL,
                                  admin2 = NULL)
  expect_equal(wh$prev, bm_no_ints[,2], tolerance=1e-6)
})

