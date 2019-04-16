context("test-equilibrium_solution.R")
library(hanojoel)

age_vector = c(0, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60)
mpl <- model_param_list_create(eta=(1/(21*365)), rP=1/20)
eqm_soln <- equilibrium_init_create(age_vector=age_vector, het_brackets=3, country = NULL, admin_unit = NULL,
                                    ft=0.4, EIR=10, model_param_list=mpl)

bm_age <- unlist(read.delim("bm_data/age.txt")[1, 2:14], use.names=FALSE)
bm_age_width <- unlist(read.delim("bm_data/age_width.txt")[1, 2:13], use.names=FALSE)
bm_age_rate <- unlist(read.delim("bm_data/age_rate.txt")[1, 2:14], use.names=FALSE)
bm_den <- unlist(read.delim("bm_data/den.txt")[1, 2:14], use.names=FALSE)

test_that("equilibrium age data", {
  expect_identical(eqm_soln$age, bm_age)
  expect_identical(eqm_soln$age_width, bm_age_width)
  expect_equal(eqm_soln$age_rate, bm_age_rate, tolerance=1e-8)
  expect_equal(eqm_soln$den, bm_den, tolerance=1e-6)
})

bm_foi_age <- unlist(read.delim("bm_data/foi_age.txt")[1, 2:14], use.names=FALSE)
bm_rel_foi <- unlist(read.delim("bm_data/rel_foi.txt")[1, 2:4], use.names=FALSE)

test_that("equilibrium foi data", {
  expect_equal(eqm_soln$foi_age, bm_foi_age, tolerance=1e-6)
  expect_equal(eqm_soln$omega, 0.73987, tolerance=1e-6)
  expect_equal(eqm_soln$rel_foi, bm_rel_foi, tolerance=1e-5)
})

bm_eir_eq <- unlist(read.delim("bm_data/eir_eq.txt")[1, 2:40], use.names=FALSE)

test_that("equilibrium eir", {
  expect_equal(as.vector(t(eqm_soln$EIR)), bm_eir_eq, tolerance=1e-6)
})

bm_x_I <- unlist(read.delim("bm_data/x_I.txt")[1, 2:14], use.names=FALSE)
bm_ib_eq <- unlist(read.delim("bm_data/ib_eq.txt")[1, 2:40], use.names=FALSE)
bm_foi_eq <- unlist(read.delim("bm_data/foi_eq.txt")[1, 2:40], use.names=FALSE)
bm_id_eq <- unlist(read.delim("bm_data/id_eq.txt")[1, 2:40], use.names=FALSE)
bm_ica_eq <- unlist(read.delim("bm_data/ica_eq.txt")[1, 2:40], use.names=FALSE)
bm_icm_eq <- unlist(read.delim("bm_data/icm_eq.txt")[1, 2:40], use.names=FALSE)
bm_icm_init_eq <- unlist(read.delim("bm_data/icm_init_eq.txt")[1, 2:4], use.names=FALSE)
bm_ca_eq <- unlist(read.delim("bm_data/ca_eq.txt")[1, 2:40], use.names=FALSE)

# Note tolerances are quite high because absolute values of some of the numbers are large.  Sometimes only have 4dp
test_that("equilibrium immunity", {
  expect_equal(eqm_soln$x_I, bm_x_I, tolerance=1e-2)
  expect_equal(as.vector(t(eqm_soln$init_IB[,,1])), bm_ib_eq, tolerance=1e-4)
  expect_equal(as.vector(t(eqm_soln$FOI)), bm_foi_eq, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$init_ID[,,1])), bm_id_eq, tolerance=1e-4)
  expect_equal(as.vector(t(eqm_soln$init_ICA[,,1])), bm_ica_eq, tolerance=1e-4)
  expect_equal(as.vector(t(eqm_soln$init_ICM[,,1])), bm_icm_eq, tolerance=1e-4)
  expect_equal(eqm_soln$ICM_init_eq, bm_icm_init_eq, tolerance=1e-4)
  expect_equal(as.vector(t(eqm_soln$cA_eq)), bm_ca_eq, tolerance=1e-7)
})

bm_foivij_eq <- unlist(read.delim("bm_data/foivij_eq.txt")[1, 2:40], use.names=FALSE)
bm_betas <- unlist(read.delim("bm_data/betaS.txt")[1, 2:40], use.names=FALSE)
bm_betaa <- unlist(read.delim("bm_data/betaA.txt")[1, 2:40], use.names=FALSE)
bm_betau <- unlist(read.delim("bm_data/betaU.txt")[1, 2:40], use.names=FALSE)

test_that("equilibrium humans", {
  expect_equal(as.vector(t(eqm_soln$FOIvij_eq)), bm_foivij_eq, tolerance=1e-8)
  expect_equal(as.vector(t(eqm_soln$betaS)), bm_betas, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$betaA)), bm_betaa, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$betaU)), bm_betau, tolerance=1e-7)
})

test_that("equilibrium mosquito", {
  expect_equal(eqm_soln$FOIv_eq, 0.00701093, tolerance=1e-7)
  expect_equal(eqm_soln$init_Iv, 0.0134728, tolerance=1e-6)
  expect_equal(eqm_soln$init_Sv, 0.949566, tolerance=1e-6)
})

test_that("equilibrium lavae", {
  expect_equal(eqm_soln$init_PL, 0.832831, tolerance=1e-4)
  expect_equal(eqm_soln$init_LL, 5.58968, tolerance=1e-4)
  expect_equal(eqm_soln$init_EL, 192.078, tolerance=1e-4)
})

s <- unlist(read.delim("bm_data/S.txt")[1, 2:40], use.names=FALSE)
t <- unlist(read.delim("bm_data/T.txt")[1, 2:40], use.names=FALSE)
d <- unlist(read.delim("bm_data/D.txt")[1, 2:40], use.names=FALSE)
a <- unlist(read.delim("bm_data/A.txt")[1, 2:40], use.names=FALSE)
u <- unlist(read.delim("bm_data/U.txt")[1, 2:40], use.names=FALSE)
p <- unlist(read.delim("bm_data/P.txt")[1, 2:40], use.names=FALSE)

test_that("equilibrium inital compartments no coverage", {
  expect_equal(as.vector(t(eqm_soln$init_S[,,1])), s, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$init_T[,,1])), t, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$init_D[,,1])), d, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$init_A[,,1])), a, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$init_U[,,1])), u, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln$init_P[,,1])), p, tolerance=1e-7)
})

cov <- 0.5
mpl_cov <- model_param_list_create(eta=(1/(21*365)), rP=1/20, itn_cov=cov, num_int=2)
eqm_soln_cov <- equilibrium_init_create(age_vector=age_vector, het_brackets=3, country = NULL, admin_unit = NULL,
                                        ft=0.4, EIR=10, model_param_list=mpl_cov)

test_that("equilibrium inital compartments no coverage", {
  expect_equal(as.vector(t(eqm_soln_cov$init_S[,,1])), s*cov, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln_cov$init_T[,,1])), t*cov, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln_cov$init_D[,,1])), d*cov, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln_cov$init_A[,,1])), a*cov, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln_cov$init_U[,,1])), u*cov, tolerance=1e-7)
  expect_equal(as.vector(t(eqm_soln_cov$init_P[,,1])), p*cov, tolerance=1e-7)
})
