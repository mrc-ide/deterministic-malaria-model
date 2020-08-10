
#' Minimal error
#'
#' @return
#' @export
#'
#' @examples
#'
minimal_error <- function(){
  set.seed(10)

  # Read in ivermectin hazards
  haz <- read.table(system.file("extdata/hazards_FINAL_020818.txt", package = "ICDMM"), header = TRUE)
  haz3_300 <- haz$d300[1:20] # 3x300

  # Define age categories, EIR and time period
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  init_EIR <- 18
  time_period <- 418

  # Generate model specific ivermectin and MDA parameters
  ivm_parms0 = ICDMM::ivm_fun(IVM_start_times = 50000,
                              time_period = time_period,
                              hazard_profile = haz3_300,
                              ivm_coverage=0.8,   # proportion of population in age range defined below that get ivm
                              ivm_min_age=4,      # mininum age of those targeted for IVM
                              ivm_max_age=150,    # maximum age of those targeted for IVM
                              bites_per_3_days=2) # assume mosquitoes feed twice every 3 days


  # this is no MDA (indicated by MDA_times = NULL)
  mda_parms0 = ICDMM::mda_fun(ncc = 2, # this is the number of groups for having correlated coverage between round - always leave as 2
                              MDA_times = NULL,   # timing of MDA rounds
                              MDA_grp_prop = 0.9,  # to do with the correlated coverage - leave as 0.9
                              MDA_cov = 0.8,       # this sets the coverage - and is proportion of population in age range defined below that get the drug each round
                              MDA_st_age = 1,      # lower age range (i.e. 1+ year olds get DP MDA)
                              MDA_en_age = 80,     # upper age range for MDA (set large to ensure all no upper limit)
                              init_age = init_age)

  m_mda_0 <- ICDMM::create_r_model(odin_model_path = system.file("extdata/old_models/odin_model_parity_GOOD_drugs_newP.R",package = "ICDMM"),
                                   het_brackets = 5,
                                   age = init_age,
                                   init_EIR = init_EIR,
                                   init_ft = 0.5,
                                   country = "Gambia",
                                   admin2 = "Upper River",
                                   mu0 = 0.1,
                                   ...=c(list(age02 = 10, # Index in age vector at which age > 2
                                              age10 = 14, # Index in age vector at which age > 10
                                              bites_Bed_new = 0.3, # No idea what these do, you need
                                              bites_Indoors_new = 0.2, # to give them values
                                              feedback = 0,  # FIX
                                              net_switch_time = 1000,  # FIX
                                              Q0_new = 0.5,  # FIX
                                              rP_DP = 0.5,# FIX
                                              rP_SP = 0.5,# FIX
                                              ccov = c(0.5, 0.5)),# FIX
                                         mda_parms0, ivm_parms0))

  path <- system.file("odin", package = "ICDMM", mustWork = TRUE)
  possible <- sub("\\.json$", "", dir(path, pattern = "\\.json$"))
  env <- asNamespace("ICDMM")
  model <- get("odin_model_parity_GOOD_drugs_newP", envir = env, mode = "function", inherits = FALSE)
  mod <- model(user = m_mda_0$state, use_dde = TRUE, unused_user_action = "message")

  ## WILL ERROR AT END OF RUN
  res <- mod$run(1:50, verbose = TRUE)
  op<- mod$transform_variables(res)
  plot(op$t ,op$clin_inc_all*1000, type="l", lwd=2, col="grey",
       xlab="Day", ylab="Clinical incidence per 1,000", ylim = c(0, max(op$clin_inc_all)*1000))
}
