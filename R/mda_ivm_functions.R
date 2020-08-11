

############################################################################################################################
# make a function that returns everthing you'd need for MDA/SMC
############################################################################################################################

mda_fun <- function(ncc = 2, MDA_times = NULL,MDA_grp_prop = 0.9, MDA_cov = 0.8, MDA_st_age = 0.5, MDA_en_age = 80, MDA_succ=0.95, MDA_drug_choice = 2){

  # set up MDA times
  L1 = as.list(MDA_times)
  if(length(MDA_times) < 9) L1 = c(L1, as.list(rep(1e5, 9 - length(MDA_times))))
  names(L1) = paste0("MDA_t",1:9)
  if(length(MDA_times) >9)print("too many MDA times")

  # translate the user defined ages in to age categories
  MDA_st_cat = which.min(abs(init_age - MDA_st_age))
  MDA_en_cat = which.min(abs(init_age - MDA_en_age))
  if(!MDA_st_age %in% init_age) print(paste("Not a pre-defined age break point.",init_age[MDA_st_cat],"used instead"))
  if(!MDA_en_age %in% init_age) print(paste("Not a pre-defined age break point.",init_age[MDA_en_cat],"used instead"))

  op = c(L1,list(
            MDA_st_cat=MDA_st_cat,
            MDA_en_cat=MDA_en_cat,
            ncc=ncc,
            MDA_grp_prop=MDA_grp_prop,
            MDA_cov = MDA_cov,
            MDA_succ = MDA_succ,
            mda_drug_choice = MDA_drug_choice))

  return(op)

}

############################################################################################################################
# make a function that returns everything you'd need for ivermectin
############################################################################################################################

ivm_fun <- function(IVM_start_times, time_period, hazard_profile, ivm_coverage=0.8, ivm_min_age=5, ivm_max_age = 200, bites_per_3_days=1){

  # function to make the ivermectin time profile
  make_IVRM = function(IVRM_st, eff_len,ttt){
    IVRM =  rep(max(ttt), length(ttt))
    for(i in 1:length(IVRM_st)){
      IVRM[ttt >= IVRM_st[i]-1 & ttt < IVRM_st[i] + eff_len -1] = IVRM_st[i]
    }
    return(IVRM)
  }

  ttt = 0:time_period
  eff_len = length(hazard_profile)
  IVRM_start = make_IVRM(IVM_start_times, eff_len, ttt)

  op = list(ttt = ttt,
            eff_len = eff_len,
            haz = hazard_profile,
            ivm_cov_par = ivm_coverage,
            ivm_min_age = ivm_min_age,
            ivm_max_age = ivm_max_age,
            B2 = bites_per_3_days,
            IVRM_start=IVRM_start
            )
  return(op)

}

############################################################################################################################
# update parameter list without
############################################################################################################################


new_params <- function(replace_params, statename){
  statename = statename[-which(names(statename) %in% names(replace_params))]
  statename_new = c(statename,replace_params)
  return(statename_new)
}


############################################################################################################################
# calculate effective mortality from ivermectin model output
############################################################################################################################

effective_mortality <- function(op){
  x1 = op$mu*(op$Sv + op$Sv_F1 + op$Sv_F2 + rowSums(op$Ev_F1) + rowSums(op$Ev_F2) + op$Iv_F1 + op$Iv_F2)

  mmat = matrix(0, nrow = length(op$Sv), ncol = ncol(op$mu_vi))

  for(j in 1:ncol(op$mu_vi)){
    mmat[,j] = op$mu_vi[1,j] * (op$Sx_F1[,j] + op$Sx_F2[,j] + op$Ix_F1[,j] + op$Ix_F2[,j] + rowSums(op$Ex_F1[,,j]) + rowSums(op$Ex_F2[,,j]))
  }
  x2 = rowSums(mmat)
  return((x1 + x2)/op$mv)
}
