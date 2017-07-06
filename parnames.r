## {'up_age_vec': array([ 0.63300008,  1.00254344,  1.50906548,  1.50040911]),
##     'up_gender_vec': array([ 0.99869263,  1.07980477]),
##     'up_HC_offer_not_prev_att': 0.197,
##     'up_QRisk_vec': array([ 0.64004603,  0.98777947,  1.5542865 ,  1.8741683 ,  2.31162219]),
##     'up_AHT_eff_age55minus_SBP': -7.312418727360097,
##     'up_HC_aht_presc_Q20plus': 0.022637501962254048,
##     'up_HC_takeup': 0.488,
##     'up_Statins_eff_extra_female': 0.87,
##     'up_Weight_eff': -1.50734176401132,
##     'up_HC_include_bp_registers': False,
##     'up_HC_statins_presc_Q20plus': 0.13952075851815965,
##     'up_Statins_eff_HDL_female': 0.04679206135008977,
##     'up_AHT_eff_age55plus_DBP_female': -3.1834320880978786,
##     'up_AHT_eff_age55plus_SBP_female': -9.17374616717331,
##     'up_AHT_eff_age55plus_SBP_male': -7.5726882673638425,
##     'up_HC_offer_prev_att': 0.197,
##     'up_HC_include_CVD_registers': False,
##     'up_bp_diagnosis_fraction': 0.05,
##     'up_Statins_comp': 0.5036726874809458,
##     'up_HC_aht_presc_Q20minus': 0.014718518984658679,
##     'up_Weight_comp': 0.30083675530591525,
##     'up_Statins_eff_male': -1.2240288710032297,
##     'up_HC_include_diabetes_registers': False,
##     'up_HC_smoker_ref': 0.03741321360811903,
##     'up_AHT_eff_age55minus_DBP': -4.011668283081843,
##     'up_Smoking_eff': 0.15521265716762533,
##     'up_Statins_dropout_rate': 0.05,
##     'up_HC_takeup_noneligible': 0.04813336062269343,
##     'up_AHT_dropout_rate': 0.05,
##     'up_HC_statins_presc_Q20minus': 0.020250955214036817,
##     'up_Statins_eff_extra_male': 0.82,
##     'up_Statins_eff_female': -1.0966198203069308,
##     'up_CVD_diagnosis_fraction': 1.0,
##     'up_HC_age_limit': [40, 74],
##     'up_AHT_comp': 0.6616448784017878,
##     'up_Smoking_relapsing_rate': 0.35,
##     'up_HC_offered': 0.197,
##     'up_physically_active': 0.6,
##     'up_eth_vec': array([ 0.99970082, 1.02063108, 0.99056265,  0.86775503]),
##     'up_AHT_eff_age55plus_DBP_male': -3.0217952479856995,
##     'up_HC_takeup_rr_not_prev_att': 1.0,
##     'up_smoker_vec': array([ 1.00018187, 0.74182106]), 'up_SES_vec': array([ 0.9902975 ,  0.98067816,  0.99360265,  1.01951527, 1.2491194 ]),
##     'up_Statins_eff_HDL_male': 0.03385041475235148,
##     'up_HC_weight_ref': 0.2758317948940889,
##     'up_HC_takeup_rr_prev_att': 1.0}

parnames <- 
c(paste0('up_age_vec',1:4),
  paste0('up_gender_vec',1:2),
  'up_HC_offer_not_prev_att',
  paste0('up_QRisk_vec',1:5),
  'up_AHT_eff_age55minus_SBP',
  'up_HC_aht_presc_Q20plus',
  'up_HC_takeup',
  'up_Statins_eff_extra_female',
  'up_Weight_eff',
  'up_HC_include_bp_registers',
  'up_HC_statins_presc_Q20plus',
  'up_Statins_eff_HDL_female',
  'up_AHT_eff_age55plus_DBP_female',
  'up_AHT_eff_age55plus_SBP_female',
  'up_AHT_eff_age55plus_SBP_male',
  'up_HC_offer_prev_att',
  'up_HC_include_CVD_registers',
  'up_bp_diagnosis_fraction',
  'up_Statins_comp',
  'up_HC_aht_presc_Q20minus',
  'up_Weight_comp',
  'up_Statins_eff_male',
  'up_HC_include_diabetes_registers',
  'up_HC_smoker_ref',
  'up_AHT_eff_age55minus_DBP',
  'up_Smoking_eff',
  'up_Statins_dropout_rate',
  'up_HC_takeup_noneligible',
  'up_AHT_dropout_rate',
  'up_HC_statins_presc_Q20minus',
  'up_Statins_eff_extra_male',
  'up_Statins_eff_female',
  'up_CVD_diagnosis_fraction',
  paste0('up_HC_age_limit',1:2),
  'up_AHT_comp',
  'up_Smoking_relapsing_rate',
  'up_HC_offered',
  'up_physically_active',
  paste0('up_eth_vec',1:4),
  'up_AHT_eff_age55plus_DBP_male',
  'up_HC_takeup_rr_not_prev_att',
  paste0('up_smoker_vec',1:2),
  paste0('up_SES_vec',1:5),
  'up_Statins_eff_HDL_male',
  'up_HC_weight_ref',
  'up_HC_takeup_rr_prev_att')

fixedpars <- c("up_HC_offer_not_prev_att",
               "up_HC_takeup",
               "up_HC_include_bp_registers",
               "up_HC_offer_prev_att", "up_HC_include_CVD_registers", 
               "up_bp_diagnosis_fraction", "up_HC_include_diabetes_registers", 
               "up_CVD_diagnosis_fraction", "up_HC_age_limit1", "up_HC_age_limit2", 
               "up_Smoking_relapsing_rate", "up_HC_offered", "up_physically_active", 
               "up_HC_takeup_rr_not_prev_att", "up_HC_takeup_rr_prev_att")

## Some parameters in this array are fixed. 
## Some of the fixed ones surely have uncertainty in them. 
## and are learnable with more research: 

## 5% yearly dropout rates for statins and AHTs (done now) 
## extra statins effects beyond chol lowering  (done now) 
## yearly rate of relapsing after quitting smoking 
## yearly proportion of CVD / high BP events diagnosed 
