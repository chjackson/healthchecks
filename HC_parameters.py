# -*- coding: utf-8 -*-
"""


HC model parameter file

containing all parameters used by the main model:

- internal parameters which are not changeable unless for a good reason (eg. values for QRISK calculation)
- uncertain parameters (ie parameters which are observed with a CI)
- assumed values (eg fraction of CVD events diagnosed)

@author: Arno Steinacher and Chris Jackson
"""



import numpy as np

P = {} # parameter dictionary: point estimates 
PU = {} # standard errors or prior parameters defining uncertainty distribution

# global HC parameters
P['HC_offered'] = 0.197 # offered HC among eligible in 2014/15
P['HC_takeup'] = 0.488 # of those offered, how many follow up HC
P['HC_takeup_rr_prev_att'] = 1.0 #  (RR for those attending at their previous invitation or those who attended despite not being eligible)
P['HC_takeup_rr_not_prev_att'] = 1.0 # (RR for those not attending at their previous invitation)
P['HC_offer_prev_att'] = 0.197
P['HC_offer_not_prev_att'] = 0.197
P['HC_age_limit'] = [40, 74] # min / max age for HC eligibility

#############################################################
# HC uptake odds ratios 
P['gender_vec'] = np.array([1, 1.08158]) # male, female
P['age_vec'] = np.array([0.6266264, 1, 1.508133, 1.653677])
P['eth_vec'] = np.array([1, 1.0288336, 0.9912068,  0.8741651])
P['SES_vec'] = np.array([1, 0.9742240, 0.9871075, 1.0066707, 1.2491194])
P['smoker_vec'] = np.array([1, 0.7436315])
P['QRisk_vec'] = np.array([0.6427019, 1, 1.5406470, 1.8957039, 2.3402212])

#############################################################
# HC uptake extra relative rates (not odds ratios) in particular groups
# to be changed in scenarios 
P['ses5_extra_uptake'] = 1 
P['sm_extra_uptake'] = 1 
P['q5_extra_uptake'] = 1 

#############################################################
# HC uptake log odds ratios and SEs on the log scale
P['gender_log'] = np.array([0.0, 0.07842257])
P['age_log'] = np.array([-0.4674048,  0.0,  0.4108643,  0.5025921 ])
P['eth_log'] = np.array([ 0.0,  0.028425776, -0.008832112, -0.134486029 ])
P['SES_log'] = np.array([ 0.0, -0.026114027, -0.012976328,  0.006648585,  0.222438803 ])
P['smoker_log'] = np.array([ 0.0, -0.2962096  ])
P['QRisk_log'] = np.array([-0.4420743,  0.0,  0.4322024,  0.6395902,  0.8502454  ])

PU['gender_log'] = {'se': np.array([ 0.004716215 , 0.004629377 ]) }
PU['age_log'] = {'se': np.array([ 0.005658008, 0.005819362, 0.006490005, 0.051974600 ]) }
PU['eth_log'] = {'se': np.array([ 0.003616423, 0.012766544, 0.014835230, 0.015422309 ]) }
PU['SES_log'] = {'se': np.array([0.007420764, 0.007464298, 0.007441611, 0.007416301, 0.007145419 ]) }
PU['smoker_log'] = {'se': np.array([ 0.003634693 , 0.006023610 ]) }
PU['QRisk_log'] = {'se': np.array([ 0.007592940, 0.007950855, 0.008972403, 0.010376484, 0.010375860])}


#############################################################
# takeup among noneligible people based on chronic condition

P['HC_takeup_noneligible'] = 0.05 # percentage of noneligible taking up HC / year
PU['HC_takeup_noneligible'] = {'l95': 0.02, 'betaa': 6.706576, 'betab': 124.044341}

# takeup among people on diabetes, bp and CVD registers
P['HC_include_diabetes_registers'] = False
P['HC_include_CVD_registers'] = False
P['HC_include_bp_registers'] = False


######################################################
# TREATMENT PARAMETERS
######################################################
# (1) who gets treated, expressed as proportions of whole population, not percent

# smoker referral to smoking cessation
P['HC_smoker_ref'] = 0.036
PU['HC_smoker_ref'] = {'l95': 0.033, 'betaa' :  511.73, 'betab' :  13681.49}

# referral to weight management
# for people with BMI >= 30
P['HC_weight_ref'] = 0.275
PU['HC_weight_ref'] = {'l95': 0.269, 'betaa' :  5810.977, 'betab' :  15280.544}

# proportion of QRisk individuals receiving statins
# for QRisk lower than 20
P['HC_statins_presc_Q20minus'] = 0.0205
PU['HC_statins_presc_Q20minus'] = {'l95': 0.0197, 'betaa' :  2423.776, 'betab' :  115762.264}

# for QRisk higher than 20
P['HC_statins_presc_Q20plus'] = 0.1423
PU['HC_statins_presc_Q20plus'] = {'l95': 0.1371, 'betaa' :  2429.835, 'betab' :  14608.207}

# proportion of people with high blood presure receiving antihypertensives
# for QRisk lower than 20
P['HC_aht_presc_Q20minus'] = 0.0154
PU['HC_aht_presc_Q20minus'] = {'l95' :  0.0146, 'betaa' :  1365.815, 'betab' :  87287.746}

# for QRisk higher than 20
P['HC_aht_presc_Q20plus'] = 0.0248
PU['HC_aht_presc_Q20plus'] =  {'l95' :  0.0205, 'betaa' :  113.7808, 'betab' :  4463.3501}


##########################################################
# (2) Compliance rates

P['Weight_comp'] = 0.5 # from Aveyard.    # big uncertainty - from belief.
PU['Weight_comp'] =  {'l95':  0.3, 'betaa':  11.25982, 'betab':  11.01507}

P['Statins_comp'] = 0.5
PU['Statins_comp'] = {'l95': 0.4, 'betaa': 47.29982, 'betab': 47.08141}

P['AHT_comp'] = 0.55 # proportion taking up AHTs
PU['AHT_comp'] = {'l95': 0.45, 'betaa': 52.55069, 'betab': 43.80485}


############################################################
# (3) Treatment effectiveness

P['Smoking_eff'] = 0.15 # proportions quitting at 1 year
PU['Smoking_eff'] = {'l95': 0.131, 'betaa': 192.5715, 'betab': 1080.4478}

# BMI change after weight management
P['Weight_eff'] = -1.5
PU['Weight_eff'] = {'std': 0.007}

# Statins effectiveness
# as measured as reduction in total cholesterol at 1 year
P['Statins_eff_male'] = -1.22
PU['Statins_eff_male'] = {'l95': -1.26}
PU['Statins_eff_male']['std'] = (P['Statins_eff_male'] - PU['Statins_eff_male']['l95']) / 1.96

P['Statins_eff_female'] = -1.16
PU['Statins_eff_female'] = {'l95': -1.23}
PU['Statins_eff_female']['std'] = (P['Statins_eff_female'] - PU['Statins_eff_female']['l95']) / 1.96

# also include effect of Statins on HDL (from CTC 2015, email from authors)
P['Statins_eff_HDL_male'] = 0.04
PU['Statins_eff_HDL_male'] = {'std': 0.006}

P['Statins_eff_HDL_female'] = 0.036
PU['Statins_eff_HDL_female'] = {'std': 0.012}

# Extra HR for effect of statins on top of totchol/HDL reduction
P['Statins_eff_extra_male'] =  0.82
P['Statins_logeff_extra_male'] =  np.log(P['Statins_eff_extra_male'])
PU['Statins_logeff_extra_male'] =  {'std': 0.03}
P['Statins_eff_extra_female'] = 0.87
P['Statins_logeff_extra_female'] =  np.log(P['Statins_eff_extra_female'])
PU['Statins_logeff_extra_female'] =  {'std': 0.08}

# Antihypertensive effectiveness

# age < 55 years, assume use ACEI
# reduction in systolic blood pressure
P['AHT_eff_age55minus_SBP'] = -6.29
PU['AHT_eff_age55minus_SBP'] = {'l95': -9.26}
PU['AHT_eff_age55minus_SBP']['std'] = (P['AHT_eff_age55minus_SBP'] - PU['AHT_eff_age55minus_SBP']['l95']) / 1.96

# reduction in diastolic blood pressure
P['AHT_eff_age55minus_DBP'] = -4.14
PU['AHT_eff_age55minus_DBP'] = {'l95': -5.81}
PU['AHT_eff_age55minus_DBP']['std'] = (P['AHT_eff_age55minus_DBP'] - PU['AHT_eff_age55minus_DBP']['l95']) / 1.96

# age > 55 years, assume use calcium blockers
P['AHT_eff_age55plus_SBP_male'] = -7.6
PU['AHT_eff_age55plus_SBP_male'] = {'l95': -7.95}
PU['AHT_eff_age55plus_SBP_male']['std'] = (P['AHT_eff_age55plus_SBP_male'] - PU['AHT_eff_age55plus_SBP_male']['l95']) / 1.96

P['AHT_eff_age55plus_SBP_female'] = -9.0
PU['AHT_eff_age55plus_SBP_female'] = {'l95': -9.32}
PU['AHT_eff_age55plus_SBP_female']['std'] = (P['AHT_eff_age55plus_SBP_female'] - PU['AHT_eff_age55plus_SBP_female']['l95']) / 1.96

P['AHT_eff_age55plus_DBP_male'] = -3.1
PU['AHT_eff_age55plus_DBP_male'] = {'l95': -3.45}
PU['AHT_eff_age55plus_DBP_male']['std'] = (P['AHT_eff_age55plus_DBP_male'] - PU['AHT_eff_age55plus_DBP_male']['l95']) / 1.96

P['AHT_eff_age55plus_DBP_female'] = -3.5
PU['AHT_eff_age55plus_DBP_female'] = {'l95': -3.82}
PU['AHT_eff_age55plus_DBP_female']['std'] = (P['AHT_eff_age55plus_DBP_female'] - PU['AHT_eff_age55plus_DBP_female']['l95']) / 1.96



########################################################################################
# OTHER PARAMETERS (not observed)
########################################################################################


P['CVD_diagnosis_fraction'] = 1.0 # yearly fraction of how many CVD events are diagnosed

P['bp_diagnosis_fraction'] = 0.05 # yearly fraction of how many instances of high blood pressure are diagnosed

P['Smoking_relapsing_rate'] = 0.35 # rate of people relapsing to smoking after quitting
P['Smoking_relapsing_rate_std'] = 0.35 # 

P['Statins_dropout_rate'] = 0.05 # yearly rate of people stopping to take statins
PU['Statins_dropout_rate'] = {'std': 0.01, 'betaa': 18.89142, 'betab': 354.12932}
#PU['Statins_dropout_rate'] = {'std': 0.005, 'betaa':82.25182 , 'betab':1553.85174 }

P['AHT_dropout_rate'] = 0.05 # yearly rate of people stopping to take AHTs
PU['AHT_dropout_rate'] = {'std': 0.01, 'betaa': 18.89142, 'betab': 354.12932}
#PU['AHT_dropout_rate'] = {'std': 0.005, 'betaa':82.25182 , 'betab':1553.85174 }

# proportion of individuals that are physically active
# Not used currently 
P['physically_active'] = 0.6

# age bounds in which CAIDE is used
P['CAIDE_age_bounds'] = [40,60]

# values used in first submission
#P['MI_sudden_death'] = 0.3
#P['Stroke_sudden_death'] = 0.3
#P['p_ihdevent_is_mi'] = 0.5 
#P['p_strokeevent_is_full'] = 0.5

P['MI_sudden_death_male'] = 0.321 # http://www.bmj.com/content/344/bmj.d8059  2010 figures 
P['MI_sudden_death_female'] = 0.299

# http://bmjopen.bmj.com/content/1/2/e000269 Fig 2
# smoothed, see stroke_mort.r
P['Stroke_sudden_death_female_80plus'] = 0.2364014
P['Stroke_sudden_death_male_80plus'] = 0.1720481
P['Stroke_sudden_death_female_below80'] = 0.06591648
P['Stroke_sudden_death_male_below80'] = 0.04727273 

# https://www.bhf.org.uk/publications/statistics/cvd-stats-2015 ch 2, Table 2.11 from Wales
P['p_ihdevent_is_mi_male'] = 0.5822785   # 4.6/(4.6+3.3) 
P['p_ihdevent_is_mi_female'] = 0.444444  # 2.4/(2.4+3.0)

# http://stroke.ahajournals.org/content/43/12/3343.long
P['p_strokeevent_is_full'] = 0.6004942 # 729 / (729 + 485)



##################################################################
#  Following are fixed data rather than uncertain parameters 
##################################################################

# QRISK2 weights - coming directly from QRISK algorithm

P['Iethrisk_male'] = np.array([    0, 0,
                                0.3567133647493443400000000,
                                0.5369559608176189800000000,
                                0.5190878419529624300000000,
                                0.2182992106490147000000000,
                                -0.3474174705898491800000000,
                                -0.3674730037922803700000000,
                                -0.3749664891426142700000000,
                                -0.1926947742531604500000000])
P['Ismoke_male'] = np.array([0,
                     0.2784649664157046200000000,
                     0.6067834395168959500000000,
                     0.7103835060989258700000000,
                     0.8626172339181202900000000])

P['Iethrisk_female'] = np.array([    0, 0,
                                  0.2671958047902151500000000,
                                  0.7147534261793343500000000,
                                  0.3702894474455115700000000,
                                  0.2073797362620235500000000,
                                  -0.1744149722741736900000000,
                                  -0.3271878654368842200000000,
                                  -0.2200617876129250500000000,
                                  -0.2090388032466696800000000])
P['Ismoke_female'] = np.array([0,
                       0.1947480856528854800000000,
                       0.6229400520450627500000000,
                       0.7405819891143352600000000,
                       0.9134392684576959600000000])





#########################################################################
# age/sex/ethnicity weights from census data

# proportion of age/sex/ethnicities in census 2011 data
# rows:     (0) age 35 - 40:
#           (1) age 40 - 45:
#           (2) age 45 - 50:
#           (3) age 50 - 55:
#           (4) age 55 - 60:
#           (5) age 60 - 65:
#           (6) age 65 - 70:
#           (7) age 70 - 75:

# columns:  (0) male_White, (1) male_Asian, (2) male_Black, (3) male_Mixed/other,
#           (4) female_White, (5) female_Asian, (6) female_Black, (7) female_Mixed/other


P['census_prop'] = np.array([[4.78,	0.76,	0.25,	0.26,	4.79,	0.76,	0.29,	0.22],
                            [5.00,	0.65,	0.25,	0.21,	5.06,	0.63,	0.28,	0.19],
                            [5.68,	0.51,	0.28,	0.18,	5.80,	0.51,	0.31,	0.16],
                            [5.86,	0.37,	0.26,	0.15,	5.97,	0.38,	0.29,	0.14],
                            [5.23,	0.33,	0.17,	0.10,	5.29,	0.34,	0.19,	0.09],
                            [4.68,	0.28,	0.09,	0.07,	4.77,	0.29,	0.11,	0.06],
                            [5.10,	0.18,	0.05,	0.05,	5.24,	0.22,	0.07,	0.05],
                            [4.00,	0.12,	0.05,	0.03,	4.22,	0.14,	0.06,	0.03],
                            [3.14,	0.12,	0.06,	0.03,	3.51,	0.12,	0.07,	0.03]])/100.0
P['census_ages'] = np.array([30,40,45,50,55,60,65,70,75])

P['DementiaLateRRs'] = 0 # dummy value to initialise 

# Annual decline in CVD incidence, as a relative risk for 1 year later
P['CVD_annual_inc_rr_male'] = 1.0
P['CVD_annual_inc_rr_female']  = 1.0

# case fatality decline (applies to both sudden and background mortality)
P['IHD_annual_cf_rr_male'] = 1.0
P['IHD_annual_cf_rr_female'] = 1.0
P['Stroke_annual_cf_rr_male'] = 1.0
P['Stroke_annual_cf_rr_female'] = 1.0

# decline in fatality after sudden events (published data for stroke is just for events)
P['Stroke_annual_sudden_rr_male'] = 1.0
P['Stroke_annual_sudden_rr_female'] = 1.0

P['CVD_extrap_horizon'] = 20

# Parameterise extra uncertainty about the baseline risk of the population
# through an extra relative QRisk centred around 1
# Error should be 0 in base case, and 0.1 in a sensitivity analysis, corresponding to 95% CI of exp(+- 2*0.1) = (0.82, 1.22)

P['QRisk_extra_logrr'] = 0.0
PU['QRisk_extra_logrr'] = {'std': 0.0}


## Utility decrements from Sullivan

P['eq_age'] = -0.0002747
P['eq_male'] = 0.0010046,	

# capture contrast in coefficients between poor/near poor (bottom fifth of US pop) and top four-fifths. apply to bottom fifth of UK pop.
P['eq_poor'] = -0.04

#icd410	acute myocardial infarct*
#icd411	oth ac ischemic hrt dis*
#icd412	old myocardial infarct
#icd413	angina pectoris*
#icd414	oth chr ischemic hrt dis*
#P['eq_'] = np.array([-0.0625727,	0.0131711])
#P['eq_'] = np.array([-0.0866826,	0.0842625])
#P['eq_'] = np.array([-0.0367975,	0.0257359])
#P['eq_'] = np.array([-0.0854255,	0.0134397])
#P['eq_'] = np.array([-0.0626527,	0.0130759])
# take one in the middle
P['eq_ihd'] = -0.0625727
	
## icd436	cva
P['eq_stroke'] = -0.1170501

## icd797	senility w/o psychosis
P['eq_dem'] = -0.2136477

# icd162	mal neo trachea/lung*
P['eq_lc'] = -0.1192427

## extra disutility for number of chronic conditions (Sullivan web table 1)
## we only need up to 4 

P['eq_2dis'] = -0.0615
P['eq_3dis'] = -0.0667
P['eq_4dis'] = -0.0433
# 5 -0.0287 0.0153
# 6 -0.0048 0.0178
# 7 0.0237 0.0189
# 8 0.0389 0.0217
# 9 0.0445 0.0251
# 10 0.1001 0.0285


## standard errors around utility decrements from Sullivan

PU['eq_age'] = {'std': 0.000165}
PU['eq_male'] = {'std': 0.0006241}
PU['eq_poor'] = {'std': 0.006}
PU['eq_ihd'] = {'std': 0.0131711}
PU['eq_stroke'] = {'std': 0.0121435}
PU['eq_dem'] = {'std': 0.0278856}
PU['eq_lc'] = {'std': 0.0430051}
PU['eq_2dis'] = {'std': 0.0081}
PU['eq_3dis'] = {'std': 0.0099}
PU['eq_4dis'] = {'std': 0.0118}
