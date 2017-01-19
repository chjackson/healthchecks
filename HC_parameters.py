# -*- coding: utf-8 -*-
"""


HC model parameter file

containing all parameters used by the main model:

- internal parameters which are not changeable unless for a good reason (eg. values for QRISK calculation)
- uncertain parameters (ie parameters which are observed with a CI)
- assumed values (eg threshold for BP to diagnose high blood pressure)

@author: arno
"""



import numpy as np

P = {} # parameter dictionary

# global HC parameters
P['HC_offered'] = 0.197 # offered HC among eligible in 2014/15
P['HC_takeup'] = 0.488 # of those offered, how many follow up HC
P['HC_takeup_prev_att'] = 0.488 #  (same, for those attending at their previous invitation or those who attended despite not being eligible)
P['HC_takeup_not_prev_att'] = 0.488 # (same, for those not attending at their previous invitation)
P['HC_offer_prev_att'] = 0.197
P['HC_offer_not_prev_att'] = 0.197
P['HC_age_limit'] = [40, 74] # min / max age for HC eligibility

#############################################################
# HC uptake odds ratios 
P['gender'] = np.array([1, 1.08158])
P['age'] = np.array([0.6266264, 1, 1.508133, 1.653677])
P['ethnicity'] = np.array([1, 1.0288336, 0.9912068,  0.8741651])
P['SES'] = np.array([1, 0.9742240, 0.9871075, 1.0066707, 1.2491194])
P['smoker'] = np.array([1, 0.7436315])
P['QR'] = np.array([0.6427019, 1, 1.5406470, 1.8957039, 2.3402212])

#############################################################
# HC uptake log odds ratios and SEs on the log scale
P['gender_log'] = np.array([0.0, 0.07842257])
P['age_log'] = np.array([-0.4674048,  0.0,  0.4108643,  0.5025921 ])
P['ethnicity_log'] = np.array([ 0.0,  0.028425776, -0.008832112, -0.134486029 ])
P['SES_log'] = np.array([ 0.0, -0.026114027, -0.012976328,  0.006648585,  0.222438803 ])
P['smoker_log'] = np.array([ 0.0, -0.2962096  ])
P['QR_log'] = np.array([-0.4420743,  0.0,  0.4322024,  0.6395902,  0.8502454  ])

P['gender_selog'] = np.array([ 0.004716215 , 0.004629377 ])
P['age_selog'] = np.array([ 0.005658008, 0.005819362, 0.006490005, 0.051974600 ])
P['ethnicity_selog'] = np.array([ 0.003616423, 0.012766544, 0.014835230, 0.015422309 ])
P['SES_selog'] = np.array([0.007420764, 0.007464298, 0.007441611, 0.007416301, 0.007145419 ])
P['smoker_selog'] = np.array([ 0.003634693 , 0.006023610 ])
P['QR_selog'] = np.array([ 0.007592940, 0.007950855, 0.008972403, 0.010376484, 0.010375860])


#############################################################
# takeup among noneligible people based on chronic condition

P['HC_takeup_noneligible'] = 0.05 # percentage of noneligible taking up HC / year
P['HC_takeup_noneligible95'] = 0.02 #  95 % belief
P['HC_takeup_noneligible_betaa'] = 6.706576
P['HC_takeup_noneligible_betab'] = 124.044341


# takeup among people on diabetes, bp and CVD registers
P['HC_include_diabetes_registers'] = False
P['HC_include_CVD_registers'] = False
P['HC_include_bp_registers'] = False


######################################################
# TREATMENT PARAMETERS
######################################################
# (1) who gets treated, expressed as proportions of whole population, not percent

# smoker referral to smoking cessation
P['HC_smoker_referral_rate'] = 0.036
P['HC_smoker_referral95'] = 0.033
P['HC_smoker_referral_betaa'] = 511.73
P['HC_smoker_referral_betab'] = 13681.49

# referral to weight management
# for people with BMI >= 30
P['HC_weight_referral_rate'] = 0.275
P['HC_weight_referral95'] = 0.269
P['HC_weight_referral_betaa'] = 5810.977
P['HC_weight_referral_betab'] = 15280.544 

# proportion of QRisk individuals receiving statins
# for QRisk lower than 20
P['HC_Q20minus_statins'] = 0.0205
P['HC_Q20minus_statins95'] = 0.0197
P['HC_Q20minus_statins_betaa'] = 2423.776
P['HC_Q20minus_statins_betab'] = 115762.264

# for QRisk higher than 20
P['HC_Q20plus_statins'] = 0.1423
P['HC_Q20plus_statins95'] = 0.1371
P['HC_Q20plus_statins_betaa'] = 2429.835
P['HC_Q20plus_statins_betab'] = 14608.207

# proportion of people with high blood presure receiving antihypertensives
# for QRisk lower than 20
P['HC_Q20minus_aht'] = 0.0154
P['HC_Q20minus_aht95'] = 0.0146
P['HC_Q20minus_aht_betaa'] = 1365.815
P['HC_Q20minus_aht_betab'] = 87287.746

# for QRisk higher than 20
P['HC_Q20plus_aht'] = 0.0248
P['HC_Q20plus_aht95'] = 0.0205
P['HC_Q20plus_aht_betaa'] = 113.7808
P['HC_Q20plus_aht_betab'] = 4463.3501


##########################################################
# (2) Compliance rates

P['Weight_compliance'] = 0.5 # from Aveyard
P['Weight_compliance95'] = 0.3 # big uncertainty - from belief.
P['Weight_compliance_betaa'] = 11.25982
P['Weight_compliance_betab'] = 11.01507

P['Statins_compliance'] = 0.5
P['Statins_compliance95'] = 0.4
P['Statins_compliance_betaa'] = 47.29982
P['Statins_compliance_betab'] = 47.08141

P['AHT_compliance'] = 0.55 # proportion taking up AHTs
P['AHT_compliance95'] = 0.45
P['AHT_compliance_betaa'] = 52.55069
P['AHT_compliance_betab'] = 43.80485


############################################################
# (3) Treatment effectiveness

P['Smoking_eff'] = 0.15 # proportions quitting at 1 year
P['Smoking_eff95'] = 0.131
P['Smoking_eff_betaa'] = 192.5715
P['Smoking_eff_betab'] = 1080.4478

# BMI change after weight management
P['Weight_eff'] = -1.5
P['Weight_eff_std'] = 0.007 

# Statins effectiveness
# as measured as reduction in total cholesterol at 1 year
P['Statins_eff_male'] = -1.22
P['Statins_eff_male95'] = -1.26
P['Statins_eff_male_std'] = (P['Statins_eff_male'] - P['Statins_eff_male95']) / 1.96

P['Statins_eff_female'] = -1.16
P['Statins_eff_female95'] = -1.23
P['Statins_eff_female_std'] = (P['Statins_eff_female'] - P['Statins_eff_female95']) / 1.96

# also include effect of Statins on HDL (from CTC 2015, email from authors)
P['Statins_eff_HDL_male'] = 0.04
P['Statins_eff_HDL_male_std'] = 0.006

P['Statins_eff_HDL_female'] = 0.036
P['Statins_eff_HDL_female_std'] = 0.012

# Extra HR for effect of statins on top of totchol/HDL reduction
P['Statins_eff_extra_male'] =  0.82
P['Statins_eff_extra_female'] = 0.87

# Antihypertensive effectiveness

# age < 55 years, assume use ACEI
# reduction in systolic blood pressure
P['AHT_eff_age55minus_SBP'] = -6.29
P['AHT_eff_age55minus_SBP95'] = -9.26
P['AHT_eff_age55minus_SBP_std'] = (P['AHT_eff_age55minus_SBP'] - P['AHT_eff_age55minus_SBP95']) / 1.96

# reduction in diastolic blood pressure
P['AHT_eff_age55minus_DBP'] = -4.14
P['AHT_eff_age55minus_DBP95'] = -5.81
P['AHT_eff_age55minus_DBP_std'] = (P['AHT_eff_age55minus_DBP'] - P['AHT_eff_age55minus_DBP95']) / 1.96

# age > 55 years, assume use calcium blockers
P['AHT_eff_age55plus_SBP_male'] = -7.6
P['AHT_eff_age55plus_SBP_male95'] = -7.95
P['AHT_eff_age55plus_SBP_male_std'] = (P['AHT_eff_age55plus_SBP_male'] - P['AHT_eff_age55plus_SBP_male95']) / 1.96

P['AHT_eff_age55plus_SBP_female'] = -9.0
P['AHT_eff_age55plus_SBP_female95'] = -9.32
P['AHT_eff_age55plus_SBP_female_std'] = (P['AHT_eff_age55plus_SBP_female'] - P['AHT_eff_age55plus_SBP_female95']) / 1.96

P['AHT_eff_age55plus_DBP_male'] = -3.1
P['AHT_eff_age55plus_DBP_male95'] = -3.45
P['AHT_eff_age55plus_DBP_male_std'] = (P['AHT_eff_age55plus_DBP_male'] - P['AHT_eff_age55plus_DBP_male95']) / 1.96

P['AHT_eff_age55plus_DBP_female'] = -3.5
P['AHT_eff_age55plus_DBP_female95'] = -3.82
P['AHT_eff_age55plus_DBP_female_std'] = (P['AHT_eff_age55plus_DBP_female'] - P['AHT_eff_age55plus_DBP_female95']) / 1.96


### OLD weight parameters for testing

#P['Weight_compliance'] = 0.58 # proportion completing treatment
#P['Weight_compliance95'] = 0.331
#P['Weight_eff_c'] = -2.0
#P['Weight_eff_c95'] = -2.02
#P['Weight_eff_nc'] = -0.7
#P['Weight_eff_nc95'] = -0.72
#P['Weight_eff_ncstd'] = 0.01



########################################################################################
# OTHER PARAMETERS (not observed)
########################################################################################


P['CVD_diagnosis_fraction'] = 1.0 # yearly fraction of how many CVD events are diagnosed

P['bp_diagnosis_fraction'] = 0.05 # yearly fraction of how many instances of high blood pressure are diagnosed

P['Smoking_relapsing_rate'] = 0.35 # rate of people relapsing to smoking after quitting

P['Statins_dropout_rate'] = 0.05 # yearly rate of people stopping to take statins

P['AHT_dropout_rate'] = 0.05 # yearly rate of people stopping to take AHTs

# proportion of individuals that are physically active
P['physically_active'] = 0.6

# age bounds in which CAIDE is used
P['CAIDE_age_bounds'] = [40,60]

P['MI_sudden_death'] = 0.3
P['Stroke_sudden_death'] = 0.3
P['p_ihdevent_is_mi'] = 0.5 
P['p_strokeevent_is_full'] = 0.5
P['CVD_background_CFR_reduction'] = True




##################################################################
#  DO NOT CHANGE ANY OF THE FOLLOWING:
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
