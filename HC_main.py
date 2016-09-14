# -*- coding: utf-8 -*-
"""

HEALTH CHECKS MODEL

Definitive version stored on github from 30 Aug 2016 onwards

Chris's changes against Arno's v 7.11 for initial github upload
* Add CalculateLY, SimulateAndResults
* Make nprocs an argument to main class
* Make randpars an argument
* Note v.712 edits with instant death rate not included

Any further changes will be documented via github

@author:
Arno Steinacher as2441@cam.ac.uk
Chris Jackson chj20@mrc-bsu.cam.ac.uk

@maintainer
Chris Jackson chj20@mrc-bsu.cam.ac.uk

"""


import os
import numpy as np
import random as rnd
import HC_additional_functions as hadd
from HC_parameters import P as parms
import time
import multiprocessing
import copy




# get rid of division through zero warnings from numpy
np.seterr(all='ignore')

class HealthChecksModel:

    def __init__(self, parent=None, population_size=1000, simulation_time=1, HealthChecks = False, randseed=0, randpars=False, nprocs=10, verbose=False):

        # all output printed during simulation
        self.verbose = verbose
        t0 = time.time()
        if self.verbose==True:
            print('initialising...')
        self.population_size = population_size  # population size in individuals
        self.simulation_time = simulation_time  # simulation time in years
        self.P = np.ndarray(())
        self.npydir = 'data/matching_npys' # directory in which npy arrays for the matching process are stored
        self.Health_Checks = HealthChecks # True if Health Checks process is simulated, False if simulated without Health Checks




        # random seeds

        np.random.seed(0)
        self.randseed = randseed # global random seed
        rnd.seed(self.randseed)



        # load parameters from HSE file
        self.LoadPopulationParameters()
        # load model parameters
        #self.parms = copy.deepcopy(parms)
        if randpars:
            self.ChangeUncertainParameters()
        else:
            self.ResetUncertainParameters()
        # load time series from ELSA data
        self.LoadELSA_bp()
        self.LoadELSA_bmi()
        self.LoadELSA_smoking()
        self.LoadELSA_chol()
        self.LoadELSA_glyhb()
        # load life table data
        self.LoadLifeTableData()
        # split QRisk and CAIDE
        self.SplitCAIDEInAnnualEvents()
        self.SplitQRiskInAnnualEvents()

        # create a new population
        self.InitialisePopulation()

        # establish age weights
        self.EstablishAgeWeights()

        # compute disease prevalences and mortality rates from LT/ONS data
        self.CalculateAllDiseaseMortalityRates()

        # measure time that this initialisaion routine takes
        self.t_init_elapsed = time.time() - t0
        hrs,mins,sec = hadd.SplitTime(self.t_init_elapsed)

        # initiate multiprocessing queue
        self.out_q = multiprocessing.Queue()
        # number of parallel processes
        self.nprocs = nprocs


       ###################################################
        # FOR CHECKING PURPOSES ONLY! DELETE BLOCK BELOW IF NOT NEEDED FURTHER
        self.bmi_entry_order = []
        self.bp_entry_order = []
        self.chol_entry_order = []
        self.smoke_entry_order = []
        self.glyhb_entry_order = []

        # under HC scenarios, are treatment effects for those offered treatment
        # applied to risk factors? Default: True, but can be set to False for checking purposes
        self.Treatment = True
        self.Treatment_Statins = True
        self.Treatment_WeightReduction = True
        self.Treatment_SmokingCessation = True
        self.Treatment_AHT = True
        self.QuitMat = np.zeros((self.population_size,self.simulation_time),dtype=bool)
        self.adjust_smoking_relapse_rates = True
        self.adjust_smoking_quit_rates = True

        # Dictionary containing information of which functions are controlled by seeded random operators
        self.RandomSeed = {'Mortality_Other': True,
                           'Mortality_Diseases': True,
                           'Disease_Incidences': True,
                           'Start_Of_Each_Timestep': True}


    def ResetUncertainParameters(self):
        '''defining parameters that are uncertain and making them changeable globally
        all numbers are given as ratios scaled to 1 (ie. 0 = 0%, 1 = 100%)'''


        # reset global parms variable
        self.parms = copy.deepcopy(parms)

        UP = {} # dictionary containing all uncertain parameters

        ######################################################
        # load parameters that are changeable within model from parms

        self.up_HC_offered = self.parms['HC_offered']; UP['up_HC_offered'] = self.up_HC_offered
        self.up_HC_takeup = self.parms['HC_takeup']; UP['up_HC_takeup'] = self.up_HC_takeup

        self.up_age_vec = self.parms['age'];         UP['up_age_vec'] = self.up_age_vec
        self.up_gender_vec = self.parms['gender'];   UP['up_gender_vec'] = self.up_gender_vec
        self.up_eth_vec = self.parms['ethnicity'];   UP['up_eth_vec'] = self.up_eth_vec
        self.up_SES_vec = self.parms['SES'];         UP['up_SES_vec'] = self.up_SES_vec
        self.up_smoker_vec = self.parms['smoker']; UP['up_smoker_vec'] = self.up_smoker_vec
        self.up_QRisk_vec = self.parms['QR'];      UP['up_QRisk_vec'] = self.up_QRisk_vec

        self.up_HC_takeup_noneligible = self.parms['HC_takeup_noneligible']
        UP['up_HC_takeup_noneligible'] = self.up_HC_takeup_noneligible

        self.up_HC_include_diabetes_registers = self.parms['HC_include_diabetes_registers']
        UP['up_HC_include_diabetes_registers'] = self.up_HC_include_diabetes_registers

        self.up_HC_include_CVD_registers = self.parms['HC_include_CVD_registers']
        UP['up_HC_include_CVD_registers'] = self.up_HC_include_CVD_registers

        self.up_HC_include_bp_registers = self.parms['HC_include_bp_registers']
        UP['up_HC_include_bp_registers'] = self.up_HC_include_bp_registers

        self.up_HC_smoker_ref = self.parms['HC_smoker_referral_rate']
        UP['up_HC_smoker_ref'] = self.up_HC_smoker_ref

        self.up_HC_weight_ref = self.parms['HC_weight_referral_rate']
        UP['up_HC_weight_ref']= self.up_HC_weight_ref

        self.up_HC_statins_presc_Q20minus = self.parms['HC_Q20minus_statins']
        UP['up_HC_statins_presc_Q20minus'] = self.up_HC_statins_presc_Q20minus

        self.up_HC_statins_presc_Q20plus = self.parms['HC_Q20plus_statins']
        UP['up_HC_statins_presc_Q20plus'] = self.up_HC_statins_presc_Q20plus

        self.up_HC_aht_presc_Q20minus = self.parms['HC_Q20minus_aht']
        UP['up_HC_aht_presc_Q20minus'] = self.up_HC_aht_presc_Q20minus

        self.up_HC_aht_presc_Q20plus = self.parms['HC_Q20plus_aht']
        UP['up_HC_aht_presc_Q20plus'] = self.up_HC_aht_presc_Q20plus

        self.up_Weight_comp = self.parms['Weight_compliance']
        UP['up_Weight_comp'] = self.up_Weight_comp

        self.up_Statins_comp = self.parms['Statins_compliance']
        UP['up_Statins_comp'] = self.up_Statins_comp

        self.up_AHT_comp = self.parms['AHT_compliance']
        UP['up_AHT_comp'] = self.up_AHT_comp

        self.up_Smoking_eff = self.parms['Smoking_eff']
        UP['up_Smoking_eff'] = self.up_Smoking_eff

        self.up_Statins_eff_male = self.parms['Statins_eff_male']
        UP['up_Statins_eff_male'] = self.up_Statins_eff_male

        self.up_Statins_eff_female = self.parms['Statins_eff_female']
        UP['up_Statins_eff_female'] = self.up_Statins_eff_female

        self.up_Statins_eff_HDL_male = self.parms['Statins_eff_HDL_male']
        UP['up_Statins_eff_HDL_male'] = self.up_Statins_eff_HDL_male

        self.up_Statins_eff_HDL_female = self.parms['Statins_eff_HDL_female']
        UP['up_Statins_eff_HDL_female'] = self.up_Statins_eff_HDL_female

        self.up_AHT_eff_age55minus_SBP = self.parms['AHT_eff_age55minus_SBP']
        UP['up_AHT_eff_age55minus_SBP'] = self.up_AHT_eff_age55minus_SBP

        self.up_AHT_eff_age55minus_DBP = self.parms['AHT_eff_age55minus_DBP']
        UP['up_AHT_eff_age55minus_DBP'] = self.up_AHT_eff_age55minus_DBP

        self.up_AHT_eff_age55plus_SBP_male = self.parms['AHT_eff_age55plus_SBP_male']
        UP['up_AHT_eff_age55plus_SBP_male'] = self.up_AHT_eff_age55plus_SBP_male

        self.up_AHT_eff_age55plus_SBP_female = self.parms['AHT_eff_age55plus_SBP_female']
        UP['up_AHT_eff_age55plus_SBP_female'] = self.up_AHT_eff_age55plus_SBP_female

        self.up_AHT_eff_age55plus_DBP_male = self.parms['AHT_eff_age55plus_DBP_male']
        UP['up_AHT_eff_age55plus_DBP_male'] = self.up_AHT_eff_age55plus_DBP_male

        self.up_AHT_eff_age55plus_DBP_female = self.parms['AHT_eff_age55plus_DBP_female']
        UP['up_AHT_eff_age55plus_DBP_female'] = self.up_AHT_eff_age55plus_DBP_female

        self.up_CVD_diagnosis_fraction = self.parms['CVD_diagnosis_fraction']
        UP['up_CVD_diagnosis_fraction'] = self.up_CVD_diagnosis_fraction

        self.up_bp_diagnosis_fraction = self.parms['bp_diagnosis_fraction']
        UP['up_bp_diagnosis_fraction'] = self.up_bp_diagnosis_fraction

        self.up_Smoking_relapsing_rate = self.parms['Smoking_relapsing_rate']
        UP['up_Smoking_relapsing_rate'] = self.up_Smoking_relapsing_rate

        self.up_Statins_dropout_rate = self.parms['Statins_dropout_rate']
        UP['up_Statins_dropout_rate'] = self.up_Statins_dropout_rate

        self.up_AHT_dropout_rate = self.parms['AHT_dropout_rate']
        UP['up_AHT_dropout_rate'] = self.up_AHT_dropout_rate

        self.up_HC_age_limit = self.parms['HC_age_limit']
        UP['up_HC_age_limit'] = self.up_HC_age_limit

        self.up_physically_active = self.parms['physically_active']
        UP['up_physically_active'] = self.up_physically_active

        ## CJ NEW STUFF
        self.up_HC_takeup_prev_att = self.parms['HC_takeup_prev_att']
        UP['up_HC_takeup_prev_att'] = self.up_HC_takeup_prev_att
        self.up_HC_takeup_not_prev_att = self.parms['HC_takeup_not_prev_att']
        UP['up_HC_takeup_not_prev_att'] = self.up_HC_takeup_not_prev_att
        self.up_HC_offer_prev_att = self.parms['HC_offer_prev_att']
        UP['up_HC_offer_prev_att'] = self.up_HC_offer_prev_att
        self.up_HC_offer_not_prev_att = self.parms['HC_offer_not_prev_att']
        UP['up_HC_offer_not_prev_att'] = self.up_HC_offer_not_prev_att
        self.up_Statins_eff_extra_male = self.parms['Statins_eff_extra_male']
        UP['up_Statins_eff_extra_male'] = self.up_Statins_eff_extra_male
        self.up_Statins_eff_extra_female = self.parms['Statins_eff_extra_female']
        UP['up_Statins_eff_extra_female'] = self.up_Statins_eff_extra_female
        self.up_cvd_sudden_death = self.parms['CVDevent_sudden_death']
        UP['up_cvd_sudden_death'] = self.up_cvd_sudden_death
        self.up_cvd_background_cfr_reduction = self.parms['CVD_background_CFR_reduction']
        UP['up_cvd_background_cfr_reduction'] = self.up_cvd_background_cfr_reduction
        self.up_Weight_eff = self.parms['Weight_eff']
        UP['up_Weight_eff'] = self.up_Weight_eff

        # make list of all uncertain parameters global
        self.UP = copy.deepcopy(UP)
        # save set of default parameters in copy list, to compare against if parameters were changed
        self.UP_R = copy.deepcopy(UP)
        # reset global parms variable




    def ChangeUncertainParameters(self):
        '''draw uncertain parameters from normal distributions around their specified mean and std'''

        # as a first step, reset all parameters to default, then draw from distributions
        self.ResetUncertainParameters()

        for i in range(4):
            self.up_age_vec[i] = np.exp(np.random.normal(self.parms['age_log'][i], self.parms['age_selog'][i]))
        for i in range(2):
            self.up_gender_vec[i] = np.exp(np.random.normal(self.parms['gender_log'][i], self.parms['gender_selog'][i]))
        for i in range(4):
            self.up_eth_vec[i] = np.exp(np.random.normal(self.parms['ethnicity_log'][i], self.parms['ethnicity_selog'][i]))
        for i in range(4):
            self.up_SES_vec[i] = np.exp(np.random.normal(self.parms['SES_log'][i], self.parms['SES_selog'][i]))
        for i in range(2):
            self.up_smoker_vec[i] = np.exp(np.random.normal(self.parms['smoker_log'][i], self.parms['smoker_selog'][i]))
        for i in range(5):
            self.up_QRisk_vec[i] = np.exp(np.random.normal(self.parms['QR_log'][i], self.parms['QR_selog'][i]))

        self.up_HC_takeup_noneligible = np.random.beta(self.parms['HC_takeup_noneligible_betaa'], self.parms['HC_takeup_noneligible_betab'])
        self.up_HC_smoker_ref = np.random.beta(self.parms['HC_smoker_referral_betaa'], self.parms['HC_smoker_referral_betab'])
        self.up_HC_weight_ref = np.random.beta(self.parms['HC_weight_referral_betaa'], self.parms['HC_weight_referral_betab'])
        self.up_HC_statins_presc_Q20minus = np.random.beta(self.parms['HC_Q20minus_statins_betaa'], self.parms['HC_Q20minus_statins_betab'])
        self.up_HC_statins_presc_Q20plus = np.random.beta(self.parms['HC_Q20plus_statins_betaa'], self.parms['HC_Q20plus_statins_betab'])
        self.up_HC_aht_presc_Q20minus = np.random.beta(self.parms['HC_Q20minus_aht_betaa'], self.parms['HC_Q20minus_aht_betab'])
        self.up_HC_aht_presc_Q20plus = np.random.beta(self.parms['HC_Q20plus_aht_betaa'], self.parms['HC_Q20plus_aht_betab'])

        self.up_Weight_comp  = np.random.beta(self.parms['Weight_compliance_betaa'], self.parms['Weight_compliance_betab'])
        self.up_Statins_comp = np.random.beta(self.parms['Statins_compliance_betaa'], self.parms['Statins_compliance_betab'])
        self.up_AHT_comp = np.random.beta(self.parms['AHT_compliance_betaa'], self.parms['AHT_compliance_betab'])

        self.up_Smoking_eff = np.random.beta(self.parms['Smoking_eff_betaa'], self.parms['Smoking_eff_betab'])
        self.up_Weight_eff = np.random.normal(self.parms['Weight_eff'], self.parms['Weight_eff_std'])

        self.up_Statins_eff_male = np.random.normal(self.parms['Statins_eff_male'], self.parms['Statins_eff_male_std'])
        self.up_Statins_eff_female = np.random.normal(self.parms['Statins_eff_female'], self.parms['Statins_eff_female_std'])
        self.up_AHT_eff_age55minus_SBP = np.random.normal(self.parms['AHT_eff_age55minus_SBP'], self.parms['AHT_eff_age55minus_SBP_std'])
        self.up_AHT_eff_age55minus_DBP = np.random.normal(self.parms['AHT_eff_age55minus_DBP'], self.parms['AHT_eff_age55minus_DBP_std'])
        self.up_AHT_eff_age55plus_SBP_male = np.random.normal(self.parms['AHT_eff_age55plus_SBP_male'], self.parms['AHT_eff_age55plus_SBP_male_std'])
        self.up_AHT_eff_age55plus_SBP_female = np.random.normal(self.parms['AHT_eff_age55plus_SBP_female'], self.parms['AHT_eff_age55plus_SBP_female_std'])
        self.up_AHT_eff_age55plus_DBP_male = np.random.normal(self.parms['AHT_eff_age55plus_DBP_male'], self.parms['AHT_eff_age55plus_DBP_male_std'])
        self.up_AHT_eff_age55plus_DBP_female = np.random.normal(self.parms['AHT_eff_age55plus_DBP_female'], self.parms['AHT_eff_age55plus_DBP_female_std'])
        self.up_Statins_eff_HDL_male = np.random.normal(self.parms['Statins_eff_HDL_male'], self.parms['Statins_eff_HDL_male_std'])
        self.up_Statins_eff_HDL_female = np.random.normal(self.parms['Statins_eff_HDL_female'], self.parms['Statins_eff_HDL_female_std'])

        UP = {} # dictionary containing all uncertain parameters

        UP['up_HC_offered'] = self.up_HC_offered
        UP['up_HC_takeup'] = self.up_HC_takeup
        UP['up_age_vec'] = self.up_age_vec
        UP['up_gender_vec'] = self.up_gender_vec
        UP['up_eth_vec'] = self.up_eth_vec
        UP['up_SES_vec'] = self.up_SES_vec
        UP['up_smoker_vec'] = self.up_smoker_vec
        UP['up_QRisk_vec'] = self.up_QRisk_vec
        UP['up_HC_takeup_noneligible'] = self.up_HC_takeup_noneligible
        UP['up_HC_include_diabetes_registers'] = self.up_HC_include_diabetes_registers
        UP['up_HC_include_CVD_registers'] = self.up_HC_include_CVD_registers
        UP['up_HC_include_bp_registers'] = self.up_HC_include_bp_registers
        UP['up_HC_smoker_ref'] = self.up_HC_smoker_ref
        UP['up_HC_weight_ref']= self.up_HC_weight_ref
        UP['up_HC_statins_presc_Q20minus'] = self.up_HC_statins_presc_Q20minus
        UP['up_HC_statins_presc_Q20plus'] = self.up_HC_statins_presc_Q20plus
        UP['up_HC_aht_presc_Q20minus'] = self.up_HC_aht_presc_Q20minus
        UP['up_HC_aht_presc_Q20plus'] = self.up_HC_aht_presc_Q20plus
        UP['up_Weight_comp'] = self.up_Weight_comp
        UP['up_Statins_comp'] = self.up_Statins_comp
        UP['up_AHT_comp'] = self.up_AHT_comp
        UP['up_Smoking_eff'] = self.up_Smoking_eff
        UP['up_Statins_eff_male'] = self.up_Statins_eff_male
        UP['up_Statins_eff_female'] = self.up_Statins_eff_female
        UP['up_Statins_eff_HDL_male'] = self.up_Statins_eff_HDL_male
        UP['up_Statins_eff_HDL_female'] = self.up_Statins_eff_HDL_female

        UP['up_AHT_eff_age55minus_SBP'] = self.up_AHT_eff_age55minus_SBP
        UP['up_AHT_eff_age55minus_DBP'] = self.up_AHT_eff_age55minus_DBP
        UP['up_AHT_eff_age55plus_SBP_male'] = self.up_AHT_eff_age55plus_SBP_male
        UP['up_AHT_eff_age55plus_SBP_female'] = self.up_AHT_eff_age55plus_SBP_female
        UP['up_AHT_eff_age55plus_DBP_male'] = self.up_AHT_eff_age55plus_DBP_male
        UP['up_AHT_eff_age55plus_DBP_female'] = self.up_AHT_eff_age55plus_DBP_female

        # other parameters,remaining unchanged
        UP['up_CVD_diagnosis_fraction'] = self.up_CVD_diagnosis_fraction
        UP['up_bp_diagnosis_fraction'] = self.up_bp_diagnosis_fraction
        UP['up_Smoking_relapsing_rate'] = self.up_Smoking_relapsing_rate
        UP['up_Statins_dropout_rate'] = self.up_Statins_dropout_rate
        UP['up_AHT_dropout_rate'] = self.up_AHT_dropout_rate
        UP['up_physically_active'] = self.up_physically_active
        UP['up_HC_age_limit'] = self.up_HC_age_limit

        # CJ new
        UP['up_HC_takeup_prev_att'] = self.up_HC_takeup_prev_att
        UP['up_HC_takeup_not_prev_att'] = self.up_HC_takeup_not_prev_att
        UP['up_HC_offer_prev_att'] = self.up_HC_offer_prev_att
        UP['up_HC_offer_not_prev_att'] = self.up_HC_offer_not_prev_att
        UP['up_Statins_eff_extra_male'] = self.up_Statins_eff_extra_male
        UP['up_Statins_eff_extra_female'] = self.up_Statins_eff_extra_female
        UP['up_cvd_sudden_death'] = self.up_cvd_sudden_death
        UP['up_cvd_background_cfr_reduction'] = self.up_cvd_background_cfr_reduction
        UP['up_Weight_eff'] = self.up_Weight_eff

        self.UP = UP



    def GetUncertainParameters(self):
        '''returns dictionary of all uncertain parameters as currently stored in the model'''
        return self.UP

    def SetUncertainParameter(self, parname, parval):
        '''updates an uncertain parameter to a new value'''

        self.UP[parname] = parval
        exec('self.%s = parval' % parname)
        try:
            print('Parameter %s updated to value %g' % (parname,eval('self.%s' % parname)) )
        except:
            pass

    def GetNumberOfCPUs(self):
        '''returns number of CPUs used for parallel tasks (matchin processes)'''
        return self.nprocs

    def SetNumberOfCPUs(self,cpus):
        '''sets number of CPUs used for parallel tasks (matchin processes)'''
        self.nprocs = cpus
        if self.verbose == True:
            print('Multiprocessing will use %d workers' % self.nprocs)


    def LoadELSA_bmi(self):
        '''loads BMI trajectory data from ELSA into an numpy array'''

        self.ELSA_bmi = np.genfromtxt('data/ELSA_bmimatch2.csv', delimiter=',', names=True, dtype=None)
        self.ELSA_bmi_numeric = np.genfromtxt('data/ELSA_bmimatch2.csv', delimiter=',', skip_header=1)


        arraynames = ['C41_bmi', 'C81_bmi', 'C41_alternative_bmi', 'C81_alternative_bmi']
        # first check whether the compressed matching arrays are already present in the data folder
        # if so, load them, if not, compute them.

        loadnew = 0

        for a in range(len(arraynames)):
            filename = '%s/%s.npy' % (self.npydir,arraynames[a])
            if os.path.isfile(filename):
                # file present, load it
                exec('self.%s = np.load(filename)' % (arraynames[a]))

            else:
                loadnew = 1
                print 'npy file %s not found, recalculating...' % arraynames[a]

        if loadnew == 1:
            print 'Routine to recalculate BMI matching template not yet implemented.'

            ###########################################################
            # ACTION: IMPLEMENT ROUTINE analogous to self.LoadELSA_bp,
            # based upon the script file 'indexing_ELSAmatch_bmi.py'



    def LoadELSA_chol(self):
        '''loads Cholesterol trajectory data from ELSA into an numpy array'''


        self.ELSA_chol = np.genfromtxt('data/ELSA_cholmatch.csv', delimiter=',', names=True, dtype=None)
        self.ELSA_chol_numeric = np.genfromtxt('data/ELSA_cholmatch.csv', delimiter=',', skip_header=1)

        arraynames = ['C1_chol', 'C1_alternative1_chol', 'C1_alternative2_chol','CH1_chol', 'CH1_alternative1_chol', 'CH1_alternative2_chol']
        # first check whether the compressed matching arrays are already present in the data folder
        # if so, load them, if not, compute them.

        loadnew = 0

        for a in range(len(arraynames)):
            filename = '%s/%s.npy' % (self.npydir,arraynames[a])
            if os.path.isfile(filename):
                # file present, load it
                exec('self.%s = np.load(filename)' % (arraynames[a]))

            else:
                loadnew = 1
                print 'npy file %s not found, recalculating...' % arraynames[a]

        if loadnew == 1:
            print 'Routine to recalculate Cholesterol matching template not yet implemented.'

            ###########################################################
            # ACTION: IMPLEMENT ROUTINE analogous to self.LoadELSA_bp,
            # based upon the script file 'indexing_ELSAmatch_chol.py'


    def LoadELSA_smoking(self):
        '''loads Smoking trajectory data from ELSA into an numpy array'''

        self.ELSA_smoke = np.genfromtxt('data/ELSA_smokematch4.csv', delimiter=',', names=True, dtype=None)
        self.ELSA_smoke_numeric = np.genfromtxt('data/ELSA_smokematch4.csv', delimiter=',', skip_header=1)


        arraynames = ['C1_smoking','C1_smoking_weights']
        # first check whether the compressed matching arrays are already present in the data folder
        # if so, load them, if not, compute them.

        loadnew = 0

        for a in range(len(arraynames)):
            filename = '%s/%s.npy' % (self.npydir,arraynames[a])
            if os.path.isfile(filename):
                # file present, load it
                exec('self.%s = np.load(filename)' % (arraynames[a]))

            else:
                loadnew = 1
                print 'npy file %s not found, recalculating...' % arraynames[a]

        if loadnew == 1:
            print 'Routine to recalculate Smoking matching template not yet implemented.'

            ###########################################################
            # ACTION: IMPLEMENT ROUTINE analogous to self.LoadELSA_bp,
            # based upon the script file 'indexing_ELSAmatch_smoking.py'



    def LoadELSA_glyhb(self):
        '''loads Blood Glucose trajectory data from ELSA into an numpy array'''


        self.ELSA_glyhb = np.genfromtxt('data/ELSA_glyhbmatch.csv', delimiter=',', names=True, dtype=None)
        self.ELSA_glyhb_numeric = np.genfromtxt('data/ELSA_glyhbmatch.csv', delimiter=',', skip_header=1)

        arraynames = ['C1_glyhb', 'C1_alternative1_glyhb', 'C1_alternative2_glyhb', 'DD1_glyhb', 'DD1_alternative1_glyhb', 'DD1_alternative2_glyhb' ]
        # first check whether the compressed matching arrays are already present in the data folder
        # if so, load them, if not, compute them.

        loadnew = 0

        for a in range(len(arraynames)):
            filename = '%s/%s.npy' % (self.npydir,arraynames[a])
            if os.path.isfile(filename):
                # file present, load it
                exec('self.%s = np.load(filename)' % (arraynames[a]))

            else:
                loadnew = 1
                print 'npy file %s not found, recalculating...' % arraynames[a]

        if loadnew == 1:

            E = self.ELSA_glyhb_numeric

            # determine bin size and values from ELSA dataset for glyhb
            age, sex, bmi, eq_nan, eq_finite, eq, diab, gly = hadd.glyhb_vectors(E)
            n_sex, n_age, n_bmi, n_eq, n_diab, n_gly = hadd.glyhb_vectorsize(E)

            # how many bin combinations for first to third matching routines are there?
            combinations = n_sex * n_age * n_bmi * n_eq * n_diab * n_gly + 1
            combinations_alternative1 =  n_bmi * n_gly * n_diab + 1
            combinations_alternative2 =  n_gly * n_diab + 1



            E_length = E.shape[0]

            # first evaluate the glyhb change over four years
            C = np.ones((combinations,E_length + 10)) * np.nan
            C_alternative1 = np.ones((combinations_alternative1,E_length + 10)) * np.nan
            C_alternative2 = np.ones((combinations_alternative2,E_length + 10)) * np.nan

            # now, also collect diabetes diagnosis
            DD = np.ones((combinations,E_length + 10)) * np.nan
            DD_alternative1 = np.ones((combinations_alternative1,E_length + 10)) * np.nan
            DD_alternative2 = np.ones((combinations_alternative2,E_length + 10)) * np.nan

            #######################################
            # first match

            for i in range(E_length):
                # evaluate combination index for row
                age_idx, sex_idx, bmi_idx, eq_idx, diab_idx, gly_idx = hadd.calculate_indices_glyhb(E, i, age, sex, bmi, eq, diab, gly)

                idx_array = hadd.calculate_Cindex_glyhb(n_sex, n_bmi, n_eq, n_diab, n_gly, age_idx, sex_idx, bmi_idx, eq_idx, diab_idx, gly_idx)

                if np.isnan(idx_array).sum() > 0:
                    print i, age_idx, sex_idx, bmi_idx, eq_idx, diab_idx, gly_idx

                # fill in bands in the idx entries
                for eq_idx2, idx in enumerate(idx_array):
                    C[idx,0] = idx
                    C[idx,1] = age[age_idx]
                    C[idx,2] = sex[sex_idx]
                    C[idx,3] = bmi[bmi_idx]
                    if idx_array.size > 1:
                        # in the case of missing data on smoking, this entry already exists
                        C[idx,4] = eq[eq_idx2]
                    else:
                        C[idx,4] = eq[eq_idx]
                    C[idx,5] = diab[diab_idx]
                    C[idx,6] = gly[gly_idx]


                    # fill in change in glyhb  for that combination
                    C[idx,i+7] = E[i,8]

                    # fill in value for diabetes diagnosis in upcoming 4 years for that combination
                    DD[idx,:7] = C[idx,:7]
                    DD[idx,i+7] = E[i,6]


            #########################################
            # second match

            for i in range(E_length):
                # evaluate combination index for row
                age_idx, sex_idx, bmi_idx, eq_idx, diab_idx, gly_idx = hadd.calculate_indices_glyhb(E, i, age, sex, bmi, eq, diab, gly)

                idx_array = hadd.calculate_Cindex_alternative1_glyhb(n_diab, n_gly, bmi_idx, diab_idx, gly_idx)

                # fill in bands in the idx entries
                idx = idx_array[0]

                C_alternative1[idx,0] = idx
                C_alternative1[idx,1] = bmi[bmi_idx]
                C_alternative1[idx,2] = diab[diab_idx]
                C_alternative1[idx,3] = gly[gly_idx]

                # fill in change in bmi  for that combination
                C_alternative1[idx,i+4] = E[i,8]


                # fill in value for diabetes diagnosis in upcoming 4 years for that combination
                DD_alternative1[idx,:4] = C_alternative1[idx,:4]
                DD_alternative1[idx,i+4] = E[i,6]

            #########################################
            # third match

            for i in range(E_length):
                # evaluate combination index for row
                age_idx, sex_idx, bmi_idx, eq_idx, diab_idx, gly_idx = hadd.calculate_indices_glyhb(E, i, age, sex, bmi, eq, diab, gly)

                idx_array = hadd.calculate_Cindex_alternative2_glyhb(n_gly, diab_idx, gly_idx)

                idx = idx_array[0]
                # fill in bands in the idx entries

                C_alternative2[idx,0] = idx
                C_alternative2[idx,1] = diab[diab_idx]
                C_alternative2[idx,2] = gly[gly_idx]

                # fill in change in bmi  for that combination
                C_alternative1[idx,i+3] = E[i,8]
                # fill in value for diabetes diagnosis in upcoming 4 years for that combination
                DD_alternative2[idx,:3] = C_alternative2[idx,:3]
                DD_alternative2[idx,i+3] = E[i,6]


            # contract C / DD
            nz = np.isfinite(C[:,7:]).sum(axis=1) # how many entries are there in the columns 7+?
            C1 = np.ones((combinations,nz.max() + 7)) * np.nan # containing delta values
            DD1 = np.ones((combinations,nz.max() + 7)) * np.nan # containing delta values


            for i in range(C1.shape[0]):
                # copy index and factors
                C1[i,:7] = C[i,:7]
                # how copy non-nan entries from the following columns
                nz = np.isfinite(C[i,7:])
                C1[i,7:(nz.sum()+7)] = C[i,7:][nz]
                # same for diabetes diagnosis values
                DD1[i,:7] = DD[i,:7]
                DD1[i,7:(nz.sum()+7)] = DD[i,7:][nz]


            # contract C_alternative1
            nz = np.isfinite(C_alternative1[:,4:]).sum(axis=1) # how many entries are there in the columns 4+?
            C1_alternative1 = np.ones((combinations_alternative1,nz.max() + 4)) * np.nan # containing delta values
            DD1_alternative1 = np.ones((combinations_alternative1,nz.max() + 4)) * np.nan


            for i in range(C1_alternative1.shape[0]):
                # copy index and factors
                C1_alternative1[i,:4] = C_alternative1[i,:4]
                # how copy non-nan entries from the following columns
                nz = np.isfinite(C_alternative1[i,4:])
                C1_alternative1[i,4:(nz.sum()+4)] = C_alternative1[i,4:][nz]
                # same for diabetes diagnosis values
                DD1_alternative1[i,:4] = DD_alternative1[i,:4]
                DD1_alternative1[i,4:(nz.sum()+4)] = DD_alternative1[i,4:][nz]


            # contract C_alternative2
            nz = np.isfinite(C_alternative2[:,3:]).sum(axis=1) # how many entries are there in the columns 3+?
            C1_alternative2 = np.ones((combinations_alternative2,nz.max() + 3)) * np.nan # containing delta values
            DD1_alternative2 = np.ones((combinations_alternative2,nz.max() + 3)) * np.nan


            for i in range(C1_alternative2.shape[0]):
                # copy index and factors
                C1_alternative2[i,:3] = C_alternative2[i,:3]
                # how copy non-nan entries from the following columns
                nz = np.isfinite(C_alternative2[i,3:])
                C1_alternative2[i,3:(nz.sum()+3)] = C_alternative2[i,3:][nz]
                # same for diabetes diagnosis values
                DD1_alternative2[i,:3] = DD_alternative2[i,:3]
                DD1_alternative2[i,3:(nz.sum()+3)] = DD_alternative2[i,3:][nz]


            self.C1_glyhb = C1.copy()
            self.C1_alternative1_glyhb = C1_alternative1.copy()
            self.C1_alternative2_glyhb = C1_alternative2.copy()

            self.DD1_glyhb = DD1.copy()
            self.DD1_alternative1_glyhb = DD1_alternative1.copy()
            self.DD1_alternative2_glyhb = DD1_alternative2.copy()




    def LoadELSA_bp(self):
        '''loads blood pressure trajectory data from ELSA into an numpy array'''


        self.ELSA_bp = np.genfromtxt('data/ELSA_bpmatch.csv', delimiter=',', names=True, dtype=None)
        self.ELSA_bp_numeric = np.genfromtxt('data/ELSA_bpmatch.csv', delimiter=',', skip_header=1)

        arraynames = ['CS1_bp', 'CD1_bp', 'CD1_alternative_bp', 'CS1_alternative_bp']
        # first check whether the compressed matching arrays are already present in the data folder
        # if so, load them, if not, compute them.

        loadnew = 0

        for a in range(len(arraynames)):
            filename = '%s/%s.npy' % (self.npydir,arraynames[a])
            if os.path.isfile(filename):
                # file present, load it
                exec('self.%s = np.load(filename)' % (arraynames[a]))

            else:
                loadnew = 1
                print 'npy file %s not found, recalculating...' % arraynames[a]

        if loadnew == 1:


            E = self.ELSA_bp_numeric
            age, sex, bmicat, bmibin, sys, dia = hadd.bp_vectors(E)
            n_sex, n_age, n_bmi, n_bmibin, n_sys, n_dia = hadd.bp_vectorsize(E)


            combinations = n_sex * n_age * n_bmi * n_sys * n_dia
            combinations_alternative = n_age * n_bmibin * n_sys # alternative matching process

            E_length = E.shape[0]

            CS = np.ones((combinations,E_length + 10)) * np.nan
            CD = np.ones((combinations,E_length + 10)) * np.nan

            CS_alternative = np.ones((combinations_alternative,E_length + 10)) * np.nan
            CD_alternative = np.ones((combinations_alternative,E_length + 10)) * np.nan

            # now loop over all entries

            for i in range(E_length):
                # evaluate combination index for row
                age_idx, sex_idx, bmi_idx, sys_idx, dia_idx = hadd.calculate_indices_bp(E,i,age,sex,bmicat,sys,dia)

                idx_array = hadd.calculate_Cindex_bp(n_sex, n_bmi, n_sys, n_dia, age_idx, sex_idx, bmi_idx, sys_idx, dia_idx)

                # fill in bands in the idx entries
                for eq_idx2, idx in enumerate(idx_array):
                    CS[idx,0] = idx
                    CS[idx,1] = age[age_idx]
                    CS[idx,2] = sex[sex_idx]
                    CS[idx,3] = bmicat[bmi_idx]
                    CS[idx,4] = sys[sys_idx]
                    CS[idx,5] = dia[dia_idx]
                    # fill in change in bp sys for that combination
                    CS[idx,i+10] = E[i,8]

                    # the same for diastolic
                    CD[idx,:6] = CS[idx,:6]
                    CD[idx,i+10] = E[i,9]

            # do the same for alternative matching
            for i in range(E_length):
                # evaluate combination index for row
                age_idx, bmi_idx, sys_idx = hadd.calculate_indices_alternative_bp(E,i,age,bmibin,sys)

                idx_array = hadd.calculate_Cindex_alternative_bp(n_bmibin, n_sys, age_idx, bmi_idx, sys_idx)

                # fill in bands in the idx entries
                for eq_idx2, idx in enumerate(idx_array):
                    CS_alternative[idx,0] = idx
                    CS_alternative[idx,1] = age[age_idx]
                    CS_alternative[idx,2] = bmibin[bmi_idx]
                    CS_alternative[idx,3] = sys[sys_idx]
                    # fill in change in bp sys for that combination
                    CS_alternative[idx,i+5] = E[i,8]

                    # the same for diastolic
                    CD_alternative[idx,:4] = CS_alternative[idx,:4]
                    CD_alternative[idx,i+5] = E[i,9]

            # now contract C in size to minimize nans
            nz = np.isfinite(CS[:,6:]).sum(axis=1) # how many entries are there in the columns 6+?

            CS1 = np.ones((combinations,nz.max() + 6)) * np.nan # containing delta sys values
            CD1 = np.ones((combinations,nz.max() + 6)) * np.nan # containing delta sys values

            for i in range(CS.shape[0]):
                # copy index and factors
                CS1[i,:6] = CS[i,:6]
                CD1[i,:6] = CD[i,:6]
                # how copy non-nan entries from the following columns
                nz = np.isfinite(CS[i,6:])
                CS1[i,6:(nz.sum()+6)] = CS[i,6:][nz]
                CD1[i,6:(nz.sum()+6)] = CD[i,6:][nz]



            # do the same for alternative match
            nz = np.isfinite(CS_alternative[:,4:]).sum(axis=1) # how many entries are there in the columns 6+?

            CS1_alternative = np.ones((combinations_alternative,nz.max() + 4)) * np.nan # containing delta sys values
            CD1_alternative = np.ones((combinations_alternative,nz.max() + 4)) * np.nan # containing delta sys values

            for i in range(CS_alternative.shape[0]):
                # copy index and factors
                CS1_alternative[i,:4] = CS_alternative[i,:4]
                CD1_alternative[i,:4] = CD_alternative[i,:4]
                # how copy non-nan entries from the following columns
                nz = np.isfinite(CS_alternative[i,4:])
                CS1_alternative[i,4:(nz.sum()+4)] = CS_alternative[i,4:][nz]
                CD1_alternative[i,4:(nz.sum()+4)] = CD_alternative[i,4:][nz]


            # save the contracted ELSA match templates as globals
            self.CS1_bp = CS1.copy()
            self.CD1_bp = CD1.copy()
            self.CS1_alternative_bp = CS1_alternative.copy()
            self.CD1_alternative_bp = CD1_alternative.copy()







    def LoadLifeTableData(self):
        '''loads life table model data for incidences and mortalities'''

        # load IHD life table data.
        # ----------------------
        # age:              col 0
        # male incidence:   col 1
        # male CFR:         col 2
        # female incidence: col 4
        # female CFR:       col 5
        self.LT_ihd = np.genfromtxt('data/lifetable_ihd.csv', delimiter=',', skip_header=3) # skip first 3 rows for ihd

        # load stroke life table data
        # ----------------------
        # age:              col 0
        # male incidence:   col 1
        # male CFR:         col 2
        # female incidence: col 4
        # female CFR:       col 5
        self.LT_stroke = np.genfromtxt('data/lifetable_stroke.csv', delimiter=',', skip_header=4) # skip first 4 rows for stroke

        # load dementia life table data
        # ----------------------
        # age:              col 0
        # male incidence:   col 1
        # male CFR:         col 2
        # female incidence: col 4
        # female CFR:       col 5
        self.LT_dementia = np.genfromtxt('data/lifetable_dementia.csv', delimiter=',', skip_header=4) # skip first 4 rows for stroke

        # load lung cancer life table data
        # -----------------------
        # age:                          col 0
        #
        # SMOKERS
        #
        # male smoker indicidence:      col 8
        # male smoker CFR:              col 9
        # female smoker incidence:      col 11
        # female smoker CFR:            col 12
        #
        # NONSMOKERS
        #
        # male nonsmoker incidence:      col 15
        # male nonsmoker CFR:            col 16
        # female nonsmoker incidence:    col 18
        # female nonsmoker CFR:          col 19

        self.LT_lungcancer = np.genfromtxt('data/lifetable_lungcancer.csv', delimiter=',', skip_header=4)


        # load ONS overall mortality data
        # ----------------------
        # age:              col 0
        # male qx:          col 2
        # female qx:        col 8
        self.LT_ONS = np.genfromtxt('data/ons_lifetable_2010-2012.csv', delimiter=',' , skip_header=7)


        # load table containing death from other causes, based on ONS data (death certificates)
        # ------------------------
        # age:              col 0
        # male mortality    col 1
        # female mortality  col 2
        self.LT_othercauses = np.genfromtxt('data/disease_mortalities_estimatingOtherCauses_2013_2013.csv', delimiter=',',skip_header=1)


        # smooth othercauses life table, which is discontinous in 5 year steps and goes only to age 90
        # this is done by a discretized diffusion equation.
        OT_male = self.LT_othercauses[:,1].copy()
        OT_female = self.LT_othercauses[:,2].copy()

        # as othercauses go only to age 95, take  ONS overall mortality for age 99 as highest other cause (age 100) (everyone is dying at that age anyway)
        OT_male[95:] = self.LT_ONS[100,2]
        OT_female[95:] = self.LT_ONS[100,8]


        # diffuse  how often?
        diffs = 10
        D = 0.2 # 'diffusion' coefficient
        OT_male_new = np.zeros((diffs+1,OT_male.size))
        OT_female_new = np.zeros((diffs+1,OT_female.size))

        OT_male_new[0] = OT_male
        OT_female_new[0] = OT_female

        # apply difference equation (= discretized diffusion equation) iteratively
        for d in range(diffs+1):
            if d>0:
                for i in range(OT_male.size):
                    if i>10 and i<(OT_male.size-1):
                        # for each point on male OT curve, define current point and point after (i+1) and before (i-1)

                        m_before = OT_male_new[d-1,i-1]
                        m_focal = OT_male_new[d-1,i]
                        m_after = OT_male_new[d-1,i+1]
                        # the same for female
                        f_before = OT_female_new[d-1,i-1]
                        f_focal = OT_female_new[d-1,i]
                        f_after = OT_female_new[d-1,i+1]

                        # apply difference operation with
                        OT_male_new[d,i] = np.exp(np.log(m_focal)+D*(np.log(m_before)+np.log(m_after)-2*np.log(m_focal)))
                        OT_female_new[d,i] = np.exp(np.log(f_focal)+D*(np.log(f_before)+np.log(f_after)-2*np.log(f_focal)))
                        #OT_female_new[d,i] = (OT_female_new[d-1,i-1] + OT_female_new[d-1,i+1])/ 2.0
                    else:
                        OT_male_new[d,i] = OT_male_new[d-1,i]
                        OT_female_new[d,i] = OT_female_new[d-1,i]

        # pick the last element in above arrays to be smoothed life table data
        self.LT_othercauses_smoothed = self.LT_othercauses.copy()
        # male other cause mortality
        self.LT_othercauses_smoothed[:,1] = OT_male_new[-1]
        # female other cause mortalitiy
        self.LT_othercauses_smoothed[:,2] = OT_female_new[-1]



    def CalculateDiseaseMortalityRate(self, male=1, disease='stroke', smoker=1):
        '''
        routine calculates probability of dying from a disease given an individual has the disease
        from CFR / prevalence / mortality data in the life tables
        This is here done for one disease, one gender only
        '''


        if male == 1:
            d_all = self.LT_ONS[:,2]
            d_other = self.LT_othercauses_smoothed[:,1]
        else:
            d_all = self.LT_ONS[:,8]
            d_other = self.LT_othercauses_smoothed[:,2]


        # load disease incidence and CFR data for given sex
        if disease == 'stroke':
            if male ==  1:
                i_disease = self.LT_stroke[:,1]
                cfr_disease = self.LT_stroke[:,2]
            else: # if female
                i_disease = self.LT_stroke[:,4]
                cfr_disease = self.LT_stroke[:,5]

        if disease == 'ihd':
            if male ==  1:
                i_disease = self.LT_ihd[:,1]
                cfr_disease = self.LT_ihd[:,2]
            else: # if female
                i_disease = self.LT_ihd[:,4]
                cfr_disease = self.LT_ihd[:,5]

        if disease == 'dementia':
            if male ==  1:
                i_disease = self.LT_dementia[:101,1]
                cfr_disease = self.LT_dementia[:101,2]
            else: # if female
                i_disease = self.LT_dementia[:101,4]
                cfr_disease = self.LT_dementia[:101,5]

        if disease == 'lungcancer':
            if male ==  1:
                if smoker == 1:
                    i_disease = self.LT_lungcancer[:,8]
                    cfr_disease = self.LT_lungcancer[:,9]
                else: # if nonsmoker
                    i_disease = self.LT_lungcancer[:,15]
                    cfr_disease = self.LT_lungcancer[:,16]
            else: # if female
                if smoker == 1:
                    i_disease = self.LT_lungcancer[:,11]
                    cfr_disease = self.LT_lungcancer[:,12]
                else: # if nonsmoker
                    i_disease = self.LT_lungcancer[:,18]
                    cfr_disease = self.LT_lungcancer[:,19]

        # define probability of death
        mort = 1 - np.exp(-cfr_disease-d_other) # probability that someone with disease dies of any cause
        mort2 = 1 - np.exp(-d_all)  # probability that someone in population dies of any cause

        # compute n_alive = number of individuals staying alive at ages
        n_alive = np.zeros(mort2.size)
        n_alive[0] = 1

        for i in range(1,mort2.size):
            n_alive[i] = n_alive[i-1] - (mort2[i-1] * n_alive[i-1])

        # now compute prevalence
        prev = np.zeros(i_disease.size)
        prev[0] = i_disease[0]

        for i in range(1,i_disease.size):
            prev[i] = (prev[i-1] + (1-prev[i-1]) * i_disease[i-1] - prev[i-1] * mort[i-1]) / (1-mort2[i-1])

        disease_mortality = cfr_disease * prev

        p_death_disease = 1 - np.exp(-disease_mortality)

        return prev, p_death_disease, cfr_disease





    def CalculateAllDiseaseMortalityRates(self):
        '''
        routine calculates probability of dying from a disease given an individual has the disease
        from CFR / prevalence / mortality data in the life tables
        This is done for all diseases
        '''

        self.prev_stroke_male, self.stroke_mort_rate_male, self.stroke_cfr_male = self.CalculateDiseaseMortalityRate(male=1, disease='stroke')
        self.prev_stroke_female, self.stroke_mort_rate_female, self.stroke_cfr_female = self.CalculateDiseaseMortalityRate(male=0, disease='stroke')

        self.prev_ihd_male, self.ihd_mort_rate_male, self.ihd_cfr_male = self.CalculateDiseaseMortalityRate(male=1, disease='ihd')
        self.prev_ihd_female, self.ihd_mort_rate_female, self.ihd_cfr_female = self.CalculateDiseaseMortalityRate(male=0, disease='ihd')

        self.prev_dementia_male, self.dementia_mort_rate_male, self.dementia_cfr_male = self.CalculateDiseaseMortalityRate(male=1, disease='dementia')
        self.prev_dementia_female, self.dementia_mort_rate_female, self.dementia_cfr_female = self.CalculateDiseaseMortalityRate(male=0, disease='dementia')

#        self.prev_lungcancer_male_s, self.lungcancer_mort_rate_male_s = self.CalculateDiseaseMortalityRate(male=1, disease='lungcancer', smoker=1)
#        self.prev_lungcancer_male_ns, self.lungcancer_mort_rate_male_ns = self.CalculateDiseaseMortalityRate(male=1, disease='lungcancer', smoker=0)
#        self.prev_lungcancer_female_s, self.lungcancer_mort_rate_female_s = self.CalculateDiseaseMortalityRate(male=0, disease='lungcancer', smoker=1)
#        self.prev_lungcancer_female_ns, self.lungcancer_mort_rate_female_ns = self.CalculateDiseaseMortalityRate(male=0, disease='lungcancer', smoker=0)






    def LoadPopulationParameters(self):
        '''loads population parameters from HSE file'''

        self.P = np.genfromtxt(
            'data/HSE_20092012.csv', delimiter=',', names=True, dtype=None)

        # Variable names that we need
        #======================================================================
        # diabetes      Diabetes register                                       0 (None), 1 (type1), 2 (type2)
        # bp1           Diagnosed with BP by doctor, not only while pregnant    can be 1, 2, -1 (yes/no/missing?)
        # compm7        Heart and circulatory system                            ?
        # male
        # age
        # ethnic16      1 "White - British"
        #               2 "White - Irish"
        #               3 "Any other white"
        #               4 "White and Black Caribbean"
        #               5 "White and Black African"
        #               6 "White and Asian"
        #               7 "Any other mixed"
        #               8 "Indian"
        #               9 "Pakistani"
        #               10 "Bangladeshi"
        #               11 "Chinese"
        #               12 "Any other Asian"
        #               13 "Caribbean"
        #               14 "African"
        #               15 "Any other Black"
        #               16 "Any other ethnic group"
        #
        # ---
        # QRisk relevant
        # ---
        # qimd          Quintile of IMD score                             1 (low deprivation) - 5 (high deprivation)
        # q_age
        # q_ethrisk     1 "White or not stated"
        #               2 "Indian"
        #               3 "Pakistani"
        #               4 "Bangladeshi"
        #               5 "Other Asian"
        #               6 "Caribbean"
        #               7 "African"
        #               8 "Chinese"
        #               9 "Other"
        # q_b_type1     Type 1 diabetes
        # q_b_type2     Type 2 diabetes
        # q_sbp         Systolic blood pressure
        # q_rati        Ratio of cholesterol
        # q_b_AF
        # q_b_ra
        # q_b_renal
        # q_b_treatedhyp Treated hypertension
        # q_fh_cvd      Fam history CVD
        # q_smoke_cat   0 "Never"
        #               1 "Ex-smoker"
        #               2 "Current <10"
        #               3 "Current 10-19"
        #               4 "Current 20+"
        # q_surv
        # q_town        Deprivation
        # q_bmi         BMI

        #======================================================================

        # check for completeness of data

        # age is complete
        # diabetes has two neg values
        ok_diab = self.P['diabetes'] >= 0
        # ditto bp
        ok_bp1 = self.P['bp1'] >= 0
        # ditto compm7
        ok_compm7 = self.P['compm7'] >= 0
        # all are male, so complete
        # q_ethrisk also complete
        # qimd also complete
        ok_q_b_type1 = self.P['q_b_type1'] >= 0
        ok_q_b_type2 = self.P['q_b_type2'] >= 0
        # sbp can be nan
        ok_q_sbp = np.isfinite(self.P['q_sbp'])
        ok_dia = np.isfinite(self.P['omdiaval'])
        # also q_rati can be nan
        ok_q_rati = np.isfinite(self.P['q_rati'])
        # q_b_AF: everyone has 0, the same for q_b_ra and q-b_renal and q_fh_cvd (how credible is that?)
        # q_b_treatedhyp is complete, also q_town
        # some q_bmi measures are also nan
        ok_q_bmi = np.isfinite(self.P['q_bmi'])
        # use only those that have complete SES data
        ok_qimd = self.P['qimd'] >= 0
        ok_gly = np.isfinite(self.P['glyhbval'])
        ok_educ = self.P['caide_educ'] >= 0
        # which data are fully complete?
        self.ok_all = ok_diab * ok_bp1 * ok_compm7 * ok_q_b_type1 * \
            ok_q_b_type2 * ok_q_sbp * ok_dia * ok_q_rati * ok_q_bmi * ok_qimd * ok_gly * ok_educ






    def InitialisePopulation(self):
        ''' initialises factors for population: '''

        # ---------------------------------------------------------------------
        # INITIALIZE ALL VECTORS AND ARRAYS

        # age of individuals
        self.age = np.zeros((self.population_size, self.simulation_time))
        self.id = np.zeros(self.population_size, dtype=int) # HSE individual's IDs
        self.gender = np.zeros(self.population_size, dtype=int)  # gender of individuals
        self.eth = np.zeros(self.population_size, dtype=int)  # ethnicity of individuals
        # Socioeconomic status of individuals = deprivation score from 1-5
        self.SES = np.zeros(self.population_size, dtype=int)
        self.alive = np.ones((self.population_size, self.simulation_time),dtype=bool)
        self.Death = np.zeros((self.population_size, self.simulation_time),dtype=bool)

        # health variables - these can change over time
        # individuals having diagnosed diabetes
        self.diabetes = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        # individuals having diagnosed CVD
        self.compm7 = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        # individuals having diagnosed hypertension
        self.bp = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.register_filter = np.zeros((self.population_size, self.simulation_time),dtype=bool)

        # DISEASE states
        self.Stroke = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.IHD = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.CVD = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.Dementia = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.LungCancer = np.zeros((self.population_size,self.simulation_time),dtype=bool)


        # cause of deaths
        self.CauseOfDeath = np.chararray(self.population_size, itemsize=12)
        self.CauseOfDeath[:] = 'None'

        # factors of individuals
        self.bmi = np.zeros((self.population_size, self.simulation_time))
        self.q_b_type1 = np.zeros((self.population_size, self.simulation_time))
        self.q_b_type2 = np.zeros((self.population_size, self.simulation_time))
        self.q_sbp = np.zeros((self.population_size, self.simulation_time)) # systolic bp
        self.dia = np.zeros((self.population_size, self.simulation_time)) # diastolic bp
        self.q_rati = np.zeros((self.population_size, self.simulation_time))
        self.chol = np.zeros((self.population_size, self.simulation_time))
        self.hdl = np.zeros((self.population_size, self.simulation_time))
        self.q_b_AF = np.zeros((self.population_size, self.simulation_time),dtype=int)
        self.q_b_ra = np.zeros((self.population_size, self.simulation_time),dtype=int)
        self.q_b_renal = np.zeros((self.population_size, self.simulation_time),dtype=int)
        self.q_b_treatedhyp = np.zeros((self.population_size, self.simulation_time),dtype=int)
        self.q_smoke_cat = np.zeros((self.population_size, self.simulation_time),dtype=int)
        self.glyhb = np.zeros((self.population_size, self.simulation_time)) # hba1c
        self.PA = np.zeros((self.population_size, self.simulation_time),dtype=int) # physical activity

        # include 'background' arrays for chol, hdl, q_rati, q_sbp, q_dia, bmi
        self.q_sbp_background = np.zeros((self.population_size, self.simulation_time)) # systolic bp
        self.dia_background = np.zeros((self.population_size, self.simulation_time)) # diastolic bp
        self.q_rati_background = np.zeros((self.population_size, self.simulation_time))
        self.chol_background = np.zeros((self.population_size, self.simulation_time))
        self.bmi_background = np.zeros((self.population_size, self.simulation_time))



        #  health variables - these stay the same
        self.q_fh_cvd = np.zeros(self.population_size)
        self.q_town = np.zeros(self.population_size)
        self.qimd = np.zeros(self.population_size)
        self.educ = np.zeros(self.population_size)




        # Health Checks variables
        # Qrisk of individuals
        self.QRisk = np.zeros((self.population_size, self.simulation_time))
        # treatments
        self.WeightReduction = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.Statins = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.SmokingCessation = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.Hypertensives = np.zeros((self.population_size, self.simulation_time),dtype=bool)

        # capture who is on treatment in a year
        self.on_statins = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.on_aht = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.on_wr = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        # time when individual last started weight management (used to implement weight regain after referral).  Set to -1 if never had weight management
        self.prev_wr_time = np.zeros((self.population_size), dtype=int) - 1


        self.WeightReduction_Offered = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.Statins_Offered = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.SmokingCessation_Offered = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.Hypertensives_Offered = np.zeros((self.population_size, self.simulation_time),dtype=bool)


        self.eligible = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.poff = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        # People offered HC, or attended anyway despite not being eligible
        self.OfferedHC = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        # time at which each individual was last offered a HC.
        # -1 if never offered
        self.prev_offer_time = np.zeros((self.population_size), dtype=int) - 1
        # Relative uptake rates
        self.RelTurnup = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        # People attending Health Checks
        self.Attending = np.zeros((self.population_size, self.simulation_time),dtype=bool)
        self.OfferedTreatment = np.zeros((self.population_size, self.simulation_time),dtype=bool)

        # Dementia scores
        self.CAIDE = np.zeros((self.population_size, self.simulation_time),dtype=int)
        self.DementiaRisk = np.zeros((self.population_size, self.simulation_time))

        # additional variables
        # within age of Health checks
        self.within_age = np.zeros((self.population_size, self.simulation_time), dtype=bool)
        self.on_register = np.zeros((self.population_size, self.simulation_time), dtype=bool)# on any register
        self.simulation_timepoints = np.ones((self.population_size, self.simulation_time), dtype=bool)
        self.current_turnup = np.zeros((self.population_size, self.simulation_time), dtype=bool)

        # Events matrices
        self.StrokeEvents = np.zeros((self.population_size, self.simulation_time), dtype=bool)
        self.IHDEvents = np.zeros((self.population_size, self.simulation_time), dtype=bool)
        self.LungCancerEvents = np.zeros((self.population_size, self.simulation_time), dtype=bool)
        self.DementiaEvents = np.zeros((self.population_size, self.simulation_time), dtype=bool)
        self.OtherEvents = np.zeros((self.population_size, self.simulation_time), dtype=bool)
        # saving scaled CVD risks
        self.CVD_risks = np.zeros((self.population_size, self.simulation_time))
        self.CVD_events = np.zeros((self.population_size, self.simulation_time), dtype=bool)

        # scaling factor for AHT treatment:
        self.AHT_prescription_scaling = np.zeros(self.simulation_time)

        # --------------------------------------------------------------------------
        # POPULATE VECTORS AND ARRAYS FROM HSE POPULATION DATA BY RANDOM SAMPLING
        # how many individuals in the dataset have complete data?

        num_ok = self.ok_all.sum()
        # sample a random selection with replacement from those
        # select = np.random.randint(0, num_ok, self.population_size)


        # alternatively, select due to population distribution

        # rows:
        #           (0) age 35 - 40: )
        #           (1) age 40 - 45:
        #           (2) age 45 - 50:
        #           (3) age 50 - 55:
        #           (4) age 55 - 60:

        #           (7) age 70 - 75:

        # columns:  (0) male_White, (1) male_Asian, (2) male_Black, (3) male_Mixed/other,
        #           (4) female_White, (5) female_Asian, (6) female_Black, (7) female_Mixed/other

        # establish how many individuals of each census combination are needed
        nr_ind = np.round(self.parms['census_prop'] * self.population_size)
        nr_ind = np.array(nr_ind,dtype = int)
        # check how much total nr_ind is off from population_size
        discrepancy = int(nr_ind.sum() - self.population_size)
        direction = np.sign(discrepancy)


        np.random.seed(0)

        # add/subtract this number from random entries, such that nr_ind sums up to population size
        for i in range(np.abs(discrepancy)):
            x = np.random.randint(0,7)
            y = np.random.randint(0,7)
            nr_ind[x,y] -= direction

        self.nr_ind = nr_ind.copy()
        self.nr_ind_original = nr_ind.copy()

        # give an unique index to every entry in nr_ind

        rows = 9
        cols  = 8
        # 7 rows and 7 cols in nr_ind
        idx = np.zeros((rows,cols),dtype=int)
        for i in range(rows):
            for j in range(cols):
                idx[i,j] = 1 + i*cols + j



        # now make an index vector for all individuals in HSE population that are ok to select from, assigning
        # the above idx to HSE pop

        self.HSE_idx = np.zeros(num_ok,dtype=int)

        self.Page = self.P['age'][self.ok_all]
        self.Psex = self.P['male'][self.ok_all]
        self.Peth = self.P['q_ethrisk'][self.ok_all]

        for i in range(num_ok):

            try:

                age_pick = self.Page[i]
                sex_pick = self.Psex[i]
                eth_pick = self.Peth[i]
            except:
                pass


            # establish if picked individual is in the right age range (30-75)
            agerange = False
            if (age_pick >= 30) * (age_pick <= 75) == 1:
                agerange = True


            # which age band is this individual in
            if agerange == True:
                bins = np.arange(30,80,5)
                age = np.array([age_pick])
                age_idx = np.digitize(age,bins)[0] - 1

                if sex_pick == 1:
                    #  male
                    if eth_pick == 1:
                        # white
                        sex_eth_idx = 0
                        # Asian
                    if (eth_pick > 1) * (eth_pick <= 5) + (eth_pick == 8) == 1:
                        sex_eth_idx = 1
                        # Black
                    if (eth_pick >= 6) * (eth_pick <= 7) == 1:
                        sex_eth_idx = 2

                    if eth_pick == 9:
                        sex_eth_idx = 3
                else:
                    #female
                    if eth_pick == 1:
                        # white
                        sex_eth_idx = 4
                        # Asian
                    if (eth_pick > 1) * (eth_pick <= 5) + (eth_pick == 8) == 1:
                        sex_eth_idx = 5
                        # Black
                    if (eth_pick >= 6) * (eth_pick <= 7) == 1:
                        sex_eth_idx = 6

                    if eth_pick == 9:
                        sex_eth_idx = 7

                # now we have for each of the selected HSE individuals with complete data (ok_all) information on which index it is in nr_ind
                self.HSE_idx[i] = 1 + age_idx*cols + sex_eth_idx


        # now go through all the category combination bins:
        #  how many individuals of each bin are needed?
        # pick out the HSE individuals of this bin, and sample the desired number of individuals into the population

        select_list = []

        # reset random number generator
        rnd.seed(0)

        # smoker status
        s = self.P['q_smoke_cat'][self.ok_all]



        for i in range(1,(rows*cols+1)):
            # establish which individuals of HSE population fall into this index/census category:

            j = np.where(self.HSE_idx == i)[0]


            # how many of this category do we need in the simulated population?
            row_nr = np.where(idx==i)[0][0]
            col_nr = np.where(idx==i)[1][0]

            num_j = self.nr_ind[row_nr,col_nr]

            # reset random number generator of numpy


            # now sample from those, add to select_list
            for jj in range(num_j):
                #print jj
                try:
                    # check where are smokers
                    sm = s[j] > 1


                    # smokers should be chosen with chance 0.6, nonsmokers with chance 0.4
                    # this comes from currently 14% smokers in ok_all HSE population, but 20% needed.
                    j_prop = np.ones(j.size)*0.35
                    j_prop[sm] = 0.65

                    # select randomly from new pool, weighted
                    r = np.random.random(j.size)
                    weighted = r<j_prop

                    j_weighted = j[weighted]

                    j_val = rnd.choice(j_weighted)


                except:
                    # pick one random element of whole num_ok if this index has no corresponding in HSE dataset
                    j_val = np.random.randint(0,num_ok)

                select_list.append(j_val)

        # in the end, transform select list into np array
        select = np.array(select_list)
        self.select = select.copy()
        # reset seed for numpy random
        np.random.seed()


        # fixed factors
        self.eth = self.P['q_ethrisk'][self.ok_all][select]
        self.SES = self.P['qimd'][self.ok_all][select]
        self.q_fh_cvd = self.P['q_fh_cvd'][self.ok_all][select]
        self.q_town = self.P['q_town'][self.ok_all][select]
        self.qimd = self.P['qimd'][self.ok_all][select]
        self.educ = self.P['caide_educ'][self.ok_all][select]

        self.age[:, 0] = self.P['age'][self.ok_all][select]
        self.gender = self.P['male'][self.ok_all][select]

        self.id = self.P['pserial'][self.ok_all][select]
        # copy in the column for age:
        for i in range(1, self.simulation_time):
            self.age[:, i] = self.age[:, i - 1] + 1

        # cap age at 99
        self.age[self.age>99]=99


        # copy in the columns of registers
        for i in range(self.simulation_time):
            self.diabetes[:, i] = self.P['diabetes'][self.ok_all][select]

            bp_vec = self.P['bp1'][self.ok_all][select]
            # set value 2 for "no hypertension" to zero as this makes more
            # sense
            bp_vec[bp_vec == 2] = 0
            self.bp[:, i] = bp_vec.copy()

            self.compm7[:, i] = self.P['compm7'][self.ok_all][select]
            self.q_b_type1[:, i] = self.P['q_b_type1'][self.ok_all][select]
            self.q_b_type2[:, i] = self.P['q_b_type2'][self.ok_all][select]


        # copy in first column for other factors
        self.bmi[:, 0] = self.P['q_bmi'][self.ok_all][select]

        self.q_sbp[:, 0] = self.P['q_sbp'][self.ok_all][select]
        self.q_rati[:, 0] = self.P['q_rati'][self.ok_all][select]
        self.q_b_AF[:, 0] = self.P['q_b_AF'][self.ok_all][select]
        self.q_b_ra[:, 0] = self.P['q_b_ra'][self.ok_all][select]
        self.q_b_renal[:, 0] = self.P['q_b_renal'][self.ok_all][select]
        self.q_b_treatedhyp[:, 0] = self.P['q_b_treatedhyp'][self.ok_all][select]
        self.q_smoke_cat[:, 0] = self.P['q_smoke_cat'][self.ok_all][select]
        self.dia[:, 0] = self.P['omdiaval'][self.ok_all][select]
        self.chol[:,0] = self.P['cholval1'][self.ok_all][select]
        self.hdl[:,0] = self.P['hdlval1'][self.ok_all][select]
        self.glyhb[:,0] = self.P['glyhbval'][self.ok_all][select]

        # also fill the first entries for 'background' risk factor values
        self.dia_background[:, 0] = self.P['omdiaval'][self.ok_all][select]
        self.chol_background[:,0] = self.P['cholval1'][self.ok_all][select]
        self.q_rati_background[:, 0] = self.P['q_rati'][self.ok_all][select]
        self.q_sbp_background[:, 0] = self.P['q_sbp'][self.ok_all][select]
        self.bmi_background[:, 0] = self.P['q_bmi'][self.ok_all][select]

        # physical activity: assume that percentage of the population are active (determined in self.parms['physically_active'])
        # this remains stable over lifetime
        np.random.seed(self.randseed)
        r = np.random.random(self.population_size)

        self.PA[(r < self.up_physically_active),:] = 1


    def Info(self, individual, timepoint=0):
        '''collects all available info for an individual at a timepoint (default timepoint=0)
        outputs this info as a dict'''
        eth =["White or not stated",
              "Indian", "Pakistani", "Bangladeshi", "Other Asian", "Caribbean",
              "African", "Chinese", "Other"]
        sex = ['female','male']
        smoke = ["Never", "Ex-smoker", "Current <10", "Current 10-19", "Current 20+"]

        i = individual
        t = timepoint
        Infodict={
        'age': self.age[i,t],
        'sex': sex[self.gender[i]],
        'eth': eth[(self.eth[i]-1)],
        'qimd': self.qimd[i],
        'BMI': self.bmi[i,t],
        'smoking': smoke[int(self.q_smoke_cat[i,t])],
        'sys BP': self.q_sbp[i,t],
        'dia BP': self.dia[i,t],
        'glyhb': self.glyhb[i,t],
        'QRisk': self.QRisk[i,t],
        'CAIDE': self.CAIDE[i,t],
        'chol': self.chol[i,t],
        'Diabetes Typ 2': self.q_b_type2[i,t]
        }

        return Infodict



    def EstablishAgeWeights(self):
        '''evaluates what proportion of males in the population are of age x at time t=0
        , compared to whole simulation:
            age_weight = number of people aged x at start of population / number of people aged x throughout simulation
        '''

        male = self.gender == 1
        female = self.gender == 0

        self.age_weight_male = np.zeros(101)
        self.age_weight_female = np.zeros(101)

        minstartage = self.age[:,0].min()
        maxstartage = self.age[:,0].max()

        minallage = self.age.min()
        maxallage = self.age.max()

        init_agespan = int(maxstartage - minstartage)

        male_start = np.histogram(self.age[male,0].flatten(),bins=(maxstartage-minstartage))
        male_all = np.histogram(self.age[male,:].flatten(),bins=(maxallage-minallage))

        # restrict bin of all ages to initial ages, as weights are only applied to initial ages anyway
        male_all_restricted = np.array(male_all[0][:init_agespan],dtype=float)
        self.age_weight_male[np.array(male_start[1][:-1],dtype=int)] = male_start[0] / male_all_restricted
        # at age 73, and 74, numbers are lumped: divide by 2
        self.age_weight_male[[73,74]] = self.age_weight_male[73]/2.0

        female_start = np.histogram(self.age[female,0].flatten(),bins=(maxstartage-minstartage))
        female_all = np.histogram(self.age[female,:].flatten(),bins=(maxallage-minallage))

        # restrict bin of all ages to initial ages, as weights are only applied to initial ages anyway
        female_all_restricted = np.array(female_all[0][:init_agespan],dtype=float)
        self.age_weight_female[np.array(female_start[1][:-1],dtype=int)] = female_start[0] / female_all_restricted
        # at age 73, and 74, numbers are lumped: divide by 2
        self.age_weight_female[[73,74]] = self.age_weight_female[73]/2.0



    def InitialiseCVDEvents(self):
        '''Based on Stroke and IHD prevalences and age weight, establish which proportion of the initial population will already have a CVD event before the simulation starts
        assign these to the individuals with highest QRisk for these age groups'''

        np.random.seed(self.randseed)

        self.stroke_init = np.zeros(self.population_size,dtype=bool)
        self.IHD_init = np.zeros(self.population_size,dtype=bool)
        # get index vector for population for later 'tracing back' highest Qrisks to individuals
        k = np.arange(self.population_size)

        for i in range(int(self.age[:,0].min()), int(self.age[:,0].max() + 1)):
            a_male = (self.age[:,0] == i) * (self.gender == 1)
            a_female = (self.age[:,0] == i) * (self.gender == 0)

            # STROKE

            # establish how many of the current age group will have disease at the beginning of the simulatin
            arisk_stroke_male = self.age_weight_male[i] * self.prev_stroke_male[i]
            arisk_stroke_female = self.age_weight_male[i] * self.prev_stroke_female[i]

            r_male = np.random.random(a_male.sum())
            male_stroke_events = (r_male < arisk_stroke_male).sum()

            r_female = np.random.random(a_female.sum())
            female_stroke_events = (r_female < arisk_stroke_female).sum()
            # if n people at this age group are selected to have had stroke before begin of simulation,
            # search for n*3 highest QRisks and assign stroke event to these
            # n*3 is chosen to relax the assumption that highest Qrisk will always cause events
            # and to account for possibility that high QRISK accounts for IHD _and_ stroke, so not just one of those

            if male_stroke_events > 0 :
                # sort by QRisk
                ind = np.argsort(self.QRisk[a_male,0])
                # of those, select the 3n highest
                ind_highest = ind[-3*male_stroke_events:]
                # of those, pick out n
                ind_select = rnd.sample(ind_highest,male_stroke_events)

                # get indices for these individuals in population vector
                pop_ind = k[a_male][ind_select]

                # assign those the illness
                self.stroke_init[pop_ind] = np.ones(pop_ind.size, dtype=bool)


            # same for females
            if female_stroke_events > 0 :
                # sort by QRisk
                ind = np.argsort(self.QRisk[a_female,0])
                # of those, select the 3n highest
                ind_highest = ind[-3*female_stroke_events:]
                # of those, pick out n
                ind_select = rnd.sample(ind_highest,female_stroke_events)
                pop_ind = k[a_female][ind_select]
                # assign those the illness
                self.stroke_init[pop_ind] = np.ones(pop_ind.size, dtype=bool)


            # IHD

            # establish how many of the current age group will have disease at the beginning of the simulatin
            arisk_IHD_male = self.age_weight_male[i] * self.prev_ihd_male[i]
            arisk_IHD_female = self.age_weight_male[i] * self.prev_ihd_female[i]

            r_male = np.random.random(a_male.sum())
            male_IHD_events = (r_male < arisk_IHD_male).sum()

            r_female = np.random.random(a_female.sum())
            female_IHD_events = (r_female < arisk_IHD_female).sum()

            # same procedure as above for stroke

            if male_IHD_events > 0 :
                # sort by QRisk
                ind = np.argsort(self.QRisk[a_male,0])
                # of those, select the 3n highest
                ind_highest = ind[-3*male_IHD_events:]
                # of those, pick out n
                ind_select = rnd.sample(ind_highest,male_IHD_events)
                # assign those the illness
                # get indices for these individuals in population vector
                pop_ind = k[a_male][ind_select]

                # assign those the illness
                self.IHD_init[pop_ind] = np.ones(pop_ind.size, dtype=bool)

            # same for females
            if female_IHD_events > 0 :
                # sort by QRisk
                ind = np.argsort(self.QRisk[a_female,0])
                # of those, select the 3n highest
                ind_highest = ind[-3*female_IHD_events:]
                # of those, pick out n
                ind_select = rnd.sample(ind_highest,female_IHD_events)
                # assign those the illness
                # get indices for these individuals in population vector
                pop_ind = k[a_female][ind_select]

                # assign those the illness
                self.IHD_init[pop_ind] = np.ones(pop_ind.size, dtype=bool)



    def SplitQRiskInAnnualEvents(self):
        '''break down QRisk into probabilities per year per event (Stroke/IHD) from life table data.
        Stores output in two arrays each for males and females:
            (1) prob_year: relative risk to develop an event over the next 10 consecutive years, given age.
                rows correspond to age, columns to 10 consecutive years
            (2) prob_event: relative risk that an event is a stroke, given that it happens and that given age.
        '''

        self.CVDprob_year_male = np.zeros((91,10))
        self.CVDprob_event_male = np.zeros(101)
        self.CVDprob_year_female = np.zeros((91,10))
        self.CVDprob_event_female = np.zeros(101)

        cvd_total_male = self.LT_ihd[:,1]+self.LT_stroke[:,1]
        cvd_total_female = self.LT_ihd[:,4]+self.LT_stroke[:,4]
        # evaluate prob_year
        for i in range(91):
            # sum up all events for the upcoming 10 years
            male_p = cvd_total_male[i:i+10]
            female_p = cvd_total_female[i:i+10]

            # convert to "q-scale", ie take into account conditional probabilities for each year of not having gotten the disease in the year before,
            # considering the 10 year interval we consider
            male_q = hadd.ConvertToQ(male_p)
            female_q = hadd.ConvertToQ(female_p)
            male_total_10yrs = male_q.sum()
            female_total_10yrs = female_q.sum()
            # now evaluate the relative risks for all cvds over this 10 year period
            with np.errstate(invalid='ignore'): # this gets rid of error messages due to dividing through zero

                self.CVDprob_year_male[i,:] = male_q / male_total_10yrs
                self.CVDprob_year_female[i,:] = female_q / female_total_10yrs

        # evaluate prob_event
        for i in range(101):
            with np.errstate(invalid='ignore'): # this gets rid of error messages due to dividing through zero
                self.CVDprob_event_male[i] = self.LT_stroke[i,1] / (self.LT_ihd[i,1]+self.LT_stroke[i,1])
                self.CVDprob_event_female[i] = self.LT_stroke[i,4] / (self.LT_ihd[i,4]+self.LT_stroke[i,4])





    def SplitCAIDEInAnnualEvents(self):
        '''break down CAIDE into probabilities per year from life table data.
        Stores output in one array each for males and females:
            (1) prob_year: relative risk to develop dementia over the next 20 consecutive years, given age.
                rows correspond to age, columns to 20 consecutive years

        '''

        self.DEMprob_year_male = np.zeros((81,20))
        self.DEMprob_year_female = np.zeros((81,20))


        # evaluate prob_year for Dementia
        for i in range(81):
            # evaluate the relative risks for all cvds over the next 20 year period
            with np.errstate(invalid='ignore'): # this gets rid of error messages due to dividing through zero
                male_q = hadd.ConvertToQ(self.LT_dementia[i:i+20,1])
                female_q = hadd.ConvertToQ(self.LT_dementia[i:i+20,4])
                self.DEMprob_year_male[i,:] = male_q / male_q.sum()
                self.DEMprob_year_female[i,:] = female_q / female_q.sum()


    def glyhbmatch_core_process(self, l, t, j_min, loopsize, out_q):
        ''' search process within GlyHb match procedure
        for loop over population to be parallelised'''

        l.acquire()
        #rnd.seed(self.randseed)

        glyhb_delta4 = np.zeros(loopsize)
        diabetes_delta4  = np.zeros(loopsize)


        #rnd.seed(self.randseed+t+j_min)

        for j_idx, j in enumerate(range(j_min,j_min + loopsize)):
            rnd.seed(self.randseed+j_idx+t)

            try:
                # enter ELSA glyhb matrix for combination index
                D = self.C1_glyhb[self.glyhb_IND[j],7:]
                D_finite = np.isfinite(D)
                d = D[D_finite]
                # same for diabetes diagnosis
                DD = self.DD1_glyhb[self.glyhb_IND[j],7:]
                DD_finite = np.isfinite(DD)
                dd = DD[DD_finite]


                # choose one value
                glyhb_delta4[j_idx] = rnd.choice(d)
                diabetes_delta4[j_idx] = rnd.choice(dd)


            except:
                # rematch 1
                try:
                    D = self.C1_alternative1_glyhb[self.glyhb_IND_alternative1[j],4:]
                    D_finite = np.isfinite(D)
                    d = D[D_finite]

                    # same for diabetes diagnosis
                    DD = self.DD1_alternative1_glyhb[self.glyhb_IND_alternative1[j],4:]
                    DD_finite = np.isfinite(DD)
                    dd = DD[DD_finite]

                    glyhb_delta4[j_idx] = rnd.choice(d)
                    diabetes_delta4[j_idx] = rnd.choice(dd)

                except:

                    # rematch 2
                    try:

                        D = self.C1_alternative2_glyhb[self.glyhb_IND_alternative2[j],3:]
                        D_finite = np.isfinite(D)
                        d = D[D_finite]

                        # same for diabetes diagnosis
                        DD = self.DD1_alternative2_glyhb[self.glyhb_IND_alternative2[j],3:]
                        DD_finite = np.isfinite(DD)
                        dd = DD[DD_finite]



                    except:

                        # assign no change
                        d = np.zeros(1)
                        if self.verbose == True:
                            print('no change assigned to glyhbmatch for individual %d' % j)

                        # no change for diabetes diagnosis means taking last value of diabetes diagnosis for individual
                        try:
                            dd = np.ones(1)  *  self.diabetes[j,t]
                        except:
                            # assume no diabetes
                            dd = np.zeros(1)



                    try:
                        glyhb_delta4[j_idx] = rnd.choice(d)
                        diabetes_delta4[j_idx] = rnd.choice(dd)
                    except:
                        #print d4,d8
                        pass


        out_q.put([j_min, glyhb_delta4, diabetes_delta4])
        l.release()



    def MatchELSA_glyhb_new(self,t=0):
        '''optimized routine to match GLYHB data from ELSA to population
        to get 4 and 1 year trajectories for each individual'''

        E = self.ELSA_glyhb_numeric

        a = self.age[:,t]
        # bin age : below 55, 55-70, above 70
        bins = np.array([0,55,70,100])
        dage = np.digitize(a,bins)

        # bin BMI
        b = self.bmi_background[:,t]
        bins = np.array([0,25,30,100])
        dbmi = np.digitize(b,bins)

        # bin deprivation
        q = self.qimd.copy()
        q[q==1] = 0
        q[q>1] = 1

        # bin glyhb
        bins = np.array([0, 5.25, 5.5, 5.75, 6, 7, 10])
        dgly = np.digitize(self.glyhb[:,t],bins)

        # get status of diabetes diagnosis / medication:
        diab_idx = self.diabetes[:,t]

        age_idx = dage-1
        bmi_idx = dbmi-1
        gly_idx = dgly-1
        sex_idx = self.gender
        eq_idx = q

        n_sex, n_age, n_bmi, n_eq, n_diab, n_gly = hadd.glyhb_vectorsize(E)

        # evaluate default index
        IND = age_idx * (n_sex * n_bmi * n_eq * n_diab * n_gly) + \
                sex_idx * (n_bmi * n_eq * n_diab * n_gly) + \
                bmi_idx * (n_eq * n_diab * n_gly) + \
                eq_idx * (n_diab * n_gly) + \
                diab_idx * (n_gly) + gly_idx

        self.glyhb_IND = np.array(IND,dtype=int)


        # evaluate alternative index 1 (second match)
        IND_alternative1 = bmi_idx * (n_diab * n_gly)  + \
                diab_idx * (n_gly) + gly_idx

        self.glyhb_IND_alternative1 = np.array(IND_alternative1,dtype=int)


        # evaluate alternative index 2 (second match)
        IND_alternative2 = diab_idx * (n_gly) + gly_idx

        self.glyhb_IND_alternative2 = np.array(IND_alternative2,dtype=int)


        # invoke change vectors for glyhb and diabetes diagnosis
        self.glyhb_delta4 =np.zeros(self.population_size) * np.nan
        self.diabetes_delta4 =np.zeros(self.population_size) * np.nan

        loopsize = self.population_size/self.nprocs
        chunksize = self.nprocs * loopsize


        t0 = time.time()
        for i in xrange(0, self.population_size, chunksize):

            # divide parallelised job into chunks loopsize or rest of array size
            proc_size = self.nprocs
            jstart_array = xrange(i,i+chunksize,loopsize)
            if jstart_array[-1] > self.population_size:
                jstart_array = xrange(i,self.population_size, loopsize)
                proc_size = len(jstart_array)

            out_q = multiprocessing.Queue()
            procs = []

            # send jobs to multiprocessing module
            for jstart in jstart_array:
                pr = multiprocessing.Process(
                            target = self.glyhbmatch_core_process,
                            args = (multiprocessing.Lock(), t, jstart, loopsize, out_q))
                procs.append(pr)
                pr.start()
                #pr.join()
                #pr.terminate()

            # retrieve results
            for m in xrange(proc_size):
                res = out_q.get()
                j_min = res[0]
                #print i, j_min
                self.glyhb_delta4[j_min:j_min+loopsize] = res[1]
                self.diabetes_delta4[j_min:j_min+loopsize] = res[2]

                self.glyhb_entry_order.append(res[0])

        t1 = time.time()
        if self.verbose == True:
            print('match time elapsed: %g sec' % (t1-t0))

        # about 0.05 percent of the population doesn't apparently match here
        # and thus gets a nan value in glyhb_delta.
        # this needs to be checked, but until then, it's safe to assume that these don't change
        delta_nans = np.isnan(self.glyhb_delta4)
        self.glyhb_delta4[delta_nans] = 0

        delta_nans = np.isnan(self.diabetes_delta4)
        self.diabetes_delta4[delta_nans] = 0

        # break down 4 year changes in change per year
        self.glyhb_delta1 = self.glyhb_delta4/4.0

        # delete internal arrays
        del delta_nans,a,b,q,dgly,dbmi,dage,age_idx,bmi_idx,gly_idx,sex_idx,eq_idx




    def cholmatch_core_process(self, l, t, j_min, loopsize, out_q):
        '''search process within cholesterin match procedure for loop
        over population to be parallelised'''


        l.acquire()
        #rnd.seed(self.randseed)
        chol_delta4 = np.zeros(loopsize)
        hdl_delta4 = np.zeros(loopsize)

        #rnd.seed(self.randseed+t*2+j_min)

        for j_idx, j in enumerate(range(j_min,j_min + loopsize)):
            rnd.seed(self.randseed+j_idx+t*2)

            try:
                D = self.C1_chol[self.chol_IND[j],6:]
                D_finite = np.isfinite(D)
                d = D[D_finite]

                # same for HDL
                DH = self.CH1_chol[self.chol_IND[j],6:]
                DH_finite = np.isfinite(DH)
                dh = DH[DH_finite]
            except:
                # no match, assign no change
                d = np.zeros(1)
                dh = np.zeros(1)
                if self.verbose == True:
                    print('no change assigned to cholmatch for individual %d' % j)

            try:
                chol_delta4[j_idx] = rnd.choice(d)
                hdl_delta4[j_idx] = rnd.choice(dh)
            except:
                try:
                    # rematch 1
                    D = self.C1_alternative1_chol[self.chol_IND_alternative1[j],4:]
                    D_finite = np.isfinite(D)
                    d = D[D_finite]
                    # same for HDL
                    DH = self.CH1_alternative1_chol[self.chol_IND_alternative1[j],4:]
                    DH_finite = np.isfinite(DH)
                    dh = DH[DH_finite]
                    try:
                        chol_delta4[j_idx] = rnd.choice(d)
                        hdl_delta4[j_idx] = rnd.choice(dh)

                    except:
                        #print d4,d8
                        pass

                except:
                    # rematch 2
                    try:
                        D = self.C1_alternative2_chol[self.chol_IND_alternative2[j],3:]
                        D_finite = np.isfinite(D)
                        d = D[D_finite]
                        # same for HDL
                        DH = self.CH1_alternative2_chol[self.chol_IND_alternative2[j],3:]
                        DH_finite = np.isfinite(DH)
                        dh = DH[DH_finite]
                    except:
                        # something went wrong in rematch 2
                        # assign zero change to that individual
                        d = np.array([0])
                        dh = np.array([0])
                        #print 'no CHOL match at time t=%d for ind %d, assigning no change' % (t,i)


                    try:
                        chol_delta4[j_idx] = rnd.choice(d)
                        hdl_delta4[j_idx] = rnd.choice(dh)
                    except:
                        #print d4,d8
                        pass


        out_q.put([j_min, chol_delta4, hdl_delta4])
        l.release()


    def MatchELSA_chol_new(self,t=0):
        '''optimized routine to match Cholesterol data from ELSA to population
        to get 4 and 1 year trajectories for each individual'''

        E = self.ELSA_chol_numeric

        a = self.age[:,t]
        # bin age : below 55, 55-70, above 70
        bins = np.array([0,55,70,100])
        dage = np.digitize(a,bins)

        # bmi category:
        b = self.bmi_background[:,t]
        bins = np.array([0,25,30,100])
        dbmi = np.digitize(b,bins)

        # categorize deprivation index
        # 2: “Fifth 1 or 2: most deprived” / 3 “Fifth 3” / 4: “Fifth 4 or 5: most advantaged”

        q = self.qimd.copy()
        q[q<3] = 2
        q[q>3] = 4

        # categorise cholesterol levels
        bins = np.array([ 0, 2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 20])
        dchol = np.digitize(self.chol_background[:,t],bins)

        age_idx = dage-1
        bmi_idx = dbmi-1
        chol_idx = dchol-2
        sex_idx = self.gender
        eq_idx = q-2

        n_sex, n_age, n_bmi, n_eq, n_chol = hadd.chol_vectorsize(E)

        IND =   age_idx * (n_sex * n_bmi * n_eq * n_chol) + \
                        sex_idx * (n_bmi * n_eq * n_chol) + \
                        bmi_idx * (n_eq * n_chol) + \
                        eq_idx * (n_chol) + chol_idx

        self.chol_IND = np.array(IND,dtype=int)

        self.chol_IND_alternative1 =  sex_idx * (n_bmi * n_chol) + \
                    bmi_idx * (n_chol) + chol_idx

        self.chol_IND_alternative2 = sex_idx * (n_chol) + chol_idx

        self.chol_delta4 = np.zeros(self.population_size) * np.nan
        self.hdl_delta4 = np.zeros(self.population_size) * np.nan

        loopsize = self.population_size/self.nprocs
        chunksize = self.nprocs * loopsize


        t0 = time.time()
        for i in xrange(0, self.population_size, chunksize):

            # divide parallelised job into chunks loopsize or rest of array size
            proc_size = self.nprocs
            jstart_array = xrange(i,i+chunksize,loopsize)
            if jstart_array[-1] > self.population_size:
                jstart_array = xrange(i,self.population_size, loopsize)
                proc_size = len(jstart_array)

            out_q = multiprocessing.Queue()
            procs = []


            # send jobs to multiprocessing module
            for jstart in jstart_array:
                pr = multiprocessing.Process(
                            target = self.cholmatch_core_process,
                            args = (multiprocessing.Lock(), t, jstart, loopsize, out_q))
                procs.append(pr)
                pr.start()
                #pr.join()
                #pr.terminate()

            # retrieve results
            for m in xrange(proc_size):
                res = out_q.get()
                j_min = res[0]
                #print i, j_min

                self.chol_delta4[j_min:j_min+loopsize] = res[1]
                self.hdl_delta4[j_min:j_min+loopsize] = res[2]
                self.chol_entry_order.append(res[0])

        t1 = time.time()
        if self.verbose == True:
            print('match time elapsed: %g sec' % (t1-t0))



        # 0.04 percent of population for some reason gets a nan in the delta file.
        # this needs checking, but since it's such a small amount of individuals
        # assign 0 change for these individuals for the moment
        delta_nans = np.isnan(self.chol_delta4)
        self.chol_delta4[delta_nans] = 0
        self.chol_delta1 = self.chol_delta4/4.0

        self.hdl_delta4[delta_nans] = 0
        self.hdl_delta1 = self.hdl_delta4/4.0


        # delete internal arrays
        del a,dage,dbmi,b,q,dchol,age_idx,bmi_idx,chol_idx,sex_idx,eq_idx




    def smokematch_core_process(self, l, t, j_min, loopsize, out_q):
        '''routine to search for matching individuals in the ELSA dataset
        to be used for the multiprocessing module'''

        l.acquire()
        #rnd.seed(self.randseed+t*3+j_min)

        nomatch_cnt = 0
        smoke_new = np.zeros(loopsize)

        for j_idx, j in enumerate(range(j_min,j_min + loopsize)):
            rnd.seed(self.randseed+j_idx + t*3)
            try:
                D = self.C1_smoking[self.smoking_IND[j],6:] # diff values
                W = self.C1_smoking_weights[self.smoking_IND[j],6:] # weight values

                D_finite = np.isfinite(D)
                W_finite = np.isfinite(W)

                d = D[D_finite]
                w_raw = W[W_finite]
                w = w_raw/w_raw.sum()


            except:
                # no match possible, assign no change
                nomatch_cnt += 1
                d = np.array(self.q_smoke_cat[j,t])
                w = np.array([1]) # in this case, weight is assumed to be 1, as no other value than current value is a likely candidate
                if self.verbose == True:
                    print('no change assigned to smoke match for individual %d' % j)

            #smoke2[i] = np.random.choice(d)
            try:
                smoke_new[j_idx] = hadd.weightedChoice(w,d)
                #smoke_new[j_idx] = rnd.choice(d)
            except:


                # assign random number between 0 and 1 if matching did not succeed
                # print 'no smoking match at time t=%d for ind %d, assigning same value as before

                smoke_new[j_idx] = self.q_smoke_cat[j,t]
                nomatch_cnt += 1



        out_q.put([j_min, smoke_new, nomatch_cnt])
        l.release()



    def MatchELSA_smoking_new(self,t=0):

        E_total = self.ELSA_smoke_numeric

        if self.ELSA_smoke_numeric.shape[1]>7:
            # if this dataset contains weights
            E = E_total[:,[0,2,3,4,5,6,7]]
        else:
            E = E_total.copy()

        # categorize age  40: <55 / 55: 55-64 / 65: 65-74 / 75: 75+
        a = self.age[:,t]

        # bin age in 8 categories (agecat8):
        bins = np.array([0,55,65,75,100])
        dage = np.digitize(a,bins)

        # categorize BMI
        # bmi category:
        b = self.bmi_background[:,t]
        bins = np.array([0,25,30,100])
        dbmi = np.digitize(b,bins)

        # categorize deprivation index
        q = self.qimd.copy()
        q[q<3] = 2
        q[q>3] = 4

        age_idx = dage - 1 # dage goes from 1 to 4
        sex_idx = self.gender
        eq_idx = q-2 # q goes from 2 to 4
        bmi_idx = dbmi-1 # dbmi goes from 1 to 3
        scat_idx = self.q_smoke_cat[:,t]

        n_sex, n_age, n_bmi, n_eq, n_scat = hadd.smoking_vectorsize(E)

        # overall indices for population
        IND =   age_idx * (n_sex * n_bmi * n_eq * n_scat) + \
                        sex_idx * (n_bmi * n_eq * n_scat) + \
                        bmi_idx * (n_eq * n_scat) + \
                        eq_idx * (n_scat) + scat_idx

        self.smoking_IND = np.array(IND,dtype=int)

        self.smoke2 = np.zeros(self.population_size)* np.nan

        loopsize = self.population_size/self.nprocs
        chunksize = self.nprocs * loopsize

        t0 = time.time()

        nomatch_cnt = 0


        for i in xrange(0, self.population_size, chunksize):

            # divide parallelised job into chunks loopsize or rest of array size
            proc_size = self.nprocs
            jstart_array = xrange(i,i+chunksize,loopsize)
            if jstart_array[-1] > self.population_size:
                jstart_array = xrange(i,self.population_size, loopsize)
                proc_size = len(jstart_array)

            out_q = multiprocessing.Queue()
            procs = []


            # send jobs to multiprocessing module
            for jstart in jstart_array:
                pr = multiprocessing.Process(
                            target = self.smokematch_core_process,
                            args = (multiprocessing.Lock(),t, jstart, loopsize, out_q))
                procs.append(pr)
                pr.start()
                #pr.join()
                #pr.terminate()

            # retrieve results
            for m in xrange(proc_size):
                res = out_q.get()
                j_min = res[0]
                #print i, j_min

                self.smoke2[j_min:j_min+loopsize] = res[1]
                self.smoke_entry_order.append(res[0])


        t1 = time.time()
        if self.verbose == True:
            print('match time elapsed: %g sec' % (t1-t0))
            if nomatch_cnt > 0:
                print('nomatch_cnt: %d' % nomatch_cnt)

        # delete internal arrays
        del a,b,q,dage,dbmi,age_idx,sex_idx,eq_idx,bmi_idx,scat_idx



    def bpmatch_core_process(self, l, t, j_min, loopsize, out_q):
        ''' search process within BP match procedure
        for loop over population to be parallelised'''

        l.acquire()
        #rnd.seed(self.randseed)

        sys_delta4 = np.zeros(loopsize)
        dia_delta4 = np.zeros(loopsize)

        #rnd.seed(self.randseed+j_min+t*5)

        for j_idx, j in enumerate(range(j_min,j_min + loopsize)):
            rnd.seed(self.randseed+j_idx + t*4)

            try:
                DS = self.CS1_bp[self.bp_IND[j],6:]
                DD = self.CD1_bp[self.bp_IND[j],6:]

                D_finite = np.isfinite(DS)

                ds = DS[D_finite]
                dd = DD[D_finite]

            except:
                # assign no change, if above goes wrong
                ds = np.zeros(1)
                dd = np.zeros(1)
                if self.verbose == True:
                    print('no change assigned to bpmatch for individual %d' % j)

            try:
                sys_delta4[j_idx] = rnd.choice(ds)
                dia_delta4[j_idx] = rnd.choice(dd)
            except:
                # rematch
                try:
                    DS = self.CS1_alternative_bp[self.bp_IND_alternative[j],4:]
                    DD = self.CD1_alternative_bp[self.bp_IND_alternative[j],4:]

                    D_finite = np.isfinite(DS)

                    ds = DS[D_finite]
                    dd = DD[D_finite]
                except:
                    # assign no change, if above goes wrong
                    ds = np.zeros(1)
                    dd = np.zeros(1)
                    if self.verbose == True:
                        print('no change assigned to bpmatch for individual %d' % j)

                try:
                    sys_delta4[j_idx] = rnd.choice(ds)
                    dia_delta4[j_idx] = rnd.choice(dd)
                except:
                    if self.verbose == True:
                        print ds,dd

        out_q.put([j_min, sys_delta4, dia_delta4])
        l.release()


    def MatchELSA_bp_new(self,t=0):
        '''optimized routine to match BP data from ELSA to population
        to get 4 year trajectories for each individual'''

        E = self.ELSA_bp_numeric
        # match age, gender, bmi, sys, dia to ELSA ELSA_bp_numeric[:,1:7] == [agecat, male, bmicat, bmibin, syscat, diacat]

        # age category:

        a = self.age[:,t]
        # bin age : below 55, 55-70, above 70
        bins = np.array([0,55,70,100])
        dage = np.digitize(a,bins)

        # bmi category:
        b = self.bmi_background[:,t]
        bins = np.array([0,25,30,100])
        dbmi = np.digitize(b,bins)
        # bmi categories are 0, 1, so just deduct
        dbmi1 = dbmi - 1
        # and in the bmibin, everyone that is not 0 is one:
        dbmibin = dbmi1.copy()
        dbmibin[dbmi1!=0] = 1

        # categorise systolic pressure
        sys = self.q_sbp_background[:,t]
        bins = np.array([0,100,110,120,130,140,150,160,170,180,190,300])
        dsys = np.digitize(sys,bins)

        # categorise diastolic pressure
        dia = self.dia_background[:,t]

        bins = np.array([0,60,70,80,90,100,200])
        ddia = np.digitize(dia,bins)

        # calculate global indices

        age_idx = dage - 1 # dage goes from 1 to 4
        sex_idx = self.gender

        bmi_idx = dbmi -1 # dbmi goes from 1 to 3
        bmibin_idx = dbmibin
        sys_idx = dsys - 1
        dia_idx = ddia - 1

        n_sex, n_age, n_bmi, n_bmibin, n_sys, n_dia = hadd.bp_vectorsize(E)
        # overall indices for population
        IND =       age_idx * (n_sex * n_bmi * n_sys * n_dia) + \
                    sex_idx * (n_bmi * n_sys * n_dia) + \
                    bmi_idx * (n_sys * n_dia) + \
                    sys_idx * (n_dia) + dia_idx

        IND = np.array(IND,dtype=int)

        IND_alternative = age_idx * (n_bmibin * n_sys) + \
                        bmibin_idx * (n_sys) + sys_idx
        IND_alternative = np.array(IND_alternative,dtype=int)


        self.bp_IND = IND
        self.bp_IND_alternative = IND_alternative

        # chunksize for internal loop
        loopsize = self.population_size/self.nprocs
        chunksize = self.nprocs * loopsize

        self.sys_delta4 = np.zeros(self.population_size) * np.nan
        self.dia_delta4 = np.zeros(self.population_size) * np.nan

        t0 = time.time()
        for i in xrange(0, self.population_size, chunksize):

            # divide parallelised job into chunks loopsize or rest of array size
            proc_size = self.nprocs
            jstart_array = xrange(i,i+chunksize,loopsize)
            if jstart_array[-1] > self.population_size:
                jstart_array = xrange(i,self.population_size, loopsize)
                proc_size = len(jstart_array)

            out_q = multiprocessing.Queue()
            procs = []


            # send jobs to multiprocessing module
            for jstart in jstart_array:
                pr = multiprocessing.Process(
                            target = self.bpmatch_core_process,
                            args = (multiprocessing.Lock(), t, jstart, loopsize, out_q))
                procs.append(pr)
                pr.start()
                #pr.join()
                #pr.terminate()

            # retrieve results
            for m in xrange(proc_size):
                res = out_q.get()
                j_min = res[0]
                #print i, j_min

                self.sys_delta4[j_min:j_min+loopsize] = res[1]
                self.dia_delta4[j_min:j_min+loopsize] = res[2]
                self.bp_entry_order.append(res[0])

        t1 = time.time()
        if self.verbose == True:
            print('match time elapsed: %g sec' % (t1-t0))

        # break down 4 year changes in change per year
        self.sys_delta1 = self.sys_delta4/4.0
        self.dia_delta1 = self.dia_delta4/4.0


        # delete internal arrays
        del a,b,sys,dia,dsys,dbmibin,dage,dbmi1,dbmi,ddia,bmi_idx,sys_idx






    def bmimatch_core_process(self, l, t, j_min, loopsize, out_q):
        '''search process within BMI match routine
        for loop over population to be parallelised

        this function contains a loop to make it computationally more 'expensive'
        than a single execution, toward optimising parallelisation, since single execution
        has too much overhead'''

        l.acquire()

        self.bmi_entry_order.append(j_min)

        nomatch_cnt = 0
        nomatch_ind = -1

        #rnd.seed(self.randseed)

        bmi_delta4 = np.zeros(loopsize)
        bmi_delta8 = np.zeros(loopsize)

        for j_idx, j in enumerate(range(j_min,j_min + loopsize)):
            rnd.seed(self.randseed+j_idx+t*5)

            try:

                D4 = self.C41_bmi[self.bmi_IND[j],5:]
                D8 = self.C81_bmi[self.bmi_IND[j],5:]

                D_finite = np.isfinite(D4)

                d4 = D4[D_finite]
                d8 = D8[D_finite]

            except:
                # no match possible - no such index in ELSA data
                # assign no change
                d4 = np.zeros(1)
                d8 = np.zeros(1)
                nomatch_cnt = 1
                nomatch_ind = j
                #print 'no BMI match at time t=%d for ind %d, assigning no change' % (t,i)

            #smoke2[i] = np.random.choice(d)
            try:

                bmi_delta4[j_idx] = rnd.choice(d4)
                bmi_delta8[j_idx] = rnd.choice(d8)
            except:
                # rematch

                try:
                    D4 = self.C41_alternative_bmi[self.bmi_IND_alternative[j],4:]
                    D8 = self.C81_alternative_bmi[self.bmi_IND_alternative[j],4:]

                    D_finite = np.isfinite(D4)

                    d4 = D4[D_finite]
                    d8 = D8[D_finite]
                except:
                # assign no change, if above goes wrong
                    d4 = np.zeros(1)
                    d8 = np.zeros(1)
                    if self.verbose == True:
                        print('no change assigned to bmimatch for individual %d' % j)

                try:

                    bmi_delta4[j_idx] = rnd.choice(d4)
                    bmi_delta8[j_idx] = rnd.choice(d8)
                except:
                    #print d4,d8
                    #print 'second BMI match failed at time t=%d for ind %d, assigning no change' % (t,i)
                    bmi_delta4[j_idx] = 0
                    bmi_delta8[j_idx] = 0
                    nomatch_cnt = 1
                    nomatch_ind = j

        out_q.put([j_min, bmi_delta4, bmi_delta8])
        l.release()




    def MatchELSA_bmi_new(self,t=0):
        '''optimized routine to match BMI data from ELSA to population
        to get 4, 8 and 1 year trajectories for each individual'''

        E = self.ELSA_bmi_numeric

        # age category:
        a = self.age[:,t]
        # bin age in 8 categories (agecat8):
        bins = np.array([0,50,55,60,65,70,75,80,100])
        dage8 = np.digitize(a,bins)

        bins = np.array([0,55,65,75,100])
        dage4 = np.digitize(a,bins)

        # bmi category:
        b = self.bmi_background[:,t]
        bins = np.arange(18,44,2)
        bins[0] = 0
        bins[-1] = 100
        dbmi = np.digitize(b,bins)

        s = self.q_smoke_cat[:,t]
        dscat = s.copy()
        dscat[s>1]=2

        age8_idx = dage8 - 1
        age4_idx = dage4 - 1

        sex_idx = self.gender
        bmi_idx = dbmi -1 # dbmi goes from 1 to 3
        scat_idx = dscat

        n_sex, n_age4, n_age8, n_bmi, n_scat = hadd.bmi_vectorsize(E)

        # overall indices for population
        IND =       age8_idx * (n_sex * n_scat * n_bmi) + \
                            sex_idx * (n_scat * n_bmi) + \
                            scat_idx * (n_bmi) + bmi_idx

        IND = np.array(IND,dtype=int)

        IND_alternative = age4_idx * (n_scat * n_bmi) + \
                        scat_idx * (n_bmi) + bmi_idx

        IND_alternative = np.array(IND_alternative,dtype=int)

        self.bmi_delta4 = np.zeros(self.population_size) * np.nan
        self.bmi_delta8 = np.zeros(self.population_size) * np.nan

        self.bmi_IND = IND
        self.bmi_IND_alternative = IND_alternative


        nomatch_cnt = 0
        nomatch_ind = []

        # chunksize for internal loop
        loopsize = self.population_size/self.nprocs
        chunksize = self.nprocs * loopsize

        t0 = time.time()
        for i in xrange(0, self.population_size, chunksize):




            # divide parallelised job into chunks loopsize or rest of array size
            proc_size = self.nprocs
            jstart_array = xrange(i,i+chunksize,loopsize)
            if jstart_array[-1] > self.population_size:
                jstart_array = xrange(i,self.population_size, loopsize)
                proc_size = len(jstart_array)

            out_q = multiprocessing.Queue()
            procs = []



            # send jobs to multiprocessing module
            for jstart in jstart_array:
                pr = multiprocessing.Process(
                            target = self.bmimatch_core_process,
                            args = (multiprocessing.Lock(), t, jstart, loopsize, out_q))
                procs.append(pr)
                pr.start()
                #pr.join()
                #pr.terminate()

            # retrieve results
            for m in xrange(proc_size):
                res = out_q.get()
                j_min = res[0]
                #print i, j_min

                self.bmi_delta4[j_min:j_min+loopsize] = res[1]
                self.bmi_delta8[j_min:j_min+loopsize] = res[2]

                self.res = res
                self.bmi_entry_order.append(res[0])
                # if nomatch appears:
#                if res[3] > 0:
#                    nomatch_cnt += 1
#                    nomatch_ind.append(res[4])


        t1 = time.time()
        if self.verbose == True:
            print('match time elapsed: %g sec' % (t1-t0))

        # if some matches did not work, give the information
        if nomatch_cnt > 0:
            nomatch_perc = (100*float(nomatch_cnt))/self.population_size
            if self.verbose == True:
                print('\tno BMI match at time t=%d for %.2f%% of population, assigned no change' % (t,nomatch_perc))
                print '\tindividuals: ', nomatch_ind

        # break down 4 year changes in change per year
        self.bmi_delta4_1 = self.bmi_delta4/4.0
        # break down 8 year changes in chage per year
        self.bmi_delta8_1 = self.bmi_delta8/8.0

        # delete internal arrays

        del a,b,s,dage4,dage8,dbmi,dscat,age8_idx,age4_idx,sex_idx,bmi_idx,scat_idx



    def RelativeUptake(self,t):
        '''evaluating update factors
        at time point t
        #==============================================================================
        # Base number: 8.1% of those eligible attend in a given year (2012/13 number)
        # Ratios to apply to scaling 8.1% (the base number) up down:
        #
        # Ratio
        # North East      1
        # North West      1.48
        # Yorkshire & Humber  1.32
        # East Midlands   1.62
        # West Midlands   1.27
        # East of England 1.61
        # London          1.48
        # South East Coast 0.80
        # South Central   1.11
        # South West      0.85
        #
        # Male        1
        # Female      1.07
        #
        # age 40-49       0.66
        # age 50-59       1
        # age 60-69       1.40
        # age 70-74       1.50
        #
        # White       1
        # South Asian 1.02
        # Black       0.99
        # Other       0.89
        #
        # Non-smoker  1
        # Smoker      0.77
        #
        # Non-drinker 1
        # Drink trivial/light 0.92
        # Drink moderate+ 1.01
        #
        # qrisk <5        1
        # qrisk 5-9.9     1.40
        # qrisk 10-14.9   1.88
        # qrisk 15-19.9   2.13
        # qrisk 20+       2.40
        #==============================================================================

        This function attempts to be a speedup to the other Uptake function by using
        numpy methods to populate the propensity vectors instead of native python loops
        CJ (31/08/2016) renamed Uptake to RelativeUptake, moved absolute uptake stuff to HEALTH CHECKS ASSESSMENT block in Simulate(), since uptake may depend on whether previous HC offer was accepted (Nick's scenario)
        '''

        # evaluate uptake factors
        # (1) age groups
        p_age = np.zeros(self.population_size)

        p40 = np.where((self.age[:,t] >= 40) * (self.age[:,t] <50))[0]
        p_age[p40] = self.up_age_vec[0]

        p50 = np.where((self.age[:,t] >= 50) * (self.age[:,t] <60))[0]
        p_age[p50] = self.up_age_vec[1]

        p60 = np.where((self.age[:,t] >= 60) * (self.age[:,t] <70))[0]
        p_age[p60] = self.up_age_vec[2]

        p70 = np.where((self.age[:,t] >= 70) * (self.age[:,t] <75))[0]
        p_age[p70] = self.up_age_vec[3]

        # (2) gender
        p_gender = np.zeros(self.population_size)

        p_m = np.where(self.gender == 1)[0]
        p_gender[p_m] = self.up_gender_vec[0]

        p_f = np.where(self.gender == 0)[0]
        p_gender[p_f] = self.up_gender_vec[1]

        # (3) ethnicity
        p_eth = np.zeros(self.population_size)

        p_e1 = np.where(self.eth == 1)[0]
        p_eth[p_e1] = self.up_eth_vec[0]

        p_e234 = np.where((self.eth == 2) + (self.eth==3) + (self.eth==4))[0]
        p_eth[p_e234] = self.up_eth_vec[1]

        p_e589 = np.where((self.eth == 5) + (self.eth == 8) + (self.eth ==9))[0]
        p_eth[p_e589] = self.up_eth_vec[3]

        p_e67 = np.where((self.eth == 6) + (self.eth ==7))[0]
        p_eth[p_e67] = self.up_eth_vec[2]

        # (4) Qrisk
        p_qrisk = np.zeros(self.population_size)
        Q = self.CalculateQrisk(t)

        p_q1 = np.where(Q<5)[0]
        p_qrisk[p_q1] = self.up_QRisk_vec[0]

        p_q2 = np.where((Q>=5)*(Q<10))[0]
        p_qrisk[p_q2] = self.up_QRisk_vec[1]

        p_q3 = np.where((Q>=10) * (Q<15))[0]
        p_qrisk[p_q3] = self.up_QRisk_vec[2]

        p_q4 = np.where((Q>=15) * (Q<20))[0]
        p_qrisk[p_q4] = self.up_QRisk_vec[3]

        p_q5 = np.where(Q>=20)[0]
        p_qrisk[p_q5] = self.up_QRisk_vec[4]

        # (5) SES
        p_ses = np.zeros(self.population_size)

        p_s1 = np.where(self.SES == 1)[0]
        p_ses[p_s1] = self.up_SES_vec[0]

        p_s2 = np.where(self.SES == 2)[0]
        p_ses[p_s2] = self.up_SES_vec[1]

        p_s3 = np.where(self.SES == 3)[0]
        p_ses[p_s3] = self.up_SES_vec[2]

        p_s4 = np.where(self.SES == 4)[0]
        p_ses[p_s4] = self.up_SES_vec[3]

        p_s5 = np.where(self.SES == 5)[0]
        p_ses[p_s5] = self.up_SES_vec[4]

        # (6) smoking
        p_smoking = np.zeros(self.population_size)

        p_nonsm = np.where(self.q_smoke_cat[:,t] <= 1)[0]
        p_smoking[p_nonsm] = self.up_smoker_vec[0]

        p_sm = np.where(self.q_smoke_cat[:,t] > 1)[0]
        p_smoking[p_sm] = self.up_smoker_vec[1]

        self.RelTurnup_woSES =  p_age * p_gender * p_eth * p_qrisk * p_smoking

        self.RelTurnup =   p_age * p_gender * p_eth * p_qrisk * p_ses * p_smoking

        # delete internal arrays
        del p_age,p_gender,p_eth,p_smoking,p_nonsm,p_sm,p_ses,p_qrisk,p_s1,p_s2,p_s3,p_s4,p_s5
        del p_q1,p_q2,p_q3,p_q4,p_q5,p_m,p_f
        del p40,p50,p60,p70



    def CalculateCAIDE(self, t=0, output=True):
        '''calculates the CAIDE Dementia Risk Score at time t
        The risk score is:

        age:    < 47 yrs        score 0
                47-53 yrs       score 3
                > 53 yrs        score 4

        education:
                >= 10 yrs       score 0
                7-9 yrs         score 2
                0-6 yrs         score 3

        sex:    female          score 0
                male            score 1

        SBP:    <= 140          score 0
                > 140           score 2

        BMI:    <= 30           score 0
                > 30            score 2

        Chol:   <= 6.5          score 0
                > 6.5           score 2

        PA:     active          score 0
                inactive        score 1


        Risk:
                score 0-5:       1.0%
                score 6-7:       1.9%
                score 8-9:       4.2%
                score 10-11:     7.4%
                score 12-15:    16.4%

        if output==True, then the CAIDE score for timepoint t will be returned
        '''

        # bin age and assign CAIDE score
        a = self.age[:,t]
        # bin age in 3 groups
        bins = np.array([0,47,53,100])
        Cage = np.digitize(a,bins)
        # replace agebins with score
        replace = {1:0, 2:3, 3:4}
        CAIDE_age = Cage.copy()
        for k, v in replace.iteritems(): CAIDE_age[Cage==k] = v

        # bin SBP
        CAIDE_bp = np.zeros(self.population_size)
        CAIDE_bp[self.q_sbp[:,t]>140] = 2

        # bin BMI
        CAIDE_bmi = np.zeros(self.population_size)
        CAIDE_bmi[self.bmi[:,t]>30] = 2

        # bin Cholesterol
        CAIDE_chol = np.zeros(self.population_size)
        CAIDE_chol[self.chol[:,t]>6.5] = 2

        CAIDE_PA = 1-self.PA[:,t]

        caide_score = self.gender + self.educ + CAIDE_age + CAIDE_bp + CAIDE_bmi + CAIDE_chol + CAIDE_PA

        self.CAIDE[:,t] = caide_score
        # risks based on caide score
        risks = np.ones(self.population_size)
        risks[caide_score >= 6] = 1.9
        risks[caide_score >= 8] = 4.2
        risks[caide_score >= 10] = 7.4
        risks[caide_score >= 12] = 16.4

        self.DementiaRisk[:,t] = risks

        # delete internal arrays
        del a,Cage,CAIDE_age, CAIDE_bp, CAIDE_bmi, CAIDE_chol, CAIDE_PA

        if output==True:
            return caide_score





    def CalculateQrisk(self, timestep):
        '''calculates the Qrisk for individuals, based on the QRisk algorithms written in C
        Q_80_model_4_0.c (females)
        Q_80_model_4_1.c (males)
        published on
        http://svn.clinrisk.co.uk/opensource/qrisk2/
        accessed last time 15/Aug/2014

        function which should give a speedup to the self.CalculateQrisk() function by eliminating for loops.
        '''

        # initiate score
        score = np.zeros(self.population_size)

        i = timestep

        extralhr_male = np.log(self.up_Statins_eff_extra_male)
        extralhr_female = np.log(self.up_Statins_eff_extra_female)

        # d = MALE ------------------------------------------------------------
        # only apply calculations for male population
        incl = self.gender == 1

        # Applying the fractional polynomial transforms
        # (which includes scaling)
        dage = self.age[incl, i]
        dage = dage / 10.0
        age_1 = pow(dage, -1)
        age_2 = pow(dage, 2)
        dbmi = self.bmi[incl, i]
        dbmi = dbmi / 10.0
        bmi_1 = pow(dbmi, -2)
        bmi_2 = pow(dbmi, -2) * np.log(dbmi)
        # centering continuous variables
        age_1 = age_1 - 0.232008963823318
        age_2 = age_2 - 18.57763671875
        bmi_1 = bmi_1 - 0.146408438682556
        bmi_2 = bmi_2 - 0.140651300549507
        rati = self.q_rati[incl, i] - 4.377167701721191
        sbp = self.q_sbp[incl, i] - 131.038314819335940
        town = self.q_town[incl] - 0.151332527399063
        surv = 0.977699398994446

        # start sum
        self.a = np.zeros(incl.sum())


        # populate the smoker risk of the current population, based on Iethrisk
        iethrisk = np.zeros(incl.sum())

        ismoke = np.zeros(incl.sum())
        for j in range(5):
            iethrisk[self.eth[incl] == j] = self.parms['Iethrisk_male'][j]
            ismoke[(self.q_smoke_cat[incl, i]) == j] = self.parms['Ismoke_male'][j]


        self.a += iethrisk
        self.a += ismoke

        # Sum from continuous values

        self.a += age_1 * -17.622554338194561
        self.a += age_2 * 0.024187318929827364
        self.a += bmi_1 * 1.7320282704272665
        self.a += bmi_2 * -7.2311754066699754
        self.a += rati * 0.17513879740122351
        self.a += sbp * 0.01016763051791969
        self.a += town * 0.029817727149672096

        # Sum from boolean values

        self.a += self.q_b_AF[incl, i] * 0.98909975261894023
        self.a += self.q_b_ra[incl, i] * 0.25418862091186112
        self.a += self.q_b_renal[incl, i] * 0.794978923043832
        self.a += self.q_b_treatedhyp[incl, i] * 0.62293594798680441
        self.a += self.q_b_type1[incl, i] * 1.333035332146393
        self.a += self.q_b_type2[incl, i] * 0.93729568281519404
        self.a += self.q_fh_cvd[incl] * 0.59233537365824229

        # sum from interaction terms
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 1) * 0.9243747443632776
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 2) * 1.9597527500081284
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 3) * 2.9993544847631153
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 4) * 5.03707352547681
        self.a += age_1 * self.q_b_AF[incl, i] * 8.2354205455482727
        self.a += age_1 * self.q_b_renal[incl, i] * -3.9747389951976779
        self.a += age_1 * self.q_b_treatedhyp[incl, i] * 7.8737743159167728
        self.a += age_1 * self.q_b_type1[incl, i] * 5.4238504414460937
        self.a += age_1 * self.q_b_type2[incl, i] * 5.0624161806530141
        self.a += age_1 * bmi_1 * 33.543752516739424
        self.a += age_1 * bmi_2 * -129.97667382572038
        self.a += age_1 * self.q_fh_cvd[incl] * 1.9279963874659789
        self.a += age_1 * sbp * 0.05234408921756202
        self.a += age_1 * town * -0.17305880749635402
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 1) * -0.0034466074038854394
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 2) * -0.0050703431499952954
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 3) * 0.00032160597999164408
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 4) * 0.0031312537144240087
        self.a += age_2 * self.q_b_AF[incl, i] * 0.0073291937255039966
        self.a += age_2 * self.q_b_renal[incl, i] * -0.026155707328653178
        self.a += age_2 * self.q_b_treatedhyp[incl, i] * 0.0085556382622618121
        self.a += age_2 * self.q_b_type1[incl, i] * 0.0020586479482670723
        self.a += age_2 * self.q_b_type2[incl, i] * -0.00023285907708541729
        self.a += age_2 * bmi_1 * 0.081184721208079499
        self.a += age_2 * bmi_2 * -0.25589190688509483
        self.a += age_2 * self.q_fh_cvd[incl] * -0.0056729073729663406
        self.a += age_2 * sbp * -0.000053658425730729933
        self.a += age_2 * town * -0.0010763305052605857
        self.a += self.on_statins[incl,i] * extralhr_male

        score[incl] = 100.0 * (1 - pow(surv, np.exp(self.a)))

        # FEMALE ------------------------------------------------------------
        # only apply calculations for female population
        incl = self.gender == 0

        # Applying the fractional polynomial transforms
        # (which includes scaling)
        dage = self.age[incl, i]
        dage = dage / 10.0
        age_1 = pow(dage, 0.5)
        age_2 = dage
        dbmi = self.bmi[incl, i]
        dbmi = dbmi / 10.0
        bmi_1 = pow(dbmi, -2)
        bmi_2 = pow(dbmi, -2) * np.log(dbmi)
        # centering continuous variables
        age_1 = age_1 - 2.099778413772583
        age_2 = age_2 - 4.409069538116455
        bmi_1 = bmi_1 - 0.154046609997749
        bmi_2 = bmi_2 - 0.144072100520134
        rati = self.q_rati[incl, i] - 3.554229259490967
        sbp = self.q_sbp[incl, i] - 125.773628234863280
        town = self.q_town[incl] - 0.032508373260498
        surv = 0.988948762416840

        # start sum
        self.a = np.zeros(incl.sum())

        # populate the smoker risk of the current population, based on Iethrisk
        iethrisk = np.zeros(incl.sum())
        ismoke = np.zeros(incl.sum())
        for j in range(5):
            iethrisk[self.eth[incl] == j] = self.parms['Iethrisk_female'][j]
            ismoke[(self.q_smoke_cat[incl, i]) == j] = self.parms['Ismoke_female'][j]

        self.a += iethrisk
        self.a += ismoke

        # Sum from continuous values

        self.a += age_1 * 3.8734583855051343
        self.a += age_2 * 0.13466343044783846
        self.a += bmi_1 * -0.15578724033330626
        self.a += bmi_2 * -3.7727795566691125
        self.a += rati * 0.15256952089196796
        self.a += sbp * 0.013216530011965356
        self.a += town * 0.064364752986401708

        # Sum from boolean values

        self.a += self.q_b_AF[incl, i] * 1.4235421148946676
        self.a += self.q_b_ra[incl, i] * 0.30214625115536481
        self.a += self.q_b_renal[incl, i] * 0.861474303972141640
        self.a += self.q_b_treatedhyp[incl, i] * 0.58893554587337038
        self.a += self.q_b_type1[incl, i] * 1.6684783657502795
        self.a += self.q_b_type2[incl, i] * 1.1350165062510138
        self.a += self.q_fh_cvd[incl] * 0.51339727757386733

        # sum from interaction terms
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 1) * 0.6891139747579299
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 2) * 0.69426328021216266
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 3) * -1.6952388644218186
        self.a += age_1 * (self.q_smoke_cat[incl, i] == 4) * -1.2150150940219255
        self.a += age_1 * self.q_b_AF[incl, i] * -3.5855215448190969
        self.a += age_1 * self.q_b_renal[incl, i] * -3.0766647922469192
        self.a += age_1 * self.q_b_treatedhyp[incl, i] * -4.0295302811880314
        self.a += age_1 * self.q_b_type1[incl, i] * -0.33441105674057786
        self.a += age_1 * self.q_b_type2[incl, i] * -3.3144806806620530
        self.a += age_1 * bmi_1 * -5.59339057972300060
        self.a += age_1 * bmi_2 * 64.363557283768898
        self.a += age_1 * self.q_fh_cvd[incl] * 0.860543376121715720
        self.a += age_1 * sbp * -0.05093211545511885900
        self.a += age_1 * town * 0.1518664540724453700
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 1) * -0.17653954858826815
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 2) * -0.23238364832785730
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 3) * 0.27343957705518263
        self.a += age_2 * (self.q_smoke_cat[incl, i] == 4) * 0.143255228745415270
        self.a += age_2 * self.q_b_AF[incl, i] * 0.49868713908070322
        self.a += age_2 * self.q_b_renal[incl, i] * 0.439303361566493860
        self.a += age_2 * self.q_b_treatedhyp[incl, i] * 0.69043857903032502
        self.a += age_2 * self.q_b_type1[incl, i] * -0.17343165660603277
        self.a += age_2 * self.q_b_type2[incl, i] * 0.48649306558679495
        self.a += age_2 * bmi_1 * 1.5223341309207974
        self.a += age_2 * bmi_2 * -12.741343620796407
        self.a += age_2 * self.q_fh_cvd[incl] * -0.27567084814151099
        self.a += age_2 * sbp * 0.0073790750039744186
        self.a += age_2 * town * -0.04874654626796409
        self.a += self.on_statins[incl,i] * extralhr_female

        score[incl] = 100.0 * (1 - pow(surv, np.exp(self.a)))

        # delete internal arrays
        del dage,age_1,age_2,dbmi,bmi_1,bmi_2,rati,sbp,town

        return score


    def GetCVDRisks(self,t):
        '''get risk for stroke and IHD for the population at time t
        expressed as probability to get either within the next year'''

        rel_risk = np.zeros(self.population_size)
        risk_stroke = np.zeros(self.population_size)

        male = self.gender == 1
        female = self.gender == 0
        current_age = np.array(self.age[:,t],dtype=int) # vector of current ages in integers

        rem = current_age - 90
        rem[rem<0] = 0
        current_age_capped = current_age.copy()
        current_age_capped[current_age>90] = 90

        # for those over 90, risks are taken from previous QRisk projections, as the risk is computed over the next 10 years
        # for men
        rel_risk[male] = self.CVDprob_year_male[current_age_capped[male],0]
        risk_stroke[male] = self.CVDprob_event_male[current_age_capped[male]]
        # for women
        rel_risk[female] = self.CVDprob_year_female[current_age_capped[female],0]
        risk_stroke[female] = self.CVDprob_event_female[current_age_capped[female]]

        # delete internal arrays
        del male,female,current_age,rem,current_age_capped

        return rel_risk, risk_stroke



    def GetDementiaRisks(self,t, age_dependent_risk = False):
        '''determining risk for dementia for the population at time ti
        expressed as probability to get dementia within the next year

        the variable 'age_dependent_risk' indicates whether CAIDE is age-adjusted.

        if age_dependent_risk=True:
            then people within the age limit for which CAIDE
        is designed (assumed to be 40-60) will be assigned a yearly risk coming from the
        mean of the 20-year risk, if age is below, then they are assigned the first relative risk element
        in the 20-year risk vector coming from CAIDE, and if above 65, they are assigned
        the x-th element in the 20-year risk vector, where x is min(age-65, 20)

        if age_dependent_risk=False:
            everyone is assigned a yearly risk coming from the
        mean of the 20-year risk
        '''

        rel_risk = np.zeros(self.population_size)
        CAIDE_risk = np.zeros(self.population_size)
        male = self.gender == 1
        female = self.gender == 0
        current_age = np.array(self.age[:,t],dtype=int) # vector of current ages in integers

        # for those over 80, risks are taken from previous years, as the risk is computed over the next 20 years
        rem = current_age - 80
        rem[rem<0] = 0
        current_age_capped = current_age.copy()
        current_age_capped[current_age>80] = 80


        if age_dependent_risk == False:
            rel_risk[male] = self.DEMprob_year_male[current_age_capped[male],:].mean()
            rel_risk[female] = self.DEMprob_year_female[current_age_capped[female],:].mean()
            CAIDE_risk = self.DementiaRisk[:,t]/100.0
        else:
            # evaluate who is at what age
            age_below = current_age_capped < self.parms['CAIDE_age_bounds'][0]
            v1 = male * age_below
            w1 = female * age_below
            rel_risk[v1] = self.DEMprob_year_male[current_age_capped[v1],0]
            rel_risk[w1] = self.DEMprob_year_female[current_age_capped[w1],0]


            age_within = (current_age_capped >= self.parms['CAIDE_age_bounds'][0]) * (current_age_capped <= self.parms['CAIDE_age_bounds'][1])
            v3 = male * age_within
            w3 = female * age_within
            #rel_risk[v3] = self.DEMprob_year_male[current_age_capped[v3],:].mean()
            #rel_risk[w3] = self.DEMprob_year_female[current_age_capped[w3],:].mean()
            rel_risk[v3] = self.DEMprob_year_male[current_age_capped[v3],0]
            rel_risk[w3] = self.DEMprob_year_female[current_age_capped[w3],0]

            ##################################
            # if age is above 60,
            # then take CAIDE from age 60, and apply lifetable based incidences
            # scaled by risk factors (ie. everyone stays on percentile)

            age_above = current_age_capped >= self.parms['CAIDE_age_bounds'][1]
#            v2 = male * age_above
#            w2 = female * age_above

            j_male = (self.age[:,t] >= self.parms['CAIDE_age_bounds'][1]) * (self.gender==1) * (self.alive[:,t]==1)
            j_female = (self.age[:,t] >= self.parms['CAIDE_age_bounds'][1]) * (self.gender==0) * (self.alive[:,t]==1)

            j_age_male = np.array(self.age[j_male,t],dtype=int)
            j_age_female = np.array(self.age[j_female,t],dtype=int)

            # get mean incidence for that age (male) from Life table
            LT_dem_male = np.zeros((100,2))
            LT_dem_female = np.zeros((100,2))
            # copy Life Table data for males
            LT_dem_male[:,0] = self.LT_dementia[:100,0]
            LT_dem_male[:,1] = self.LT_dementia[:100,1]
            # and for females
            LT_dem_female[:,0] = self.LT_dementia[:100,0]
            LT_dem_female[:,1] = self.LT_dementia[:100,4]

            L_male = LT_dem_male[j_age_male,1]
            L_female = LT_dem_female[j_age_female,1]

            # from those who are above 60, where in the simulation were they at age 60?
            age_above_male = self.parms['CAIDE_age_bounds'][1] - j_age_male + t
            age_above_female = self.parms['CAIDE_age_bounds'][1] - j_age_female + t

            # if age60 is negative, then this means that individual was already over 60 at start of simulation
            # so give them first CAIDE value
            age_above_male[age_above_male<0]=0
            age_above_female[age_above_female<0]=0

            # what is the CAIDE derived risk for those people?
            C_select_male = self.DementiaRisk[j_male,age_above_male]/100.0
            C_select_female = self.DementiaRisk[j_female,age_above_female]/100.0

            # at a given age, how many are there with one of the 5 risks?
            nC_male = np.zeros((100,5))
            nC_female = np.zeros((100,5))
            C = np.unique(self.DementiaRisk[self.DementiaRisk>0])/100.0

            for a in range(100):
                for i in range(5):
                    nC_male[a,i] = (C_select_male[j_age_male==a]==C[i]).sum()
                    nC_female[a,i] = (C_select_female[j_age_female==a]==C[i]).sum()
                #nC[i] = (C_all==C[i]).sum()

            nC_select_male = nC_male[j_age_male]
            nC_select_female = nC_female[j_age_female]

            # (1) scaling factor for C
            try:
                with np.errstate(invalid='ignore'):
                    x_male = (L_male[:,np.newaxis] * nC_select_male.sum(axis=1)[:,np.newaxis])/(C.sum() * nC_select_male)
                    x_female = (L_female[:,np.newaxis] * nC_select_female.sum(axis=1)[:,np.newaxis])/(C.sum() * nC_select_female)

                    # which Caide risks are there actually for each individual, ie. which column in the C vector
                    c_index_male = np.digitize(C_select_male,C)-1
                    c_index_female = np.digitize(C_select_female,C)-1

                    x_select_male = x_male[np.arange(c_index_male.size),c_index_male]
                    x_select_female = x_female[np.arange(c_index_female.size),c_index_female]

                    rel_risk[j_male] = x_select_male
                    rel_risk[j_female] = x_select_female
            except:
                print i,C_select_male,C_select_female

            # CAIDE-based risks are taken from current CAIDE
            CAIDE_risk = self.DementiaRisk[:,t]/100.0
            # however, if person is above the upper bounds of CAIDE (eg 60 years)
            # this will be based on the risk at age 60 or at beginning of simulation, if
            # individual was over 60 at the start
            # this accounts for the assumption that no intervention can change Dementia risk after age 60
            CAIDE_risk[j_male] = C_select_male
            CAIDE_risk[j_female] = C_select_female

        # delete internal arrays
        del male,female,j_male,j_female,j_age_male,j_age_female,current_age,rem,current_age_capped

        return rel_risk, CAIDE_risk



    def Simulate(self,bmimatchyrs=4,newroutines=True):
        '''runs a simulation over a few years
        bmimatchyrs: is BMI matched every 4 or 8 years? (default:4)
        '''

        # measure time
        t0 = time.time()

        h = ['without Health Check', 'with Health Checks']

        for i in range(0, self.simulation_time):

            # define if random seed is globally applied at beginning of each timestep:
            if self.RandomSeed['Start_Of_Each_Timestep'] == True:
                np.random.seed(i+10101)
                rnd.seed(self.randseed)
            else:
                np.random.seed()
                rnd.seed(self.randseed)

            if self.verbose == True:
                print('\n-------------------------\n** simulating time t=%d (%s)...' % (i,h[self.Health_Checks]))

            ############################
            #   MATCH TRAJECTORIES
            ############################
            # at beginning and at every 4 years, rematch population

            # (1) for BMI
            # (2) for BP
            # (3) for smoking
            if bmimatchyrs == 4:
                if np.mod(i,4) == 0:
                    if self.verbose == True:
                        print('\tmatching BMI...')
                    if newroutines==False:
                        self.MatchELSA_bmi(t=i)
                    else:
                        self.MatchELSA_bmi_new(t=i)
            else:
                if np.mod(i,8) == 0:
                    if self.verbose == True:
                        print('\tmatching BMI...')
                    if newroutines==False:
                        self.MatchELSA_bmi(t=i)
                    else:
                        self.MatchELSA_bmi_new(t=i)

            # bmi: apply the yearly delta

            if i<(self.simulation_time-1):
                if bmimatchyrs == 4:
                    self.bmi_background[:,i+1] = self.bmi_background[:,i] + self.bmi_delta4_1
                else:
                    self.bmi_background[:,i+1] = self.bmi_background[:,i] + self.bmi_delta8_1
                # cap BMI at 50 should values go over this threshold
                # likewise at 10, should values go under this threshold
                bmimin = 10
                bmimax = 50
                bmi_toolow = self.bmi_background[:,i+1] < bmimin
                self.bmi_background[bmi_toolow,i+1] = bmimin
                bmi_toohigh = self.bmi_background[:,i+1] > bmimax
                self.bmi_background[bmi_toohigh,i+1] = bmimax


            ##################################################################
            # WEIGHT MANAGEMENT TREATMENT - 'shifted trajectory
            ##################################################################
            bmi_change = np.zeros(self.population_size)

            if self.Treatment_WeightReduction == True:
                wr = self.on_wr[:,i] == 1
                bmi_change[wr] = self.up_Weight_eff
                ## Maximum weight loss after one year, but this is regained over 5 years
                wr_year = i - self.prev_wr_time
                bmi_change[wr * (wr_year == 2)] = self.up_Weight_eff * 3/5 # ref Amy Ahern, unpublished.
                bmi_change[wr * (wr_year == 3)] = self.up_Weight_eff * 2/5 # linearly interpolate...
                bmi_change[wr * (wr_year == 4)] = self.up_Weight_eff * 1/5
                bmi_change[wr * (wr_year >= 5)] = 0

            # then apply changes to level for this year
            self.bmi[:,i] = self.bmi_background[:,i]  +  bmi_change


            # bp: apply the yearly delta
            if np.mod(i,4) == 0:
                if self.verbose == True:
                    print('\tmatching BP...')
                if newroutines==False:
                    self.MatchELSA_bp(t=i)
                else:
                    self.MatchELSA_bp_new(t=i)
            if i<(self.simulation_time-1):


                self.dia_background[:,i+1] = self.dia_background[:,i] + self.dia_delta1
                self.q_sbp_background[:,i+1] = self.q_sbp_background[:,i] + self.sys_delta1

                # restrict bp values to HSE bounds, should they be outside as a result of matching process
                sys_toolow = self.q_sbp_background[:,i+1] < self.q_sbp_background[:,0].min()
                self.q_sbp_background[sys_toolow,i+1] = self.q_sbp_background[:,0].min()
                sys_toohigh = self.q_sbp_background[:,i+1] > self.q_sbp_background[:,0].max()
                self.q_sbp_background[sys_toohigh,i+1] = self.q_sbp_background[:,0].max()

                dia_toolow = self.dia_background[:,i+1] < self.dia_background[:,0].min()
                self.dia_background[dia_toolow,i+1] = self.dia_background[:,0].min()
                dia_toohigh = self.dia_background[:,i+1] > self.dia_background[:,0].max()
                self.dia_background[dia_toohigh,i+1] = self.dia_background[:,0].max()



            ##################################################################
            # AHT TREATMENT - 'shifted trajectory
            ##################################################################
            # apply ATH treatment directly at same time point

            sbp_change = np.zeros(self.population_size)
            dia_change = np.zeros(self.population_size)


            # who is currently on Hypertensives treatment
            aht = self.on_aht[:,i] == 1

            # reduction of SBP/DBP by gender and age
            age55minus = self.age[:,i] <= 55
            age55plus = self.age[:,i] > 55

            male = self.gender == 1
            female = self.gender == 0

            if self.Treatment_AHT == True:
                # find changes
                sbp_change[(aht * age55minus)] = self.up_AHT_eff_age55minus_SBP
                dia_change[(aht * age55minus)] = self.up_AHT_eff_age55minus_DBP
                sbp_change[(aht * age55plus * male)] = self.up_AHT_eff_age55plus_SBP_male
                sbp_change[(aht * age55plus * female)] = self.up_AHT_eff_age55plus_SBP_female
                dia_change[(aht * age55plus * male)] = self.up_AHT_eff_age55plus_DBP_male
                dia_change[(aht * age55plus * female)] = self.up_AHT_eff_age55plus_DBP_female

            self.q_sbp[:,i] = self.q_sbp_background[:,i] + sbp_change
            self.dia[:,i] = self.dia_background[:,i] +  dia_change


            # restrict bp values to HSE bounds, should they be outside as a result of matching process
            sys_toolow = self.q_sbp[:,i] < self.q_sbp[:,0].min()
            self.q_sbp[sys_toolow,i] = self.q_sbp[:,0].min()
            sys_toohigh = self.q_sbp[:,i] > self.q_sbp[:,0].max()
            self.q_sbp[sys_toohigh,i] = self.q_sbp[:,0].max()

            dia_toolow = self.dia[:,i] < self.dia[:,0].min()
            self.dia[dia_toolow,i] = self.dia[:,0].min()
            dia_toohigh = self.dia[:,i] > self.dia[:,0].max()
            self.dia[dia_toohigh,i] = self.dia[:,0].max()


            # cholesterol
            if np.mod(i,4) == 0:
                if self.verbose == True:
                    print('\tmatching Cholesterol...')
                if newroutines==False:
                    self.MatchELSA_chol(t=i)
                else:
                    self.MatchELSA_chol_new(t=i)

            if i<(self.simulation_time-1):




                # background chol levels don't change with treatments
                self.chol_background[:,i+1] = self.chol_background[:,i] + self.chol_delta1

                # same for background
                chol_toolow = self.chol_background[:,i+1] < self.chol_background[:,0].min()
                self.chol_background[chol_toolow,i+1] = self.chol_background[:,0].min()
                chol_toohigh = self.chol_background[:,i+1] > self.chol_background[:,0].max()
                self.chol_background[chol_toohigh,i+1] = self.chol_background[:,0].max()




            ##################################################################
            # STATIN TREATMENT - 'shifted trajectory
            ##################################################################
            # apply statin treatment directly at same time point

            # first check who is currently on statin treatment
            stat_male = self.on_statins[:,i] * (self.gender==1)
            stat_female = self.on_statins[:,i] * (self.gender==0)

            chol_change = np.zeros(self.population_size)

            if self.Treatment_Statins == True:
                chol_change[stat_male] = self.up_Statins_eff_male
                chol_change[stat_female] = self.up_Statins_eff_female

            # then apply changes to cholesterol level for this year
            self.chol[:,i] = self.chol_background[:,i] +  chol_change

            # restrict those cholesterol values that are outside the HSE bounds for cholesterol into the bounds
            chol_toolow = self.chol[:,i] < self.chol[:,0].min()
            self.chol[chol_toolow,i] = self.chol[:,0].min()
            chol_toohigh = self.chol[:,i] > self.chol[:,0].max()
            self.chol[chol_toohigh,i] = self.chol[:,0].max()

            # if no treatment, hdl is remaining the same as at start
            self.hdl[:,i] = self.hdl[:,0]
            if self.Treatment_Statins == True:
                # if on statins treatment, hdl is increased

                self.hdl[stat_male,i] = self.hdl[stat_male,0] +  self.up_Statins_eff_HDL_male
                self.hdl[stat_female,i] = self.hdl[stat_female,0] +  self.up_Statins_eff_HDL_female

            self.q_rati[:,i] = self.chol[:,i] / self.hdl[:,i]





            # smoking:
            if np.mod(i,2) == 0:
                if self.verbose == True:
                    print('\tmatching smoking...')
                if newroutines==False:
                    self.MatchELSA_smoking(t=i)
                else:
                    self.MatchELSA_smoking_new(t=i)



            if i>2:
                if self.adjust_smoking_relapse_rates == True:
                    ##### ADJUST RELAPSE RATES
                    # this is done looking back
                    ## currently ~2% relapse within 5 years. add 35% to this, and assume that the relapse happens within 2 years of quitting
                    # (Hawkings 2010 show that 75% of relapses happen within 2 years of quitting, so this is a simplified assumption of that finding)

                    # check who quitted 2 years ago
                    quitted = (self.q_smoke_cat[:,i-3] > 1) * (self.q_smoke_cat[:,i-2] == 1)
                    # of those, 30% take up smoking again within next 2 years
                    r = np.random.random(self.population_size)

                    relapsing = quitted * (r<self.up_Smoking_relapsing_rate)
                    self.smoke2[relapsing] = 3


            if i<(self.simulation_time-1):

                self.q_smoke_cat[:,i+1] = self.smoke2

                #### ADJUST SMOKING QUIT RATES
                # should be in the ~5% range
                if self.adjust_smoking_quit_rates == True:
                    this_year_smoker = self.q_smoke_cat[:,i] > 1
                    next_year_exsmoker = self.q_smoke_cat[:,i+1] == 1
                    next_year_quitting = this_year_smoker * next_year_exsmoker


                    # reduce the amount of quitters by 25% to get towards the 5% range as at the moment they are about 11%
                    r = np.random.random(self.population_size)

                    quit_reverse = next_year_quitting * (r<0.25)

                    self.q_smoke_cat[quit_reverse,i+1] = 3
                    self.smoke2[quit_reverse] = 3


                    # update smoking quit rates after weighting
                    this_year_smoker = self.q_smoke_cat[:,i] > 1
                    next_year_exsmoker = self.q_smoke_cat[:,i+1] == 1
                    next_year_quitting = this_year_smoker * next_year_exsmoker

                    self.background_quit_rate = next_year_quitting.sum()/float(this_year_smoker.sum())
                    #print old_quit_rate, self.background_quit_rate


                    # should someone being an ex-smoker revert back to 0 in the matching process, override with current value
                    revert_to_zero = (self.q_smoke_cat[:,i]>0) * (self.q_smoke_cat[:,i+1]==0)
                    self.q_smoke_cat[revert_to_zero,i+1] = self.q_smoke_cat[revert_to_zero,i]
                else:
                    self.background_quit_rate = 0.05
            else:
                self.background_quit_rate = 0.05



            # GlyHb
            if np.mod(i,4) == 0:
                if self.verbose == True:
                    print('\tmatching GlyHb...')
                if newroutines==False:
                    self.MatchELSA_glyhb(t=i)
                else:
                    self.MatchELSA_glyhb_new(t=i)

                # assign


            if i<(self.simulation_time-1):
                self.glyhb[:,i+1] = self.glyhb[:,i] + self.glyhb_delta1

                # cap glyhb at values above or below min / max of baseline population

                glyhb_toolow = self.glyhb[:,i+1] < self.glyhb[:,0].min()
                self.glyhb[glyhb_toolow,i+1] = self.glyhb[:,0].min()

                glyhb_toohigh = self.glyhb[:,i+1] > self.glyhb[:,0].max()
                self.glyhb[glyhb_toohigh,i+1] = self.glyhb[:,0].max()



            # if any of the risk factors are outside margins of HSE dataset, reset them to margins






#



            #########################################
            # CALCULATE QRISK
            #########################################


            self.QRisk[:,i] = self.CalculateQrisk(i)

            age = np.array(self.age[:,i], dtype=int)
            male = self.gender == 1
            female = self.gender == 0



            # reset back seeds to timestep setting

            if self.RandomSeed['Start_Of_Each_Timestep'] == True:
                np.random.seed(i+10101)
                rnd.seed(self.randseed)
            else:
                np.random.seed()
                rnd.seed(self.randseed)


            #########################################
            # HEALTH CHECKS ASSESSMENT
            #########################################

            if self.Health_Checks == True:

                # evaluate relative uptake at current time point
                self.RelativeUptake(i)
                rnd.seed(self.randseed)

                # -----------------------------------------------------------------
                # check Health Checks eligibility
                self.on_register[:, i] = (self.bp[:, i] > 0) + (self.compm7[:, i] > 0) + (self.diabetes[:, i] > 0)
                self.within_age[:, i] = (self.age[:, i] >= self.up_HC_age_limit[0]) * (self.age[:, i] <= self.up_HC_age_limit[1])

                # depending on whether people on diabetes, blood pressure and CVD registers are allowed to take part, assign eligibility

                if self.up_HC_include_diabetes_registers == False:
                    db_reg = self.diabetes[:, i]
                else:
                    db_reg = np.zeros(self.population_size, dtype=bool)

                if self.up_HC_include_CVD_registers == False:
                    cvd_reg = self.compm7[:, i]
                else:
                    cvd_reg = np.zeros(self.population_size, dtype=bool)

                if self.up_HC_include_bp_registers == False:
                    bp_reg = self.bp[:, i]
                else:
                    bp_reg = np.zeros(self.population_size, dtype=bool)

                # depending on above values, determine overall limiting of eligibility by applying combination of bp, CVD and bp register settings
                self.register_filter[:,i] = (db_reg + cvd_reg + bp_reg) == False


                # include individuals on CVD and diabetes registers, but not on bp registers
                self.eligible[:, i] = self.register_filter[:,i] * (self.within_age[:, i] > 0) * (self.alive[:, i] == 1)


                # some people get a HC despite being non-eligible
                not_eligible = np.where((self.eligible[:,i]==0) * (self.alive[:, i] == 1))[0]
                n_eligible_other = int(not_eligible.size * self.up_HC_takeup_noneligible)
                self.eligible_other_idx = rnd.sample(not_eligible,n_eligible_other)
                self.eligible[self.eligible_other_idx,i] = 1
                # determine how many are attending it based on the switching
                # probabilities

                r = np.random.random(self.population_size)
                # for those who have already attended a health check within last 4 years, set propensity for attending again to 5% of original value
                rel_uptake = self.RelTurnup.copy()
                red_turnup_factor = 0.01
                if i == 1:
                    # for year 1 in simulation, if HC was last year
                    rel_uptake[self.Attending[:,0] == 1] *= red_turnup_factor
                if i == 2:
                    # for year 2 in simulation, if HC was last year
                    rel_uptake[self.Attending[:,1] == 1] *= red_turnup_factor
                if i == 3:
                    # for year 3 in simulation, if HC was last year
                    rel_uptake[self.Attending[:,2] == 1] *= red_turnup_factor
                if i == 4:
                    # for year 4 in simulation, if HC was last year
                    rel_uptake[self.Attending[:,3] == 1] *= red_turnup_factor
                if i > 4:
                    # for all other years in simulation, if HC was last 4 year
                    rel_uptake[self.Attending[:,i-4:i].sum(axis=1) >= 1] *= red_turnup_factor

                r = np.random.random(self.population_size)

                ever_offered = self.prev_offer_time >= 0
                ## did they attend at their last offer time?
                prev_attended = self.Attending[range(self.population_size) , self.prev_offer_time]
                prev_accepted = ever_offered * prev_attended
                prev_declined = ever_offered * (1 - prev_attended)

                poff = np.zeros((self.population_size)) + self.up_HC_offered  # Prob offered HC can depend on previous attendance
                poff[prev_accepted==1] = self.up_HC_offer_prev_att
                poff[prev_declined==1] = self.up_HC_offer_not_prev_att

                self.poff[:,i] = poff 
                offered = r < (self.alive[:,i] * self.eligible[:, i] * poff)
                self.OfferedHC[offered, i] = 1
                ## update last offer time to include people offered HC this time
                self.prev_offer_time[self.OfferedHC[:,i] == 1] = i

                pup = np.zeros((self.population_size)) + self.up_HC_takeup  # Baseline prob of attendance if offered
                pup[prev_accepted==1] = self.up_HC_takeup_prev_att
                pup[prev_declined==1] = self.up_HC_takeup_not_prev_att
                pup *= rel_uptake
                pup[pup > 1] = 1

                att_prob = self.alive[:,i] * self.eligible[:, i] * poff * pup
                attending = r < att_prob

#                self.pup = pup
#                self.att_prob = att_prob

                self.Attending[attending, i] = 1


                #############################
                # TREATMENT
                #############################


                # (1) Statins
                #np.random.seed(self.randseed)
                # if QRisk is lower than 20, apply rate HC_statins_presc_Q20minus
                statins_q20minus = np.array(self.Attending[:,i] * (self.QRisk[:,i]<20) * (np.random.random(self.population_size) < self.up_HC_statins_presc_Q20minus), dtype=bool)
                # if QRisk is higher than 20, apply rate HC_statins_presc_Q20plus
                statins_q20plus = np.array(self.Attending[:,i] * (self.QRisk[:,i]>=20) * (np.random.random(self.population_size) < self.up_HC_statins_presc_Q20plus), dtype=bool)
                self.Statins_Offered[statins_q20minus,i] = 1
                self.Statins_Offered[statins_q20plus,i] = 1
                # statins compliance

                statins_compliance = np.zeros(self.population_size, dtype=bool)
                r = np.random.random(self.population_size)
                statins_compliance[r < self.up_Statins_comp] = 1

                # those who are prescribed statins and who comply, put them onto statins treatment effect
                self.Statins[self.Statins_Offered[:,i] * statins_compliance,i] = 1
                # put individuals on
                if i<(self.simulation_time-1):
                        # who is put on statins treatment
                        self.on_statins[self.Statins_Offered[:,i] * statins_compliance,(i+1):] = 1
                        # who drops out from statins treatment
                        r = np.random.random(self.population_size)
                        dropout = r < self.up_Statins_dropout_rate
                        self.on_statins[dropout,(i+1):] = 0




                # (2) Weight Reduction

                # ~40% of people with BMI>30 referred to exercise
                r = np.random.random(self.population_size)
                pw_comp = r < self.up_Weight_comp
                r = np.random.random(self.population_size)
                pw_ref = np.array(self.Attending[:,i] * (self.bmi[:,i]>=30) * (r < self.up_HC_weight_ref), dtype=bool)
                self.WeightReduction_Offered[pw_ref,i] = 1

                # compliance (about 50%)
                wred_comp = np.array(self.WeightReduction_Offered[:,i] * pw_comp , dtype=bool)
                self.on_wr[wred_comp,(i+1):] = 1  # labelled as being on weight reduction for ever, though treatment effect diminishes over time
                self.prev_wr_time[wred_comp] = i # store the last time they started weight reduction, to implement weight regain
                self.WeightReduction[wred_comp,i] = 1


                # (3) Smoking Cessation

                # up_HC_smoker_ref % of smokers referred to cessation
                # up_Smoking_eff % of smokers quit at 1 year

                r = np.random.random(self.population_size)
                p_cess = r < self.up_HC_smoker_ref
                cess = np.array(self.Attending[:,i] * (self.q_smoke_cat[:,i]>1) * p_cess, dtype=bool)
                self.SmokingCessation_Offered[cess,i] = 1

                # among those on Cessation program, some quit

                r = np.random.random(self.population_size)

                # Smoking_eff gives *TOTAL* quit rate in health checks scenario,
                # so don't apply quitting for those who already quit irrespective of HC ('background' quitters)

                quitting = cess * (r < (self.up_Smoking_eff - self.background_quit_rate))
                self.SmokingCessation[quitting,i] = 1
                # for checking purposes, include array of quitters in matrix QuitMat
                self.QuitMat[quitting,i] = True
                if i<(self.simulation_time-1):
                    if self.Treatment_SmokingCessation == True:
                        self.q_smoke_cat[quitting,i+1]= 1# ex-smoker
                        # update smokematch: individuals that quit smoking stay quitted for next year
                        #self.smoke2[quitting] = 1



                # (4) Anti-Hypertensives - given to people

                # if QRisk is lower than 20, apply rate HC_aht_presc_Q20minus
                # check who is hypertensive and (sbp > 140) scale prescription rate to those
                i_BP_high = self.q_sbp[:,i] > 140

                i_att = (self.Attending[:,i]==1)
                # how many are hypertensive vs not among attending?
                hc_hypert = (i_att * i_BP_high).sum()
                BP_high_proportion = hc_hypert / float(i_att.sum())
                self.AHT_prescription_scaling[i] = 1./BP_high_proportion

                p_Q20minus = np.random.random(self.population_size) < (self.up_HC_aht_presc_Q20minus * self.AHT_prescription_scaling[i] )
                p_Q20plus = np.random.random(self.population_size) < (self.up_HC_aht_presc_Q20plus * self.AHT_prescription_scaling[i])

                aht_q20minus = np.array(self.Attending[:,i] * (self.QRisk[:,i]<20) * p_Q20minus * i_BP_high, dtype=bool)
                # if QRisk is higher than 20, apply rate HC_aht_presc_Q20plus
                aht_q20plus = np.array(self.Attending[:,i] * (self.QRisk[:,i]>=20) * p_Q20plus * i_BP_high, dtype=bool)

                self.Hypertensives_Offered[aht_q20minus,i] = 1
                self.Hypertensives_Offered[aht_q20plus,i] = 1

                # AHT compliance
                aht_compliance = np.zeros(self.population_size, dtype=bool)
                r = np.random.random(self.population_size)
                aht_compliance[r < self.up_AHT_comp] = 1

                # those who are prescribed AH and who comply, put them onto AHT treatment effect
                self.Hypertensives[self.Hypertensives_Offered[:,i] * aht_compliance,i] = 1

                # also, those who are prescribed AHT are diagnosed for high blood pressure and hence put on bp register
                self.bp[self.Hypertensives[:,i]==1,i] = 1




                if i<(self.simulation_time-1):
                        # who is put on statins treatment
                        self.on_aht[self.Hypertensives_Offered[:,i] * aht_compliance,(i+1):] = 1
                        # who drops out from statins treatment
                        r = np.random.random(self.population_size)
                        dropout = r < self.up_AHT_dropout_rate
                        self.on_aht[dropout,(i+1):] = 0



                # all those who are offered either Statins, AHT, Weight Management, or Smoking Cessation will be categorised as being "Offered Treatment"
                self.OfferedTreatment[self.Statins_Offered[:,i]==1,i] = 1
                self.OfferedTreatment[self.Hypertensives_Offered[:,i]==1,i] = 1
                self.OfferedTreatment[self.WeightReduction_Offered[:,i]==1,i] = 1
                self.OfferedTreatment[self.SmokingCessation_Offered[:,i]==1,i] = 1






            #################################################
            # DISEASES
            #################################################

            if self.RandomSeed['Disease_Incidences'] == True:
                np.random.seed(self.randseed+i)
            else:
                np.random.seed()


            # initialise diseases according to prevalence for all ages

            if i==0:
                self.InitialiseCVDEvents()
                self.Stroke[self.stroke_init,i:] = 1
                self.IHD[self.IHD_init,i:] = 1

                self.CVD[self.stroke_init,i:] = 1
                self.CVD[self.IHD_init,i:] = 1


            # CVD: STROKE AND IHD
            rel_risk, stroke_risk = self.GetCVDRisks(i)
            # draw random numbers - whoever is below the relative risk, develops an event next year,
            # based on QRisk

            cvd_risk = (self.QRisk[:,i]/100.0)*rel_risk # divide by 100 as QRisk is in percentages
            r = np.random.random(self.population_size)
            cvd_event = r<cvd_risk

            # save cvd_risk vector
            self.CVD_risks[:,i] = cvd_risk



            # now for those who develop an event, determine whether it's stroke or IHD
            r = np.random.random(self.population_size)
            stroke_event = r<=stroke_risk
            ihd_event = r>stroke_risk

            # now, fill in Stroke/IHD disease states
             # apply events and diseases only to those still alive
            curr_alive = self.alive[:,i] == 1

            self.Stroke[cvd_event*stroke_event * curr_alive,i:] = 1
            self.StrokeEvents[cvd_event * stroke_event * curr_alive,i] = 1
            self.IHD[cvd_event*ihd_event*curr_alive,i:] = 1
            self.IHDEvents[cvd_event * ihd_event * curr_alive,i] = 1

            # enter CVD events
            self.CVD_events[:,i] = cvd_event * curr_alive

            # enter CVD status
            self.CVD[self.CVD_events[:,i] == 1,i:] = 1


            # if CVD event happened, put diagnosed cases on compm7 register
            r = np.random.random(self.population_size)

            diagnosed_CVD_cases = (self.CVD[:,i] > 0) * (r < self.up_CVD_diagnosis_fraction) * (self.alive[:,i] == 1)
            self.compm7[diagnosed_CVD_cases,i:] = 1



            # DEMENTIA
            # compute CAIDE score
            self.CalculateCAIDE(t=i, output=False)
            rel_risk,CAIDE_risk = self.GetDementiaRisks(i,age_dependent_risk=True)
            # draw random numbers - whoever is below the relative risk, develops an event next year,
            # based on CAIDE risk score
            dem_risk = rel_risk*CAIDE_risk
            r = np.random.random(self.population_size)
            with np.errstate(invalid='ignore'):
                dem_event = r<dem_risk

            self.Dementia[dem_event,i:] = 1
            self.DementiaEvents[dem_event * (self.alive[:,i]==1),i] = 1


            ###################################################
            # DIABETES
            # assign values of diabetes diagnosis here
            diabetes_diagnosed = self.diabetes_delta4 == 1

            # apply diagnose rate self.diabetes_delta4, which is change within 4 years, for single years

            sy = np.mod(i,4)


            diab_quarter_diagnosis_subpopulation = np.zeros(self.population_size,dtype=bool)
            diab_quarter_diagnosis_subpopulation[sy::4] = 1

            if i<(self.simulation_time-1):
                self.diabetes[diabetes_diagnosed * diab_quarter_diagnosis_subpopulation,(i+1):] = 1
                # assume that all newly diagnosed diabetes individuals are diabetes type 2
                self.q_b_type2[diabetes_diagnosed,(i+1):] = 1


            ###################################################
            # BLOOD PRESSURE REGISTER
            # other registered conditions
            # high blood pressure, defined as sbp > 150
            bp_high = self.q_sbp[:,i] > 150
            bp_diagnosed = np.random.random(self.population_size) < self.up_bp_diagnosis_fraction
            self.bp[bp_high * bp_diagnosed,i:] = 1





            # LUNG CANCER
            # evaluate incidences based on lung cancer life table for age, sex, smoking category


            self.LC_incidence = np.zeros(self.population_size)


            LC_male_s_inc = self.LT_lungcancer[age,8]
            LC_female_s_inc = self.LT_lungcancer[age,11]
            LC_male_ns_inc = self.LT_lungcancer[age,15]
            LC_female_ns_inc = self.LT_lungcancer[age,18]

            # as a 'strict' rule, everyone that has ever smoked is considered a smoker
            #nsmoker = self.q_smoke_cat.max(axis=1) <= 1
            #smoker = self.q_smoke_cat.max(axis=1) > 0
            # this would be a more relaxed rule, in which the current smoking status only determines whether
            # one is a smoker or not.
            nsmoker = self.q_smoke_cat[:,i]<=1
            smoker = self.q_smoke_cat[:,i]>1

            male_s = male * smoker
            male_ns = male * nsmoker
            female_s = female * smoker
            female_ns = female * nsmoker

            self.LC_incidence[male_s] = LC_male_s_inc[male_s]
            self.LC_incidence[female_s] = LC_female_s_inc[female_s]
            self.LC_incidence[male_ns] = LC_male_ns_inc[male_ns]
            self.LC_incidence[female_ns] = LC_female_ns_inc[female_ns]


            r = np.random.random(self.population_size)
            LC_event = r<self.LC_incidence
            self.LungCancer[LC_event,i:] = 1
            self.LungCancerEvents[LC_event * (self.alive[:,i]==1),i] = 1




            #############################
            # MORTALITY
            #############################
            #np.random.seed(self.randseed+i*10)
            # those who are 99, automatically die
            if i>0:

                # all those who are 99 and were 98 last year, die this year
                # this is to distinguish them from age staying at 99 for indexing purpose at LT lookups

                thisyear99 = self.age[:,i] == 99
                lastyear98 = self.age[:,i-1] == 98

                exits = (thisyear99 * lastyear98 * (self.alive[:,i]==1)) == 1

                nodeath = self.Death.sum(axis=1) == 0

                self.Death[exits*nodeath,i] = 1
                self.CauseOfDeath[exits*nodeath] = 'Other'
                self.alive[exits*nodeath,i:] = 0

            # death by other causes

            if self.RandomSeed['Mortality_Other'] == False:
                np.random.seed()
            else:
                np.random.seed(i+101)


            other_male = self.LT_othercauses_smoothed[age,1]
            other_female = self.LT_othercauses_smoothed[age,2]
            self.other_mortality = np.zeros(self.population_size)
            self.other_mortality[male] = other_male[male]
            self.other_mortality[female] = other_female[female]
            r = np.random.random(self.population_size)
            OTdies = (r<self.other_mortality)
            self.alive[OTdies,i:] = 0
            nodeath = self.Death.sum(axis=1) == 0
            self.Death[OTdies*nodeath,i] = 1
            self.CauseOfDeath[OTdies*nodeath] = 'Other'


            if self.RandomSeed['Mortality_Diseases'] == False:
                np.random.seed()
            else:
                np.random.seed(i+102)


            # mortality by Stroke

            self.stroke_mortality = np.zeros(self.population_size)
            stroke_male = self.LT_stroke[age,2]
            stroke_female = self.LT_stroke[age,5]

            self.stroke_mortality[male] = stroke_male[male]
            self.stroke_mortality[female] = stroke_female[female]
            self.stroke_mortality *= self.up_cvd_background_cfr_reduction

            # convert rate to probability of death following exponential distribution
            self.stroke_mortality = 1 -np.exp(-self.stroke_mortality)

#            if (self.up_cvd_sudden_death > 0): 
#            p_stroke = cvd_risk*stroke_risk 
#            self.stroke.mortality = (self.stroke_mortality - self.up_cvd_sudden_death * p_stroke) / (1 - p_stroke)
            
            r = np.random.random(self.population_size)
            STdies = (r<self.stroke_mortality) * (self.Stroke[:,i]==1)
            self.alive[STdies,i:] = 0
            nodeath = self.Death.sum(axis=1) == 0
            self.Death[STdies*nodeath,i] = 1
            self.CauseOfDeath[STdies*nodeath] = 'Stroke'



            # mortality by IHD

            self.ihd_mortality = np.zeros(self.population_size)
            ihd_male = self.LT_ihd[age,2]
            ihd_female = self.LT_ihd[age,5]

            self.ihd_mortality[male] = ihd_male[male]
            self.ihd_mortality[female] = ihd_female[female]
            self.ihd_mortality *= self.up_cvd_background_cfr_reduction

            # convert rate to probability of death following exponential distribution
            self.ihd_mortality = 1 -np.exp(-self.ihd_mortality)

#            if (self.up_cvd_sudden_death > 0):
#                p_ihd = cvd_risk*ihd_risk 
#                self.ihd.mortality = (self.ihd_mortality - self.up_cvd_sudden_death * p_ihd) / (1 - p_ihd)

            
            r = np.random.random(self.population_size)
            Idies = (r<self.ihd_mortality) * (self.IHD[:,i]==1) # * (self.CVD_events[:,i] != 1)
            self.alive[Idies,i:] = 0
            nodeath = self.Death.sum(axis=1) == 0
            self.Death[Idies*nodeath,i] = 1
            self.CauseOfDeath[Idies*nodeath] = 'IHD'

            ######
            # mortality by sudden death following a CVD event
            r = np.random.random(self.population_size)
            sudden_death_risk = self.up_cvd_sudden_death
            CVDdeath = (r<sudden_death_risk) * (self.CVD_events[:,i] == 1)
            nodeath = self.Death.sum(axis=1) == 0
            self.Death[CVDdeath*nodeath,i] = 1
            self.alive[CVDdeath*nodeath,i:] = 0
            self.CauseOfDeath[CVDdeath*nodeath*(self.StrokeEvents[:,i] == 1)] = 'Stroke'
            self.CauseOfDeath[CVDdeath*nodeath*(self.IHDEvents[:,i] == 1)] = 'IHD'

            # mortality by Dementia

            self.dem_mortality = np.zeros(self.population_size)
            dem_male = self.LT_dementia[age,2]
            dem_female = self.LT_dementia[age,5]
            self.dem_mortality[male] = dem_male[male]
            self.dem_mortality[female] = dem_female[female]
            # convert CFR to probability of death within next year
            self.dem_mortality = 1 -np.exp(-self.dem_mortality)

            r = np.random.random(self.population_size)
            Ddies = (r<self.dem_mortality) * (self.Dementia[:,i]==1)
            self.alive[Ddies,i:] = 0

            nodeath = self.Death.sum(axis=1) == 0
            self.Death[Ddies*nodeath,i] = 1
            self.CauseOfDeath[Ddies*nodeath] = 'Dementia'



            # mortality by Lung Cancer

            self.LC_mortality = np.zeros(self.population_size)
            LC_male_s_mort = self.LT_lungcancer[age,9]
            LC_female_s_mort = self.LT_lungcancer[age,12]
            LC_male_ns_mort = self.LT_lungcancer[age,16]
            LC_female_ns_mort = self.LT_lungcancer[age,19]

            # as a 'strict' rule, everyone that has ever smoked is considered a smoker
            #nsmoker = self.q_smoke_cat.max(axis=1) <= 1
            #smoker = self.q_smoke_cat.max(axis=1) > 1
            # this would be a more relaxed rule, in which the current smoking status only determines whether
            # one is a smoker or not.
            nsmoker = self.q_smoke_cat[:,i]<=1
            smoker = self.q_smoke_cat[:,i]>1

            male_s = male * smoker
            male_ns = male * nsmoker
            female_s = female * smoker
            female_ns = female * nsmoker

            self.LC_mortality[male_s] = LC_male_s_mort[male_s]
            self.LC_mortality[female_s] = LC_female_s_mort[female_s]
            self.LC_mortality[male_ns] = LC_male_ns_mort[male_ns]
            self.LC_mortality[female_ns] = LC_female_ns_mort[female_ns]

            # convert CFR to probability of death within next year
            self.LC_mortality = 1 - np.exp(-self.LC_mortality)


            r = np.random.random(self.population_size)
            Ldies = (r<self.LC_mortality) * (self.LungCancer[:,i]==1)
            self.alive[Ldies,i:] = 0
            nodeath = self.Death.sum(axis=1) == 0
            self.Death[Ldies*nodeath,i] = 1
            self.CauseOfDeath[Ldies*nodeath] = 'Lung Cancer'











            # delete internal arrays
            del r, female_ns, female_s, male_s, male_ns, smoker,bp_high,bp_diagnosed
            del rel_risk,CAIDE_risk, dem_risk
            del stroke_event,ihd_event, cvd_risk,cvd_event

            try:
                del p_Q20minus, p_Q20plus
                del aht_q20minus, aht_q20plus
            except:
                pass

        ################################################################################
        # SIMULATION COMPUTATIONAL SUMMARY STATISTICS
        # output of initialisation and simulation time

        hrs1,mins1,sec1 = hadd.SplitTime(self.t_init_elapsed)
        if self.verbose == True:
            print('initialisation time %02d:%02d:%02d (hrs:min:sec)' % (hrs1, mins1, sec1) )


        self.t_elapsed = time.time() - t0

        hrs2,mins2,sec2 = hadd.SplitTime(self.t_elapsed)
        if self.verbose == True:
            print('simulation time %02d:%02d:%02d (hrs:min:sec)' % (hrs2, mins2, sec2) )

        hrs3,mins3,sec3 = hadd.SplitTime(self.t_elapsed+self.t_init_elapsed)

        print('total execution time %02d:%02d:%02d (hrs:min:sec)' % (hrs3, mins3, sec3) )

        return self.t_elapsed




    ###########################################################
    # CALCULATING SIMULATION OUTPUTS
    ###########################################################

    def CalculateQALY_CVD(self):
        '''calculates QALY for stroke and IHD

        if people die within simulation, all their life years during simulation are counted
        if people still live with age x at the end of simulation, their prospective age of death is estimated
        using life expectancy at age x from ONS life tables, assuming their disease state at age x remains unchanged for the rest of their life

        '''
        self.QALY_CVD = np.zeros(self.population_size)

#        # load life expectancy at age from life table data
#        LE = np.genfromtxt('data/ons_lifetable_2010-2012.csv',delimiter=',', skip_header=7)
#        LE_male = LE[:,5]
#        LE_female = LE[:,11]

        ihd_weight = 0.63
        stroke_weight = 0.518


        self.qaly_cvd_matrix = np.zeros((self.population_size,self.simulation_time))

        strokes = (self.alive==1) * (self.Stroke==1)
        IHD = (self.alive==1) * (self.IHD==1)
        healthy = (self.alive==1) * (self.IHD==0) * (self.Stroke==0)

        self.qaly_cvd_matrix[IHD] = ihd_weight
        self.qaly_cvd_matrix[strokes] = stroke_weight
        self.qaly_cvd_matrix[healthy] = 1

        self.QALY_CVD = self.qaly_cvd_matrix.sum(axis=1)


    def CalculateQALY(self):
        '''calculates QALY for all diseases: stroke, IHD, Lung Cancer, Dementia

        '''
        self.QALY = np.zeros(self.population_size)

#        # load life expectancy at age from life table data
#        LE = np.genfromtxt('data/ons_lifetable_2010-2012.csv',delimiter=',', skip_header=7)
#        LE_male = LE[:,5]
#        LE_female = LE[:,11]
        ihd_weight = 0.63
        lc_weight = 0.561
        stroke_weight = 0.518
        dem_weight = 0.4

        self.qaly_matrix = np.zeros((self.population_size,self.simulation_time))

        strokes = (self.alive==1) * (self.Stroke==1)
        IHD = (self.alive==1) * (self.IHD==1)
        lc = (self.alive == 1) * (self.LungCancer==1)
        dem = (self.alive==1) * (self.Dementia==1)
        healthy = (self.alive==1) * (self.IHD==0) * (self.Stroke==0) * (self.LungCancer==0) * (self.Dementia==0)

        self.qaly_matrix[healthy] = 1
        self.qaly_matrix[IHD] = ihd_weight
        self.qaly_matrix[lc] = lc_weight
        self.qaly_matrix[strokes] = stroke_weight
        self.qaly_matrix[dem] = dem_weight


        self.QALY = self.qaly_matrix.sum(axis=1)


    def CalculateMeanLifeExpectancy(self):
        '''calculates mean life expectancy for males and females'''


        male = (self.gender == 1)[:,np.newaxis]
        female = (self.gender == 0)[:,np.newaxis]

        self.life_expectancy_male = self.age[self.Death*male].mean()
        self.life_expectancy_female = self.age[self.Death*female].mean()





    def CalculateSmokers_baseline(self, uptoyear=1):
        '''evaluate who smoked and who quitted in first  years of simulation
        up to year specified (and then take mean to get average for year


        compares smoking status in time point x with situation with time point y
            assuming that x<=y
            we want to know:

            at time point x:
                -(a) does someone quit?
            from comparing time point y back to x:
                -(b) is someone still smoking at y when s/he was smoking at x?
                -(c) is someone still non-smoker at y when s/he was nonsmoking at x?
                -(d) is someone still ex-smoker at y when s/he was exsmoking at x?
                -(e) did someone quit between x and y?
                -(f) did someone relapse between x and y?

            quitting means negative number (change to exsmoking), so -1
            relapsing means positive number (change to smoking), so +1
            if still smoking/non-smoking/ex-smoker (no change), we have 0

            Thus, with information from smoking status 'self.q_smoke_cat' at time point x,
            this gives information on what changed to time point y

            returns index in a x/y matrix of size (simulation_size,simulation_size)
            '''


        SM = np.zeros((self.population_size,self.simulation_time*self.simulation_time),dtype=int)

        for i in range(self.simulation_time):
            for j in range(i,self.simulation_time):

                si = self.q_smoke_cat[:,i]
                sj = self.q_smoke_cat[:,j]

                # no change::,
                nochange = si==sj

                # quitting:
                quitting = (si>1) * (sj==1)

                # relapsing:
                relapsing = (si==1) * (sj>1)

                # calculating index from i and j
                idx = i * self.simulation_time + j

                SM[nochange, idx] = 0
                SM[quitting, idx] = -1
                SM[relapsing, idx] = 1

        # now compute numbers of people and relevant measures for years of interest

        # group a: smokers attending and not-attending overall
        a1 = 0
        a2 = 0
        a3 = 0
        a4 = 0
        a5 = 0
        a6 = 0

        # group b: smokers, attending and referred
        b1 = 0
        b2 = 0
        b3 = 0
        b4 = 0
        b5 = 0

        # group c: smokers, attending and not-referred
        c1 = 0
        c2 = 0
        c3 = 0
        c4 = 0
        c5 = 0

        # group d: smokers, not-attending
        d1 = 0
        d2 = 0
        d3 = 0
        d4 = 0

        # initiate lists containing QRisk for groups a at year 0, 1, and 5
        Q_intervals = [0,1,5]

        aQ2 = [[],[],[]]
        aQ3 = [[],[],[]]
        aQ5 = [[],[],[]]
        aQ6 = [[],[],[]]

        bQ1 = [[],[],[]]
        bQ2 = [[],[],[]]
        bQ3 = [[],[],[]]
        bQ4 = [[],[],[]]
        bQ5 = [[],[],[]]

        cQ1 = [[],[],[]]
        cQ2 = [[],[],[]]
        cQ3 = [[],[],[]]
        cQ4 = [[],[],[]]
        cQ5 = [[],[],[]]

        dQ1 = [[],[],[]]
        dQ2 = [[],[],[]]
        dQ3 = [[],[],[]]
        dQ4 = [[],[],[]]

        for j in range(0,uptoyear):
            # how many are going to HC in year 1?
            att = self.Attending[:,j] == 1
            nonatt = self.Attending[:,j] == 0

            # of which smokers and nonsmokers
            att_s = (self.q_smoke_cat[:,j] > 1) * att
            att_ns = (self.q_smoke_cat[:,j] <= 1) * att

            nonatt_s = (self.q_smoke_cat[:,j] > 1) * nonatt
            nonatt_ns = (self.q_smoke_cat[:,j] <= 1) * nonatt

            # group a: smokers attending and not-attending overall
            a1 += att.sum()
            a2 += att_s.sum()
            a3 += att_ns.sum()
            a4 += nonatt.sum()
            a5 += nonatt_s.sum()
            a6 += nonatt_ns.sum()

            # compute QRISKs for year 0, 1 and 3 for the 'a' group:

            for xidx,x in enumerate(Q_intervals):
                aQ2[xidx].extend(self.QRisk[att_s, j+x].tolist())
                aQ3[xidx].extend(self.QRisk[att_ns, j+x].tolist())
                aQ5[xidx].extend(self.QRisk[nonatt_s, j+x].tolist())
                aQ6[xidx].extend(self.QRisk[nonatt_ns, j+x].tolist())

            # get indices from SM matrix for comparison with 1st / 5th year afterwards
            idx_1yr = j * self.simulation_time + j+1
            idx_5yr = j * self.simulation_time + j+5
            idx_relapse5yr = (j+1) * self.simulation_time + j+5

            quitting_after_1yr = SM[:,idx_1yr] == -1
            quitting_after_5yr = SM[:,idx_5yr] == -1
            still_smoking_after_5yr = (SM[:,idx_5yr] == 0) * (self.q_smoke_cat[:,j]>1)
            relapsed_within_5yr = SM[:,idx_relapse5yr] == 1

            # group b: smokers, attending and referred
            referred = self.SmokingCessation_Offered[:,j] == 1

            b1 += referred.sum()
            b2 += (referred * quitting_after_1yr).sum()
            b3 += (referred * still_smoking_after_5yr).sum()
            b4 += (referred * relapsed_within_5yr).sum()
            b5 += (referred * quitting_after_5yr).sum()

            for xidx,x in enumerate(Q_intervals):
                bQ1[xidx].extend(self.QRisk[referred, j+x].tolist())
                bQ2[xidx].extend(self.QRisk[(referred * quitting_after_1yr), j+x].tolist())
                bQ3[xidx].extend(self.QRisk[(referred * still_smoking_after_5yr), j+x].tolist())
                bQ4[xidx].extend(self.QRisk[(referred * relapsed_within_5yr), j+x].tolist())
                bQ5[xidx].extend(self.QRisk[(referred * quitting_after_5yr), j+x].tolist())


            # group c: smokers, attending and not referred
            not_referred = (self.SmokingCessation_Offered[:,j] == 0) * att_s

            c1 += not_referred.sum()
            c2 += (not_referred * quitting_after_1yr).sum()
            c3 += (not_referred * still_smoking_after_5yr).sum()
            c4 += (not_referred * relapsed_within_5yr).sum()
            c5 += (not_referred * quitting_after_5yr).sum()

            for xidx,x in enumerate(Q_intervals):
                cQ1[xidx].extend(self.QRisk[not_referred, j+x].tolist())
                cQ2[xidx].extend(self.QRisk[(not_referred * quitting_after_1yr), j+x].tolist())
                cQ3[xidx].extend(self.QRisk[(not_referred * still_smoking_after_5yr), j+x].tolist())
                cQ4[xidx].extend(self.QRisk[(not_referred * relapsed_within_5yr), j+x].tolist())
                cQ5[xidx].extend(self.QRisk[(not_referred * quitting_after_5yr), j+x].tolist())

            # group d: smokers, not attending
            d1 += (nonatt_s * quitting_after_1yr).sum()
            d2 += (nonatt_s * still_smoking_after_5yr).sum()
            d3 += (nonatt_s * relapsed_within_5yr).sum()
            d4 += (nonatt_s * quitting_after_5yr).sum()

            for xidx,x in enumerate(Q_intervals):
                dQ1[xidx].extend(self.QRisk[(nonatt_s * quitting_after_1yr), j+x].tolist())
                dQ2[xidx].extend(self.QRisk[(nonatt_s * still_smoking_after_5yr), j+x].tolist())
                dQ3[xidx].extend(self.QRisk[(nonatt_s * relapsed_within_5yr), j+x].tolist())
                dQ4[xidx].extend(self.QRisk[(nonatt_s * quitting_after_5yr), j+x].tolist())





        Smoker_Stats={}
        Smoker_Stats['a.0 Years covered'] = uptoyear
        Smoker_Stats['a.1 Number of attending individuals in one year'] = a1 / float(uptoyear)
        Smoker_Stats['a.2 Attending smokers'] = a2 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['a.2.QRISK%d Attending smokers' % Q_intervals[r]] = np.mean(aQ2[r])
        Smoker_Stats['a.3 Attending non-smokers'] = a3 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['a.3.QRISK%d Attending non-smokers' % Q_intervals[r]] = np.mean(aQ3[r])
        Smoker_Stats['a.4 Number of nonattending individuals in one year'] = a4 / float(uptoyear)
        Smoker_Stats['a.5 Non-Attending smokers'] = a5 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['a.5.QRISK%d Non-Attending smokers' % Q_intervals[r]] = np.mean(aQ5[r])
        Smoker_Stats['a.6 Non-Attending non-smokers'] = a6 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['a.6.QRISK%d Non-Attending non-smokers' % Q_intervals[r]] = np.mean(aQ6[r])


        Smoker_Stats['b.1 Referred to SC'] = b1 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['b.1.QRISK%d Referred to SC' % Q_intervals[r]] = np.mean(bQ1[r])

        Smoker_Stats['b.2 Referred to SC - quitting after 1 yr'] = b2 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['b.2.QRISK%d Referred to SC - quitting after 1 yr' % Q_intervals[r]] = np.mean(bQ2[r])

        Smoker_Stats['b.3 Referred to SC - still smoking after 5 yrs'] = b3 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['b.3.QRISK%d Referred to SC - still smoking after 5 yrs' % Q_intervals[r]] = np.mean(bQ3[r])

        Smoker_Stats['b.4 Referred to SC - relapse within 5 yrs'] = b4 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['b.4.QRISK%d Referred to SC - relapse within 5 yrs' % Q_intervals[r]] = np.mean(bQ4[r])

        Smoker_Stats['b.5 Referred to SC - ex-smokers after 5 yrs'] = b5 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['b.5.QRISK%d Referred to SC - ex-smokers after 5 yrs' % Q_intervals[r]] = np.mean(bQ5[r])


        Smoker_Stats['c.1 Not referred to SC'] = c1 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['c.1.QRISK%d Not referred to SC' % Q_intervals[r]] = np.mean(cQ1[r])

        Smoker_Stats['c.2 Not referred to SC - quitting after 1 yr'] = c2 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['c.2.QRISK%d Not referred to SC - quitting after 1 yr' % Q_intervals[r]] = np.mean(cQ2[r])

        Smoker_Stats['c.3 Not referred to SC - still smoking after 5 yrs'] = c3 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['c.3.QRISK%d Not referred to SC - still smoking after 5 yrs' % Q_intervals[r]] = np.mean(cQ3[r])

        Smoker_Stats['c.4 Not referred to SC - relapse within 5 yrs'] = c4 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['c.4.QRISK%d Not referred to SC - relapse within 5 yrs' % Q_intervals[r]] = np.mean(cQ4[r])

        Smoker_Stats['c.5 Not referred to SC - ex-smokers after 5 yrs'] = c5 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['c.5.QRISK%d Not referred to SC - ex-smokers after 5 yrs' % Q_intervals[r]] = np.mean(cQ5[r])


        Smoker_Stats['d.1 Not attending - quitting after 1 yr'] = d1 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['d.1.QRISK%d Not attending - quitting after 1 yr' % Q_intervals[r]] = np.mean(dQ1[r])

        Smoker_Stats['d.2 Not attending - still smoking after 5 yrs'] = d2 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['d.2.QRISK%d Not attending - still smoking after 5 yrs' % Q_intervals[r]] = np.mean(dQ2[r])

        Smoker_Stats['d.3 Not attending - relapse within 5 yrs'] = d3 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['d.3.QRISK%d Not attending - relapse within 5 yrs' % Q_intervals[r]] = np.mean(dQ3[r])

        Smoker_Stats['d.4 Not attending - ex-smokers after 5 yrs'] = d4 / float(uptoyear)
        for r in range(3):
            Smoker_Stats['d.4.QRISK%d Not attending - ex-smokers after 5 yrs' % Q_intervals[r]] = np.mean(dQ4[r])


        self.SStats_y1 = Smoker_Stats






    def CalculateSmokers_baseline_old(self,uptoyear=1):
        '''evaluate who smoked and who quitted in first  years of simulation
        up to year specified (and then take mean to get average for year'''

        a1 = 0
        a2 = 0
        a3 = 0


        b1 = 0
        b2 = 0
        b3 = 0
        c1 = 0
        c2 = 0
        c3 = 0
        c4 = 0
        d1 = 0
        d2 = 0
        d3 = 0
        d4 = 0
        e1 = 0
        e2 = 0
        e3 = 0
        e4 = 0
        e5 = 0

        q1 = 0
        q2 = 0
        q3 = 0
        q4 = 0
        q5 = 0
        q6 = 0
        q7 = 0
        q8 = 0
        q9 = 0
        q10 = 0
        q11 = 0
        q12 = 0


        for j in range(0,uptoyear):
            # how many are going to HC in year 1?
            att = self.Attending[:,j] == 1
            nonatt = self.Attending[:,j] == 0

            # of which smokers and nonsmokers
            att_s = (self.q_smoke_cat[:,j] > 1) * att
            att_ns = (self.q_smoke_cat[:,j] <= 1) * att

            nonatt_s = (self.q_smoke_cat[:,j] > 1) * nonatt
            nonatt_ns = (self.q_smoke_cat[:,j] <= 1) * nonatt


            s_offered = 0
            s_offered_quit = 0
            s_offered_notquit = 0
            s_offered_quit5yrs = 0
            s_offered_relapse = 0

            q_offered = 0
            q_offered_quit = 0
            q_offered_quit5yrs = 0
            q_offered_relapse = 0

            s_notoffered = 0
            s_notoffered_quit = 0
            s_notoffered_notquit = 0
            s_notoffered_quit5yrs = 0
            s_notoffered_relapse = 0

            q_notoffered = 0
            q_notoffered_quit = 0
            q_notoffered_quit5yrs = 0
            q_notoffered_relapse = 0

            s_notattending = 0
            s_notattending_quit = 0
            s_notattending_notquit = 0
            s_notattending_quit5yrs = 0
            s_notattending_relapse = 0

            q_notattending = 0
            q_notattending_quit = 0
            q_notattending_quit5yrs = 0
            q_notattending_relapse = 0

            # array containing number of years to be devided in average
            uptoyear_array = np.ones(12) * uptoyear
            ns_notattending = 0


            # of attending smokers, how many quitted and not, how many relapsed within 5 years
            for i in range(self.population_size):
                if self.Attending[i,j] == 1:
                    #  attending

                    if self.SmokingCessation_Offered[i,j] == True:
                        s_offered += 1
                        q_offered += self.QRisk[i,j]
                        if self.q_smoke_cat[i,j+1] == 1:
                            s_offered_quit += 1
                            q_offered_quit += self.QRisk[i,j+1]
                            if self.q_smoke_cat[i,j+5] > 1:
                                s_offered_relapse += 1
                                q_offered_relapse += self.QRisk[i,j+5]

                        else:
                            s_offered_notquit += 1
                        if self.q_smoke_cat[i,j+5] == 1:
                            s_offered_quit5yrs += 1
                            q_offered_quit5yrs += self.QRisk[i,j+5]

                    else:
                        if self.q_smoke_cat[i,j] > 1:
                            s_notoffered += 1
                            q_notoffered += self.QRisk[i,j]
                            if self.q_smoke_cat[i,j+1] == 1:
                                s_notoffered_quit += 1
                                q_notoffered_quit += self.QRisk[i,j+1]
                                if self.q_smoke_cat[i,j+5] > 1:
                                    s_notoffered_relapse += 1
                                    q_notoffered_relapse += self.QRisk[i,j+5]
                            else:
                                s_notoffered_notquit += 1
                            if self.q_smoke_cat[i,j+5] == 1:
                                s_notoffered_quit5yrs += 1
                                q_notoffered_quit5yrs += self.QRisk[i,j+5]



                else:
                    # if currently smoker:
                    if self.q_smoke_cat[i,j] > 1:
                    # nonattending smokers
                        s_notattending += 1
                        q_notattending += self.QRisk[i,j]
                        if self.q_smoke_cat[i,j+1] == 1:
                            s_notattending_quit += 1
                            q_notattending_quit += self.QRisk[i,j+1]

                            if self.q_smoke_cat[i,j+5] > 1:
                                s_notattending_relapse += 1
                                q_notattending_relapse += self.QRisk[i,j+5]
                        else:
                            s_notattending_notquit += 1

                        if self.q_smoke_cat[i,j+5] == 1:
                            s_notattending_quit5yrs += 1
                            q_notattending_quit5yrs += self.QRisk[i,j+5]

                    else:
                        # nonattending nonsmokers:
                        ns_notattending += 1



            a1 += att.sum()
            a2 += att_s.sum()
            a3 += att_ns.sum()

            b1 += nonatt.sum()
            b2 += nonatt_s.sum()
            b3 += nonatt_ns.sum()


            c1 += s_offered
            c2 += s_offered_quit
            c3 += s_offered_quit5yrs
            c4 += s_offered_relapse

            d1 += s_notoffered
            d2 += s_notoffered_quit
            d3 += s_notoffered_quit5yrs
            d4 += s_notoffered_relapse

            e1 += s_notattending
            e2 += s_notattending_quit
            e3 += s_notattending_quit5yrs
            e4 += s_notattending_relapse
            e5 += ns_notattending


            # now add all QRisks.
            # if any of the qroups has zero elements, then reduce amount of years to be divided later by 1

            try:
                q1 += q_offered / float(s_offered)
            except:
                uptoyear_array[0] -= 1
            try:
                q2 += q_offered_quit / float(s_offered_quit)
            except:
                uptoyear_array[1] -= 1
            try:
                q3 += q_offered_quit5yrs / float(s_offered_quit5yrs)
            except:
                uptoyear_array[2] -= 1
            try:
                q4 += q_offered_relapse / float(s_offered_relapse)
            except:
                uptoyear_array[3] -= 1

            try:
                q5 += q_notoffered / float(s_notoffered)
            except:
                uptoyear_array[4] -= 1
            try:
                q6 += q_notoffered_quit / float(s_notoffered_quit)
            except:
                uptoyear_array[5] -= 1
            try:
                q7 += q_notoffered_quit5yrs / float(s_notoffered_quit5yrs)
            except:
                uptoyear_array[6] -= 1
            try:
                q8 += q_notoffered_relapse / float(s_notoffered_relapse)
            except:
                uptoyear_array[7] -= 1

            try:
                q9 += q_notattending / float(s_notattending)
            except:
                uptoyear_array[8] -= 1
            try:
                q10 += q_notattending_quit / float(s_notattending_quit)
            except:
                uptoyear_array[9] -= 1
            try:
                q11 += q_notattending_quit5yrs / float(s_notattending_quit5yrs)
            except:
                uptoyear_array[10] -= 1
            try:
                q12 += q_notattending_relapse / float(s_notattending_relapse)
            except:
                uptoyear_array[11] -= 1

        self.uptoyear_array = uptoyear_array

        Smoker_Stats={}
        Smoker_Stats['a.0 Years covered'] = uptoyear
        Smoker_Stats['a.1 Number of attending individuals in one year'] = a1 / float(uptoyear)
        Smoker_Stats['a.2 Attending smokers'] = a2 / float(uptoyear)
        Smoker_Stats['a.3 Attending non-smokers'] = a3 / float(uptoyear)

        Smoker_Stats['b.1 Number of nonattending individuals in one year'] = b1 / float(uptoyear)
        Smoker_Stats['b.2 Non-Attending smokers'] = b2 / float(uptoyear)
        Smoker_Stats['b.3 Non-Attending non-smokers'] = b3 / float(uptoyear)

        Smoker_Stats['c.1 Smokers offered CS'] = c1 / float(uptoyear)
        Smoker_Stats['c.2 Smokers offered CS, quitting, and ex-smokers after 1 year'] = c2 / float(uptoyear)
        Smoker_Stats['c.3 Smokers offered CS, quitting, and ex-smokers after 5 years'] = c3 / float(uptoyear)
        Smoker_Stats['c.4 Smokers offered CS, quitting, and relapsing'] = c4 / float(uptoyear)
        Smoker_Stats['cQ.1 mean QRISK of smokers offered CS'] = q1 / np.float64(uptoyear_array[0])
        Smoker_Stats['cQ.2 mean QRISK of smokers offered CS, quitting, and ex-smokers after 1 year'] = q2 / np.float64(uptoyear_array[1])
        Smoker_Stats['cQ.3 mean QRISK of smokers offered CS, quitting, and ex-smokers after 5 years'] = q3 / np.float64(uptoyear_array[2])
        Smoker_Stats['cQ.4 mean QRISK of smokers offered CS, quitting, and relapsing'] = q4 / np.float64(uptoyear_array[3])



        Smoker_Stats['d.1 Smokers not offered CS'] = d1 / float(uptoyear)
        Smoker_Stats['d.2 Smokers not offered CS, quitting, and ex-smokers after 1 year'] = d2 / float(uptoyear)
        Smoker_Stats['d.3 Smokers not offered CS, quitting, and ex-smokers after 5 years'] = d3  / float(uptoyear)
        Smoker_Stats['d.4 Smokers not offered CS, quitting, and relapsing'] = d4 / float(uptoyear)
        Smoker_Stats['dQ.1 mean QRISK of smokers not offered CS'] = q5 / np.float64(uptoyear_array[4])
        Smoker_Stats['dQ.2 mean QRISK of smokers not offered CS, quitting, and ex-smokers after 1 year'] = q6 / np.float64(uptoyear_array[5])
        Smoker_Stats['dQ.3 mean QRISK of smokers not offered CS, quitting, and ex-smokers after 5 years'] = q7 / np.float64(uptoyear_array[6])
        Smoker_Stats['dQ.4 mean QRISK of smokers not offered CS, quitting, and relapsing'] = q8 / np.float64(uptoyear_array[7])


        Smoker_Stats['e.1 Smokers not attending HC'] = e1   / float(uptoyear)
        Smoker_Stats['e.2 Smokers not attending HC, quitting, and ex-smokers after 1 year'] = e2 / float(uptoyear)
        Smoker_Stats['e.3 Smokers not attending HC, quitting, and ex-smokers after 5 years'] = e3  / float(uptoyear)
        Smoker_Stats['e.4 Smokers not attending HC, quitting, and relapsing'] = e4 / float(uptoyear)
        Smoker_Stats['eQ.1 mean QRISK of smokers not attending HC'] = q9   / np.float64(uptoyear_array[8])
        Smoker_Stats['eQ.2 mean QRISK of smokers not attending HC, quitting, and ex-smokers after 1 year'] = q10 / np.float64(uptoyear_array[9])
        Smoker_Stats['eQ.3 mean QRISK of smokers not attending HC, quitting, and ex-smokers after 5 years'] = q11  / np.float64(uptoyear_array[10])
        Smoker_Stats['eQ.4 mean QRISK of smokers not attending HC, quitting, and relapsing'] = q12 / np.float64(uptoyear_array[11])

        Smoker_Stats['e.5 Non-Smokers not attending HC'] = e5 / float(uptoyear)

        self.SStats_y1 = Smoker_Stats



    def CalculateBPStats(self,uptoyear=1):
        '''calculate statistics for Blood pressure for intenal flow diagrams
        for the first uptoyear years'''


        a1 = 0
        a2 = 0

        b1 = 0
        b2 = 0
        b3 = 0
        b4 = 0
        b5 = 0
        b6 = 0
        b7 = 0
        b8 = 0

        c1 = 0
        c2 = 0
        c3 = 0

        d1 = 0
        d2 = 0
        d3 = 0
        d4 = 0
        d5 = 0


        for j in range(0,uptoyear):
            # how many are going to HC in year 1?
            att = self.Attending[:,j] == 1
            nonatt = self.Attending[:,j] == 0

            # how many of those are above/below QRISK 20?
            # and what's their mean BP?

            att_q20plus = (self.QRisk[:,j] >=20) * att
            att_q20 = (self.QRisk[:,j] < 20) * att
            nonatt_q20plus = (self.QRisk[:,j] >=20) * nonatt
            nonatt_q20 = (self.QRisk[:,j] < 20) * nonatt

            att_q20plus_sbp = self.q_sbp[att_q20plus,j].mean()
            att_q20_sbp = self.q_sbp[att_q20,j].mean()
            nonatt_q20plus_sbp = self.q_sbp[nonatt_q20plus,j].mean()
            nonatt_q20_sbp = self.q_sbp[nonatt_q20,j].mean()

            # who gets treatment offered, who takes it up?
            aht_offered = (self.Hypertensives_Offered[:,j]>0) * att_q20plus
            aht_notoffered = (self.Hypertensives_Offered[:,j]==0) * att_q20plus
            aht_takeup = (self.Hypertensives[:,j]>0) * att_q20plus



            #  mean change in AHT over next 4 years
            sbp_change_aht_takeup = (self.q_sbp[aht_takeup,j+4] - self.q_sbp[aht_takeup,j]).mean()
            sbp_change_aht_notoffered = (self.q_sbp[aht_notoffered,j+4] - self.q_sbp[aht_notoffered,j]).mean()
            sbp_change_att20 = (self.q_sbp[att_q20,j+4] - self.q_sbp[att_q20,j]).mean()
            sbp_change_nonatt20plus = (self.q_sbp[nonatt_q20plus,j+4] - self.q_sbp[nonatt_q20plus,j]).mean()
            sbp_change_nonatt20 = (self.q_sbp[nonatt_q20,j+4] - self.q_sbp[nonatt_q20,j]).mean()



            a1 += att.sum()
            a2 += nonatt.sum()

            b1 += att_q20plus.sum()
            b2 += att_q20plus_sbp.sum()
            b3 += att_q20.sum()
            b4 += att_q20_sbp.sum()
            b5 += nonatt_q20plus.sum()
            b6 += nonatt_q20plus_sbp.sum()
            b7 += nonatt_q20.sum()
            b8 += nonatt_q20_sbp.sum()

            c1 += aht_offered.sum()
            c2 += aht_takeup.sum()
            c3 += aht_notoffered.sum()

            d1 += sbp_change_aht_takeup.sum()
            d2 += sbp_change_aht_notoffered.sum()
            d3 += sbp_change_att20.sum()
            d4 += sbp_change_nonatt20plus.sum()
            d5 += sbp_change_nonatt20.sum()


        BP_Stats={}
        BP_Stats['a.0 Years covered'] = uptoyear
        BP_Stats['a.1 Number of attending individuals in one year'] = a1 / float(uptoyear)
        BP_Stats['a.2 Number of nonattending individuals in one year'] = a2 / float(uptoyear)

        BP_Stats['b.1 Number of attending individuals with QRisk >20'] = b1 / float(uptoyear)
        BP_Stats['b.2 their mean BP'] = b2 / float(uptoyear)
        BP_Stats['b.3 Number of attending individuals with QRisk <20'] = b3 / float(uptoyear)
        BP_Stats['b.4 their mean BP'] = b4 / float(uptoyear)
        BP_Stats['b.5 Number of nonattending individuals with QRisk >20'] = b5 / float(uptoyear)
        BP_Stats['b.6 their mean BP'] = b6 / float(uptoyear)
        BP_Stats['b.7 Number of nonattending individuals with QRisk <20'] = b7 / float(uptoyear)
        BP_Stats['b.8 their mean BP'] = b8 / float(uptoyear)

        BP_Stats['c.1 AHT offered'] = c1 / float(uptoyear)
        BP_Stats['c.2 AHT takeup'] = c2 / float(uptoyear)
        BP_Stats['c.3 AHT notoffered'] = c3 / float(uptoyear)


        BP_Stats['d.1 mean BP change for AHT takeup group'] = d1 / float(uptoyear)
        BP_Stats['d.2 mean BP change for AHT notoffered group'] = d2 / float(uptoyear)
        BP_Stats['d.3 mean BP change for attending group with QRISK <20'] = d3 / float(uptoyear)
        BP_Stats['d.4 mean BP change for nonattending group with QRISK >20'] = d4 / float(uptoyear)
        BP_Stats['d.5 mean BP change for nonattending group with QRISK <20'] = d5 / float(uptoyear)



        self.BPstats = BP_Stats








    def CalculateBMIStats(self,uptoyear=1):
        '''calculate statistics for Blood pressure for intenal flow diagrams
        for the first uptoyear years'''


        a1 = 0
        a2 = 0

        b1 = 0
        b2 = 0
        b3 = 0
        b4 = 0
        b5 = 0
        b6 = 0
        b7 = 0
        b8 = 0

        c1 = 0
        c2 = 0
        c3 = 0

        d1 = 0
        d2 = 0
        d3 = 0
        d4 = 0
        d5 = 0


        for j in range(0,uptoyear):
            # how many are going to HC in year 1?
            att = self.Attending[:,j] == 1
            nonatt = self.Attending[:,j] == 0

            # how many of those are above/below QRISK 20?
            # and what's their mean BP?

            att_q20plus = (self.QRisk[:,j] >=20) * att
            att_q20 = (self.QRisk[:,j] < 20) * att
            nonatt_q20plus = (self.QRisk[:,j] >=20) * nonatt
            nonatt_q20 = (self.QRisk[:,j] < 20) * nonatt

            att_q20plus_bmi = self.bmi[att_q20plus,j].mean()
            att_q20_bmi = self.bmi[att_q20,j].mean()
            nonatt_q20plus_bmi = self.bmi[nonatt_q20plus,j].mean()
            nonatt_q20_bmi = self.bmi[nonatt_q20,j].mean()

            # who gets treatment offered, who takes it up?
            WR_offered = (self.WeightReduction_Offered[:,j]>0) * att_q20plus
            WR_notoffered = (self.WeightReduction_Offered[:,j]==0) * att_q20plus
            WR_takeup = (self.WeightReduction[:,j]>0) * att_q20plus



            #  mean change in AHT over next 4 years
            bmi_change_WR_takeup = (self.bmi[WR_takeup,j+4] - self.bmi[WR_takeup,j]).mean()
            bmi_change_WR_notoffered = (self.bmi[WR_notoffered,j+4] - self.bmi[WR_notoffered,j]).mean()
            bmi_change_att20 = (self.bmi[att_q20,j+4] - self.bmi[att_q20,j]).mean()
            bmi_change_nonatt20plus = (self.bmi[nonatt_q20plus,j+4] - self.bmi[nonatt_q20plus,j]).mean()
            bmi_change_nonatt20 = (self.bmi[nonatt_q20,j+4] - self.bmi[nonatt_q20,j]).mean()



            a1 += att.sum()
            a2 += nonatt.sum()

            b1 += att_q20plus.sum()
            b2 += att_q20plus_bmi.sum()
            b3 += att_q20.sum()
            b4 += att_q20_bmi.sum()
            b5 += nonatt_q20plus.sum()
            b6 += nonatt_q20plus_bmi.sum()
            b7 += nonatt_q20.sum()
            b8 += nonatt_q20_bmi.sum()

            c1 += WR_offered.sum()
            c2 += WR_takeup.sum()
            c3 += WR_notoffered.sum()

            d1 += bmi_change_WR_takeup.sum()
            d2 += bmi_change_WR_notoffered.sum()
            d3 += bmi_change_att20.sum()
            d4 += bmi_change_nonatt20plus.sum()
            d5 += bmi_change_nonatt20.sum()


        BMI_Stats={}
        BMI_Stats['a.0 Years covered'] = uptoyear
        BMI_Stats['a.1 Number of attending individuals in one year'] = a1 / float(uptoyear)
        BMI_Stats['a.2 Number of nonattending individuals in one year'] = a2 / float(uptoyear)

        BMI_Stats['b.1 Number of attending individuals with QRisk >20'] = b1 / float(uptoyear)
        BMI_Stats['b.2 their mean BMI level'] = b2 / float(uptoyear)
        BMI_Stats['b.3 Number of attending individuals with QRisk <20'] = b3 / float(uptoyear)
        BMI_Stats['b.4 their mean BMI level'] = b4 / float(uptoyear)
        BMI_Stats['b.5 Number of nonattending individuals with QRisk >20'] = b5 / float(uptoyear)
        BMI_Stats['b.6 their mean BMI level'] = b6 / float(uptoyear)
        BMI_Stats['b.7 Number of nonattending individuals with QRisk <20'] = b7 / float(uptoyear)
        BMI_Stats['b.8 their mean BMI level'] = b8 / float(uptoyear)

        BMI_Stats['c.1 Weight Reduction offered'] = c1 / float(uptoyear)
        BMI_Stats['c.2 Weight Reduction takeup'] = c2 / float(uptoyear)
        BMI_Stats['c.3 Weight Reduction notoffered'] = c3 / float(uptoyear)


        BMI_Stats['d.1 mean BMI change for Weight Reduction takeup group'] = d1 / float(uptoyear)
        BMI_Stats['d.2 mean BMI change for Weight Reduction notoffered group'] = d2 / float(uptoyear)
        BMI_Stats['d.3 mean BMI change for attending group with QRISK <20'] = d3 / float(uptoyear)
        BMI_Stats['d.4 mean BMI change for nonattending group with QRISK >20'] = d4 / float(uptoyear)
        BMI_Stats['d.5 mean BMI change for nonattending group with QRISK <20'] = d5 / float(uptoyear)



        self.BMIstats = BMI_Stats


    def CalculateCholStats(self,uptoyear=1):
        '''calculate statistics for Blood pressure for intenal flow diagrams
        for the first uptoyear years'''


        a1 = 0
        a2 = 0

        b1 = 0
        b2 = 0
        b3 = 0
        b4 = 0
        b5 = 0
        b6 = 0
        b7 = 0
        b8 = 0

        c1 = 0
        c2 = 0
        c3 = 0

        d1 = 0
        d2 = 0
        d3 = 0
        d4 = 0
        d5 = 0


        for j in range(0,uptoyear):
            # how many are going to HC in year 1?
            att = self.Attending[:,j] == 1
            nonatt = self.Attending[:,j] == 0

            # how many of those are above/below QRISK 20?
            # and what's their mean BP?

            att_q20plus = (self.QRisk[:,j] >=20) * att
            att_q20 = (self.QRisk[:,j] < 20) * att
            nonatt_q20plus = (self.QRisk[:,j] >=20) * nonatt
            nonatt_q20 = (self.QRisk[:,j] < 20) * nonatt

            att_q20plus_chol = self.chol[att_q20plus,j].mean()
            att_q20_chol = self.chol[att_q20,j].mean()
            nonatt_q20plus_chol = self.chol[nonatt_q20plus,j].mean()
            nonatt_q20_chol = self.chol[nonatt_q20,j].mean()

            # who gets treatment offered, who takes it up?
            stat_offered = (self.Statins_Offered[:,j]>0) * att_q20plus
            stat_notoffered = (self.Statins_Offered[:,j]==0) * att_q20plus
            stat_takeup = (self.Statins[:,j]>0) * att_q20plus



            #  mean change in AHT over next 4 years
            chol_change_stat_takeup = (self.chol[stat_takeup,j+4] - self.chol[stat_takeup,j]).mean()
            chol_change_stat_notoffered = (self.chol[stat_notoffered,j+4] - self.chol[stat_notoffered,j]).mean()
            chol_change_att20 = (self.chol[att_q20,j+4] - self.chol[att_q20,j]).mean()
            chol_change_nonatt20plus = (self.chol[nonatt_q20plus,j+4] - self.chol[nonatt_q20plus,j]).mean()
            chol_change_nonatt20 = (self.chol[nonatt_q20,j+4] - self.chol[nonatt_q20,j]).mean()



            a1 += att.sum()
            a2 += nonatt.sum()

            b1 += att_q20plus.sum()
            b2 += att_q20plus_chol.sum()
            b3 += att_q20.sum()
            b4 += att_q20_chol.sum()
            b5 += nonatt_q20plus.sum()
            b6 += nonatt_q20plus_chol.sum()
            b7 += nonatt_q20.sum()
            b8 += nonatt_q20_chol.sum()

            c1 += stat_offered.sum()
            c2 += stat_takeup.sum()
            c3 += stat_notoffered.sum()

            d1 += chol_change_stat_takeup.sum()
            d2 += chol_change_stat_notoffered.sum()
            d3 += chol_change_att20.sum()
            d4 += chol_change_nonatt20plus.sum()
            d5 += chol_change_nonatt20.sum()


        C_Stats={}
        C_Stats['a.0 Years covered'] = uptoyear
        C_Stats['a.1 Number of attending individuals in one year'] = a1 / float(uptoyear)
        C_Stats['a.2 Number of nonattending individuals in one year'] = a2 / float(uptoyear)

        C_Stats['b.1 Number of attending individuals with QRisk >20'] = b1 / float(uptoyear)
        C_Stats['b.2 their mean Chol level'] = b2 / float(uptoyear)
        C_Stats['b.3 Number of attending individuals with QRisk <20'] = b3 / float(uptoyear)
        C_Stats['b.4 their mean Chol level'] = b4 / float(uptoyear)
        C_Stats['b.5 Number of nonattending individuals with QRisk >20'] = b5 / float(uptoyear)
        C_Stats['b.6 their mean Chol level'] = b6 / float(uptoyear)
        C_Stats['b.7 Number of nonattending individuals with QRisk <20'] = b7 / float(uptoyear)
        C_Stats['b.8 their mean Chol level'] = b8 / float(uptoyear)

        C_Stats['c.1 Statins offered'] = c1 / float(uptoyear)
        C_Stats['c.2 Statins takeup'] = c2 / float(uptoyear)
        C_Stats['c.3 Statins notoffered'] = c3 / float(uptoyear)


        C_Stats['d.1 mean Chol change for Statins takeup group'] = d1 / float(uptoyear)
        C_Stats['d.2 mean Chol change for Statins notoffered group'] = d2 / float(uptoyear)
        C_Stats['d.3 mean Chol change for attending group with QRISK <20'] = d3 / float(uptoyear)
        C_Stats['d.4 mean Chol change for nonattending group with QRISK >20'] = d4 / float(uptoyear)
        C_Stats['d.5 mean Chol change for nonattending group with QRISK <20'] = d5 / float(uptoyear)



        self.Cstats = C_Stats



    def CalculateSmokers(self):
        '''evaluate who smoked and who quitted'''

        att = self.Attending.sum(axis=1)>0
        nonatt = self.Attending.sum(axis=1)==0
        # create matrix of smoking
        self.smokers= np.zeros((self.population_size,self.simulation_time),dtype=bool)

        # identify quitters:
        # define a quitter as someone who stops smoking (getting from cat 2+ down to 1)
        # create matrix of quitting
        self.quitters= np.zeros((self.population_size,self.simulation_time),dtype=bool)
        for i in range(self.population_size):
            for j in range(self.simulation_time):
                if self.q_smoke_cat[i,j]>1:
                    self.smokers[i,j] = True
                if j>0:
                    if self.q_smoke_cat[i,j-1]>1:
                        if self.q_smoke_cat[i,j]==1:
                            self.quitters[i,j] = True



        s=0
        ns=0

        for i in range(self.population_size):
            if att[i] == True:
                idx = np.where(self.Attending[i] == True)[0]

                if self.q_smoke_cat[i,idx[0]]>1:
                    s+=1
                else:
                    ns+=1


        # how many of those smokers are offered smoking cessation and who of those quit within a year
        # quitrates for 1 year
        s_offered_quit=0
        s_notoffered_quit=0
        s_notattending_quit = 0
        # quitrates for 5 years after quitting
        s_offered_quit5yrs = 0
        s_notoffered_quit5yrs = 0
        s_notattending_quit5yrs = 0

        s_offered_relapse5yrs = 0
        s_notoffered_relapse5yrs = 0
        s_notattending_relapse5yrs = 0

        s_offered_relapse10yrs = 0
        s_notoffered_relapse10yrs = 0
        s_notattending_relapse10yrs = 0

        s_offered_nonquit=0
        s_notoffered_notquit=0
        s_notattending_notquit = 0
        s_offered_notquit5yrs = 0
        s_notoffered_notquit5yrs = 0
        s_notattending_quit5yrs = 0
        s_notattending_notquit5yrs = 0

        s_notattending = 0
        ns_notattending = 0


        for i in range(self.population_size):

            if att[i] == True:
                # if person is attending
                idx = np.where(self.Attending[i]==True)[0]
                # if it is a smoker
                if self.q_smoke_cat[i,idx[0]]>1:
                    if idx[0]<(self.simulation_time-7):
                        # if smoking cessation is offered
                        if self.SmokingCessation_Offered[i,idx[0]]==1:
                            if self.q_smoke_cat[i,idx[0]+1]==1:
                                s_offered_quit+=1
                                # check relapse
                                if self.q_smoke_cat[i,idx[0]+5]!=1:
                                    s_offered_relapse5yrs += 1
                            else:
                                s_offered_nonquit+=1
                            # check whether still smoking  5 years post quitting
                            if self.q_smoke_cat[i,idx[0]+5]==1:
                                s_offered_quit5yrs+=1
                            else:
                                s_offered_notquit5yrs+=1
                        # if smoking cessation is not offered
                        else:
                            if self.q_smoke_cat[i,idx[0]+1]==1:
                                s_notoffered_quit+=1
                                if self.q_smoke_cat[i,idx[0]+5]!=1:
                                    s_notoffered_relapse5yrs += 1
                            else:
                                s_notoffered_notquit+=1
                            # check whether still smoking  5 years post quitting
                            if self.q_smoke_cat[i,idx[0]+5]==1:
                                s_notoffered_quit5yrs+=1
                            else:
                                s_notoffered_notquit5yrs+=1
            else:
                # if person is not attending

                if self.smokers[i].sum()>0:
                    # check when they are smokers
                    s_notattending += 1
                    idx1 = np.where(self.smokers[i]==True)[0]
                    if idx1[0]<(self.simulation_time-7):

                        if self.quitters[i,idx1[0]+1]==True:
                            # if quitted the year afterward

                            s_notattending_quit += 1

                            if self.q_smoke_cat[i,idx1[0]+5]>1:
                                s_notattending_relapse5yrs += 1
                        else:
                            s_notattending_notquit += 1


                        # if this person after 5 years exsmoker
                        if self.q_smoke_cat[i,idx1[0] + 5] == 1:
                            s_notattending_quit5yrs += 1
                        else:
                            s_notattending_notquit5yrs += 1



                else:
                    # nonattending nonsmoker
                    ns_notattending += 1




        Smoker_Stats={}
        Smoker_Stats['S offered SC and (a) quitting within 1 year'] = s_offered_quit
        Smoker_Stats['S offered SC and (b) not quitting within 1 year'] = s_offered_nonquit
        Smoker_Stats['S offered SC and (c) ex-smoker after 5 years'] = s_offered_quit5yrs
        Smoker_Stats['S offered SC and (d) not ex-smoker after 5 years'] = s_offered_notquit5yrs
        Smoker_Stats['S offered SC and (e) quitting and relapsing within 5 years'] = s_offered_relapse5yrs


        Smoker_Stats['S offered SC, total'] = s_offered_quit + s_offered_nonquit

        Smoker_Stats['S not offered SC and (a) quitting within 1 year'] = s_notoffered_quit
        Smoker_Stats['S not offered SC and (b) not quitting within 1 year'] = s_notoffered_notquit
        Smoker_Stats['S not offered SC and (c) ex-smoker after 5 years'] = s_notoffered_quit5yrs
        Smoker_Stats['S not offered SC and (d) not ex-smoker after 5 years'] = s_notoffered_notquit5yrs
        Smoker_Stats['S not offered SC and (e) quitting and relapsing within 5 years'] = s_notoffered_relapse5yrs

        Smoker_Stats['S not offered SC, total'] = s_notoffered_quit + s_notoffered_notquit

        Smoker_Stats['non-attending smokers'] = s_notattending
        Smoker_Stats['non-attending smokers and (a) quitting within 1 year'] = s_notattending_quit
        Smoker_Stats['non-attending smokers and (b) not quitting within 1 year'] = s_notattending_notquit
        Smoker_Stats['non-attending smokers and (c) ex-smoker after 5 years'] = s_notattending_quit5yrs
        Smoker_Stats['non-attending smokers and (d) not ex-smoker after 5 years'] = s_notattending_notquit5yrs
        Smoker_Stats['non-attending smokers and (e) quitting and relapsing within 5 years'] = s_notattending_relapse5yrs

        Smoker_Stats['non-attending nonsmokers']  =ns_notattending
        Smoker_Stats['non-attending total']  = att.sum()
        Smoker_Stats['non-attending total']  = nonatt.sum()

        Smoker_Stats['total smokers at time of first HC'] = s
        Smoker_Stats['total nonsmokers at time of first HC'] = ns

        self.SmokerStats = Smoker_Stats



    def LE_per_extra_years_smoked(self):
        '''for smokers, calculate how much life expectancy changes per extra years smoked'''

        # compare smokers with quitters
        ID = np.unique(self.id)
        ID_num = np.zeros(ID.size,dtype = int)
        ID_smoker = np.zeros(ID.size,dtype = int)


        for i in range(ID.size):
            idx = np.where(self.id == ID[i])[0]
            # smoker or not?
            ID_smoker[i] = self.q_smoke_cat[idx[0],0] > 1
            ID_num[i] = idx.size



        # which of these are interesting?
        # we look only at those ID's that have more than 4 entries and are smokers

        check = (ID_smoker == 1) * (ID_num > 4)

        check_ID = ID[check]

        dLEdY = []

        for j in range(check.sum()):
            z = np.where(self.id == check_ID[j])[0]
            ys = np.zeros(z.size)
            ad = np.zeros(z.size)

            try:
                for i_idx, i in enumerate(z):

                    death = np.where(self.Death[i]==True)[0][0]
                    age_death = self.age[i,death]
                    years_smoked = ((self.q_smoke_cat[i] > 1)*(self.alive[i])).sum()
                    ys[i_idx] = years_smoked
                    ad[i_idx] = age_death

                    # now sort both in ascending order
                srt = np.argsort(ys)
                ys_sorted = ys[srt]
                ad_sorted = ad[srt]

                # now check what is the change of life expectancy wrt change in extra years smoked

                deltaLE = []

                for r in range(z.size):
                    for q in range(r):
                        #print r,q
                        year_diff = ys_sorted[r] - ys_sorted[q]
                        LE_diff = ad_sorted[r] - ad_sorted[q]
                        #print year_diff, LE_diff
                        dLE = LE_diff/year_diff
                        if np.isfinite(dLE) == True:
                            deltaLE.append(dLE)

                deltaLE = np.array(deltaLE)

            except:
                deltaLE=np.ones(10) * np.nan

            dLEdY.append(deltaLE.mean())


        dLEdY = np.array(dLEdY)

        nonnan = np.isfinite(dLEdY)
        dLEdY_nonnan = dLEdY[nonnan]

        return dLEdY_nonnan



    def CalculateDiseasePrevalences(self):

        '''calculates disease prevalences for males and females'''

        self.prev_dementia = np.zeros((100,2))
        self.prev_lungcancer = np.zeros((100,2))
        self.prev_ihd = np.zeros((100,2))
        self.prev_stroke = np.zeros((100,2))
        self.prev_diabetes = np.zeros((100,2))

        for i in range(20,100):

            age = self.age == i

            alive = self.alive * age # people alive at this age
            justdied = self.Death * age # people dying in this year of life

            incl = alive + justdied # all individuals at a certain age that are either alive or just die at this age

            # male

            incl_male = incl * (self.gender==1)[:,np.newaxis]
            agegroup_male = float(incl_male.sum())

            self.prev_dementia[i,0] = self.Dementia[incl_male].sum()/agegroup_male
            self.prev_lungcancer[i,0] = self.LungCancer[incl_male].sum()/agegroup_male
            self.prev_ihd[i,0] = self.IHD[incl_male].sum()/agegroup_male
            self.prev_stroke[i,0] = self.Stroke[incl_male].sum()/agegroup_male
            # for diabetes, assume 6.5 as the cutoff point
            high_glyhb = self.glyhb>6.5
            self.prev_diabetes[i,0] = (high_glyhb * incl_male).sum()/agegroup_male

            # female

            incl_female = incl * (self.gender==0)[:,np.newaxis]
            agegroup_female = float(incl_female.sum())

            self.prev_dementia[i,1] = self.Dementia[incl_female].sum()/agegroup_female
            self.prev_lungcancer[i,1] = self.LungCancer[incl_female].sum()/agegroup_female
            self.prev_ihd[i,1] = self.IHD[incl_female].sum()/agegroup_female
            self.prev_stroke[i,1] = self.Stroke[incl_female].sum()/agegroup_female
            # for diabetes, assume 6.5 as the cutoff point
            high_glyhb = self.glyhb>6.5
            self.prev_diabetes[i,1] = (high_glyhb * incl_female).sum()/agegroup_female


    def SaveParameters(self):
        '''saves the current set of parameters and model settings in a separate csv file'''

        t = time.gmtime()

        # create folder with calender day containing parameter files
        f='parameter_logfiles_%d-%02d-%02d' % (t.tm_year, t.tm_mon, t.tm_mday)

        if not os.path.exists(f):
            os.makedirs(f)


        # which parameter was changed?
        # get parameter index
        s = copy.deepcopy(self.UP_R)
        s1 = copy.deepcopy(self.UP)

        par_idx = 0
        par_changed = 0
        par_value = 0

        items = s.items()
        items.sort()

        for keys, values in items:
            # check whether the item is float, list or array
            if type(values) == float or type(values) == bool:
                if s1[keys] != values:
                    par_changed = par_idx
                    par_value = s1[keys]
                    break
                par_idx += 1

            else:
                for i in range(len(values)):
                    if s1[keys][i] != values[i]:
                        par_changed = par_idx
                        par_value = s1[keys][i]
                        break
                    par_idx += 1

        if par_changed == 0:
            parval_str = '000000'
        else:
            parval_str = str(par_value * 1e6)

        # now a unique parameter file can be created, using time, parameter changed, and parameter value to ensure uniqueness of file name
        self.par_filename = '%s/parameters_t%04d%02d%02d_%02d%02d%02d_p%02d_v%s' % (f,t.tm_year, t.tm_mon, t.tm_mday,t.tm_hour,t.tm_min,t.tm_sec,par_changed,parval_str[:4])


        outfile = open('%s.log' % (self.par_filename), 'w' )
        outfile.write('Health Check Simulation Run, %04d/%02d/%02d at %02d:%02d:%02d\n\n' % (t.tm_year, t.tm_mon, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec))
        outfile.write('population_size\t%d\n' % self.population_size)
        outfile.write('simulation_time\t%d\n' % self.simulation_time)
        outfile.write('Health Check\t%s\n\n' % self.Health_Checks)

        outfile.write('MODEL PARAMETERS')
        if par_changed == 0:
            outfile.write(' (unchanged):\n\n')
        else:
            outfile.write(' (changed):\n\n')

        for key, value in items:
            outfile.write( str(key) + '\t' + str(value) + '\n' )


        return t




    def SaveResults(self):
        '''saves all relevant results in pytables format with unique specifier'''

        day = time.asctime()

        self.savepath='results_data/%s_%s_%s%'  % (day[-4:], day[4:7], day[8:10])
        if not os.path.exists(self.savepath):
            os.makedirs(self.savepath)

        self.resultfilename = 'results_popsize_%d_simtime_%d_HC%d_ID%d.h5' % (self.population_size, self.simulation_time, self.Health_Checks, int(time.time()))


    def ReleaseMemory(self):
        '''this deletes all numpy arrays of current class and releases their memory
        don't do this if you want to keep data for this simulation!!'''

        attributes = []
        for attr, value in self.__dict__.iteritems():
            # check if attribute is numpy array

            isnumpy = eval('type(%s.%s).__module__ == np.__name__' % ('self',attr))

            if isnumpy == True:
                attributes.append(attr)


        # then going through list and deleting objects
        for a in attributes:
            exec('del %s.%s' % ('self',a))


    def CalculateLY(self):
        self.LY = self.age[self.Death]


    def Run(self):
        self.Simulate()
        self.CalculateQALY_CVD()
        self.CalculateQALY()
        self.CalculateLY()
