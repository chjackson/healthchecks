# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 08:36:15 2015

additional function files for the Health Checks Microsimulation model

@author: arno
"""



import numpy as np


def weightedChoice(weights, objects):
    """Return a random item from objects, with the weighting defined by weights 
    (which must sum to 1)."""
    cs = np.cumsum(weights) #An array of the weights, cumulatively summed.
    idx = np.sum(cs < np.random.random()) #Find the index of the first weight over a random value.
    return objects[idx]
    
    

def inc_to_prev(inc):
    ''' 
    converts incidences array to prevalences array
    '''
    
    prev = np.zeros(inc.size)
    prev[0] = inc[0]
    
    for i in range(1,inc.size):
        prev[i] = prev[i-1] + (1-prev[i-1]) * inc[i]
    
    return prev



def inc_to_prev_new(inc,mort):
    ''' 
    converts incidences array to prevalences array, considering only alives
    '''
    
    prev = np.zeros(inc.size)
    prev[0] = inc[0]
    
    for i in range(1,inc.size):
        prev[i] = prev[i-1] + (1-prev[i-1]) * inc[i] * (1-mort[i-1])
    
    return prev    


def mort_to_alive(mort):
    ''' 
    converts ONS mortality rates to vector of alive population
    '''
    
    alive = np.ones(mort.size)
    
    
    for i in range(1,mort.size):
        alive[i] = alive[i-1] - alive[i-1] * mort[i-1]
    
    return alive    
    
    

def inc_to_prev_new2(inc,mort):
    ''' 
    converts ONS mortality rates to vector of alive population
    and then uses the alive vector for converting incidences to prevalences
    
    '''
    
    alive = mort_to_alive(mort)

    prev = np.zeros(inc.size)
    prev[0] = inc[0]
    
    for i in range(1,inc.size):
        prev[i] = prev[i-1] + (1-prev[i-1]) * inc[i] * alive[i]
    
    return prev    


def inc_to_prev_corrected(inc, mortality_all, mortality_other, CFR):
    ''' 
    converts ONS mortality rates to vector of alive population
    and then uses the alive vector for converting incidences to prevalences
    
    '''
    
    # first, establish number of people alive at ages    
    n_alive = mort_to_alive(mortality_all)

    prev = np.zeros(inc.size)
    prev[0] = inc[0]
    
    for i in range(1,inc.size):
        prev[i] = prev[i-1] + (1-prev[i-1]) * inc[i] * alive[i]
    
    return prev    



    

def ConvertToQ(p):
    '''converts a vector p containing annual probabilities (ie incidences on Dismod life tables)
    to a vector q containing annual probabilities within the period for which the
    risk score holds (ie 10 years in QRisk, 20 years for CAIDE)
    
    ' converting from the p to the q scale'
    '''
    
    q = np.zeros(p.shape)
    
    p_1 = 1-p
    for i in range(p.size):
        if i == 0:
            q[i] =  p[i]
        elif i == 1:
            q[i] = p[i] * p_1[0]
        else:
            q[i] =  p[i] * np.prod(p_1[0:i])
    
    return q


def ConvertToP(q):
    '''opposite of ConvertToQ(p)
    
    ' converting back from the q to the p scale'
    '''
    
    p = np.zeros(q.shape)
    
    
    for i in range(q.size):
        if i == 0:
            p[i] =  q[i]
        elif i == 1:
            p[i] = q[i] / (1-q[0])
        else:
            conv_prob = 1
            for j in range(i):
                conv_prob *= (1-q[j])
            p[i] = q[i] / conv_prob
    
    return p





def SplitTime(reftime):
    '''splits the time in seconds in a more human readable format
    hrs:min:sec'''
    
    hrs = np.floor(reftime/3600.0)
    mins = np.floor(np.remainder(reftime/60.0,60))
    sec = np.floor(np.remainder(reftime,60))
    
    return hrs,mins,sec

# determine vector sizes for indexing during ELSA matching process


def bp_vectors(E):
    
    age = np.unique(E[:,1])
    sex = np.unique(E[:,2])
    bmicat = np.unique(E[:,3])
    bmibin = np.unique(E[:,4])
    sys = np.unique(E[:,5])
    dia = np.unique(E[:,6])
    
    return age, sex, bmicat, bmibin, sys, dia


def bp_vectorsize(E):
    
    age, sex, bmicat, bmibin, sys, dia = bp_vectors(E)
    
    # how many category bins are there?
    n_sex = sex.size
    n_age = age.size
    n_bmi = bmicat.size
    n_bmibin = bmibin.size
    n_sys = sys.size
    n_dia = dia.size    
    
    return n_sex, n_age, n_bmi, n_bmibin, n_sys, n_dia
    


def glyhb_vectors(E):
    
    age = np.unique(E[:,1])
    sex = np.unique(E[:,2])
    bmi = np.unique(E[:,3])
    # many in the eq category have missing values
    eq_nan = np.isnan(E[:,4])
    eq_finite = np.isfinite(E[:,4])
    eq = np.unique(E[eq_finite,4])
    diab = np.unique(E[eq_finite,5])
    gly = np.unique(E[:,7])
    
    return age, sex, bmi, eq_nan, eq_finite, eq, diab, gly


def glyhb_vectorsize(E):
    
    age, sex, bmi, eq_nan, eq_finite, eq, diab, gly = glyhb_vectors(E)
    
    n_sex = sex.size
    n_age = age.size
    n_bmi = bmi.size
    n_eq = eq.size
    n_diab = diab.size
    n_gly = gly.size
    
    return n_sex, n_age, n_bmi, n_eq, n_diab, n_gly
    


def chol_vectors(E):
    
    age = np.unique(E[:,1])
    sex = np.unique(E[:,2])
    bmi = np.unique(E[:,3])
    # many in the eq category have missing values
    eq_nan = np.isnan(E[:,4])
    eq_finite = np.isfinite(E[:,4])
    eq = np.unique(E[eq_finite,4])
    chol = np.unique(E[:,5])    
    
    return age, sex, bmi, eq_nan, eq_finite, eq, chol
    
    
def chol_vectorsize(E):
    
    age, sex, bmi, eq_nan, eq_finite, eq, chol = chol_vectors(E)
    
    # how many category bins are there?
    n_sex = sex.size
    n_age = age.size
    n_bmi = bmi.size
    n_eq = eq.size
    n_chol = chol.size
    
    return n_sex, n_age, n_bmi, n_eq, n_chol


def smoking_vectors(E):
    
    age = np.unique(E[:,1])
    sex = np.unique(E[:,2])
    bmi = np.unique(E[:,3])
    # many in the eq category have missing values
    eq_nan = np.isnan(E[:,4])
    eq_finite = np.isfinite(E[:,4])
    eq = np.unique(E[eq_finite,4])
    scat = np.unique(E[:,5]) # smoking category
    
    return age, sex, bmi, eq_nan, eq_finite, eq, scat

    
def smoking_vectorsize(E):

    age, sex, bmi, eq_nan, eq_finite, eq, scat = smoking_vectors(E)
    
    n_sex = sex.size
    n_age = age.size
    n_bmi = bmi.size
    n_eq = eq.size
    n_scat = scat.size
    
    return n_sex, n_age, n_bmi, n_eq, n_scat


def bmi_vectors(E):
    
    age8 = np.unique(E[:,1])
    age4 = np.unique(E[:,2])
    sex = np.unique(E[:,3])
    s_finite = np.isfinite(E[:,4])
    scat = np.unique(E[s_finite,4])
    bmi = np.unique(E[:,5])
    
    return age8, age4, sex, s_finite, scat, bmi


def bmi_vectorsize(E):

    age8, age4, sex, s_finite, scat, bmi = bmi_vectors(E)

    n_sex = sex.size
    n_age4 = age4.size
    n_age8 = age8.size
    n_bmi = bmi.size
    n_scat = scat.size     
    
    return n_sex, n_age4, n_age8, n_bmi, n_scat

###############################################################################
# MATCHING BP
###############################################################################

def calculate_indices_bp(E,row,age,sex,bmicat,sys,dia):
    '''computes individual indices for compressing ELSA data for use in matching
    for blood pressure'''

    age_idx = np.where(E[row,1]==age)[0][0]
    sex_idx = np.where(E[row,2]==sex)[0][0]
    bmi_idx = np.where(E[row,3]==bmicat)[0][0]
    
    sys_idx = np.where(E[row,5]==sys)[0][0]
    dia_idx = np.where(E[row,6]==dia)[0][0]
    
    return age_idx, sex_idx, bmi_idx, sys_idx, dia_idx   
    
    
def calculate_indices_alternative_bp(E,row,age,bmibin,sys):
    '''computes individual indices for compressing ELSA data for use in matching
    for blood pressure, 2nd matching'''
    
    age_idx = np.where(E[row,1]==age)[0][0]
    bmi_idx = np.where(E[row,4]==bmibin)[0][0]    
    sys_idx = np.where(E[row,5]==sys)[0][0]
    
    return age_idx, bmi_idx, sys_idx        



def calculate_Cindex_bp(n_sex, n_bmi, n_sys, n_dia, age_idx, sex_idx, bmi_idx, sys_idx, dia_idx):
    '''computes overall combination index from individual incides
    for blood pressure'''
   
    idx_array = np.zeros(1,dtype=int)

    idx =   age_idx * (n_sex * n_bmi * n_sys * n_dia) + \
            sex_idx * (n_bmi * n_sys * n_dia) + \
            bmi_idx * (n_sys * n_dia) + \
            sys_idx * (n_dia) + dia_idx
    
    idx_array[0] = idx
        
    return  idx_array
    

def calculate_Cindex_alternative_bp(n_bmibin, n_sys, age_idx, bmi_idx, sys_idx):
    '''computes overall combination index from individual incides
    for blood pressure, 2nd matching'''
    
    idx_array = np.zeros(1,dtype=int)

    idx =   age_idx * (n_bmibin * n_sys) + \
            bmi_idx * (n_sys) + sys_idx
    
    idx_array[0] = idx
        
    return  idx_array




###############################################################################
# MATCHING CHOLESTEROL
###############################################################################


def calculate_indices_chol(E, row, age, sex, bmi, eq, chol):
    '''computes individual indices for all factors, given a row in the E matrix'''
    
    
    
    age_idx = np.where(E[row,1]==age)[0][0]
    sex_idx = np.where(E[row,2]==sex)[0][0]
    bmi_idx = np.where(E[row,3]==bmi)[0][0]
    if np.isnan(E[row,4]):
        eq_idx = np.nan
    else:
        eq_idx = np.where(E[row,4]==eq)[0][0]
    chol_idx = np.where(E[row,5]==chol)[0][0]
    
    return age_idx, sex_idx, bmi_idx, eq_idx, chol_idx




def calculate_Cindex_chol(E,age_idx, sex_idx, bmi_idx, eq_idx, chol_idx):
    '''computes overall combination index from individual incides'''
    
    n_sex, n_age, n_bmi, n_eq, n_chol = chol_vectorsize(E)
    if np.isfinite(eq_idx):
        
        idx_array = np.zeros(1,dtype=int)
    
        idx =   age_idx * (n_sex * n_bmi * n_eq * n_chol) + \
                sex_idx * (n_bmi * n_eq * n_chol) + \
                bmi_idx * (n_eq * n_chol) + \
                eq_idx * (n_chol) + chol_idx
        
        idx_array[0] = idx
        
    
    else:
        # deprivation index is missing - three overall indices assigned
        # as missing data means matching all potentiall bands for this cateogory is possible
    
        idx_array = np.zeros(3, dtype=int)
        
        for eq_idx in range(n_eq):
            
            idx =   age_idx * (n_sex * n_bmi * n_eq * n_chol) + \
                    sex_idx * (n_bmi * n_eq * n_chol) + \
                    bmi_idx * (n_eq * n_chol) + \
                    eq_idx * (n_chol) + chol_idx
            
            idx_array[eq_idx] = idx        
    
    return  idx_array





def calculate_Cindex_alternative1_chol(sex_idx, bmi_idx, chol_idx):
    '''computes overall combination index from individual incides'''
       
    idx_array = np.zeros(1,dtype=int)

    idx =   sex_idx * (n_bmi * n_chol) + \
            bmi_idx * (n_chol) + chol_idx
    
    idx_array[0] = idx
        
    return  idx_array
    


def calculate_Cindex_alternative2_chol(sex_idx, chol_idx):
    '''computes overall combination index from individual incides'''
       
    idx_array = np.zeros(1,dtype=int)

    idx =   sex_idx * (n_chol) + chol_idx
    
    idx_array[0] = idx
        
    return  idx_array
    
    


###############################################################################
# MATCHING SMOKING
###############################################################################



###############################################################################
# MATCHING BLOOD GLUCOSE
###############################################################################


def calculate_indices_glyhb(E,row, age, sex, bmi, eq, diab, gly):
    '''computes individual indices for all factors, given a row in the E matrix'''
    
    
    age_idx = np.where(E[row,1]==age)[0][0]
    sex_idx = np.where(E[row,2]==sex)[0][0]
    bmi_idx = np.where(E[row,3]==bmi)[0][0]
    if np.isnan(E[row,4]):
        eq_idx = np.nan
    else:
        eq_idx = np.where(E[row,4]==eq)[0][0]
    diab_idx = np.where(E[row,5]==diab)[0][0]
    gly_idx = np.where(E[row,7]==gly)[0][0]
    
    return age_idx, sex_idx, bmi_idx, eq_idx, diab_idx, gly_idx





def calculate_Cindex_glyhb(n_sex, n_bmi, n_eq, n_diab, n_gly, age_idx, sex_idx, bmi_idx, eq_idx, diab_idx, gly_idx):
    '''computes overall combination index from individual incides'''
    
    if np.isfinite(eq_idx):
        
        idx_array = np.zeros(1,dtype=int)
    
#        idx =   age_idx * (n_sex * n_bmi * n_eq * n_gly) + \
#                sex_idx * (n_bmi * n_eq * n_gly) + \
#                bmi_idx * (n_eq * n_gly) + \
#                eq_idx * (n_gly) + gly_idx

        idx =   age_idx * (n_sex * n_bmi * n_eq * n_diab * n_gly) + \
                sex_idx * (n_bmi * n_eq * n_diab * n_gly) + \
                bmi_idx * (n_eq * n_diab * n_gly) + \
                eq_idx * (n_diab * n_gly) + \
                diab_idx * (n_gly) + gly_idx

        
        idx_array[0] = idx
        
    
    else:
        # deprivation index is missing - three overall indices assigned
        # as missing data means matching all potentiall bands for this cateogory is possible
    
        idx_array = np.zeros(2, dtype=int)
        
        for eq_idx in range(n_eq):
            
#            idx =   age_idx * (n_sex * n_bmi * n_eq * n_gly) + \
#                    sex_idx * (n_bmi * n_eq * n_gly) + \
#                    bmi_idx * (n_eq * n_gly) + \
#                    eq_idx * (n_gly) + gly_idx
                    
            idx =   age_idx * (n_sex * n_bmi * n_eq * n_diab * n_gly) + \
            sex_idx * (n_bmi * n_eq * n_diab * n_gly) + \
            bmi_idx * (n_eq * n_diab * n_gly) + \
            eq_idx * (n_diab * n_gly) + \
            diab_idx * (n_gly) + gly_idx
            
            idx_array[eq_idx] = idx        
    
    return  idx_array
    
    
        
def calculate_Cindex_alternative1_glyhb(n_diab, n_gly, bmi_idx, diab_idx, gly_idx):
    '''computes overall combination index from individual incides'''
       
    idx_array = np.zeros(1,dtype=int)

    idx =   bmi_idx * (n_diab * n_gly)  + \
    diab_idx * (n_gly) + gly_idx
    
    idx_array[0] = idx
        
    return  idx_array
    


def calculate_Cindex_alternative2_glyhb(n_gly, diab_idx, gly_idx):
    '''computes overall combination index from individual incides'''
       
    idx_array = np.zeros(1,dtype=int)

    idx = diab_idx * (n_gly) + gly_idx
    
    idx_array[0] = idx
        
    return  idx_array



###############################################################################
# MATCHING BMI
###############################################################################
