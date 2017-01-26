# Script to produce results for all scenarios
# Produces CSVs which can then be manipulated in scenarios.r to get graphs or tables
## This version uses the subset boosting trick

## Restrict age range of baseline population, e.g. 35-40


import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 200000
n_cpus = 4

prefix = "results/paper/scen"
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

nsc = 15 # number of scenarios
npops = 5 # number of subpopulations to calculate mean outcome for 
npopsn = 12 # number of subpopulations to calculate size of
nouts = 38 # number of outputs (LY, QALY etc).  8 + 3*(number of event outcomes=5)*(number of ages=2)
baseage_min = 40 
baseage_max = 45 

M = np.zeros((nsc, npops, nouts))
S = np.zeros((nsc, npops, nouts))
N = np.zeros((nsc, npopsn), dtype=int) # subpopulation sizes and other count data in main run
N[...,0] = ps

def SaveResults(M, S, N, nsc, npops, nouts, prefix):
    np.savetxt("%s_mean_%s.csv" % (prefix,run), M.reshape((nsc,npops*nouts)) , delimiter=",") # rows are populations (eligible, treated...), cols are outputs (LY, QALY...)
    np.savetxt("%s_SD_%s.csv" % (prefix,run), S.reshape((nsc,npops*nouts)), delimiter=",") #
    np.savetxt("%s_N_%s.csv" % (prefix,run), N, fmt="%d", delimiter=",") #
    
## Without HC
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, baseage_min=baseage_min, baseage_max=baseage_max)
H.Run()
## Basic HC model

def initmodels():
    H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, baseage_min=baseage_min, baseage_max=baseage_max)
    return H1

def runmodels(H1, H2, M, S, N, r):
    H2.Run()
    M[r,],S[r,],N[r,] = gr.GetResults_paper(H1, H2, M[r,], S[r,], N[r,])
    SaveResults(M, S, N, nsc, npops, nouts, prefix)
    r += 1
    return M, S, N, r

r = 0

## Base case scenario: relative to no HCs 

H1 = initmodels()
M, S, N, r = runmodels(H, H1, M, S, N, r)

## All other scenarios: relative to base case

H2 = initmodels()
H2.SetUncertainParameter('up_HC_include_bp_registers', 1)
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [50, 74.9])
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [40, 80.9])
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [50, 80.9])
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_takeup', 0.85)
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
ses = H2.GetUncertainParameters()['up_SES_vec']
H2.SetUncertainParameter('up_SES_vec', [ses[0], ses[1], ses[2], ses[3]*1.1, ses[4]*1.2])
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_smoker_vec', [1, 1])
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
qr = H2.GetUncertainParameters()['up_QRisk_vec']
H2.SetUncertainParameter('up_QRisk_vec', [qr[0], qr[1], qr[2], qr[3]*1.2, qr[4]*1.2])
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_offer_not_prev_att', 0.4)
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('up_HC_aht_presc_Q20plus', 0.06)  # was 0.0248
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_smoker_ref', 0.09) # was 0.036
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_weight_ref', 0.6875) # was 0.275
M, S, N, r = runmodels(H1, H2, M, S, N, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
H2.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('up_HC_aht_presc_Q20plus', 0.06)  # was 0.0248
H2.SetUncertainParameter('up_HC_smoker_ref', 0.09) # was 0.036
H2.SetUncertainParameter('up_HC_weight_ref', 0.6875) # was 0.275
M, S, N, r = runmodels(H1, H2, M, S, N, r)

assert r == nsc
