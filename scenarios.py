# Script to produce results for all scenarios
# Produces CSVs which are then converted to results tables in tables.r 

import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 25000
n_cpus = 4

prefix = "results/paperfeb/scen"
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

nsc = 18 # number of scenarios
npops = 6 # number of subpopulations to calculate mean outcome for 
npopsn = 19 # number of subpopulations to calculate size of
nouts = 38 # number of outputs (LY, QALY etc).  8 + 3*(number of event outcomes=5)*(number of ages=2)
noutsp = 4 # number of post-HC outputs 
baseage_min = 40  # Restrict age range of baseline population
baseage_max = 45 

M = np.zeros((nsc, npops, nouts))
S = np.zeros((nsc, npops, nouts))
N = np.zeros((nsc, npopsn), dtype=int) # subpopulation sizes and other count data in main run
P = np.zeros((nsc, noutsp, 2)) # mean and SD together 
N[...,0] = ps

def SaveResults(M, S, N, P, nsc, npops, nouts, prefix):
    np.savetxt("%s_mean_%s.csv" % (prefix,run), M.reshape((nsc,npops*nouts)) , delimiter=",") # rows are populations (eligible, treated...), cols are outputs (LY, QALY...)
    np.savetxt("%s_SD_%s.csv"   % (prefix,run), S.reshape((nsc,npops*nouts)), delimiter=",") #
    np.savetxt("%s_N_%s.csv"    % (prefix,run), N, fmt="%d", delimiter=",") #
    np.savetxt("%s_post_%s.csv" % (prefix,run), P.reshape((nsc,noutsp*2)), fmt="%s", delimiter=",")
    
## Without HC
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, baseage_min=baseage_min, baseage_max=baseage_max)
H.Run()
## Basic HC model

def initmodels():
    H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, baseage_min=baseage_min, baseage_max=baseage_max)
    return H1

def runmodels(H1, H2, M, S, N, P, r):
    H2.Run()
    M[r,],S[r,],N[r,],P[r,] = gr.GetResults_paper(H1, H2, M[r,], S[r,], N[r,], P[r,])
    SaveResults(M, S, N, P, nsc, npops, nouts, prefix)
    r += 1
    return M, S, N, P, r

r = 0

## Base case scenario: relative to no HCs 

H1 = initmodels()
M, S, N, P, r = runmodels(H, H1, M, S, N, P, r)

## All other scenarios: relative to base case

H2 = initmodels()
H2.SetUncertainParameter('up_HC_include_bp_registers', 1)
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [50, 74.9])
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [40, 80.9])
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [50, 80.9])
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
takeup = H2.GetUncertainParameters()['up_HC_takeup']
H2.SetUncertainParameter('up_HC_takeup', takeup*1.3)
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

# See odds.r for calculation of these odds ratios, which correspond to 30% increases in
# probability of attendance for the corresponding category.
H2 = initmodels()
H2.SetUncertainParameter('ses5_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('sm_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('q5_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
patt = H2.GetUncertainParameters()['up_HC_offer_not_prev_att']
H2.SetUncertainParameter('up_HC_offer_not_prev_att', patt*1.3)
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

# attenders keep attending.  Compare with no HC, not with current practice
Hatt = initmodels()
Hatt.SetUncertainParameter('up_HC_takeup_rr_prev_att', 0.7/takeup)
Hatt.SetUncertainParameter('up_HC_takeup_rr_not_prev_att', 0.3/takeup)
M, S, N, P, r = runmodels(H, Hatt, M, S, N, P, r)

# attenders keep attending, and target non-attenders.  Compare with new base 
H2 = initmodels()
H2.SetUncertainParameter('up_HC_takeup_rr_prev_att', 0.7/takeup)
H2.SetUncertainParameter('up_HC_takeup_rr_not_prev_att', 0.3/takeup)
H2.SetUncertainParameter('up_HC_offer_not_prev_att', patt*1.3)
M, S, N, P, r = runmodels(Hatt, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('up_HC_aht_presc_Q20plus', 0.06)  # was 0.0248
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_smoker_ref', 0.09) # was 0.036
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_weight_ref', 0.6875) # was 0.275
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
H2.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('up_HC_aht_presc_Q20plus', 0.06)  # was 0.0248
H2.SetUncertainParameter('up_HC_smoker_ref', 0.09) # was 0.036
H2.SetUncertainParameter('up_HC_weight_ref', 0.6875) # was 0.275
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)

# high treatment, uptake, eligibility
H2 = initmodels()
# eligibility
H2.SetUncertainParameter('up_HC_include_bp_registers', 1)
# uptake
H2.SetUncertainParameter('up_HC_takeup', takeup*1.3)
# treatment
H2.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
H2.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('up_HC_aht_presc_Q20plus', 0.06)  # was 0.0248
H2.SetUncertainParameter('up_HC_smoker_ref', 0.09) # was 0.036
H2.SetUncertainParameter('up_HC_weight_ref', 0.6875) # was 0.275
M, S, N, P, r = runmodels(H1, H2, M, S, N, P, r)


assert r == nsc
