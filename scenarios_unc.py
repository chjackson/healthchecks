import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 25000
n_cpus = 4
prefix = "results/papermar/unc"
randpars = True

## used for running a large population in a set of batches
## this program is called by a bash script with the batch number as a command line argument
## results are saved to CSV files with the run/batch number in the file name
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

if (len(sys.argv) > 2):
    n_cpus = int(sys.argv[2])

nsc = 18
npops = 6 # number of subpopulations
npopsn = 19 # number of subpopulations to calculate size of
nouts = 38 # number of outputs (LY, QALY etc)
noutsp = 4  # number of post-HC outputs 
baseage_min = 40  # Restrict age range of baseline population
baseage_max = 45 

M = np.zeros((nsc, npops, nouts))
S = np.zeros((nsc, npops, nouts))
N = np.zeros((nsc, npopsn), dtype=int) # subpopulation sizes and other count data in main run
P = np.zeros((nsc, noutsp, 2))
N[...,0] = ps

if (run==1):
    try:
        os.remove("%s_mean.csv" % (prefix))
        os.remove("%s_SD.csv" % (prefix))
        os.remove("%s_N.csv" % (prefix))
        os.remove("%s_post.csv" % (prefix))
        os.remove("%s_pars.csv" % (prefix))
    except OSError:
        pass

def SaveResults(M, S, N, P, nsc, npops, nouts, prefix):
    # arrange results as one row per scenario
    # block of results from current run appended to previous runs
    runid = np.array([run]*nsc) # vector of run IDs, first column of result block
    Msave = np.column_stack((runid, M.reshape((nsc,npops*nouts))))
    Ssave = np.column_stack((runid, S.reshape((nsc,npops*nouts))))
    Nsave = np.column_stack((runid, N))
    Psave = np.column_stack((runid, P.reshape((nsc,noutsp*2))))

    with open("%s_mean.csv" % (prefix), 'a') as f:
        np.savetxt(f, Msave, delimiter=",")
    with open("%s_SD.csv" % (prefix), 'a') as f:
        np.savetxt(f, Ssave, delimiter=",")
    with open("%s_N.csv" % (prefix), 'a') as f:
        np.savetxt(f, Nsave, delimiter=",", fmt="%d")
    with open("%s_post.csv" % (prefix), 'a') as f:
        np.savetxt(f, Psave, delimiter=",")
    
## Without HC
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, randpars=randpars, baseage_min=baseage_min, baseage_max=baseage_max)
H.Run()

## Basic HC model

def initmodels():
    if (randpars):
        UP = H.GetUncertainParameters()
    else: UP = None
    H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, pars=UP, baseage_min=baseage_min, baseage_max=baseage_max)
    return H1

def runmodels(H1, H2, H, M, S, N, P, r):
    H2.Run()
    M[r,],S[r,],N[r,],P[r,] = gr.GetResults_paper(H1, H2, H, M[r,], S[r,], N[r,], P[r,])
    r += 1
    return M, S, N, P, r

r = 0

H1 = initmodels()
M, S, N, P, r = runmodels(H, H1, H, M, S, N, P, r)

## All other scenarios: relative to base case

H2 = initmodels()
H2.SetUncertainParameter('up_HC_include_bp_registers', 1)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [50, 74.9])
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [40, 80.9])
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_age_limit', [50, 80.9])
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
takeup = H2.GetUncertainParameters()['up_HC_takeup']
H2.SetUncertainParameter('up_HC_takeup', takeup*1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('ses5_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('sm_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('q5_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
patt = H2.GetUncertainParameters()['up_HC_offer_not_prev_att']
H2.SetUncertainParameter('up_HC_offer_not_prev_att', patt*1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

# attenders keep attending.  Compare with no HC, not with current practice
Hatt = initmodels()
Hatt.SetUncertainParameter('up_HC_takeup_rr_prev_att', 0.7/takeup)
Hatt.SetUncertainParameter('up_HC_takeup_rr_not_prev_att', 0.3/takeup)
M, S, N, P, r = runmodels(H, Hatt, H, M, S, N, P, r)

# attenders keep attending, and target non-attenders.  Compare with new base 
H2 = initmodels()
H2.SetUncertainParameter('up_HC_takeup_rr_prev_att', 0.7/takeup)
H2.SetUncertainParameter('up_HC_takeup_rr_not_prev_att', 0.3/takeup)
H2.SetUncertainParameter('up_HC_offer_not_prev_att', patt*1.3)
M, S, N, P, r = runmodels(Hatt, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('up_HC_aht_presc_Q20plus', 0.06)  # was 0.0248
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_smoker_ref', 0.09) # was 0.036
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_weight_ref', 0.6875) # was 0.275
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H2 = initmodels()
H2.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
H2.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('up_HC_aht_presc_Q20plus', 0.06)  # was 0.0248
H2.SetUncertainParameter('up_HC_smoker_ref', 0.09) # was 0.036
H2.SetUncertainParameter('up_HC_weight_ref', 0.6875) # was 0.275
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

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
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

SaveResults(M, S, N, P, nsc, npops, nouts, prefix)

# save parameter values for current iteration
UP = H.GetUncertainParameters()
par_vals = [ v for v in UP.values() ]
with open("%s_pars.csv" % (prefix), 'a') as f:
    pars = np.hstack(par_vals)
    pars = pars.reshape((1, pars.size))
    np.savetxt(f, pars, delimiter=",")

assert r == nsc
