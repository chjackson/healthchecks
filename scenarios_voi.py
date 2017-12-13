import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 200000
n_cpus = 4
prefix = os.path.expanduser('~/scratch/hc/healthchecks/results/paperrev/voi')
randpars = True

## used for running a large population in a set of batches
## this program is called by a bash script with the batch number as a command line argument
## results are saved to CSV files with the run/batch number in the file name
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

if (len(sys.argv) > 2):
    n_cpus = int(sys.argv[2])

nsc = 2 # number of scenarios (was 18 in first submission)

npops = 6 # number of subpopulations
npopsn = 19 # number of subpopulations to calculate size of
nouts = 46 # number of outputs (LY, QALY etc) (was 38 in first submission)
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
        os.remove("%s_senspars.csv" % (prefix))
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

def initmodels(HealthChecks = True):
    if (randpars):
        UP = H.GetUncertainParameters()
    else: UP = None
    H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=HealthChecks, nprocs=n_cpus, randseed=run, pars=UP, baseage_min=baseage_min, baseage_max=baseage_max)
    return H1

def runmodels(H1, H2, H, M, S, N, P, r):
    H2.Run()
    M[r,],S[r,],N[r,],P[r,] = gr.GetResults_paper(H1, H2, H, M[r,], S[r,], N[r,], P[r,])
    r += 1
    return M, S, N, P, r

r = 0

H1 = initmodels()
M, S, N, P, r = runmodels(H, H1, H, M, S, N, P, r)

# high treatment, uptake, eligibility
H2 = initmodels()
# eligibility
H2.SetUncertainParameter('HC_include_bp_registers', 1)
# uptake
H2.SetUncertainParameter('HC_takeup', H2.up_HC_takeup*1.3)
# treatment
H2.SetUncertainParameter('HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('HC_statins_presc_Q20plus', 0.36)  # was 14%
H2.SetUncertainParameter('HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('HC_aht_presc_Q20plus', 0.06)  # was 0.0248
H2.SetUncertainParameter('HC_smoker_ref', 0.09) # was 0.036
H2.SetUncertainParameter('HC_weight_ref', 0.6875) # was 0.275
H2.SetUncertainParameter('HC_age_limit', [40, 79.9]) # was 40,74
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)

H1.ReleaseMemory()
H2.ReleaseMemory()

SaveResults(M, S, N, P, nsc, npops, nouts, prefix)

# save parameter values for current iteration (base case) 

UP = H.GetUncertainParameters()
del UP['DementiaLateRRs']
del UP['CVD_risks']
par_names = np.array(UP.keys(), dtype=str).reshape(1, len(UP))
par_vals = [ v for v in UP.values() ]
par_sizes = [ np.array(v).size for v in UP.values() ]
# create vector of parameter names for the CSV header
namsrep = np.repeat(par_names, par_sizes)
ind = np.hstack( [ v for v in map(range, par_sizes) ] ) + 1
namsrep1 = []
for i in xrange(len(namsrep)):
    namsrep1.append(namsrep[i] + str(ind[i]))
namsrep1 = np.array(namsrep1, dtype=str).reshape(1, len(namsrep1))

with open("%s_pars.csv" % (prefix), 'a') as f:
    if (run==1):
        np.savetxt(f, namsrep1, delimiter=",", fmt="%s")
    pars = np.hstack(par_vals)
    pars = pars.reshape((1, pars.size))
    np.savetxt(f, pars, delimiter=",")

H.ReleaseMemory()
