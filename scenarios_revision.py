import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 200000
n_cpus = 4
prefix = os.path.expanduser('~/scratch/hc/healthchecks/results/paperrev/scen')
randpars = True

## used for running a large population in a set of batches
## this program is called by a bash script with the batch number as a command line argument
## results are saved to CSV files with the run/batch number in the file name
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

if (len(sys.argv) > 2):
    n_cpus = int(sys.argv[2])

nsc = 22 # number of scenarios (was 18 in first submission)

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


## All other scenarios: relative to base case

H2 = initmodels()
H2.SetUncertainParameter('HC_include_bp_registers', 1)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_age_limit', [50, 74.9])
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_age_limit', [40, 79.9])
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_age_limit', [50, 79.9])
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
takeup = H2.up_HC_takeup
H2.SetUncertainParameter('HC_takeup', takeup*1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('ses5_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('sm_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('q5_extra_uptake', 1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
patt = H2.up_HC_offer_not_prev_att
H2.SetUncertainParameter('HC_offer_not_prev_att', patt*1.3)
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

# attenders keep attending.  Compare with no HC, not with current practice
Hsens = initmodels()
takeup = Hsens.up_HC_takeup
Hsens.SetUncertainParameter('HC_takeup_rr_prev_att', 0.7/takeup)
Hsens.SetUncertainParameter('HC_takeup_rr_not_prev_att', 0.3/takeup)
M, S, N, P, r = runmodels(H, Hsens, H, M, S, N, P, r)

# attenders keep attending, and target non-attenders.  Compare with new base 
H2 = initmodels()
H2.SetUncertainParameter('HC_takeup_rr_prev_att', 0.7/takeup)
H2.SetUncertainParameter('HC_takeup_rr_not_prev_att', 0.3/takeup)
H2.SetUncertainParameter('HC_offer_not_prev_att', patt*1.3)
M, S, N, P, r = runmodels(Hsens, H2, H, M, S, N, P, r)
H2.ReleaseMemory()
Hsens.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('HC_statins_presc_Q20plus', 0.36)  # was 14%
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('HC_aht_presc_Q20plus', 0.06)  # was 0.0248
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_smoker_ref', 0.09) # was 0.036
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_weight_ref', 0.6875) # was 0.275
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

H2 = initmodels()
H2.SetUncertainParameter('HC_statins_presc_Q20minus', 0.05) # was 2%
H2.SetUncertainParameter('HC_statins_presc_Q20plus', 0.36)  # was 14%
H2.SetUncertainParameter('HC_aht_presc_Q20minus', 0.04) # was 0.0154
H2.SetUncertainParameter('HC_aht_presc_Q20plus', 0.06)  # was 0.0248
H2.SetUncertainParameter('HC_smoker_ref', 0.09) # was 0.036
H2.SetUncertainParameter('HC_weight_ref', 0.6875) # was 0.275
M, S, N, P, r = runmodels(H1, H2, H, M, S, N, P, r)
H2.ReleaseMemory()

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
H2.ReleaseMemory()


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
H1.ReleaseMemory()


    
# Sensitivity analysis with CVD incidence and case fatality declining in the future 
# Compare no HC with current practice

Hsens0 = initmodels(HealthChecks=False)
Hsens0.SetUncertainParameter('CVD_annual_inc_rr_male', 0.9567661) # from geom ave of IHD and stroke figures 
Hsens0.SetUncertainParameter('CVD_annual_inc_rr_female', 0.9582724)  # from geom ave of IHD and stroke figures
Hsens0.SetUncertainParameter('IHD_annual_cf_rr_male', 0.964)
Hsens0.SetUncertainParameter('IHD_annual_cf_rr_female', 0.958)
Hsens0.SetUncertainParameter('Stroke_annual_sudden_rr_male', 0.9397142) # (0.12/0.21)^{1/9}
Hsens0.SetUncertainParameter('Stroke_annual_sudden_rr_female', 0.9397142)
Hsens0.Run()

Hsens1 = initmodels()
Hsens1.SetUncertainParameter('CVD_annual_inc_rr_male', 0.9567661) # from geom ave of IHD and stroke figures 
Hsens1.SetUncertainParameter('CVD_annual_inc_rr_female', 0.9582724)  # from geom ave of IHD and stroke figures
Hsens1.SetUncertainParameter('IHD_annual_cf_rr_male', 0.964)
Hsens1.SetUncertainParameter('IHD_annual_cf_rr_female', 0.958)
Hsens1.SetUncertainParameter('Stroke_annual_sudden_rr_male', 0.9397142) # (0.12/0.21)^{1/9}
Hsens1.SetUncertainParameter('Stroke_annual_sudden_rr_female', 0.9397142)
M, S, N, P, r = runmodels(Hsens0, Hsens1, Hsens0, M, S, N, P, r)

# Sensitivity analysis with CVD incidence and case fatality declining in the future 
# Compare current practice with extra eligibility, uptake, treatment

Hsens2 = initmodels()
Hsens2.SetUncertainParameter('CVD_annual_inc_rr_male', 0.9567661) # from geom ave of IHD and stroke figures 
Hsens2.SetUncertainParameter('CVD_annual_inc_rr_female', 0.9582724)  # from geom ave of IHD and stroke figures
Hsens2.SetUncertainParameter('IHD_annual_cf_rr_male', 0.964)
Hsens2.SetUncertainParameter('IHD_annual_cf_rr_female', 0.958)
Hsens2.SetUncertainParameter('Stroke_annual_sudden_rr_male', 0.9397142) # (0.12/0.21)^{1/9}
Hsens2.SetUncertainParameter('Stroke_annual_sudden_rr_female', 0.9397142)
# eligibility
Hsens2.SetUncertainParameter('HC_include_bp_registers', 1)
# uptake
Hsens2.SetUncertainParameter('HC_takeup', Hsens2.up_HC_takeup*1.3)
# treatment
Hsens2.SetUncertainParameter('HC_statins_presc_Q20minus', 0.05) # was 2%
Hsens2.SetUncertainParameter('HC_statins_presc_Q20plus', 0.36)  # was 14%
Hsens2.SetUncertainParameter('HC_aht_presc_Q20minus', 0.04) # was 0.0154
Hsens2.SetUncertainParameter('HC_aht_presc_Q20plus', 0.06)  # was 0.0248
Hsens2.SetUncertainParameter('HC_smoker_ref', 0.09) # was 0.036
Hsens2.SetUncertainParameter('HC_weight_ref', 0.6875) # was 0.275
Hsens2.SetUncertainParameter('HC_age_limit', [40, 79.9]) # was 40,74
M, S, N, P, r = runmodels(Hsens1, Hsens2, Hsens0, M, S, N, P, r)

Hsens0.ReleaseMemory()
Hsens1.ReleaseMemory()
Hsens2.ReleaseMemory()



# Sensitivity analysis with extra uncertainty around baseline risk of population
# Compare no HC vs HC

Hsens0 = initmodels(HealthChecks=False)
Hsens0.SetUncertainHyperparameter('QRisk_extra_logrr', 'std', 0.1)
Hsens0.ChangeUncertainParameters()
Hsens0.Run()

Hsens1 = initmodels()
Hsens1.SetUncertainHyperparameter('QRisk_extra_logrr', 'std', 0.1)
Hsens1.ChangeUncertainParameters()
M, S, N, P, r = runmodels(Hsens0, Hsens1, Hsens0, M, S, N, P, r)

# Compare HC base case vs increased eligibility, uptake, treatment

Hsens2 = initmodels()
Hsens2.SetUncertainHyperparameter('QRisk_extra_logrr', 'std', 0.1)
Hsens2.ChangeUncertainParameters()
# eligibility
Hsens2.SetUncertainParameter('HC_include_bp_registers', 1)
# uptake
Hsens2.SetUncertainParameter('HC_takeup', Hsens2.up_HC_takeup*1.3)
# treatment
Hsens2.SetUncertainParameter('HC_statins_presc_Q20minus', 0.05) # was 2%
Hsens2.SetUncertainParameter('HC_statins_presc_Q20plus', 0.36)  # was 14%
Hsens2.SetUncertainParameter('HC_aht_presc_Q20minus', 0.04) # was 0.0154
Hsens2.SetUncertainParameter('HC_aht_presc_Q20plus', 0.06)  # was 0.0248
Hsens2.SetUncertainParameter('HC_smoker_ref', 0.09) # was 0.036
Hsens2.SetUncertainParameter('HC_weight_ref', 0.6875) # was 0.275
Hsens2.SetUncertainParameter('HC_age_limit', [40, 79.9]) # was 40,74
M, S, N, P, r = runmodels(Hsens1, Hsens2, Hsens0, M, S, N, P, r)

Hsens1.ReleaseMemory()
Hsens2.ReleaseMemory()

SaveResults(M, S, N, P, nsc, npops, nouts, prefix)

# assert r == nsc


# save parameter values for current iteration (sensitivity analysis with extra uncertainty on QRisk) 

UP = Hsens.GetUncertainParameters()
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

with open("%s_senspars.csv" % (prefix), 'a') as f:
    if (run==1):
        np.savetxt(f, namsrep1, delimiter=",", fmt="%s")
    pars = np.hstack(par_vals)
    pars = pars.reshape((1, pars.size))
    np.savetxt(f, pars, delimiter=",")

Hsens.ReleaseMemory()
