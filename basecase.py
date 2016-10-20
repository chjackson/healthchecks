import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 25000
n_cpus = 4
prefix = "results/unc/"
randpars = True

## used for running a large population in a set of batches
## this program is called by a bash script with the batch number as a command line argument
## results are saved to CSV files with the run/batch number in the file name
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

if (len(sys.argv) > 2):
    n_cpus = int(sys.argv[2])

nsc = 1 # number of scenarios
npops = 8 # number of subpopulations
nouts = 6 # number of outputs (LY, QALY etc)

M = np.zeros((nsc, npops, nouts))
S = np.zeros((nsc, npops, nouts))
N = np.zeros((nsc, 9), dtype=int) # subpopulation sizes and other count data in main run
N[...,0] = ps
NT = np.zeros((nsc, 4), dtype=int) # subpop sizes in artificial treated subset runs.
P = np.zeros((nsc, 36), dtype=object)
ST = np.zeros((nsc, 16, 41), dtype=object)

def SaveResults(M, S, N, NT, P, ST, nsc, npops, nouts, prefix):
    np.savetxt("%s_mean_%s.csv" % (prefix,run), M.reshape((nsc,npops*nouts)) , delimiter=",") # rows are populations (eligible, treated...), cols are outputs (LY, QALY...)
    np.savetxt("%s_SD_%s.csv" % (prefix,run), S.reshape((nsc,npops*nouts)), delimiter=",") #
    np.savetxt("%s_N_%s.csv" % (prefix,run), N, fmt="%d", delimiter=",") #
    np.savetxt("%s_NT_%s.csv" % (prefix,run), NT, fmt="%d", delimiter=",") #
    np.savetxt("%s_process_%s.csv" % (prefix,run), P, delimiter=",") #
    np.savetxt("%s_short_%s.csv" % (prefix,run), ST.reshape((nsc*16,41)), fmt="%s", delimiter=",") #
    
## Without HC
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, randpars=randpars)
H.Run()
## Basic HC model

if (randpars):
    UP = H.GetUncertainParameters()
else: UP = None

def initmodels():
    H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, pars=UP)

    ## Different models to boost size of treated subsets
    ## Increases the precision of the estimates of post-treatment
    ## outcomes without biasing them, since treatment allocation given
    ## risk group is random

    HS = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, pars=UP)
    HS.SetUncertainParameter('up_HC_statins_presc_Q20plus', 1.0)
    HS.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.0205 / 0.1423)
    HS.SetUncertainParameter('up_Statins_comp', 1.0)
    HA = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, pars=UP)
    HA.SetUncertainParameter('up_HC_aht_presc_Q20plus', 1.0)
    HA.SetUncertainParameter('up_HC_aht_presc_Q20minus', 0.0154 / 0.0248)
    HA.SetUncertainParameter('up_AHT_comp', 1.0)
    HW = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, pars=UP)
    HW.SetUncertainParameter('up_HC_weight_ref', 1.0)
    HW.SetUncertainParameter('up_Weight_comp', 1.0)
    HC = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, pars=UP)
    HC.SetUncertainParameter('up_HC_smoker_ref', 1.0)
    return H1, HS, HA, HW, HC

r = 0

## Base case scenario
H1, HS, HA, HW, HC = initmodels()
for h in [H1, HS, HA, HW, HC]:
    h.Run()
M[r,],S[r,],N[r,],NT[r,],P[r,],ST[r,] = gr.GetResults_allSub(H, H1, HS, HA, HW, HC, M[r,], S[r,], N[r,], NT[r,], P[r,], ST[r,])
SaveResults(M, S, N, NT, P, ST, nsc, npops, nouts, prefix)
r += 1
