import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 25000
n_cpus = 4
prefix = "results/paperfeb/unc"
randpars = True

## used for running a large population in a set of batches
## this program is called by a bash script with the batch number as a command line argument
## results are saved to CSV files with the run/batch number in the file name
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

if (len(sys.argv) > 2):
    n_cpus = int(sys.argv[2])

nsc = 1
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
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, baseage_min=baseage_min, baseage_max=baseage_max)
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
    SaveResults(M, S, N, P, nsc, npops, nouts, prefix)
    r += 1
    return M, S, N, P, r

r = 0

H1 = initmodels()
M, S, N, P, r = runmodels(H, H1, H, M, S, N, P, r)
