import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 250000
n_cpus = 4
prefix = "results/paper/base"
randpars = False

## used for running a large population in a set of batches
## this program is called by a bash script with the batch number as a command line argument
## results are saved to CSV files with the run/batch number in the file name
if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

if (len(sys.argv) > 2):
    n_cpus = int(sys.argv[2])

npops = 5 # number of subpopulations
npopsn = 12 # number of subpopulations to calculate size of
nouts = 38 # number of outputs (LY, QALY etc)

M = np.zeros((npops, nouts))
S = np.zeros((npops, nouts))
N = np.zeros(npopsn, dtype=int) # subpopulation sizes and other count data in main run
N[0] = ps

def SaveResults(M, S, N, npops, nouts, prefix):
    np.savetxt("%s_mean_%s.csv" % (prefix,run), M.reshape((1,npops*nouts)) , delimiter=",") # rows are populations (eligible, treated...), cols are outputs (LY, QALY...)
    np.savetxt("%s_SD_%s.csv" % (prefix,run), S.reshape((1,npops*nouts)), delimiter=",") #
    np.savetxt("%s_N_%s.csv" % (prefix,run), N, fmt="%d", delimiter=",") #
    
## Without HC
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, randpars=randpars, baseage_min=40, baseage_max=45)
H.Run()

if (randpars):
    UP = H.GetUncertainParameters()
else: UP = None

## With HC, current practice
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, pars=UP, baseage_min=40, baseage_max=45)
H1.Run()

### TABLE 2 baseline demographics
## Do this with one run of the model, as big a population as will fit in memory.
## Allows calculation of IQRs, which is tricky with combining aggregate outcomes from multiple runs. 
## Assume estimates of percentages/means precise enough with ps=250000 

B = gr.BaselineChars(H, H1)
np.savetxt("results/baseline.csv", B, delimiter=",")
M, S, N = gr.GetResults_paper(H, H1, M, S, N)
SaveResults(M, S, N, npops, nouts, prefix)
    
