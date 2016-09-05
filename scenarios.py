# Script to produce results for all scenarios
# Produces CSVs which can then be manipulated in scenarios.r to get graphs or tables

# run time per scenario
# 3 minutes on fast desktop, 4 cores, pop of 250000
# 20 minutes on cluster for pop of 1000000

### TODO
### Run and present graph for QALYs gained per HC
### Effects in particular treated populations 

import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 250000
n_cpus = 4
prefix = "scen"    

nsc = 20 # number of scenarios
npops = 4 # number of subpopulations
nouts = 6 # number of outputs (LY, QALY etc)

M = np.zeros((nsc, npops, nouts))
S = np.zeros((nsc, npops, nouts))
N = np.zeros((nsc, 9), dtype=int) # subpopulation sizes and other count data in main run
N[...,0] = ps
# NT = np.zeros((nsc, 4), dtype=int) # subpop sizes in artificial treated subset runs.
# treated subset results not presented for the moment

def SaveResults(M, S, N, nsc, npops, nouts, prefix):
    np.savetxt("results/%s_mean.csv" % prefix, M.reshape((nsc,npops*nouts)) , delimiter=",") # rows are populations (eligible, treated...), cols are outputs (LY, QALY...)
    np.savetxt("results/%s_SD.csv" % prefix, S.reshape((nsc,npops*nouts)), delimiter=",") #
    np.savetxt("results/%s_N.csv" % prefix, N, fmt="%d", delimiter=",") #

    
## Without HC
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus)
H.Run()

r = 0

## With HC
## Include people with CVD and on diabetes registers
## First one here is the base case for all scenarios - don't include any people on registers. r=0
comb = 8 # number of combinations for settings (2^3)
reg_array = np.zeros((comb,3),dtype=bool)
c = 0
for i in [False,True]:
    for j in [False,True]:
        for k in [False,True]:
            reg_array[c,0] = i
            reg_array[c,1] = j
            reg_array[c,2] = k
            c+=1
for i in range(comb):
    H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
    H1.SetUncertainParameter('up_HC_include_diabetes_registers', reg_array[i,0])
    H1.SetUncertainParameter('up_HC_include_CVD_registers', reg_array[i,1])
    H1.SetUncertainParameter('up_HC_include_bp_registers', reg_array[i,2])
    H1.Run()
    M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
    SaveResults(M, S, N, nsc, npops, nouts, prefix)
    r += 1

## Change age threshold for eligibility
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
H1.SetUncertainParameter('up_HC_age_limit', [50, 74])
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

## Overall uptake rate increased to 'best practice'
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
H1.SetUncertainParameter('up_HC_takeup', 0.85)
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

## Increased uptake of people in deprived areas
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
ses = H1.GetUncertainParameters()['up_SES_vec']
H1.SetUncertainParameter('up_SES_vec', [ses[0], ses[1], ses[2], ses[3]*1.1, ses[4]*1.2])
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

## Increased uptake of smokers
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
H1.SetUncertainParameter('up_smoker_vec', [1, 1])
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

## Increased uptake of high risk individuals
## TODO 20% or 30% here?
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
qr = H1.GetUncertainParameters()['up_QRisk_vec']
H1.SetUncertainParameter('up_QRisk_vec', [qr[0], qr[1], qr[2], qr[3]*1.2, qr[4]*1.2])
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

## Increased treatment rates to 'best practice'
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
H1.SetUncertainParameter('up_HC_statins_presc_Q20minus', 0.05) # was 2%
H1.SetUncertainParameter('up_HC_statins_presc_Q20plus', 0.36)  # was 14%
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

# Subsequent attendance more likely if attended before, less likely if not
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
H1.SetUncertainParameter('up_HC_takeup_prev_att', 0.7)
H1.SetUncertainParameter('up_HC_takeup_not_prev_att', 0.3)
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

# Subsequent attendance less likely if attended before, target never-attenders
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
H1.SetUncertainParameter('up_HC_takeup_prev_att', 0.3)
H1.SetUncertainParameter('up_HC_takeup_not_prev_att', 0.7)
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

# No extra effect of statins
H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
H1.SetUncertainParameter('up_Statins_eff_extra_male', 1)
H1.SetUncertainParameter('up_Statins_eff_extra_female', 1)
H1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(H, H1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

# Sudden death from CVD events, no background mortality reduction
HS0 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus)
HS0.SetUncertainParameter('up_cvd_sudden_death', 0.1)
HS0.Run()
HS1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
HS1.SetUncertainParameter('up_cvd_sudden_death', 0.1)
HS1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(HS0, HS1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

# Sudden death from CVD events, with no background mortality reduction
HS0 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus)
HS0.SetUncertainParameter('up_cvd_sudden_death', 0.1)
HS0.SetUncertainParameter('up_cvd_background_cfr_reduction', 0.3)
HS0.Run()
HS1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
HS1.SetUncertainParameter('up_cvd_sudden_death', 0.1)
HS1.SetUncertainParameter('up_cvd_background_cfr_reduction', 0.3)
HS1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(HS0, HS1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1

# Sudden death from CVD events, with no background mortality reduction
HS0 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus)
HS0.SetUncertainParameter('up_cvd_sudden_death', 0.4)
HS0.Run()
HS1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus)
HS1.SetUncertainParameter('up_cvd_sudden_death', 0.4)
HS1.Run()
M[r,],S[r,],N[r,] = gr.GetResults_all(HS0, HS1, M[r,], S[r,], N[r,])
SaveResults(M, S, N, nsc, npops, nouts, prefix)
r += 1


assert r == nsc

## These CSVs are manipulated in scenarios.r to get graphs or tables
