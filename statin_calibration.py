## Calibrate the additional effect of statins on top of cholesterol reduction to achieve the CVD event rate reduction observed in trial data 
## Cholesterol Treatment Trialists' (CTT) Collaboration, Lancet 2015, "Efficacy and safety of LDL-lowering therapy among men and women: meta-analysis of individual data from 174000 participants in 27 randomised trials"
# "The proportional reductions per 1.0 mmol/L reduction in LDL cholesterol in major vascular events were similar overall for women (rate ratio [RR] 0.84, 99% CI 0.78-0.91) and men (RR 0.78, 99% CI 0.75-0.81)"

from __future__ import division # do float instead of integer division
import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 250000
n_cpus = 4
randpars = False
baseage_min = 30 # trial populations weren't restricted to narrow age range
baseage_max = 75

if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0

H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, randpars=randpars, baseage_min=baseage_min, baseage_max=baseage_max)
H.Run()

H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, randpars=randpars, baseage_min=baseage_min, baseage_max=baseage_max)
## boost treated population for efficiency: assume full prescription and initial adherence, but assume existing dropout rate same as trials
H1.SetUncertainParameter('HC_statins_presc_Q20minus', 0.0205 / 0.1423) 
H1.SetUncertainParameter('HC_statins_presc_Q20plus', 1.0)
H1.SetUncertainParameter('Statins_comp', 1.0)
# Parameters to calibrate 
H1.SetUncertainParameter('Statins_eff_extra_male', 0.87)
H1.SetUncertainParameter('Statins_eff_extra_female', 1.02)
H1.Run()

allpop = np.ones(H1.population_size, dtype=bool)
HC1 = gr.FirstHC(H1) # first HC time for all
stat_firsthc = (H1.Statins[allpop,HC1]>0) # are statins given at first HC

for i in (1,0):   # men, women
    stat = stat_firsthc * (H1.gender==i) # men or women who get statins at first HC
    HC1i = HC1[stat] # first HC time for this subgroup
    cvdevents_nohc = cvdevents_hc = 0
    for j in range(5):  # sum CVD events over first five years after HC
        cvdevents_nohc += H.CVD_events[stat, HC1i + j].sum()
        cvdevents_hc += H1.CVD_events[stat, HC1i + j].sum()
    rr = cvdevents_hc / cvdevents_nohc
    print stat.sum(), cvdevents_nohc, cvdevents_hc, rr

#Extra HRs to add to QRisk: stored in HC_parameters.py
#P['Statins_eff_extra_male'] =  0.82
#P['Statins_eff_extra_female'] = 0.87
# lead to: 
# RR 0.84, 0.78 for men, women.

## Uncertainty interval, to match the CIs in trials 0.78 (0.75-0.81), 0.84 (0.78-0.91)
# pars 0.78, 0.8  gets RRs 0.75, 0.78
# pars 0.87, 1.02 gets 0.81, 0.91
# giving SE 0.03, 0.08 for log HRs, using upper limit 
