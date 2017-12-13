## Validate CVD mortality from the model against ONS data 2011

import os
import sys
import HC_main as hc
import numpy as np
import getresults as gr

st = 70
ps = 2000
n_cpus = 4

if (len(sys.argv) > 1):
    run = int(sys.argv[1])
else: run = 0
    
baseage_min = 40  # Restrict age range of baseline population
baseage_max = 45 
    
## Without HC
H = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=False, nprocs=n_cpus, randseed=run, baseage_min=baseage_min, baseage_max=baseage_max)
#H.SetUncertainParameter('MI_sudden_death', 0.2)
#H.SetUncertainParameter('Stroke_sudden_death', 0.2)
H.Run()

H1 = hc.HealthChecksModel(population_size=ps, simulation_time=st, HealthChecks=True, nprocs=n_cpus, randseed=run, baseage_min=baseage_min, baseage_max=baseage_max)
#H1.SetUncertainParameter('MI_sudden_death', 0.2)
#H1.SetUncertainParameter('Stroke_sudden_death', 0.2)
H1.Run()

print (H1.QALY.mean() - H.QALY.mean()) * 365.25

# 3.7434243427 with 0.2 
# 4.00322678021 with 0.3, ie bigger impacts

maxage = 45 + H.simulation_time
atrisk_male = np.zeros((H.simulation_time, maxage))
atrisk_female = np.zeros((H.simulation_time, maxage))
died_ihd_male  = np.zeros((H.simulation_time, maxage))
died_stroke_male  = np.zeros((H.simulation_time, maxage))
died_ihd_female  = np.zeros((H.simulation_time, maxage))
died_stroke_female  = np.zeros((H.simulation_time, maxage))

for i in range(0, H.simulation_time-1):  # loop over times 
    print i
    for j in range(maxage):  # loop over ages (in practice starts at 40, but keep general)
        agegrp = (H.age[:,i]==j)
        male = (H.gender[agegrp] == 1)
        female = (H.gender[agegrp] == 0)
        # number of people of age j still alive at start of time 0
        if (i==0): 
            atrisk_male[i, j] += sum(male)
            atrisk_female[i, j] += sum(female)
        # number of people of age j still alive at end of time i
        else:
            atrisk_male[i+1, j] += sum(H.alive[agegrp, i] * male)
            atrisk_female[i+1, j] += sum(H.alive[agegrp, i] * female)
        # number of people of age j dying from cause at end of time i
        died_ihd_male[i, j] += sum(H.Death[agegrp, i] * (H.CauseOfDeath[agegrp] == 'IHD') * male)
        died_stroke_male[i, j] += sum(H.Death[agegrp, i] * (H.CauseOfDeath[agegrp] == 'Stroke') * male)
        died_ihd_female[i, j] += sum(H.Death[agegrp, i] * (H.CauseOfDeath[agegrp] == 'IHD') * female)
        died_stroke_female[i, j] += sum(H.Death[agegrp, i] * (H.CauseOfDeath[agegrp] == 'Stroke') * female)
    
np.savetxt("mort_valid_curr/atrisk_male_%s_%s.csv" % run, atrisk_male, delimiter=",")
np.savetxt("mort_valid_curr/atrisk_female_%s.csv" % run, atrisk_female, delimiter=",")
np.savetxt("mort_valid_curr/died_ihd_male_%s.csv" % run, died_ihd_male, delimiter=",")
np.savetxt("mort_valid_curr/died_stroke_male_%s.csv" % run, died_stroke_male, delimiter=",")
np.savetxt("mort_valid_curr/died_ihd_female_%s.csv" % run, died_ihd_female, delimiter=",")
np.savetxt("mort_valid_curr/died_stroke_female_%s.csv" % run, died_stroke_female, delimiter=",")

## fishy.
## should be lower mortality for stroke compared to 0.2 scenario
## currently getting results close to 0.3 ones 
