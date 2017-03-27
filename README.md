# Health Checks Microsimulation Model

## Requirements

* Python 2.x

* Python libraries:
	- numpy (scientific computing library)
	- multiprocessing (for parallelisation)
	- random

## Simulating a population

* Main model file: [HC_main.py](HC_main.py)

```python
import HC_main as hc
simulation_time = 60 # number of years to simulate
n_cpus = 8           # number of CPU cores for parallelisation
H = hc.HealthChecksModel(population_size=250000,
                         simulation_time=simulation_time,
						 HealthChecks=False,
						 nprocs=n_cpus)
H.Run()
```

The python object `H` contains all simulated outcomes.

Examples of extracting aggregate results from the simulated population are given in `getresults.py`.

The memory of a typical desktop computer permits a maximum of about 250000 people per model run.

Shell scripts are therefore used to simulate larger populations in batches.

## Scripts to reproduce the analyses in the paper 

###  Base case only 

* Shell script [basecase.sh](basecase.sh), which calls the Python script [basecase.py](basecase.py).

###  All scenarios

* Shell script [scenarios.sh](scenarios.sh), which calls the Python script [scenarios.py](scenarios.py).

###  All scenarios, with statistical uncertainty

* Shell script [unc.sh](unc.sh), which calls the Python script [scenarios_unc.py](scenarios_unc.py).  Requires the SLURM high performance computing facility.

###  Key constants 

* `nruns` in the shell scripts specifies the number of batches
* `ps` in the Python scripts gives the number of people in each batch.
* The total population size is then `nruns*ps`.

###  Formatting results 

* The Python library [getresults.py](getresults.py) specifies which aggregate results are extracted from the simulated population. 

* The Python scripts output aggregate results from the simulated population in text files.

* The R code in the file [tables.r](tables.r) formats the information in these output files to produce the tables in the paper.
