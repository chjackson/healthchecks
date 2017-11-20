#!/bin/bash
#
#SBATCH --job-name=scenarios_revision
#SBATCH --output=scenarios_revision.out
#SBATCH -p mrc-bsu-sand
#SBATCH -A MRC-BSU-SL2
#SBATCH --cpus-per-task=4
#SBATCH --nodes=5
#SBATCH --ntasks=20
#SBATCH --time=30:00:00

start=`date +%s`

npsa=100 # number of iterations of probabilistic sensitivity analysis 
n_cpus=4 # number of CPUs to use for within-run 
nruns=20 # number of model runs to do in parallel 
nbatch=$((npsa/nruns)) # number of parallel batches to run serially 

DATE=`date +%Y_%m_%d-%H_%M`
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

for j in `seq 1 $nbatch`;
do
    for i in `seq 1 $nruns`;
    do
	k=$(((j-1)*nruns + i))
	echo "Starting run $k"
	srun --exclusive -n 1 python scenarios_revision.py $k $n_cpus &
    done
    wait
done
## Remainder in case npsa is not divisible by n runs
nrem=$((npsa % nruns))
for i in `seq 1 $nrem`;
do
    k=$((nbatch*nruns + i))
    echo "Starting run $k"
    srun --exclusive -n 1 python scenarios_revision.py $k $n_cpus &
done
wait

end=`date +%s`
runtime=$((end-start))
echo "RUN TIME: $runtime SECONDS"



############################################################

### 16 cores is the max 
### how many nodes is the max? 
# On Darwin, SL1 and SL2 users are limited to 1024 cores on Sandy Bridge/Westmere in use at any one time and a maximum wallclock runtime of 36 hours per job. On Wilkes, SL1 and SL2 are limited to 288 cores (48 GPUs). 
# mrc-bsu-sand is MRC-BSU-SL2, so 
# 1024/16 = 64 nodes? 

# Flag cpus-per-task is for setting number of CPUs used in running one instance of the model. In effect, this is the variable n_cpus in the HC_main_xyz.py (and HC_results_bashloop.py) model file

# Flag nodes specifies how many nodes are used on the cluster. 

# Flag ntasks specifies how many tasks there are to run. In this case: How often will the script  HC_results_bashloop.py be run? Since each node has 16 cores, ntasks = nodes * (16/n_cpus).  This allows multiple tasks to get sent to the same node, if the CPUs are available.

# nruns and n_cpus are just redefinitions of nodes and ntasks. It would be nice to have ntasks and cpus_per_task be defined dynamically, but then the script would not work, as all #SBATCH definitions  must come first in the script before any other ‘proper’ bash script line is executed, because every #SBATCH after executed script lines are ignored.

# The --exclusive flag to srun prevents the job steps (i.e. the tasks launched by different sruns) sharing cpus. The & at the end of the line is essential to allow the steps to run concurrently rather than sequentially, but then you need the wait command at the end to stop the script, and therefore the whole job, from exiting before all tasks have returned.
