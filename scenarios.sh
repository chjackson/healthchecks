#!/bin/bash
#
#SBATCH --job-name=HC_scenarios
#SBATCH --output=scenarios.out
#SBATCH -p mrc-bsu-sand
#SBATCH -A MRC-BSU-SL2
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00

start=`date +%s`

nruns=50

for i in `seq 1 $nruns`;
do
    echo "STARTING RUN $i..."
    python scenariosub.py $i
done

end=`date +%s`
runtime=$((end-start))
echo "RUN TIME: $runtime SECONDS"
