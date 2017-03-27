#!/bin/bash

# For single-server use, not cluster.  For cluster see scenarios-hpc.sh

start=`date +%s`

## number of population batches of size 25000 to generate
nruns=40

for i in `seq 1 $nruns`;
do
    echo "STARTING RUN $i..."
    python scenarios.py $i
done

end=`date +%s`
runtime=$((end-start))
echo "RUN TIME: $runtime SECONDS"
