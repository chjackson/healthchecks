#!/bin/bash

# For single-server use, not cluster.  For cluster see scenarios-hpc.sh

start=`date +%s`

nruns=5

for i in `seq 1 $nruns`;
do
    echo "STARTING RUN $i..."
    python scenariosub_fixage.py $i
done

end=`date +%s`
runtime=$((end-start))
echo "RUN TIME: $runtime SECONDS"
