#!/bin/bash

start=`date +%s`

nruns=10

for i in `seq 1 $nruns`;
do
    echo "STARTING RUN $i..."
    python basecase.py $i
done

end=`date +%s`
runtime=$((end-start))
echo "RUN TIME: $runtime SECONDS"
