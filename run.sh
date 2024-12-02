#!/bin/bash -l
#$ -e ./logs/
#$ -o ./logs/

set -eux

user="dengel"

while IFS=$'\t' read -r v1 v2 v3 v4 v5 v6; do

    if [[ "$v1" == "quant.p" ]]; then
        continue
    fi

    # Check if the variables contain tab characters and remove them if found
    if [[ "$v1" == *$'\t'* || "$v2" == *$'\t'* || "$v3" == *$'\t'* || "$v4" == *$'\t'* ]]; then
        v1=$(echo "$v1" | tr -d '\t')
        v2=$(echo "$v2" | tr -d '\t')
        v3=$(echo "$v3" | tr -d '\t')
        v4=$(echo "$v4" | tr -d '\t')
        v5=$(echo "$v5" | tr -d '\t')
        v6=$(echo "$v6" | tr -d '\t')
    fi

    # Submit jobs based on $v3 value
    if [[ "$v3" == "raw" ]]; then
        for i in {1..5}; do
            echo "Submitting qsub for Raw data (run $i): v1=$v1, v2=$v2, v3=$v3, v4=$v4 v5=${v5//,/} v6=$v6"
            #qsub -N "bigscale2_DEV3_${v1}_${v2}_${v3}_${v5//,/}_run${i}_$(date +%m-%d-%y_%I-%M-%S%p)" bigscale2_DEVer4.qsub $v1 $v2 $v3 $v4 $v5 $v6 $i
        done
    else
        echo "Submitting qsub for Norm data: v1=$v1, v2=$v2, v3=$v3, v4=$v4 v5=${v5//,/} v6=$v6"
        #qsub -N "bigscale2_DEV3_${v1}_${v2}_${v3}_${v5//,/}_$(date +%m-%d-%y_%I-%M-%S%p)" bigscale2_DEVer4.qsub $v1 $v2 $v3 $v4 $v5 $v6
    fi

done <config.tsv

sleep 10

# Wait for all jobs to complete
echo "Waiting for jobs to complete..."
while [[ $(qstat -u $user | grep -c "bigscale2_") -gt 0 ]]; do
    echo "Jobs still in the queue. Sleeping for 60 seconds..."
    qstat -u $user  # Log the output of qstat to see the current job list
    sleep 10  # Check every 60 seconds
done
echo "All jobs completed."

bash run_aggr.sh



