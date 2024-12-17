#!/bin/bash -l
#$ -e ./logs/
#$ -o ./logs/

set -eux

user="dengel"
calc_model="true"

mainoutputdir="/projectnb/wax-dk/david/AllBS2Output"
echo "All Data is saved in '$mainoutputdir'."

mkdir -p "$mainoutputdir" # make general output directory.

if [ $? -eq 0 ]; then
  echo "Folder '$mainoutputdir' created successfully."

elif [ -d "$mainoutputdir"]; then
  echo "Directory '$mainoutputdir' already exists."

else
  echo "Failed to create folder '$mainoutputdir'."
  exit 1  # Terminate the script with an error code
fi

Home_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script is running from: $Home_dir"

models="$Home_dir/model"
echo "Models will be saved in: '$models'"

date_run=$(date +%m-%d-%y_%I-%M-%S%p)
echo "Run started at $date_run"

if [[ "$calc_model" = true ]]; then
  mkdir -p "$models"


  if [ $? -eq 0 ]; then
    echo "Folder '$models' created successfully."

  elif [ -d "$models"]; then
    echo "Directory '$models' already exists."

  else
    echo "Failed to create folder '$models'."
  exit 1  # Terminate the script with an error code
  fi
    echo "calc.model is true. Running R script to calculate network model..."
    bash run_calcmodel.sh $Home_dir $models $date_run
    
    # Wait for the calcmodel job to complete
    while qstat -u $user | grep -q "calcmodel."; do
      job_name=$(qstat -u $user | grep "calcmodel." | awk '{print $3}')
      echo "Waiting for calcmodel job ($job_name) to finish..."
      sleep 300  # Check every 5 minutes
    done

    echo "calcmodel job completed. Proceeding with the network creation."
    
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
            echo $i
            echo "Submitting qsub for Raw data (run $i): v1=$v1, v2=$v2, v3=$v3, v4=$v4 v5=${v5//,/} v6=$v6"
            qsub -N "bigscale2_DE_${v1}_${v2}_${v3}_${v5//,/}_run${i}_$date_run" bigscale2_DEVer5.qsub $v1 $v2 $v3 $v4 $v5 $v6 $i $mainoutputdir $models $date_run
        done
    else
        echo "Submitting qsub for Norm data: v1=$v1, v2=$v2, v3=$v3, v4=$v4 v5=${v5//,/} v6=$v6"
        qsub -N "bigscale2_DE_${v1}_${v2}_${v3}_${v5//,/}_$date_run" bigscale2_DEVer5.qsub $v1 $v2 $v3 $v4 $v5 $v6 $mainoutputdir $models $date_run
    fi

    done <config.tsv
    
else
    echo "calc.model is false. Network model will not be precalculated."
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
            echo $i
            echo "Submitting qsub for Raw data (run $i): v1=$v1, v2=$v2, v3=$v3, v4=$v4 v5=${v5//,/} v6=$v6"
            qsub -N "bigscale2_DE_${v1}_${v2}_${v3}_${v5//,/}_run${i}_$date_run" bigscale2_DEVer4.qsub $v1 $v2 $v3 $v4 $v5 $v6 $i $mainoutputdir $date_run
        done
    else
        echo "Submitting qsub for Norm data: v1=$v1, v2=$v2, v3=$v3, v4=$v4 v5=${v5//,/} v6=$v6"
        qsub -N "bigscale2_DE_${v1}_${v2}_${v3}_${v5//,/}_$date_run" bigscale2_DEVer4.qsub $v1 $v2 $v3 $v4 $v5 $v6 $mainoutputdir $date_run
    fi

    done <config.tsv
fi





# Wait for all jobs to complete
echo "Waiting for jobs to complete..."
while [[ $(qstat -u $user | grep -c "bigscale2_") -gt 0 ]]; do
    echo "Jobs still in the queue. Sleeping for 60 seconds..."
    qstat -u $user  # Log the output of qstat to see the current job list
    sleep 600  # Check every 60 seconds
done
echo "All jobs completed."

bash run_aggr.sh $date_run $mainoutputdir



