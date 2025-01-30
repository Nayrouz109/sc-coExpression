


#!/bin/bash
# Submit SLURM jobs sequentially with dependency

PREV_JOBID=""
while read -r FILE; do
    if [[ -z "$PREV_JOBID" ]]; then
        # Submit the first job without dependency
        JOBID=$(sbatch /home/nelazzabi/rewiring/scripts/analysis/runCor-single.sh "$FILE" | awk '{print $4}')
    else
        # Submit the next job with dependency on the previous
        JOBID=$(sbatch --dependency=afterok:$PREV_JOBID /home/nelazzabi/rewiring/scripts/analysis/runCor-single.sh "$FILE" | awk '{print $4}')
    fi
    PREV_JOBID=$JOBID
done < file_list.txt

