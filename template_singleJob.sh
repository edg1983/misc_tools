#!/bin/bash
# author Edoardo Giacopuzzi
# template single job submission

#$ -q queue [long.qc/short.qc]
#$ -cwd -j y
#$ -o [log_file/log_folder]
#$ -N job_name
#$ -pe shmem N_threads

echo "started on `hostname` at `date`"

#Commands here

echo "Stop " $(date)
exit 0;