#!/bin/bash
# author Edoardo Giacopuzzi
# template array job submission

#$ -q queue [long.qc/short.qc]
#$ -cwd -j y
#$ -o [log_file/log_folder]
#$ -t start-stop:by [1-10:1]
#$ -N job_name
#$ -pe shmem N_threads

#Example to process files from a list
inputfile=$(sed -n "${SGE_TASK_ID}p" ${filelist})
basename=${inputfile##*/}
prefix=${basename%.bam}

echo "started on `hostname` at `date`"

#Commands here
# $SGE_TASK_ID = array job number

echo "Stop " $(date)
exit 0;