#!/bin/bash

# Specify a job name
#$ -N ${jobname}

# --- Parameters for the Queue Master ---
# Project name and target queue
#$ -P brc.prjc
#$ -q ${queue}
#$ -cwd -j y
#$ -pe shmem ${cores}
# Print some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Job ID: $SGE_JOB_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"
module load Anaconda3
source /apps/eb/ivybridge/software/Anaconda3/5.3.0/etc/profile.d/conda.sh
conda activate biochem

python ${py_filename} ${py_args}
exit $?