#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N extract-PE		#give name to the job

cd $PBS_O_WORKDIR
source $HOME/khmerEnv/bin/activate
module load GNU/4.7.1

extract-paired-reads.py $input

qstat -f ${PBS_JOBID}

