#!/bin/bash -login
#PBS -l walltime=72:00:00,nodes=1:ppn=8,mem=124Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N filter_abund		#give name to the job


cd $PBS_O_WORKDIR
source $HOME/khmerEnv/bin/activate
module load GNU/4.7.1

filter-abund.py --threads 8 -V $input $files

qstat -f ${PBS_JOBID}

