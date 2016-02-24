#!/bin/bash -login
#PBS -l walltime=72:00:00,nodes=1:ppn=8,mem=124Gb
#mdiag -A ged	
#PBS -m abe			#send email to myself
#PBS -N load-into-counting		#give name to the job

cd $PBS_O_WORKDIR
source $HOME/khmerEnv/bin/activate
module load GNU/4.7.1

load-into-counting.py -k 20 -N 4 -x 30e9 --threads 8 $graph $input_files

qstat -f ${PBS_JOBID}


