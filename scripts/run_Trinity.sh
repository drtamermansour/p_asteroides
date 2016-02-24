#!/bin/bash
#PBS -l walltime=7:00:00:00,nodes=1:ppn=8,mem=360Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N trinity		#give name to the job

module load trinity/6.0.2
module load Java/1.7.0_51
module load SAMTools/0.1.19


cd $PBS_O_WORKDIR

Trinity --seqType fq --max_memory 350G --CPU 7 \
--left $(echo ${lf[*]} | tr ' ' ',') \
--right $(echo ${rt[*]} | tr ' ' ',')


qstat -f ${PBS_JOBID}

