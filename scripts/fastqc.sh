#!/bin/bash -login
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=1Gb

#PBS -m abe		#send email to myself
#PBS -N T_FastQC	#give name to the job

module load FastQC/0.11.3

cd $PBS_O_WORKDIR

fastqc -f fastq -noextract ${INPUT_FILE}

qstat -f ${PBS_JOBID}

