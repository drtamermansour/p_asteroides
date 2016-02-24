#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=5,mem=20Gb
#mdiag -A ged	
#PBS -m abe			#send email to myself
#PBS -N Ad_remove		#give name to the job


module load Trimmomatic/0.33 

output_dir=${output_dir}
cd ${output_dir}

temp=$(basename "$R1_INPUT")
new_dir=${temp%_R1_*}
						
mkdir ${new_dir}
cd ${new_dir}

java -jar $TRIM/trimmomatic PE -threads 4 -phred33 ${R1_INPUT} ${R2_INPUT} s1_pe.fq s1_se.fq s2_pe.fq s2_se.fq ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:25


qstat -f ${PBS_JOBID}

