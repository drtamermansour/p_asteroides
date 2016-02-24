#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=24Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N blastn		#give name to the job


module load BLAST+/2.2.30
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB

cd $PBS_O_WORKDIR

blastx -query $input \
       -db $DB \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp stitle" \
       -num_threads 1 \
       -max_target_seqs 10  -out $input.bx ## -perc_identity 90 -qcov_hsp_perc 50

qstat -f ${PBS_JOBID}

