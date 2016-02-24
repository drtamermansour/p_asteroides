#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=8,mem=64Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N blastn		#give name to the job


module load BLAST+/2.2.30
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB

cd $PBS_O_WORKDIR

blastx -query $input \
       -db nr \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp stitle salltitles staxids sskingdoms sscinames" \
       -num_threads 8 \
       -max_target_seqs 4  -out $input.bx ## -perc_identity 90 -qcov_hsp_perc 50

qstat -f ${PBS_JOBID}

