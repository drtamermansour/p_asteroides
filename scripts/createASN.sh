#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=64Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N create_ASN		#give name to the job


module load tbl2asn/20150331


cd $PBS_O_WORKDIR

tbl2asn -t template.sbt -p. -Y assembly.cmt -M t -j "[organism=Porites astreoides] [moltype=transcribed_RNA] [tech=TSA]"
#tbl2asn -t template.sbt -p. -M t -n "Porites astreoides"

qstat -f ${PBS_JOBID}

