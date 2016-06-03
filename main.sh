#!/bin/sh

## construction of the basic diretory structure
git clone https://github.com/drtamermansour/p_asteroides.git
cd p_asteroides
p_asteroides=$(pwd)

## create working directory and define paths for raw data and scripts
mkdir -p $p_asteroides/{data,QC_raw,b_adap_remove}
script_path=$p_asteroides/scripts
data_path=$p_asteroides/data
QC_raw=$p_asteroides/QC_raw
trimmed_data=${p_asteroides}/b_adap_remove


##Data Striping on Scratch
lfs setstripe --count -1 ${data_path}
## copy the data to MSU-HPC
scp tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/*.* ${data_path}/. ## keep copy at home dir
chmod u-w ${data_path}/*.*

## calculation of input read counts
cd ${data_path}
for f in *.gz; do zcat $f | wc -l;done > rawInput_readCounts
for f in *.gz; do echo $f;done >> rawInput_readCounts
head -n16 rawInput_readCounts | awk '{ sum+=$1} END {print sum,sum/4,sum/8}' ## 393,511,591

## QC check for the original data files 
cd ${QC_raw}
for f in ${data_path}/*.gz; do qsub -v INPUT_FILE="$f" ${script_path}/fastqc.sh; done
mv ${data_path}/{*.zip,*.html} ${QC_raw}/.

#============================================================================================
## adaptor removal only 
cd ${trimmed_data}
bash ${script_path}/run_adapter_trimmer.sh ${data_path} ${trimmed_data} ${script_path}/adapter_removal.sh
#java -jar $TRIM/trimmomatic PE -threads 4 -phred33 ${R1_INPUT} ${R2_INPUT} s1_pe.fq s1_se.fq s2_pe.fq s2_se.fq ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:25

## combine single reads 
for f in ${trimmed_data}/*; do if [ -d "${f}" ]; then echo "${f}"; cat "${f}"/s1_se.fq "${f}"/s2_se.fq > "${f}"/$(basename "$f").s_se.fq; fi; done;

## QC check
#for dir in ${trimmed_data}/*; do cd $dir; for f in *_pe.fq *.s_se.fq; do qsub -v INPUT_FILE="$f" ${script_path}/fastqc.sh; done; done;
#============================================================================================
## prepare the files for subsequent analysis
for dir in ${trimmed_data}/*; do if [ -d "${dir}" ]; then echo $dir; cd $dir; qsub -v output=$(basename "$dir").s_pe.fq ${script_path}/interleave.sh; fi; done;

## filter abundance
mkdir ${p_asteroides}/c_abundFilter
cd ${p_asteroides}/c_abundFilter
mv ${trimmed_data}/*/*.s_[ps]e.fq .
input_files=()
for f in *.fq; do input_files+=($f); done
qsub -v graph=$"countgraph_k20.kt",input_files="${input_files[*]}" ${script_path}/load-into-counting.sh
#load-into-counting.py -k 20 -N 4 -x 30e9 --threads 8 $graph $input_files
qsub -v input=$"countgraph_k20.kt",files="${input_files[*]}" ${script_path}/filter_abund.sh
#filter-abund.py -V $input $files
mkdir trimmedData
mv *.s_[ps]e.fq trimmedData/.
for f in filter_abund.e*;do grep -B 2 "^output in" $f; done > filter_abund.summary;
## break out the orphaned and still-paired reads and rename files (this step ends with .s_pe.fq and .s_se.fq for each sample)
#for i in *.s_pe.*.abundfilt; do extract-paired-reads.py $i; done
for i in *.s_pe.*.abundfilt; do qsub -v input=$i $script_path/extract-paired-reads.sh; done
##  combine the orphaned reads into a single file & rename pe files
for i in *.s_se.fq.abundfilt; do
   pe_orphans=$(basename $i .s_se.fq.abundfilt).s_pe.fq.abundfilt.se
   cat $i $pe_orphans > $(basename $i .abundfilt)
done
#rm *.abundfilt.se
for i in *.abundfilt.pe; do mv $i $(basename $i .abundfilt.pe); done
## split interleaved files & ## merge the single reads to the end of left reads (this step reform the data into .s_pe1_se.fq & .s_pe2.fq for each sample )
for i in *.s_pe.fq; do qsub -v input=$i $script_path/split-paired-reads.sh; done
#for i in *.s_pe.fq; do split-paired-reads.py $i; done
for f in *.s_pe.fq; do
    echo $f  >> check_file_ends;
    base=$(basename $f);
    tail -n 8 $f | head -n 4 >> check_file_ends;
    tail -n 4 $base.1 >> check_file_ends;
    tail -n 4 $f >> check_file_ends;
    tail -n 4 $base.2 >> check_file_ends;
done
for f in *.s_pe.fq.1; do echo $f; base=$(basename $f .s_pe.fq.1); cat $f $base.s_se.fq > $base.s_pe1_se.fq; done
#rm *.s_pe.fq.1
for f in *.s_pe.fq.2; do mv $f $(basename $f .s_pe.fq.2).s_pe2.fq; done

## calculation of filtered read counts
for f in *.s_pe1_se.fq ; do wc -l $f;done > filtered_readCounts
cat filtered_readCounts | awk '{ sum+=$1} END {print sum,sum/4}' ## 393,511,591

## Trinity
lf_files=()
for f in *.s_pe1_se.fq; do lf_files+=($f); done;
rt_files=()
for f in *.s_pe2.fq; do rt_files+=($f); done;
qsub -v lf="${lf_files[*]}",rt="${rt_files[*]}" ${script_path}/run_Trinity.sh
echo $(grep "^>" Trinity.fasta | wc -l) ## 881402

## SeqClean (trim polyA tails and remove dust "low complexity seq")
module load SeqClean/20130718
seqclean_dir=$p_asteroides/c_abundFilter/trinity_out_dir/seqclean
mkdir $seqclean_dir
cd $seqclean_dir
seqclean $p_asteroides/c_abundFilter/trinity_out_dir/Trinity.fasta
## no of transcripts
echo $(grep "^>" Trinity.fasta.clean | wc -l) ## 881362

## remove small transcripts
perl ${script_path}/removesmalls.pl 201 Trinity.fasta.clean > Trinity.fasta.clean.201
trinity_transcriptome=$seqclean_dir/Trinity.fasta.clean.201
## get the longest isoform from my asselbly (edit the Ctr+v,tab at coomand line)
cat Trinity.fasta.clean.201 | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | sed 's/_i/_i\t/;s/ len=/\t/;s/ path=/\tpath=/' | sort -t $'\t' -k1,1 -k3,3nr | sort -t $'\t' -k1,1 -u -s | cut -f 1,2,5 | sed 's/Ctr+v,tab/./' | tr "\t" "\n"  | fold -w 60 > Trinity.fasta.clean.201.longest
echo $(($(grep "^>" Trinity.fasta.clean | wc -l) - $(grep "^>" Trinity.fasta.clean.201 | wc -l))) ## 1077
###########
## Functional annotation of whole transcriptome assembly
## a) generation of gene map
module load trinity/6.0.2
cd $seqclean_dir
perl ${script_path}/trinity_util/get_Trinity_gene_to_trans_map.pl $trinity_transcriptome > gene_trans_map
gene_transcript_map=$seqclean_dir/gene_trans_map

## b) Generation of possible ORFs
module load TransDecoder/2.0.1
TransDecoder.LongOrfs -t $trinity_transcriptome
LongOrfs=$seqclean_dir/$(basename $trinity_transcriptome).transdecoder_dir/longest_orfs.pep

## c) Run blast to search Trinity transcripts & Transdecoder-predicted proteins
module load BLAST+/2.2.30
## header of blast output= qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
mkdir ${p_asteroides}/uniprot
cd ${p_asteroides}/uniprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot

#wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
#gunzip uniref90.fasta.gz
#makeblastdb -in uniref90.fasta -dbtype prot

mkdir $seqclean_dir/blastx_dir
cd $seqclean_dir/blastx_dir
cp $trinity_transcriptome trans_truncated.fa
perl ${script_path}/splitFasta.pl trans_truncated.fa 300  ## 500 for uniprot_uniref90

# blasting all chunks to the targrt DB database
for f in subset*_trans_truncated.fa; do
  qsub -v input=$f,DB=${p_asteroides}/uniprot/uniprot_sprot.fasta,label="uniprot_sprot" ${script_path}/blastx_targetDB.sh;
done
cat subset*_trans_truncated.fa.bx > ../uniprot_sprot.blastx.outfmt6             ## 5614133

#for f in subset*_trans_truncated.fa; do
#  qsub -v input=$f,DB=${p_asteroides}/uniprot/uniref90.fasta,label="uniprot_uniref90" ${script_path}/blastx_targetDB.sh;
#done
#cat subset*_trans_truncated.fa.bx > ../uniprot_uniref90.blastx.outfmt6
#rm subset*.fa subset*.blastx

cd ../
cat uniprot_sprot.blastx.outfmt6 | awk '$11 <= 1e-3' > uniprot_sprot.blastx.outfmt6.sig     ## 2498562
sort -k1,1 -k12,12nr -k11,11n  uniprot_sprot.blastx.outfmt6.sig | sort -u -k1,1 --merge > uniprot_sprot.blastx.outfmt6.sig.best             ## 206780
header='qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\tqcovhsp\tstitle'
sed -i -e 1i"${header}" uniprot_sprot.blastx.outfmt6.sig.best
###########
## Abundance estimation
cd $seqclean_dir
qsub -v index="salmon_index",transcriptome="$trinity_transcriptome" ${script_path}/salmonIndex.sh

cd $p_asteroides/c_abundFilter
for f in $p_asteroides/c_abundFilter/*.s_pe.fq.1; do if [ -f $f ]; then
 identifier=$(basename ${f%_L00*.s_pe.fq.1}); echo $identifier;
fi;done | uniq > identifiers.txt
identifiers=$p_asteroides/c_abundFilter/identifiers.txt

while read identifier;do
 ls ${identifier}_L00*.s_pe.fq.1 ${identifier}_L00*.s_pe2.fq ${identifier}_L00*.s_se.fq
 qsub -v index="$seqclean_dir/salmon_index",identifier=$identifier ${script_path}/salmonQuant_PE.sh
 qsub -v index="$seqclean_dir/salmon_index",identifier=$identifier ${script_path}/salmonQuant_SE.sh
#  wl=$(cat ${identifier}_L00*.s_pe.fq.1 ${identifier}_L00*.s_se.fq | wc -l);
#  c=$(($wl/4)); echo $identifier $c >> readCounts
done < $identifiers
find ./*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
python $script_path/gather-counts.py -i "$(pwd)"
echo "transcript"$'\t'"length" > transcripts.lengthes
sf=$(find ./*.quant -name \*.sf | head -n1)
cat $sf | grep -v "^#" | awk -F "\t" -v OFS='\t' '{print $1,$2}' >> transcripts.lengthes

module load R/3.0.1
while read identifier;do
  echo $(pwd) $identifier
  Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "$gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

## exclude the unexpressed transcripts from the annoatation files
cd $seqclean_dir
comm -12 <(cat $seqclean_dir/uniprot_sprot.blastx.outfmt6.sig.best | tail -n+2 | awk '{print $1}' |sort) <(cat $p_asteroides/c_abundFilter/unexp_isoformTPM | awk '{print $1}' | sort) > unexpIDs  ## 2657
grep -v -w -F -f unexpIDs uniprot_sprot.blastx.outfmt6.sig.best >  uniprot_sprot.blastx.outfmt6.sig.best.exp                      ## 204124

head -n1 uniprot_sprot.blastx.outfmt6.sig.best > uniprot_sprot.blastx.outfmt6.sig.best.unexp
grep -w -F -f unexpIDs uniprot_sprot.blastx.outfmt6.sig.best >>  uniprot_sprot.blastx.outfmt6.sig.best.unexp                    ## 2657
## exclude the unexpressed transcripts from the transcriptome
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $trinity_transcriptome --output_fasta_fp Trinity.clean.201.exp.fasta --seq_id_fp $p_asteroides/c_abundFilter/unexp_isoformTPM --negate
exp_transcriptome=$seqclean_dir/Trinity.clean.201.exp.fasta        ## 868905
echo $(($(grep "^>" Trinity.fasta.clean.201 | wc -l) - $(grep "^>" Trinity.clean.201.exp.fasta | wc -l))) ## 11380

tail -n+2 $p_asteroides/c_abundFilter/unexp_isoformTPM | awk '{print ">"$1"|"}' | grep -F -f - $LongOrfs | sed 's/>//' > LongOrfs.key
filter_fasta.py --input_fasta_fp $LongOrfs --output_fasta_fp $LongOrfs.exp --seq_id_fp LongOrfs.key --negate
###################
## transfer the annotation to the fasta files
## the structure of UniprotKB header (stitle in Blast): http://www.uniprot.org/help/fasta-headers
## db ['sp' for UniProtKB/Swiss-Prot and 'tr' for TrEMBL]|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]PE=ProteinExistence SV=SequenceVersion
cat $seqclean_dir/uniprot_sprot.blastx.outfmt6.sig.best.exp | awk -F '\t' '{print $1,$17}' > $seqclean_dir/uniprot_sprot.blastx.outfmt6.sig.best.exp.key
blastx_key=$seqclean_dir/uniprot_sprot.blastx.outfmt6.sig.best.exp.key
while read line; do
  if [[ $line == \>* ]]; then
    newline=$(echo $line | awk -F '[> ]' '{print $2}' | grep -w -f - $blastx_key)
    if [ "$newline" != "" ];then
      echo ">"$newline;
    else echo $line; fi
  else echo $line; fi
done < $exp_transcriptome > ${exp_transcriptome%.fasta}.ann.fasta
ann_exp_transcriptome=${exp_transcriptome%.fasta}.ann.fasta

## change the Fasta header to remove special characters interfering with blast server
sed 's/|/_/g' $ann_exp_transcriptome > ${ann_exp_transcriptome%.ann.fasta}.ann2.fasta
ann_exp_transcriptome2=${ann_exp_transcriptome%.ann.fasta}.ann2.fasta

## change the Fasta header to remove len and path phrases added by trinity to allow ASN conversion
sed 's/len=.*$//g' $ann_exp_transcriptome2 > ${ann_exp_transcriptome2%.ann2.fasta}.ann3.fasta
ann_exp_transcriptome3=${ann_exp_transcriptome2%.ann2.fasta}.ann3.fasta

## remove univec contaminants
## The UniVec Database: http://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/
## Contamination in Sequence Databases: http://www.ncbi.nlm.nih.gov/tools/vecscreen/contam/
## About VecScreen (parameters and categories): http://www.ncbi.nlm.nih.gov/tools/vecscreen/about/
## Interpretation of VecScreen Results: http://www.ncbi.nlm.nih.gov/tools/vecscreen/interpretation/
mkdir -p $p_asteroides/UniVec/db
cd $p_asteroides/UniVec/db
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
mv UniVec UniVec.fasta
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
#mv UniVec_Core UniVec_Core.fasta
module load BLAST+/2.2.30
makeblastdb -in UniVec.fasta -dbtype nucl
#makeblastdb -in UniVec_Core.fasta -dbtype nucl
cd $p_asteroides/UniVec
input=$ann_exp_transcriptome3
DB=$p_asteroides/UniVec/db/UniVec.fasta
blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query $input -db $DB -num_threads 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore score" -out univec_blastn.out ## 442213
cat univec_blastn.out | awk '{print $1}' | sort | uniq | wc -l ## 183041

## choose all hit with score more than the weak hit for terminal match
cat univec_blastn.out | awk 'BEGIN{OFS="\t";} {if($14>=16) print $1,$7,$8,$9,$14;}' | sort -k1,1 -k2,2n -k3,3nr -k5,5nr > univec_blastn.out.sig ## 442213

## remove contained matches and merge overlapping or neighbering matches(assign highest score)
b_qseqid=""
while read qseqid qstart qend qlen score;do
if [ "$qseqid" == "$b_qseqid" ];then
  max=$(( b_score > score ? b_score : score ))
  if [ "$qend" -le "$b_qend" ]; then b_score=$max; ## keep first with highest score
  elif [ "$qstart" -le "$((b_qend+24))" ];then
#    echo -e "$b_qseqid\t$b_qstart\t$b_qend\t$b_qlen\t$max"; ## echo first with highest score
    b_qend=$qend;b_score=$max;## keep second with  highest score
  else
    echo -e "$b_qseqid\t$b_qstart\t$b_qend\t$b_qlen\t$b_score"; ## echo first
    b_qstart=$qstart;b_qend=$qend;b_score=$score;## keep second
  fi
elif [ "$b_qseqid" != "" ];then
  echo -e "$b_qseqid\t$b_qstart\t$b_qend\t$b_qlen\t$b_score"; ## echo first
  b_qseqid=$qseqid;b_qstart=$qstart;b_qend=$qend;b_qlen=$qlen;b_score=$score; ## keep second
else b_qseqid=$qseqid;b_qstart=$qstart;b_qend=$qend;b_qlen=$qlen;b_score=$score; ## keep first
fi; done < univec_blastn.out.sig > univec_blastn.out.sig2 ## 233312

## select moderate and strong hits
cat univec_blastn.out.sig2 | awk '(((($2<=25)||(($4-$3)<=25))&&$5>=19)||(($2>25)&&(($4-$3)>25)&&($5>=25)))' | sort -k1,1 -k2,2n -k3,3nr -k5,5nr > univec_blastn.out.sig3 ## 1182
#confirmation#cat univec_blastn.out.sig3 | sort -u -k1,1 -k2,2 -k3,3 --merge | wc -l ## 1182
awk '{print $1,$2,$3,$4;}' univec_blastn.out.sig3 > cont.report.len

#cat report_FASTA.txt | awk -F '[ .]' '{print $4,$5,$7}' > cont.report
#sort cont.report | awk '{print $1}' | uniq -c | awk '$1==1 {print $2}' | grep -F -w -f - cont.report > cont.report.uniq
#sort cont.report | awk '{print $1}' | uniq -c | awk '$1>1 {print $2}' | grep -F -w -f - cont.report > cont.report.dup  ## I made the decision manually by selecting the line with more deep coordinates
#sed -n '2p' cont.report.dup >> cont.report.uniq
#grep "^>" $exp_transcriptome | awk -F '[ >=]' '{print $2,$4}' | sed 's/|/_/g' > trans_len
#Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=F,row.names=NULL); head(data1); data2=read.table(args[2],header=F,row.names=NULL); head(data2); dataMerge=merge(data1,data2,by="V1",all.x=T);write.table(dataMerge,args[3], sep=" ", quote=F, row.names=F, col.names=F);' cont.report.uniq trans_len cont.report.len

cat cont.report.len | awk '($4-$3) < 201 && ($2-1) < 201 {print $1}' > excludeIDs  ## 75
## remove excludeIDs (this would remove their duplicates as well if exist)
cat excludeIDs | awk '{print $1}' | grep -v -F -w -f - cont.report.len > cont.report.len.keep ## 1106 (i.e. we removed 75 transcript, one of them has dup matches)
#sort cont.report.len.keep | awk '{print $1}' | uniq -c | awk '($1>1) {print $2}' | grep -F -w -f - cont.report.len | sort > cont.report.len.dup ## 2
cat cont.report.len.keep | awk '($4-$3) <= ($2-1) {print $1,$2}' | sort -k1,1 -k2,2n | sort -u -k1,1 --merge > suffix_cut ## 571 ## cut from $2-len
cat cont.report.len.keep | awk '($4-$3) > ($2-1) {print $1,$3}' | sort -k1,1 -k2,2nr | sort -u -k1,1 --merge > prefix_cut ## 535 ## cut from 1-$3

module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $ann_exp_transcriptome3 --output_fasta_fp suffix_cut.fa --seq_id_fp suffix_cut
while read line; do
  if [[ $line == \>* ]]; then
    suffix=$(echo $line | awk -F '[> ]' '{print $2}' | grep -w -f - suffix_cut | awk '{print $2}')
    echo $line;
  else echo ${line:0:$suffix-1}; fi
done < suffix_cut.fa > suffix_cut.fa.fixed

while read contig;do mark=$(echo $contig | cut -d" " -f 1); grep $mark suffix_cut.fa.fixed | sed 's/>//';done < prefix_cut > suffix_cut_pass
filter_fasta.py --input_fasta_fp $ann_exp_transcriptome3 --output_fasta_fp temp_prefix_cut.fa --seq_id_fp suffix_cut --negate
filter_fasta.py --input_fasta_fp temp_prefix_cut.fa --output_fasta_fp prefix_cut.1.fa --seq_id_fp prefix_cut
filter_fasta.py --input_fasta_fp suffix_cut.fa.fixed --output_fasta_fp prefix_cut.2.fa --seq_id_fp suffix_cut_pass
cat prefix_cut.1.fa prefix_cut.2.fa > prefix_cut.fa
while read line; do
  if [[ $line == \>* ]]; then
    prefix=$(echo $line | awk -F '[> ]' '{print $2}' | grep -w -f - prefix_cut | awk '{print $2}')
    echo $line;
  else echo ${line:$prefix}; fi
done < prefix_cut.fa > prefix_cut.fa.fixed
filter_fasta.py --input_fasta_fp $ann_exp_transcriptome3 --output_fasta_fp ${ann_exp_transcriptome3%.ann3.fasta}.UniVec.fasta --seq_id_fp cont.report.len --negate
univec_ann_exp_tran=${ann_exp_transcriptome3%.ann3.fasta}.UniVec.fasta

filter_fasta.py --input_fasta_fp suffix_cut.fa.fixed --output_fasta_fp suffix_cut.fa.fixed_final --seq_id_fp suffix_cut_pass --negate
cat suffix_cut.fa.fixed_final >> $univec_ann_exp_tran
cat prefix_cut.fa.fixed >> $univec_ann_exp_tran
## trimmed 1106 & excluded 75 from 868905 ==> 868830

## exclude the univec transcripts from the annoatation files
cd $seqclean_dir
while read line;do echo $line | sed 's/|/_/'; done < uniprot_sprot.blastx.outfmt6.sig.best.exp > uniprot_sprot.blastx.outfmt6.sig.best.exp2
grep -v -w -F -f $p_asteroides/UniVec/excludeIDs uniprot_sprot.blastx.outfmt6.sig.best.exp2 >  uniprot_sprot.blastx.outfmt6.sig.best.exp2.univec    ## 204110

sed 's/TR\([0-9]*\)|c\([0-9]*\)_g/TR\1_c\2_g/g' $LongOrfs.exp > $LongOrfs.exp2
cat $p_asteroides/UniVec/excludeIDs | awk '{print ">"$1"|"}' | grep -F -f - $LongOrfs.exp2 | sed 's/>//' > LongOrfs.exp2.key
filter_fasta.py --input_fasta_fp $LongOrfs.exp2 --output_fasta_fp $LongOrfs.exp2.univec --seq_id_fp LongOrfs.exp2.key --negate
##################
## remove the Forgin contamination sequences (FCS). Based on GenBank report
mkdir -p $p_asteroides/FCS
cd $p_asteroides/FCS
# download the NCBI_FCS report
# creat Trim file from the trim section & creat exclude file from the exclude section & dup from the duplicated sequences section
cat FCS_report_dup | sed 's/lcl|//g' | awk '{print $1}' > FCS_report_dup_keep
cat FCS_report_dup | sed 's/lcl|//g' | awk '{i=1; while(i < NF-1) {print $i;i++;}}' | grep -v -F -w -f FCS_report_dup_keep  > FCS_report_dup_exclude

sed -i 's/\.\./\t/g' FCS_report_trim
cat FCS_report_trim | awk -F '\t' '{if(NF > 5)print $1}' > FCS_report_trim_exclude
cat FCS_report_trim | awk -F '\t' '(NF == 5) && ($2-$4) < 201 && ($3-1) < 201 {print $1}' >> FCS_report_trim_exclude 
grep -v -F -w -f FCS_report_trim_exclude FCS_report_trim > FCS_report_trim.keep 
cat FCS_report_trim.keep | awk '($2-$4) <= ($3-1) {print $1,$3}' | sort -k1,1 -k2,2n | sort -u -k1,1 --merge > suffix_cut ## 571 ## cut from $3-len
cat FCS_report_trim.keep | awk '($2-$4) > ($3-1) {print $1,$4}' | sort -k1,1 -k2,2nr | sort -u -k1,1 --merge > prefix_cut ## 535 ## cut from 1-$4

module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp suffix_cut.fa --seq_id_fp suffix_cut
while read line; do
  if [[ $line == \>* ]]; then
    suffix=$(echo $line | awk -F '[> ]' '{print $2}' | grep -w -f - suffix_cut | awk '{print $2}')
    echo $line;
  else echo ${line:0:$suffix-1}; fi
done < suffix_cut.fa > suffix_cut.fa.fixed

while read contig;do mark=$(echo $contig | cut -d" " -f 1); grep $mark suffix_cut.fa.fixed | sed 's/>//';done < prefix_cut > suffix_cut_pass
filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp temp_prefix_cut.fa --seq_id_fp suffix_cut --negate
filter_fasta.py --input_fasta_fp temp_prefix_cut.fa --output_fasta_fp prefix_cut.1.fa --seq_id_fp prefix_cut
filter_fasta.py --input_fasta_fp suffix_cut.fa.fixed --output_fasta_fp prefix_cut.2.fa --seq_id_fp suffix_cut_pass
cat prefix_cut.1.fa prefix_cut.2.fa > prefix_cut.fa
while read line; do
  if [[ $line == \>* ]]; then
    prefix=$(echo $line | awk -F '[> ]' '{print $2}' | grep -w -f - prefix_cut | awk '{print $2}')
    echo $line;
  else echo ${line:$prefix}; fi
done < prefix_cut.fa > prefix_cut.fa.fixed

cat FCS_report_exclude FCS_report_trim FCS_report_dup_exclude > FCS_report_all
filter_fasta.py --input_fasta_fp $univec_ann_exp_tran --output_fasta_fp ${univec_ann_exp_tran%.UniVec.fasta}.FCS.fasta --seq_id_fp FCS_report_all --negate
FCS_ann_exp_tran=${univec_ann_exp_tran%.UniVec.fasta}.FCS.fasta

filter_fasta.py --input_fasta_fp suffix_cut.fa.fixed --output_fasta_fp suffix_cut.fa.fixed_final --seq_id_fp suffix_cut_pass --negate
## remove the annotation from trimmed transcripts
cat suffix_cut.fa.fixed_final | awk '{print $1}' > suffix_cut.fa.fixed.noAnn
cat prefix_cut.fa.fixed | awk '{print $1}' > prefix_cut.fa.fixed.noAnn
cat suffix_cut.fa.fixed.noAnn >> $FCS_ann_exp_tran
cat prefix_cut.fa.fixed.noAnn >> $FCS_ann_exp_tran
## trimmed 26 & excluded 1568 from 868830 ==> 867262

## exclude the FCS transcripts from the annoatation files
cd $seqclean_dir
grep -v -w -F -f $p_asteroides/FCS/FCS_report_all uniprot_sprot.blastx.outfmt6.sig.best.exp2.univec >  uniprot_sprot.blastx.outfmt6.sig.best.exp2.FCS    ## 204071

cat $p_asteroides/FCS/FCS_report_all | awk '{print ">"$1"|"}' | grep -F -f - $LongOrfs.exp2.univec | sed 's/>//' > LongOrfs.exp2.key2
filter_fasta.py --input_fasta_fp $LongOrfs.exp2.univec --output_fasta_fp $LongOrfs.exp2.FCS --seq_id_fp LongOrfs.exp2.key2 --negate   ## 578372
##################
## Assessement of the transcriptome
cd $seqclean_dir
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl $trinity_transcriptome >  $trinity_transcriptome.MatzStat
perl ${script_path}/seq_stats.pl $exp_transcriptome > $exp_transcriptome.MatzStat
perl ${script_path}/seq_stats.pl $univec_ann_exp_tran > $univec_ann_exp_tran.MatzStat
perl ${script_path}/seq_stats.pl $FCS_ann_exp_tran > $FCS_ann_exp_tran.MatzStat

module load trinity/6.0.2
TrinityStats.pl $trinity_transcriptome >  $trinity_transcriptome.TrinityStat
TrinityStats.pl $exp_transcriptome > $exp_transcriptome.TrinityStat
TrinityStats.pl $univec_ann_exp_tran > $univec_ann_exp_tran.TrinityStat
TrinityStats.pl $FCS_ann_exp_tran > $FCS_ann_exp_tran.TrinityStat

## calc the the no of Complete ORFs
#grep "type:complete" $LongOrfs | wc -l  ##225219
#grep -A1 "type:complete" $LongOrfs | grep -v "^--" > $LongOrfs.complete ## 225219
#grep "^>" $LongOrfs.complete | awk -F '[>|]' '{print $2"|"$3}' | sort | uniq | wc -l ## 124271
grep "type:complete" $LongOrfs.exp2.univec | wc -l  ## 223844
grep -A1 "type:complete" $LongOrfs.exp2.univec | grep -v "^--" > $LongOrfs.exp2.univec.complete
grep "^>" $LongOrfs.exp2.univec.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 123425

grep "type:complete" $LongOrfs.exp2.FCS | wc -l  ## 223741
grep -A1 "type:complete" $LongOrfs.exp2.FCS | grep -v "^--" > $LongOrfs.exp2.FCS.complete 
grep "^>" $LongOrfs.exp2.FCS.complete | awk -F '[>|]' '{print $2}' | sort | uniq | wc -l ## 123339

## instaling and running assemblathon2
#cd ${script_path}
#git clone https://github.com/ucdavis-bioinformatics/assemblathon2-analysis.git
#cd assemblathon2-analysis ## you have to be in this folder so that the file can use the FAlite.pm script
#perl ./assemblathon_stats.pl ${seqclean_dir}/Trinity.fasta.clean > ${seqclean_dir}/trinity_assemblathon2.stat
#perl ./assemblathon_stats.pl ${p_asteroides}/data/Porites.Astreoides.uniprot2013.fa > ${p_asteroides}/data/trinity_assemblathon2.stat

## http://deweylab.biostat.wisc.edu/detonate/vignette.html
#cd $ass_dir
#module load DETONATE/1.8.1
#lf="${lf_files[*]}"
#rt="${rt_files[*]}"
#rsem-eval-calculate-score --paired-end $(echo ${lf[*]} | tr ' ' ',') $(echo ${rt[*]} | tr ' ' ',') \
#			  trinity_out_dir/seqclean/Trinity.fasta.clean \
#			  Trinity.fasta.clean.detonate \
#			  400 \
#			  --transcript-length-parameters rsem-eval/true_transcript_length_distribution/mouse.txt \
#			  -p 16

#####################
## TSA Submission Guide
## http://www.ncbi.nlm.nih.gov/genbank/tsaguide
## create ASN file
## https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2#tbl
module load tbl2asn/20150331
mkdir ${p_asteroides}/ASN
cd ${p_asteroides}/ASN
## download the template.sbt & assembly.cmt
cat $FCS_ann_exp_tran | awk '{print $1}' > allTrans.fsa
#cd ~/bin
#wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
#gunzip linux64.tbl2asn.gz
#mv linux64.tbl2asn tbl2asn
#chmod 755 tbl2asn
tbl2asn -t template.sbt -p. -Y assembly.cmt -M t -j "[organism=Porites astreoides] [moltype=transcribed_RNA] [tech=TSA]"
#qsub $script_path/createASN.sh
##############
## copy the assemblies to the other server
cd ${p_asteroides}/ASN
scp allTrans.fsa tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/NCBI_submission/allTrans.fsa
scp allTrans.sqn tmansour@loretta.hpcf.upr.edu:/storage/prcen/coral/NCBI_submission/allTrans.sqn
##############
## compare the new Assembly with the older assembly
#bash $script_path/compareVSoldTrans.sh
##################
