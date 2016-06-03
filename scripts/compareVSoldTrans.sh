## compare the new Assembly with the older assembly
## convert the old assembly to trinity format 
cd $p_asteroides/data
cat Porites.Astreoides.uniprot2013.fa | awk '/^>/ {if(N>0) printf("\n"); printf("%s ",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | awk '{print ">",$(NF-1),$(NF-4),$NF}' | sed 's/> tr=/>TR1|c1_g/;s/ /_i/;s/_of_[0-9]*_in_tr[0-9]*//' | tr " " "\n"  | fold -w 60 > Porites.Astreoides.uniprot2013.trinformat.fasta
## run SeqClean on the old Assembly
seqclean Porites.Astreoides.uniprot2013.trinformat.fasta
## I will ignore fitering the old assembly to remove short ests
## get the longest isoform from the old assembly
cat Porites.Astreoides.uniprot2013.trinformat.fasta.clean | awk '/^>/ {if(N>0) printf("\n"); printf("%s ",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | sed 's/_i/_i /' | awk '{printf("%s %d\n",$0,length($3));}' | sort -k1,1 -k4,4nr | sort -k1,1 -u -s | sed 's/ /./' | cut -d" " -f 1,2 | tr " " "\n" | fold -w 60 > Porites.Astreoides.uniprot2013.trinformat.fasta.clean.longest
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl $data_path/Porites.Astreoides.uniprot2013.trinformat.fasta.clean > $data_path/Porites.Astreoides.uniprot2013.trinformat.fasta.clean.MatzStat
module load trinity/6.0.2
TrinityStats.pl $data_path/Porites.Astreoides.uniprot2013.trinformat.fasta.clean > $data_path/Porites.Astreoides.uniprot2013.trinformat.fasta.clean.TrinityStat

## compare the whole assemblies
cd $seqclean_dir
mkdir  TrinityAll_BlastDB && cd TrinityAll_BlastDB
cp ../Trinity.fasta.clean.201 .
module load BLAST+/2.2.30
makeblastdb -in Trinity.fasta.clean.201 -input_type fasta -dbtype nucl
blastn -query $data_path/Porites.Astreoides.uniprot2013.trinformat.fasta.clean \
-db Trinity.fasta.clean.201 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
-max_target_seqs 10  -out oldAssVsnewAss ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k12,12nr -k11,11n  oldAssVsnewAss | sort -u -k1,1 --merge > oldAssVsnewAss_best
wc -l oldAssVsnewAss_best ## 313803 (out of 316231)
cat oldAssVsnewAss_best | awk '$11 <= 1e-5' | wc -l        ## 313803
cat oldAssVsnewAss_best | awk '$13 < $14' > oldAssVsnewAss_best_better
wc -l oldAssVsnewAss_best_better ## 226134
cat oldAssVsnewAss_best | awk '$13 == $14' > oldAssVsnewAss_best_equal
wc -l oldAssVsnewAss_best_equal ## 5774
cat oldAssVsnewAss_best | awk '$13 > $14' > oldAssVsnewAss_best_less
wc -l oldAssVsnewAss_best_less ## 81895
cat oldAssVsnewAss_best_better | awk '{ sum+=$14} END {print sum}' ## 506894995
cat oldAssVsnewAss_best_better | awk '{ sum+=$13} END {print sum}' ## 290066489 (i.e. difference of 216828506 =~217Mb)
cat oldAssVsnewAss_best_less | awk '{ sum+=$14} END {print sum}' ## 132389312
cat oldAssVsnewAss_best_less | awk '{ sum+=$13} END {print sum}' ## 181193999 (i.e. difference of 48804687 =~49Mb)

## compare the longest isoforms
cd $seqclean_dir
mkdir  Trinity_BlastDB && cd Trinity_BlastDB
cp ../Trinity.fasta.clean.201.longest .
module load BLAST+/2.2.30
makeblastdb -in Trinity.fasta.clean.201.longest -input_type fasta -dbtype nucl
blastn -query $data_path/Porites.Astreoides.uniprot2013.trinformat.fasta.clean.longest \
-db Trinity.fasta.clean.201.longest \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 -evalue 1e-5 \
-max_target_seqs 10  -out oldAssVsnewAss ## -perc_identity 90 -qcov_hsp_perc 50
sort -k1,1 -k12,12nr -k11,11n  oldAssVsnewAss | sort -u -k1,1 --merge > oldAssVsnewAss_best
wc -l oldAssVsnewAss_best ## 115860 (out of 120592)  ## 117304
cat oldAssVsnewAss_best | awk '$11 <= 1e-5' | wc -l        ## 115860  ## 117304
cat oldAssVsnewAss_best | awk '$13 < $14' > oldAssVsnewAss_best_better
wc -l oldAssVsnewAss_best_better ## 90581  ## 91655
cat oldAssVsnewAss_best | awk '$13 == $14' > oldAssVsnewAss_best_equal
wc -l oldAssVsnewAss_best_equal ## 4682 ## 4677
cat oldAssVsnewAss_best | awk '$13 > $14' > oldAssVsnewAss_best_less
wc -l oldAssVsnewAss_best_less ## 20597  ## 20972
cat oldAssVsnewAss_best_better | awk '{ sum+=$14} END {print sum}' ## 140843586 ## 138694936
cat oldAssVsnewAss_best_better | awk '{ sum+=$13} END {print sum}' ## 80890520 (i.e. difference of 59953066 =~60Mb) ## 79737533 (difference 58957403)
cat oldAssVsnewAss_best_less | awk '{ sum+=$14} END {print sum}' ## 26434723 ## 27276896
cat oldAssVsnewAss_best_less | awk '{ sum+=$13} END {print sum}' ## 32603280 (i.e. difference of 6168557 =~6Mb) ## 34212511 (difference 6935615)

