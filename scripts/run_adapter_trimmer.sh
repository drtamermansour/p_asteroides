#!/bin/sh

if [ $# -lt 3 ]
then
    printf "\nUsage run_adapter_trimmer.sh [input dir] [output dir] [script]\n"
    exit 0
fi

input_dir="$1"
output_dir="$2"
script="$3"
		
cd $input_dir

for f in *_R1_*.gz		# loop for multiple data files in each directory
do
    input_one=$input_dir"/""$f"
    input_two=$input_dir"/"$(echo "$f" | sed s/_R1_/_R2_/)
    qsub -v R1_INPUT="$input_one",R2_INPUT="$input_two",output_dir="$output_dir" "$script" 

done

 