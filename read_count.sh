#!/bin/bash

cd $1
main_dir=$(pwd)
printf "Sample Name,Read Count\n" >> $main_dir/read_count.csv

sample_dirs=(*/)
for dir in ${sample_dirs[@]}; do
	
  printf "Directory: $dir\n"
  cd $dir
  samples=(*.fastq)
  for sam in ${samples[@]}; do

    printf "  Sample: $sam\n"
    line_count=$(wc -l | cut -d ' ' -f1)
    #printf "    DEBUG: $line_count\n"
    read_count=$(( $line_count / 4 ))
    printf "    Read count: $read_count\n"
    printf "$sam,$read_count\n" >> $main_dir/read_count.csv
	  printf "    Sample DONE! Logged in read_count.csv.\n"

  done
  cd ..

done
