#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 1
#$ -M jggatter@mit.edu
#############################################

#NXN_alignment_count
#Used in tandem with NXGN_alignment_count to extract alignment statistics from .txt files outputted by bowtie2 and compile them in a .csv
#Please change the ROUS header to include your own username!
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#July 3rd, 2018

if [ $# -ne 3 ]; then 
  printf "Please format as ./NXN_alignment_count.sh [Path to txts dir] [Lesion] [Write flag]"
  printf "2nd arg: default for all NXN lesions or input an individual NXN lesion surrounded by single quotation marks\n"
  printf "3rd arg: a for appending to an existing .csv file, w for writing a new .csv file.\n"
  printf "Format example: ./NXN_alignment_count.sh ANA_txts 'ANA' w\n" 
fi

declare -A NXN_sample_groups=([0]='2660' [1]='2661' [2]='2663' [3]='2664' [4]='2665' [5]='2667' 
                              [6]='2659' [7]='2668' [8]='2669' [9]='2671' [10]='2672' [11]='2673' 
                              [12]='2675' [13]='2676' [14]='2677' [15]='2679' [16]='2680' [17]='2681')

lesion_parameter=$2
if [ lesion_parameter = "default" ]; then
  declare -A NXN_lesions=([0]='ANA' [1]='ANC' [2]='ANG' [3]='ANT'
                 [4]='CNA' [5]='CNC' [6]='CNG' [7]='CNT'
                 [8]='GNA' [9]='GNC' [10]='GNG' [11]='GNT'
                 [12]='TNA' [13]='TNC' [14]='TNG' [15]='TNT')
else declare -A NXN_lesions=([0]=$lesion_parameter)
fi

main_dir=$(pwd)                                           #Only run this script in the directory which is the parent of the txts_dir
txts_dir=$1                                               #First argument should be the directory of txt files
write_argument=$3                                         #Only valid option is 'w', which stands for write 
                                        
if [ write_argument = "w" ]; then                         #If the write argument 'w' was entered, set up the .csv file
  if [ -f $main_dir/alignments.csv ]; then                #If the .csv already exists, it is removed
    rm $main_dir/alignments.csv
    printf "Replaced alignments.csv\n"
  fi
    printf "Sample Group,Total Aligned Reads for M13mp7,Alignment Percentage for M13mp7,Total Aligned reads for 1k,Alignment Percentage for 1k\n" >> $main_dir/alignments.csv
fi

cd $main_dir/$txts_dir
for NXN_group in ${NXN_sample_groups[@]}; do
  
  total_read_count_M13=0
  total_read_count_1k=0
  alignment_total_M13=0
  alignment_total_1k=0
  printf "Reading files for ${NXN_group}...\n"

  for NXN_lesion in ${NXN_lesions[@]}; do

    M13_txt=(${NXN_group}_M13mp7_${NXN_lesion}.txt)                    #Read specific M13mp7 .txt file
    line_count=1
    while IFS= read -r line; do

      if [ $line_count -eq 1 ]; then                                    #Make certain .txt files do not contain warning/error messages at the top or the algorithm will not work
        line1_num=$(echo $line | cut -d ' ' -f1)                        #Extracts total read count from line 1 and adds to total read count variable            
        (( total_read_count_M13 = $total_read_count_M13 + $line1_num ))  
      fi
      if [ $line_count -eq 4 ] || [ $line_count -eq 5 ]; then 
        line45_num=$(echo ${line:4} | cut -d ' ' -f1)                   #Extracts aligned counts from lines 4 and 5 and adds to aligned total variable
        (( alignment_total_M13 = $alignment_total_M13 + $line45_num ))  
      fi
      if [ $line_count -gt 5 ]; then break; fi
      (( line_count = $line_count + 1 ))

    done < "$M13_txt"

    txt_1k=(${NXN_group}_1k_${NXN_lesion}.txt)                         #Read specific 1k .txt file
    line_count=1
    while IFS= read -r line; do
        
      if [ $line_count -eq 1 ]; then
        line1_num=$(echo $line | cut -d ' ' -f1)
        (( total_read_count_1k = $total_read_count_1k + $line1_num ))  
      fi
      if [ $line_count -eq 4 ] || [ $line_count -eq 5 ] ; then
        line45_num=$(echo ${line:4} | cut -d ' ' -f1)
        (( alignment_total_1k = $alignment_total_1k + $line45_num ))
      fi
      if [ $line_count -gt 5 ]; then break; fi
      (( line_count = $line_count + 1 ))

    done < "$txt_1k"
    printf "  Done!\n"
  done

  #Following block of code calculates the percentages with two digits following the decimal point
  #Provided that the result is not 0, integer division guarantees at least 3 digits
  (( percentage_reads_M13 = (${alignment_total_M13} * 10000) / ${total_read_count_M13} ))
  (( percentage_reads_1k = (${alignment_total_1k} * 10000) / ${total_read_count_1k} ))
  #Provided that the result is not 0, the two digits after the decimal point are the last two digits in the string/number...
  #...while the digits before the point are the string/number with the digits after the point deleted.
  # EX: 85 aligned out of 100 total reads:
  # 85 * 10000 / 100 = 8500
  # Two digits after = "00", Digits before = "85"
  # Percentage = "85.00"
  if [ ${#percentage_reads_1k} -gt 1 ]; then percentage_reads_1k=${percentage_reads_1k%${percentage_reads_1k:(-2)}}.${percentage_reads_1k:(-2)}; fi
  if [ ${#percentage_reads_M13} -gt 1 ]; then percentage_reads_M13=${percentage_reads_M13%${percentage_reads_M13:(-2)}}.${percentage_reads_M13:(-2)}; fi
  
  printf "D18-${NXN_group},${alignment_total_M13},${percentage_reads_M13},${alignment_total_1k},${percentage_reads_1k}\n" >> $main_dir/alignments.csv

done