#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 1
#$ -M jggatter@mit.edu
#############################################

#fasta_converter
#Converts the sequence from .seq files into an appropriately formatted .fasta reference file.
#Please change the ROUS header to include your own username!
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#July 3rd, 2018

if [ $# -ne 1 ]; then
  printf "Invalid number of arguments. Please enter only the .txt files directory without a '/'.\n"
  exit 1
fi

main_dir=$(pwd) #New FASTA files will be stored in main_dir/fastas
if [ -d $main_dir/fastas ]; then    #If the .fasta already exists, it is removed
      rm -rf $main_dir/fastas
      printf "Replacing fastas dir\n"
fi
mkdir fastas
txt_files_dir=$1

#Array of lesions that will substitute for the 2nd discontiguous N region
declare -A refs=([1]='ANA' [2]='ANC' [3]='ANG' [4]='ANT'
                 [5]='CNA' [6]='CNC' [7]='CNG' [8]='CNT'
                 [9]='GNA' [10]='GNC' [11]='GNG' [12]='GNT'
                 [13]='TNA' [14]='TNC' [15]='TNG' [16]='TNT'
                 [17]='ANGA' [18]='ANGC' [19]='ANGG' [20]='ANGT'
                 [21]='CNGA' [22]='CNGC' [23]='CNGG' [24]='CNGT'
                 [25]='GNGA' [26]='GNGC' [27]='GNGG' [28]='GNGT'
                 [29]='TNGA' [30]='TNGC' [31]='TNGG' [32]='TNGT')

cd $txt_files_dir
txt_files=*_ref.txt

for txt in ${txt_files[@]}; do 

  printf "Reading $txt...\n"

  while IFS= read -r line || [[ -n "$line" ]]; do   #The sequence is the first line of the .txt file.
		sequence=$line
	done < "$txt"

  search_sequence=$sequence    #Preserves sequence. Copies sequence into search_sequence, which will be repeatedly truncated
  i=0
  NNN_second_index=0           #The index of the 2nd discontiguous instance of N 
  while (( i < 2 )); do        #Repeats only once
    
    NNN_index=$(expr index "$search_sequence" NNN)            #Search for first NNN
    NNN_second_index=$(( $NNN_second_index + $NNN_index ))    #Target index is the sum of the first index of NNN and the substring between the 1st NNN and 2nd discontiguous instance of N  
    search_sequence=${search_sequence:$(( $NNN_index + 2 ))}  #Truncate sequence to search for next discontigious instance of N
    i=$(( i+ 1 ))

  done

  #Extract first half and second half excluding the 2nd discontiguous N
	sequence_first_half=${sequence:0:$(( $NNN_second_index + 1 ))}
	sequence_second_half=${sequence:$(( $NNN_second_index + 4 ))}

  txt_truncated=${txt%_*}   #Truncates "_ref.txt" from .txt file name for .fasta naming purposes.
  for ref in ${refs[@]}; do

    printf ">${txt_truncated}_$ref\n" >> $main_dir/fastas/${txt_truncated}_$ref.fasta #Fasta header
    fasta_sequence=$sequence_first_half$ref$sequence_second_half
    line_num=$(( (${#fasta_sequence} / 74) + 1 )) #Calculates number of 74-character-long lines that the sequence will be printed into
    for (( i = 0; i < line_num; i++ )); do
      printf "${fasta_sequence:0:74}\n" >> $main_dir/fastas/${txt_truncated}_$ref.fasta    #Extract and print first 74 characters of string
      fasta_sequence=${fasta_sequence:74}     #Extract substring from index 74 onwards for next line
    done
    printf "  ${txt_truncated}_$ref.fasta has been created in $main_dir/fastas\n" 

  done

done