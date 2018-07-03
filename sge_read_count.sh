#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 7
#$ -M jggatter@mit.edu
#############################################

#sge_read_count.sh
#Reads .fastq files and counts the number of reads for each sample and percentage for each barcode in a .csv file.
#Please change the ROUS header to include your own username!
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#July 3rd, 2018

if [ $# -ne 1 ]; then
  printf "Invalid number of arguments. Please enter only the parent directory without a '/'.\n"
  exit 1
fi

parent_dir=$1                         #The parent directory must be specified as the first argument
csv_dir=$(pwd)                        #Where the .csv will be created and updated
if [ -f sge_read_count.csv ]; then    #If the .csv already exists, it is removed
  rm $csv_dir/sge_read_count.csv
  printf "Replaced sge_read_count.csv\n"
fi
printf "Sample Group,Filename,Read Count\n" >> $csv_dir/sge_read_count.csv

cd $parent_dir/03 
sample_groups=(*)   #Step into dir 03 to get all group names 
cd ..

declare -a barcode_read_counts=(0 0 0 0 0 0 0 0)    #Keep tally of all barcode-related counts (sam3 counts)
for group in ${sample_groups[@]}; do
	
	printf "Sample group: $group\n"
	printf "$group,,\n" >> $csv_dir/sge_read_count.csv

  if [ "$(ls -A 01)" ]; then
      
  	cd 01    #Read counts for the post-trimmomatic .fastq files
  	posttrim_samples=(*$group*P.fastq)
  	for sam in ${posttrim_samples[@]}; do

      line_count=$(wc -l $sam | cut -d' ' -f1)  #cut cleaves the filename from the output
      read_count=$(( $line_count / 4 ))
      printf ",$sam,$read_count\n" >> $csv_dir/sge_read_count.csv

  	done
  	cd ..

  else
    printf "  WARNING: No fastq files found in dir 01, continuing to 02...\n"
  fi


  if [ "$(ls -A 02)" ]; then
  
  	cd 02   #Read counts for the post-pear .fastq files
    pear_samples=(*$group*assembled*.fastq)
  	for sam2 in ${pear_samples[@]}; do
  	  
      line_count=$(wc -l $sam2 | cut -d' ' -f1)   #cut cleaves the filename from the output
      read_count=$(( $line_count / 4 ))
      printf ",$sam2,$read_count\n" >> $csv_dir/sge_read_count.csv

  	done
    cd ..
  else
    printf "  WARNING: No fastq files found in dir 02, continuing to 03/$group...\n"
  fi
	
	cd 03/$group                                 #Steps into the group folder and stores all barcode directories
	barcode_dirs=(*[^.fastq])
	for (( i = 0; i < ${#barcode_dirs[@]}; i++ )); do # A for-each loop is not used as information is stored in different arrays by index.
    
    if ! [ "$(ls -A ${barcode_dirs[$i]})" ]; then
      printf "  WARNING: No fastq files found in 03/$group/${barcode_dirs[$i]}, continuing to next barcode dir...\n"
      continue
    fi

    cd ${barcode_dirs[$i]}                        #Steps into the barcode folder and stores all .fastq file
    matchbarcode_samples=(*.fastq)
    for sam3 in ${matchbarcode_samples[@]}; do    #Read counts for the post-match_barcode.sh script processing

      line_count=$(wc -l $sam3 | cut -d' ' -f1)   #cut cleaves the filename from the output
      read_count=$(( $line_count / 4 ))
      printf ",$sam3,$read_count\n" >> $csv_dir/sge_read_count.csv
      (( barcode_read_counts[$i] = ${barcode_read_counts[$i]} + $read_count ))   #Update barcode read count
      printf "    DONE! Logged in sge_read_count.csv.\n"

    done
    cd ..

	done
  cd ../..
	
done

printf "Assembling percentages...\n"
printf "\nBarcode,Total Read Count,Percentage\n" >> $csv_dir/sge_read_count.csv     #Begins the second part of the .csv, which is barcode-related information
barcodes_total_count=0
for (( i = 0; i < ${#barcode_dirs[@]} ; i++ )); do
  let "barcodes_total_count = $barcodes_total_count + ${barcode_read_counts[$i]}"   #Sum of all barcode counts    
done

declare -a barcode_percentages=(0 0 0 0 0 0 0 0)
for (( i = 0; i < ${#barcode_dirs[@]} ; i++ )); do
  #barcode_percentages[$i]=$(bc -l <<< "scale = 2; 100 * (${barcode_read_counts[$i]} / ${barcodes_total_count})")    #Calculating barcode's read count percentage. Rous doesn't like this method.
  barcode_percentages[$i]=$(( (100 * ${barcode_read_counts[$i]}) / ${barcodes_total_count} ))    #Calculating barcode's read count percentage as an integer. Rous likes this method.
  printf "${barcode_dirs[$i]},${barcode_read_counts[$i]},${barcode_percentages[$i]}\n" >> $csv_dir/sge_read_count.csv   #Storing all barcode-related information in .csv
done

printf "Script finished.\n"