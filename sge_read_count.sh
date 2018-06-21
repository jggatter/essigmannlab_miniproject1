#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 7
#$ -M jggatter@mit.edu
#############################################

if [ $# -ne 1 ]
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
	
	printf "DEBUG: group $group\n"
	printf "$group,,\n" >> $csv_dir/sge_read_count.csv

  if [ "$(ls -A 01)" ]; then
      
  	cd 01    #Read counts for the post-trimmomatic .fastq files
  	posttrim_samples=(*$group*P.fastq)
  	for sam in ${posttrim_samples[@]}; do
  		
      printf "  DEBUG: sam $sam\n"
      printf "    DEBUG: Counting reads...\n"
      line_count=$(wc -l $sam | cut -d' ' -f1)  #cut cleaves the filename from the output
      #printf "      DEBUG: Line count $line_count\n"
      read_count=$(( $line_count / 4 ))
      printf "      DEBUG: Read count $read_count\n"
      printf ",$sam,$read_count\n" >> $csv_dir/sge_read_count.csv

  	done
  	cd ..

  else
    printf "  DEBUG: No fastq files found in dir 01, continuing to 02...\n"
  fi


  if [ "$(ls -A 02)" ]; then
  
  	cd 02   #Read counts for the post-pear .fastq files
    pear_samples=(*$group*assembled*.fastq)
  	for sam2 in ${pear_samples[@]}; do
  	  
      printf "  DEBUG: sam2 $sam2\n"
      printf "    DEBUG: Counting reads...\n"
      line_count=$(wc -l $sam2 | cut -d' ' -f1)   #cut cleaves the filename from the output
      #printf "      DEBUG: Line count $line_count\n"
      read_count=$(( $line_count / 4 ))
      printf "      DEBUG: Read count $read_count\n"
      printf ",$sam2,$read_count\n" >> $csv_dir/sge_read_count.csv

  	done
    cd ..
  else
    printf "  DEBUG: No fastq files found in dir 02, continuing to 03/$group...\n"
  fi
	
	cd 03/$group                                 #Steps into the group folder and stores all barcode directories
	barcode_dirs=(*[^.fastq])
	for (( i = 0; i < ${#barcode_dirs[@]}; i++ )); do # A for-each loop is not used as information is stored in different arrays by index.
    
    if ! [ "$(ls -A ${barcode_dirs[$i]})" ]; then
      printf "  DEBUG: No fastq files found in 03/$group/${barcode_dirs[$i]}, continuing to next barcode dir...\n"
      continue
    fi

    cd ${barcode_dirs[$i]}                        #Steps into the barcode folder and stores all .fastq files
    printf "  DEBUG: barcode ${barcode_dirs[i]}\n"
    
    matchbarcode_samples=(*.fastq)
    for sam3 in ${matchbarcode_samples[@]}; do    #Read counts for the post-match_barcode.sh script processing

      printf "    DEBUG: sam3 $sam3\n"
      printf "    DEBUG: Counting reads...\n"
      line_count=$(wc -l $sam3 | cut -d' ' -f1)   #cut cleaves the filename from the output
      #printf "      DEBUG: Line count $line_count\n"
      read_count=$(( $line_count / 4 ))
      printf "      DEBUG: Read count $read_count\n"
      printf ",$sam3,$read_count\n" >> $csv_dir/sge_read_count.csv
      let "barcode_read_counts[$i] = ${barcode_read_counts[$i]} + $read_count"   #Update barcode read count
      printf "      DEBUG: ${barcode_dirs[i]} count: ${barcode_read_counts[i]}\n"
      printf "      DEBUG: DONE! Logged in sge_read_count.csv.\n"

    done
    cd ..

	done
  cd ../..
	
done
printf "DEBUG: done with part 1\n"

printf "\nBarcode,Total Read Count,Percentage\n" >> $csv_dir/sge_read_count.csv     #Begins the second part of the .csv, which is barcode-related information
barcodes_total_count=0
for (( i = 0; i < ${#barcode_dirs[@]} ; i++ )); do
  let "barcodes_total_count = $barcodes_total_count + ${barcode_read_counts[$i]}"   #Sum of all barcode counts    
done
printf "Total count: $barcodes_total_count\n"

declare -a barcode_percentages=(0 0 0 0 0 0 0 0)
for (( i = 0; i < ${#barcode_dirs[@]} ; i++ )); do
  #barcode_percentages[$i]=$(bc -l <<< "scale = 2; 100 * (${barcode_read_counts[$i]} / ${barcodes_total_count})")    #Calculating barcode's read count percentage. Rous doesn't like this method.
  barcode_percentages[$i]=$(( (100 * ${barcode_read_counts[$i]}) / ${barcodes_total_count} ))    #Calculating barcode's read count percentage as an integer. Rous likes this method.
  printf "DEBUG: barcode ${barcode_dirs[$i]}, read_counts ${barcode_read_counts[$i]}, ${barcode_percentages[$i]} percent\n"
  printf "${barcode_dirs[$i]},${barcode_read_counts[$i]},${barcode_percentages[$i]}\n" >> $csv_dir/sge_read_count.csv   #Storing all barcode-related information in .csv
done