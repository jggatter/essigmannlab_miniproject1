#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 7
#$ -M jggatter@mit.edu
#############################################

#targeted_bowtie2_aligner
#This script requires a lot of disk space, processing power, and time.
#Uses bowtie2 to align .fasta reference files to .fastq paired-end reads (forward and reverse) and store output in .txt files.
#Please change the ROUS header to include your own username!
#Version 1.0
#Author: James Gatter, jggatter [at] mit.edu
#July 3rd, 2018

export PATH=/net/bmc-pub7/data0/essigmannlab/dependencies/anaconda3/bin:$PATH
export PATH=/bin:$PATH

if [ $# -ne 3 ]; then 
  printf "Please format as ./targeted_bowtie2_aligner.sh [Path to parent dir] [NXN lesion] [NXGN lesion]\n"
  printf "Formatted example: ./targeted_bowtie2_aligner.sh parent_dir ANA ANGA\n" 
  exit 1
fi

main_dir=$(pwd)
parent_dir=$1
NXN_lesion=$2
NXGN_lesion=$3

cd $main_dir/fastas
ref_fastas=(*.fasta)
for ref in ${ref_fastas[@]}; do bowtie2-build -f $ref ${ref%".fasta"}; done #Index files are built for each reference .fasta file.

declare -A NXGN_sample_groups=([1]='2662' [2]='2666' [3]='2670' [4]='2674' [5]='2678' [6]='2682')
declare -A NXN_sample_groups=([1]='2660' [2]='2661' [3]='2663' [4]='2664' [5]='2665' [6]='2667' 
                              [7]='2659' [8]='2668' [9]='2669' [10]='2671' [11]='2672' [12]='2673' 
                              [13]='2675' [14]='2676' [15]='2677' [16]='2679' [17]='2680' [18]='2681')

#The following code is run for the NXN lesion.

if [ ! -d $main_dir/${NXN_lesion}_txts ]; then    #If the dir doesn't exist, it is made
      printf "Making ${NXN_lesion}_txts dir\n"
      mkdir $main_dir/${NXN_lesion}_txts
fi
for NXN_group in ${NXN_sample_groups[@]}; do

    #The 1k and M13mp7 .fasta files are located and their extensions are deleted for bowtie2 semantics
    cd $main_dir/fastas
    NXN_1k_fasta=(1k_$NXN_lesion.fasta)
    NXN_M13_fasta=(M13mp7_$NXN_lesion.fasta)
    NXN_1k_fasta=${NXN_1k_fasta%".fasta"}
    NXN_M13_fasta=${NXN_M13_fasta%".fasta"}
    
    #The .fastq files are located in the [parent]/01 directory.
    cd $main_dir/$parent_dir/01
		NXN_onep=(*${NXN_group}_1P.fastq)
		NXN_twop=(*${NXN_group}_2P.fastq)
    
    #The alignment command is run for both the 1k and M13mp7 .fasta references.
    #The fasta is aligned to the paired end read .fastq files, an output .sam file is created, and the alignment statistics are stored in a .txt file
    printf "NXN Aligning $NXN_onep and $NXN_twop to $NXN_ref\n"
    bowtie2 -x $main_dir/fastas/$NXN_1k_fasta -1 $main_dir/$parent_dir/01/${NXN_onep} -2 $main_dir/$parent_dir/01/${NXN_twop} -S $main_dir/fastas/${NXN_group}_${NXN_1k_fasta}.sam &>> $main_dir/${NXN_lesion}_txts/${NXN_group}_${NXN_1k_fasta}.txt
    printf "NXN Aligning $NXN_onep and $NXN_twop to $NXN_ref\n"
    bowtie2 -x $main_dir/fastas/$NXN_M13_fasta -1 $main_dir/$parent_dir/01/${NXN_onep} -2 $main_dir/$parent_dir/01/${NXN_twop} -S $main_dir/fastas/${NXN_group}_${NXN_M13_fasta}.sam &>> $main_dir/${NXN_lesion}_txts/${NXN_group}_${NXN_M13_fasta}.txt

done

#The same code is now run for the NXGN lesion.

if [ ! -d $main_dir/${NXGN_lesion}_txts ]; then
      printf "Making ${NXGN_lesion}_txts dir\n"
      mkdir $main_dir/${NXGN_lesion}_txts
fi
for NXGN_group in ${NXGN_sample_groups[@]}; do

		cd $main_dir/fastas
    NXGN_1k_fasta=(1k_${NXGN_lesion}.fasta)
    NXGN_M13_fasta=(M13mp7_${NXGN_lesion}.fasta)
    NXGN_1k_fasta=${NXGN_1k_fasta%".fasta"}
    NXGN_M13_fasta=${NXGN_M13_fasta%".fasta"}

    cd $main_dir/$parent_dir/01
    NXGN_onep=(*${NXGN_group}_1P.fastq)
    NXGN_twop=(*${NXGN_group}_2P.fastq)

    printf "NXGN Aligning $NXGN_onep and $NXGN_twop to $NXGN_ref\n"
    bowtie2 -x $main_dir/fastas/$NXGN_1k_fasta -1 $main_dir/$parent_dir/01/${NXGN_onep} -2 $main_dir/$parent_dir/01/${NXGN_twop} -S $main_dir/fastas/${NXGN_group}_${NXGN_1k_fasta}.sam &>> $main_dir/${NXGN_lesion}_txts/${NXGN_group}_${NXGN_1k_fasta}.txt
    printf "NXGN Aligning $NXGN_onep and $NXGN_twop to $NXGN_ref\n"
    bowtie2 -x $main_dir/fastas/$NXGN_M13_fasta -1 $main_dir/$parent_dir/01/${NXGN_onep} -2 $main_dir/$parent_dir/01/${NXGN_twop} -S $main_dir/fastas/${NXGN_group}_${NXGN_M13_fasta}.sam &>> $main_dir/${NXGN_lesion}_txts/${NXGN_group}_${NXGN_M13_fasta}.txt

done