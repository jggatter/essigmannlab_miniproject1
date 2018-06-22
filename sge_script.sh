#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -pe whole_nodes 7
#$ -M klkim@mit.edu
###############################

:'
NGS MiSeq Analysis Pipeline
(Submission script for Sun Grid Engine job scheduler)
Version 1.0
Author: Lina Kim, klkim [at] mit.edu
July 5, 2017
'


module add jre/1.7.0-25
module add pear/0.9.5
module add bowtie2/2.2.6
module add bwa/0.7.12
module add samtools/0.1.19

# Step 0: Setup with necessary variables and directories
clear

###----- PARAMETERS TO EDIT -----
ID=''
SAM_END=''
IN_PATH=data/$ID/
PATH_TO_TRIM=''
declare -a SAMPLES=('');
declare -a BARCODES=('');
###----- END EDIT -----

declare -A refs=([1]='AXA' [2]='AXC' [3]='AXG' [4]='AXT'
                 [5]='CXA' [6]='CXC' [7]='CXG' [8]='CXT'
                 [9]='GXA' [10]='GXC' [11]='GXG' [12]='GXT'
                 [13]='TXA' [14]='TXC' [15]='TXG' [16]='TXT')

mkdir logs
mkdir $ID

# Step 1: Adapter trimming with Trimmomatic
mkdir $ID/01

echo 'Trimming adapters...'

for SAM in "${SAMPLES[@]}"
do
  echo $SAM
  logFile=logs/$ID\_$SAM.log.txt
  touch $logFile
  PATH_TO=$IN_PATH$SAM$SAM_END/$ID\_$SAM
  java -jar $PATH_TO_TRIM/trimmomatic-0.36.jar PE $PATH_TO\_1\_sequence.fastq $PATH_TO\_2\_sequence.fastq -baseout $ID/01/$ID\_$SAM.fastq ILLUMINACLIP:$PATH_TO_TRIM/adapters/NexteraPE-PE.fa:2:30:10:1:true 1>> $logFile
done

# Step 2: Concatenating pair-end reads with PEAR
echo 'Concatenating paired-end reads...'

mkdir $ID/02
for SAM in "${SAMPLES[@]}"
do
  echo $SAM
  logFile=logs/$ID\_$SAM.log.txt
  FILE_NAME=$ID/01/$ID\_$SAM
  pear -f $FILE_NAME\_1P.fastq -r $FILE_NAME\_2P.fastq -o $ID/02/$ID\_$SAM -n 10 -b 33 1>> $logFile
done

# Step 3: Selecting reads matching 13-mers around barcode, replacing low-quality (Q<=10) bases with N
echo 'Selecting matching reads and performing quality check...'

mkdir $ID/03
bc_list=''
for i in ${BARCODES[@]}; do
  bc_list+=${i},
done
bc_list=${bc_list:0:${#bc_list}-1}

for SAM in "${SAMPLES[@]}"
do
  echo $SAM
  FILE_NAME=$ID/02/$ID\_$SAM.assembled.fastq,$ID/02/$ID\_$SAM.unassembled.forward.fastq,$ID/02/$ID\_$SAM.unassembled.reverse.fastq
  mkdir $ID/03/$SAM
  for BC in ${BARCODES[@]}; do
    mkdir $ID/03/$SAM/$BC
    logFile=logs/$ID\_$SAM.log.txt
  done
  python scripts/match_barcode.py -i $FILE_NAME -r $ID -s $SAM -b $bc_list 1>> $logFile
done

# Step 4: Align to reference genome
echo 'Aligning to reference M13 genome...'

mkdir $ID/04

# Build 16 reference sequences
# Build bowtie2 index
cd data/ref
bash write_context_refs.sh
# and also index for samtools
samtools faidx M13mp7_ref.fasta
bowtie2-build -f M13mp7_ref.fasta M13mp7_ref
for (( i=1; i<17; i++ )); do
  samtools faidx ref$i.fasta
  bowtie2-build -f ref$i.fasta ref$i
done

cd -

# Make alignments, convert files to BAM
for SAM in "${SAMPLES[@]}"
do
  echo $SAM
  mkdir $ID/04/$SAM
  for BC in "${BARCODES[@]}"; do
    mkdir $ID/04/$SAM/$BC
    for (( i=1; i<17; i++ )); do
      ref=${refs[$i]}
      logFile=logs/$ID\_$SAM\_$BC\_$ref.log.txt
      touch $logFile
      FILE_NAME=$ID/03/$SAM/$BC/$SAM\_$BC\_$ref.matched.fastq
      bowtie2 -x data/ref/ref$i $FILE_NAME -S $ID/04/$SAM/$BC/$SAM\_$BC\_$ref.sam 1>> $logFile
      samtools view -bS $ID/04/$SAM/$BC/$SAM\_$BC\_$ref.sam | samtools sort - $ID/04/$SAM/$BC/$SAM\_$BC\_$ref.sort 1>> $logFile
      samtools index $ID/04/$SAM/$BC/$SAM\_$BC\_$ref.sort.bam
    done
  done
done

echo 'Alignments complete.'

# Step 5: Pileup mapping results for mutation data
echo 'Creating pileup mapping results...'

mkdir $ID/05

# Create position list file
touch pos_list.txt
for i in {6244..6280}; do
  echo -e 'M13mp7_ref\t'$i >> pos_list.txt
done

inp_list=''

for SAM in "${SAMPLES[@]}"
do
  echo $SAM
  mkdir $ID/05/$SAM
  for BC in "${BARCODES[@]}"; do
    mkdir $ID/05/$SAM/$BC
    for (( i=1; i<17; i++ )); do
      ref=${refs[$i]}
      inp_list+=$ID/05/$SAM/$BC/$SAM\_$BC\_$ref.pileup,
      samtools mpileup -AB -d1000000 -f data/ref/ref$i.fasta -l pos_list.txt $ID/04/$SAM/$BC/$SAM\_$BC\_$ref.sort.bam > $ID/05/$SAM/$BC/$SAM\_$BC\_$ref.pileup
    done
  done
done
echo 'Pileup mapping results complete.'

inp_list=${inp_list:0:${#inp_list}-1}

# Step 6: Analyze mutation data from pileup files
echo 'Analyzing mutation data from pileup files...'
python scripts/from_pileup.py -i $inp_list -o analysis_01.csv