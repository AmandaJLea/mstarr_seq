#!/bin/bash

# do in interactive session first
# export PATH="$HOME/miniconda/bin:$PATH"
# source $HOME/miniconda/bin/activate
# conda install -c bioconda cutadapt

# working directory
cd /scratch/tmp/alea/mSTARR/

# for f in `cat 14Nov2019_Novaseq_rerun.txt`; do cat process2.sh | sed -e s/FILEINFO/$f/g > process2.${f}.sh; done
# rm commands2.sh; touch commands2.sh
# for f in `cat 14Nov2019_Novaseq_rerun.txt` ; do echo 'sh process2.'$f'.sh' >> commands2.sh; done
# sbatch -a 1-24%24 array2.sh

# modules
module load samtools
module load bedtools

# set paths
f2=FILEINFO
R1=/Genomics/ayroleslab2/alea/archive_raw_fastq/from_RAJ_mSTARR/14Nov2019_Novaseq_rerun/${f2}_R1_001.fastq.gz
R2=/Genomics/ayroleslab2/alea/archive_raw_fastq/from_RAJ_mSTARR/14Nov2019_Novaseq_rerun/${f2}_R2_001.fastq.gz
out_DIR=/scratch/tmp/alea/mSTARR/bams

f=run2_FILEINFO
cutadapt --nextseq-trim 20 -e 0.05 --overlap 2 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length=20 --trim-n -o $out_DIR/${f}.trimmed_R1.fastq.gz -p $out_DIR/${f}.trimmed_R2.fastq.gz $R1 $R2

echo 'trimming done'

bwa mem -t 8 /Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_all_chr.fa $out_DIR/${f}.trimmed_R1.fastq.gz $out_DIR/${f}.trimmed_R2.fastq.gz > $out_DIR/${f}.sam

echo 'mapping done'

rm $out_DIR/${f}.trimmed_R1.fastq.gz
rm $out_DIR/${f}.trimmed_R2.fastq.gz
# rm $out_DIR/${f}.sam

# bam to bed
samtools view -bSq 1 $out_DIR/${f}.sam > $out_DIR/${f}.bam
samtools sort -n -o $out_DIR/${f}.temp.bam $out_DIR/${f}.bam 
samtools fixmate $out_DIR/${f}.temp.bam $out_DIR/${f}.temp2.bam
samtools sort -n -o $out_DIR/${f}.temp3.bam $out_DIR/${f}.temp2.bam 
bedtools bamtobed -i $out_DIR/${f}.temp3.bam -bedpe | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6}'> ${f}.pe.bed

rm $out_DIR/${f}.temp3.bam
rm $out_DIR/${f}.temp2.bam
rm $out_DIR/${f}.temp.bam

# get pileups
# bedtools makewindows -g /Genomics/ayroleslab2/alea/ref_genomes/hg38/chrNameLength.txt -w 200 > /Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_200bp_windows.bed
intervals=/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_200bp_windows.bed
awk '$6-$2<2000' ${f}.pe.bed | awk '$6-$2>0' | awk '$1==$4' | awk '{OFS="\t";print $1,$2,$6}' > ${f}.filt_pe.bed
bedtools coverage -a $intervals -b ${f}.filt_pe.bed | awk '$4>0'> ${f}_200kb_cov.bed

rm ${f}.filt_pe.bed

echo 'pileups done'
