#########
# combine files from R1 and R2
#########

#!/bin/bash

cat run1_L31223_S1_L001.pe.bed run2_L31223_S1_L001.pe.bed> L31223.combined.bed
cat run1_L31224_S2_L001.pe.bed run2_L31224_S2_L001.pe.bed> L31224.combined.bed
cat run1_L31225_S3_L001.pe.bed run2_L31225_S3_L001.pe.bed> L31225.combined.bed
cat run1_L31232_S10_L001.pe.bed run2_L31232_S10_L001.pe.bed> L31232.combined.bed
cat run1_L31233_S11_L001.pe.bed run2_L31233_S11_L001.pe.bed> L31233.combined.bed
cat run1_L31234_S12_L001.pe.bed run2_L31234_S12_L001.pe.bed> L31234.combined.bed
cat run1_L31241_S19_L001.pe.bed run2_L31241_S19_L001.pe.bed> L31241.combined.bed
cat run1_L31242_S20_L001.pe.bed run2_L31242_S20_L001.pe.bed> L31242.combined.bed
cat run1_L31243_S21_L001.pe.bed run2_L31243_S21_L001.pe.bed> L31243.combined.bed
cat run1_L31250_S28_L001.pe.bed run2_L31250_S28_L001.pe.bed> L31250.combined.bed
cat run1_L31251_S29_L001.pe.bed run2_L31251_S29_L001.pe.bed> L31251.combined.bed
cat run1_L31252_S30_L001.pe.bed run2_L31252_S30_L001.pe.bed> L31252.combined.bed
cat run1_L31259_S45_L002.pe.bed run2_L31259_S37_L002.pe.bed> L31259.combined.bed
cat run1_L31260_S46_L002.pe.bed run2_L31260_S38_L002.pe.bed> L31260.combined.bed
cat run1_L31261_S47_L002.pe.bed run2_L31261_S39_L002.pe.bed> L31261.combined.bed
cat run1_L31268_S54_L002.pe.bed run2_L31268_S46_L002.pe.bed> L31268.combined.bed
cat run1_L31269_S55_L002.pe.bed run2_L31269_S47_L002.pe.bed> L31269.combined.bed
cat run1_L31270_S56_L002.pe.bed run2_L31270_S48_L002.pe.bed> L31270.combined.bed
cat run1_L31277_S63_L002.pe.bed run2_L31277_S55_L002.pe.bed> L31277.combined.bed
cat run1_L31278_S64_L002.pe.bed run2_L31278_S56_L002.pe.bed> L31278.combined.bed
cat run1_L31279_S65_L002.pe.bed run2_L31279_S57_L002.pe.bed> L31279.combined.bed
cat run1_L31286_S72_L002.pe.bed run2_L31286_S64_L002.pe.bed> L31286.combined.bed
cat run1_L31287_S73_L002.pe.bed run2_L31287_S65_L002.pe.bed> L31287.combined.bed
cat run1_L31288_S74_L002.pe.bed run2_L31288_S66_L002.pe.bed> L31288.combined.bed

#########
# get counts per 200bp window from combined files
#########

#!/bin/bash

# working directory
cd /scratch/tmp/alea/mSTARR/bed

# for f in `cat LIDs.txt`; do cat process1.sh | sed -e s/FILEINFO/$f/g > process1.${f}.sh; done
# rm commands1.sh; touch commands1.sh
# for f in `cat LIDs.txt` ; do echo 'sh process1.'$f'.sh' >> commands1.sh; done
# sbatch -a 1-24%24 array1.sh

module load bedtools

intervals=/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_200bp_windows.bed

f=FILEINFO
awk '$6-$2<2000' ${f}.combined.bed | awk '$6-$2>0' | awk '$1==$4' | awk '{OFS="\t";print $1,$2,$6}' > ${f}.filt_pe.bed

bedtools coverage -a $intervals -b ${f}.filt_pe.bed > ${f}_200kb_cov.bed

rm ${f}.filt_pe.bed

