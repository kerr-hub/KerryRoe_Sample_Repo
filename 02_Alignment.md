# Alignment

I aligned sequences to a reference genome using bwa mem to produce bam files for each sample.

Download reference genome from NCBI and unzip/rename to ref_genome
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/401/065/GCA_013401065.1_ASM1340106v1/GCA_013401065.1_ASM1340106v1_genomic.fna.gz

gunzip GCA_013401065.1_ASM1340106v1_genomic.fna.gz

mv GCA_013401065.1_ASM1340106v1_genomic.fna.gz ref_genome.fna
```

Index reference genome 
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=index_refgenome
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 20G
#SBATCH --time 2:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load bwa

bwa index ref_genome.fna
```
Create directory and make sample lists to run an array job
1) Unalias ls
2) List trimmed fastqc files to a txt
3) Remove every other row
4) Remove the uneeded characters from start and end of names 
5) Move file lists to 02_alignment
6) Split file list into groups of 3 samples
7) Change name of file lists
```
mkdir 02_alignment
cd 01_adaptertrimming
unalias ls
ls -1 *.fastq.gz > tr_filelist.txt
sed -i -e n\;d tr_filelist.txt 
sed -i 's/\(.*\)............/\1/' tr_filelist.txt
mv tr_filelist.txt ../02_alignment
cd ../02_alignment/
split tr_filelist.txt -l 3
mv xaa tr_fl_1.txt
mv xab tr_fl_2.txt
mv xac tr_fl_3.txt
...
```
Run alignment using BWA mem and samblaster to remove PCR duplicates. Outputs binary BAMfile instead of SAM file for each sample. 
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --array=1-6
#SBATCH --job-name=alignment_bwa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 20G
#SBATCH -c 8
#SBATCH --time 24:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load StdEnv/2020
module load bwa
module load samblaster
module load samtools

# $SLURM_ARRAY_TASK_ID is the id of the current script copy

SAMPLELIST="tr_fl_${SLURM_ARRAY_TASK_ID}.txt"

for SAMPLE in `cat $SAMPLELIST`; do

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

bwa mem -M -t 8 ./ref_genome.fna ../01_adaptertrimming/${SAMPLE}_R1.fastq.gz ../01_adaptertrimming/${SAMPLE}_R2.fastq.gz | samblaster --removeDups | samtools view -h -b -@8 -o ./${SAMPLE}_aligned.bam

done

# remove PCR duplicates with samblaster
# | samtools view -h -b -@8 -o ./${SAMPLE}_aligned.bam
# samtools view: Converts the deduplicated SAM stream into a BAM file (binary and compressed).
# -h: Include header in output.
# -b: Output BAM format (instead of default SAM).
# -@8: Use 8 threads.
# -o ./${SAMPLE}_aligned.bam: Save the final aligned BAM file with the sample name.
```
