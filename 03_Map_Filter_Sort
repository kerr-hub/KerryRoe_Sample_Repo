# Map, filter, and sort bam files

BAMs were mapped and sorted using samtools. Alignments with a MAPQ score of less than 20 were removed. Batch script also produces depth files for each sample. 
I ended up running variant alling on all my samples, and then checking depths in the vcfs, removing bad samples, and rerunning variant calling. 

Set up for batch script
1) Make directory for this step
2) Copy in sample name files for batching
3) Make directories for output
```
mkdir 03_map_filter_sort
cp ./02_alignment/tr_fl_*.txt ./03_map_filter_sort
cd 03_map_filter_sort
mkdir mapped
mkdir sorted
mkdir depth_reps
```
Run batch script 
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --array=1-10
#SBATCH --job-name=map_sort_minq20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 40G
#SBATCH -c 16
#SBATCH --time 12:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load samtools

SAMPLELIST="tr_fl_${SLURM_ARRAY_TASK_ID}.txt"

for SAMPLE in `cat $SAMPLELIST`; do

## Convert to bam file for storage (including all the mapped reads)
  # -F 4 means retain only mapped reads
        samtools view -bS -F 4 ../02_alignment/${SAMPLE}_aligned.bam  > ./mapped/${SAMPLE}_mapped.bam

## Filter the mapped reads (to only retain reads with high mapping quality)
  # Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)

        samtools view -h -q 20 ./mapped/${SAMPLE}_mapped.bam | samtools view -buS - | samtools sort -o ./sorted/${SAMPLE}_sorted_minq20.bam
        samtools index ./sorted/${SAMPLE}_sorted_minq20.bam
```

# samtools depth -aa absolutely all sites even empty ones
# could also do -f and use a txt with list of bams want to use
# will produce a column in each file with the depth (cut -f 3 only keeps this colum)
