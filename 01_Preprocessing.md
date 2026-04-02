# Pre-processing

Whole genomce sequences were generated two batches with up to 58 samples using paired-end 150 base pair reads on an Illumina NovaSeq at Genome Quebec, Montreal, QC. 

## Fastqc

To examine the qaulity of the raw reads I ran fastqc, and exported each factqc report in multiqc to view them all together

```
#!/bin/bash
#SBATCH --job-name=fastqc_reports
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --account=def-vlf
#SBATCH --mem 48G
#SBATCH -c 1
#SBATCH --time 10:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load fastqc

mkdir fastqc_reports

for file in *_R1.fastq.gz; do
    fastqc "$file" "${file/_R1/_R2}" -o fastqc_reports;
done

```
MultiQC 
```
# Download multiqc -  I had to download rust and cargo dependancies before downloading multiqc

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# restart setion to load cargo and rust to path
# download multiqc

pip install multiqc

# enter directory with fastqc reports in it and run multiqc
cd fastqc_reports

# run multiqc
multiqc .
```

Download html file on desktop to view

## Fastp Trimming

Make multiple sample lists to run multiple samples at once:
1) Remove any aliasing that will alter the output of ls
2) Get file names
3) Remove every second entry
4) Remove file extension (.fastq.gz) and last three characters upstream (to get just file name, not _R1/_R2)
5) Split file into chunks of 5 lines (30/5 = 10 files)
```
ls *.gz -1 > file_list1.txt
cp file_list1.txt file_list2.txt
sed -i -e n\;d file_list2.txt
sed -i 's/\(.*\)............/\1/' file_list2.txt
split file_list2.txt -a 2 -l 5
# -a 2 → Specifies that the suffix for output files should be 2 characters long (e.g., xaa, xab, xac, etc.).
# -l 5 → Splits the file into chunks of 5 lines each.

```
Re name split files
```
mv xaa file_split1.txt
mv xab file_split2.txt
...
```
Move split files into 01_adaptertrimming
```
mv file_split* ../01_adaptertrimming/
```

Remove adaptors and trim low qaulity sequences. Illumina Noveseq data can have low-quality 3’ tails with false instances of G (Lou & Therkildsen, 2022). I used a sliding window approach to trim the 3’ end of reads, where reads were trimmed after 4 consecutive bases with a quality score of less than 20. 

```
#!/bin/bash

#SBATCH --job-name=fastp
#SBATCH --account=def-vlf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --array=1-7  #1-12 for batch 1
#SBATCH --mem 40G
#SBATCH -c 8
#SBATCH --time 6:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load fastp

mkdir fastp_reps

# one way to pass sample list = first passed parameter on command line
# SAMPLELIST=$1

# or, just set the file name directly
SAMPLELIST="file_split${SLURM_ARRAY_TASK_ID}.txt"

# loop through sample list
while IFS=$'\t' read -r -a array; do
fastp --detect_adapter_for_pe --dedup -g -r --cut_right_window_size 4 --cut_right_mean_quality 20 -w 8 -i ../00_raw/${array[0]}_R1.fastq.gz -I ../00_raw/${array[0]}_R2.fastq.gz \
-o ./${array[0]}_trimmed_R1.fastq.gz -O ./${array[0]}_trimmed_R2.fastq.gz \
--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-h ./fastp_reps/${array[0]}.html
done < $SAMPLELIST

# --cut_right_window_size 4 --cut_right_mean_quality 20 will cut right end after 4 bases with less than 20 quality
# -w is worker thread number
# -g is cut polyg tail
# -r, --cut_right - move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
# -W, --cut_window_size  the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
# Length filtering is enabled by default, but you can disable it by -L
# -p overrepresentation analysis

```
