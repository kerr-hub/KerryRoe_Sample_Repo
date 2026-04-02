# Variant calling with bcftools

I chose to perform variant calling with bcftools, which has been shown to outperform other genotyping softwares when working with 
low-coverage, non-model species WGS (Lefouili & Nam, 2022). Variant sites or single nucleotide polymorphisms (SNPs) were called by 
passing BAM files through bcftools mpileup and bcftools call, where mpileup produces genotype likelihood from BAMs and directly passes 
(‘pipes’) the output to bcftools call to a generate variant call format (VCF) for each individual (Danecek et al., 2021). 

I had issues with the vcf coming out truncated and not properly bgzipped, but this is the order I got it to work with. 

```
#!/bin/bash

#SBATCH --account=def-vlf
#SBATCH --job-name=mpileup_nobad
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 120G
#SBATCH -c 16
#SBATCH --time 72:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load samtools
module load bcftools

# BAM list file: one BAM per line
BAMLIST="./filelist_nobad.txt"

# joint varaint calling acroos all BAMS
bcftools mpileup -a AD,ADF,DP,ADR,SP \
    -Ou -f ./ref_genome.fna \
    -b $BAMLIST | \
bcftools call --threads 16 -mv -f GQ -Ob -o ../nobad_mpileup/batch1_nobad.bcf

# Then convert to bgzipped VCF + index
bcftools view -Oz ../nobad_mpileup/batch1_nobad.bcf -o ../nobad_mpileup/batch1_nobad.vcf.gz
bcftools index -t ../nobad_mpileup/batch1_nobad.vcf.gz

```
