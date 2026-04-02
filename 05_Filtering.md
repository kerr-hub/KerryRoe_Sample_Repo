
# Filtering VCF

To check for poor samples and determine appropriate filtering parameters I generated the following statisitcs and plotted in R to visualize.
As outlined by the following tutorial https://speciationgenomics.github.io/
The tutorial subsets the data first, but I wanted to be able to see averages for all variant sites. 

1) Minor allele frequency for each variant site
2) Average depth per individual
3) Average depth per variant site
4) Average quality per site
5) Missingness per variant site
6) Missingness per individual
7) Heterzygosity/inbreeding coefficeint for each individual 


Generate stats with vcftools
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=qc_full
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 80G
#SBATCH -c 16
#SBATCH --time 12:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

module load vcftools

VCF=batch1_joint_all.vcf.gz
OUT=./qc

module load vcftools
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
vcftools --gzvcf $VCF --depth --out $OUT
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
vcftools --gzvcf $VCF --site-quality --out $OUT
vcftools --gzvcf $VCF --missing-indv --out $OUT
vcftools --gzvcf $VCF --missing-site --out $OUT
vcftools --gzvcf $VCF --het --out $OUT
```
## Graphing in R

Quality
```
library(tidyverse)
#quality
var_qual <- read_delim("./qc.lqual", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
ggsave("var_qual_density.png", plot = a, width = 6, height = 4, dpi = 300)
```

Variant mean depth
```
var_depth <- read_delim("./qc.ldepth.mean", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
a + theme_light() + xlim(0, 100)
# Save zoomed plot (0–100)
a_zoom <- a + xlim(0, 100)
ggsave("var_depth_density_0to100.png", plot = a_zoom, width = 6, height = 4, dpi = 300)
summary(var_depth$mean_depth)
quantile(var_depth$mean_depth, probs = c(0.05, 0.95), na.rm = TRUE)
```
Variant missingness
```
var_miss <- read_delim("./qc.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
# without X11 capabilities
a <- ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

ggsave("missingness_density.png", plot = a, width = 6, height = 4, dpi = 300)

```
Minor allele frequency
```
var_freq <- read_delim("./qc.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_freq$maf)

ggsave("maf.png", plot = a, width = 6, height = 4, dpi = 300)
```
Mean depth per individual
```
depth <- read_table("qc.idepth")
ggplot(depth, aes(x = reorder(INDV, MEAN_DEPTH), y = MEAN_DEPTH)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Mean sequencing depth",
    title = "Mean depth per sample (before filtering)"
  ) +
  theme_minimal(base_size = 12)
ggsave("inddepth.png", plot = a, width = 6, height = 4, dpi = 300)
print(ind_depth)
```
Missing data per individual
```
imiss <- read_table("qc.imiss")
ggplot(imiss, aes(x = reorder(INDV, F_MISS), y = F_MISS)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Sample", y = "Missing data fraction", 
       title = "Missing data per individual") +
  theme_minimal()

ggsave("indmissing.png", plot = a, width = 15, height = 10, dpi = 300)
summary(var_miss$fmiss)

print(ind_miss)
```
Heterozygosity and inbreeding per individual 
```
het <- read_table("qc_before_filter.het")
a <- ggplot(het, aes(x = reorder(INDV, F), y = F)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    x = "Sample",
    y = "Inbreeding coefficient (F)",
    title = "Per-sample heterozygosity (PLINK --het)"
  )
ggsave("ind_het_histogram.png", plot = a, width = 15, height = 10, dpi = 300)
```
Individuals with an average sequencing depth less than 3 were removed, or missingness over 20% were removed form the dataset. 

### Sex chromosomes
Next I identified Contigs that are likely sex chromosomes were removed (identified using minimap2, aligning to sex chromosomes from the Atlantic Puffin, assuming 10% divergence between the species. 
```
#!/bin/bash
#SBATCH --job-name=sex_chroms_10
#SBATCH --mem 66G
#SBATCH -c 1
#SBATCH --time 24:01:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca


module load  minimap2/2.28

#NOTE: minimap2 command must be run in directory that either contains the reference genomes themselves or else symlinks to the reference genomes
#NOTE: individual query contigs are split up by minimap2 for alignment; current filtering strategy throws out contig if one segment aligns to a sex chromosome
#NOTE: Using a conservative mapq score of 13 (corresponding to a 5% chance of incorrect alignment) in order to say particular contig segments do indeed align to a sex chromosome


#all of the entries in a given position in these lists correspond to a single scaffolded reference genome (e.g. every position 0 in these lists is a set that corresponds to the limnodromus scolopaceus)
chromosomal_refGenomeList=("GCA_947846985.1_FratArc_v3_genomic.fna")
chromosomal_speciesNames=("fratercula_arctica" "fratercula_arctica")
z_chromosome_names=("CANTUZ010000026.1")
w_chromosome_names=("CANTUZ010000025.1")

scaffold_refGenomeList=("ref_genome.fna")
scaffold_speciesNames=("cepphus_grylle")

minimap2 -ax asm10 "${chromosomal_refGenomeList[$SLURM_ARRAY_TASK_ID]}"  "${scaffold_refGenomeList[$SLURM_ARRAY_TASK_ID]}" > "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_alignedto_"${chromosomal_speciesNames[$SLURM_ARRAY_TASK_ID]}".sam

awk -v Z="${z_chromosome_names[$SLURM_ARRAY_TASK_ID]}" '($3 == Z && $5 >= 20) {print $1}' "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_alignedto_"${chromosomal_speciesNames[$SLURM_ARRAY_TASK_ID]}".sam | sort -u > "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_z_contigs.txt

awk -v W="${w_chromosome_names[$SLURM_ARRAY_TASK_ID]}" '($3 == W && $5 >= 20) {print $1}' "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_alignedto_"${chromosomal_speciesNames[$SLURM_ARRAY_TASK_ID]}".sam | sort -u > "${scaffold_speciesNames[$SLURM_ARRAY_TASK_ID]}"_w_contigs.txt
```
```
cat *_w_contigs.txt *_z_contigs.txt > all_sex_linked_contigs.txt

bcftools query -f '%CHROM\t%END\n' removed_nofilter.vcf.gz > contig_lengths.txt

# Keep only contigs to remove
grep -Ff sex_contigs_clean_final.txt contig_lengths.txt > sex_contigs.bed

bcftools view -T ^sex_contigs.bed -Oz -o nosex_nofilter.vcf.gz removed_nofilter.vcf.gz

```
### NOTE: I ran A PCA to check for batch effects and outliers before finishing final filtering. I found pairs of siblings which I removed. Code in 06_PCA. 

Before filtering by quality and missingness, I masked genotypes with a depth of lower than 4, converting them to NA. 
```
bcftools +setGT final_samples.vcf.gz \
  -Oz -o DP4_masked.vcf.gz \
  -- -t q -n . -i 'FMT/DP<4'
```
Next I identified loci that only had up to 20% missing across all populations. This way bigger populations with 6-8 samples could have 1 missing, but smaller ones could not.
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=pop_missing_filter_vcf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 86G
#SBATCH --time 12:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e


module load vcftools

for pop in pop*.txt
do
  name=${pop}

  vcftools --gzvcf DP4_masked.vcf.gz \
  --keep ${pop} \
  --max-missing 0.8 \
  --recode --out ${name}
done
```
Now find common loci across all 11 populations that have no more than 20% missingness per population
```
# amke pos file for each vcf

zgrep -v "^#" file.vcf.gz | cut -f1,2 > file.pos # for each pop vcf


# remove these sites from the vcf
bcftools view \
  -T shared_sorted.pos \
  -Oz -o REDO_20miss_dp4.vcf.gz \
  final_samples_nofilter.vcf.gz

```
Now filter with vcftools, removing low quality snps, snps with maf <1%
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=filter_vcf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 20G
#SBATCH --time 2:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

# set filters
MAF=0.01
QUAL=30
MIN_DEPTH=6
MAX_DEPTH=10
VCF_OUT=../analysis/FINAL.vcf.gz
VCF_IN=./20miss_dp4.vcf.gz

module load vcftools

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--recode --stdout | gzip -c > \
$VCF_OUT
```
Make a dataset with no maf filtering for PCA
```
#!/bin/bash
#SBATCH --account=def-vlf
#SBATCH --job-name=filter_vcf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=19kr21@queensu.ca
#SBATCH --mem 20G
#SBATCH --time 2:00:00
#SBATCH -o %x-%j.o
#SBATCH -e %x-%j.e

# set filters
QUAL=30
MIN_DEPTH=6
MAX_DEPTH=10
VCF_OUT=../population_analysis/nomaf_20miss_dp4.vcf.gz
VCF_IN=./20miss_dp4.vcf.gz

module load vcftools

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--recode --stdout | gzip -c > \
$VCF_OUT
```
