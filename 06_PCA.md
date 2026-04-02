### PCA to look for outliers / batch effects

First need to make ld decay pruned dataset from the nomaf dataset

NOTE: code below is from final PCA with final dataset
```
# in scratch/poulation_analysis

module load plink
plink --vcf nomaf_miss90_nodp.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out ./nomaf_ld

# prune and create pca
plink --vcf nomaf_miss90_nodp --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract nomaf_ld.prune.in \
--make-bed --pca --out nomaf_pruned

```
Plot PCA 
```
# read in data
library(tidyverse)

pca <- read_table("nomaf_pruned.eigenvec", col_names = FALSE)
eigenval <- scan("nomaf_pruned.eigenval")

# remove nuisance column
pca <- pca[,-2]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# location - batch1
Location <- rep(NA, length(pca$ind))
Location[grep("SF", pca$ind)] <- "Sud Fulgoya, Norway"
Location[grep("SOD", pca$ind)] <- "Finland"
Location[grep("118", pca$ind)] <- "Akpaqarvik, NU"
Location[grep("NWT18", pca$ind)] <- "Akpaqarvik, NU"
Location[grep("621", pca$ind)] <- "Green Isl., NU"
Location[grep("816", pca$ind)] <- "Saglek, NFL"
Location[grep("501", pca$ind)] <- "Saglek, NFL"
Location[grep("170", pca$ind)] <- "Middle Lawn Isl., NFL"
Location[grep("171401923", pca$ind)] <- "Middle Lawn Isl., NFL"
Location[grep("684", pca$ind)] <- "Kent Isl., NB"
Location[grep("BAC", pca$ind)] <- "Baccalieu Isl., NFL"
Location[grep("KI", pca$ind)] <- "Kelly's Isl, NFL"
Location[grep("LBI", pca$ind)] <- "Little Bell Isl., NFL"
Location[grep("KIW1", pca$ind)] <- "Kelly's Isl., NFL"
Location[grep("N1", pca$ind)] <- "Nain, NFL"
Location[grep("N5", pca$ind)] <- "Nain, NFL"
Location[grep("NWT01", pca$ind)] <- "Coats Isl., NU"
Location[grep("NWT02", pca$ind)] <- "Coats Isl., NU"
Location[grep("NWT03", pca$ind)] <- "Coats Isl., NU"
Location[grep("STL", pca$ind)] <- "St. Lawrence"
Location[grep("EBI", pca$ind)] <- "Mitivik, NU"
Location[grep("QIK", pca$ind)] <- "Qikiqtarjuaq, NU"
Location[grep("429", pca$ind)] <- "Finland"
Location[grep("649", pca$ind)] <- "Finland"
Location[grep("616", pca$ind)] <- "Finland"
Location[grep("SVAL", pca$ind)] <- "Svalbard"
Location[grep("ICE", pca$ind)] <- "Flatey, Iceland"
Location[grep("66", pca$ind)] <- "Country Isl., NS"
Location[grep("6603", pca$ind)] <- "Kent Isl., NB"
Location[grep("MN", pca$ind)] <- "Great Duck Isl., MN"
Location[grep("M7", pca$ind)] <- "Makkovik, NFL"
# combine - if you want to plot each in different colours
spp_loc <- paste0(Location) 

# remake data.frame
pca <- as_tibble(data.frame(pca,Location))


# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
#0074D9", "#B10DC9", "#85144b", "#FF4136", "#FF851B"
"#CAB2D6", "#6A3D9A", "#33A02C",
    "#E31A1C", "#FF7F00", "#A6CEE3", "#B15928", "#000000", "blue

########### plot 1-2

b <- ggplot(pca, aes(PC1, PC2, col = Location, shape = Location)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
    "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A",
    "#FB9A99", "#FDBF6F", "#0074D9", "#B10DC9", "#85144b",
    "#FF4136", "#FF851B"
  )) +
  scale_shape_manual(values = 1:length(unique(pca$Location))) +
  labs(
    x = "PC1 (7.4%)",
    y = "PC2 (6.0%)"
  ) +
  theme_classic() +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 14)
  )

ggsave("pca_individuals_final.png",
       plot = b,
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)

```
There seems to be a few pairs of siblings (individuals that group together outsideof their populations/colony cluster)

Remove sibling that has more missing data if there is a difference

NOTE: Code below is from first PCA run with all samples
```
cd nosibs

plink \
  --bfile ../nomaf_pruned \
  --remove remove_samples.txt \
  --allow-extra-chr \
  --make-bed \
  --out sibsrem_nomaf_pruned

plink --bfile sibsrem_nomaf_pruned --allow-extra-chr --pca --out nosib_pca

```
Ran PCA with final dataset using first bit of code 

### For PCA with population allele freqs

```
# in rorqual /scratch/population_analysis

sed 's/:/\t/' nomaf_ld.prune.out > remove_sites.txt

bcftools view -T ^remove_sites.txt ./REDO_nomaf_30miss.vcf.gz -Oz -o pruned_nomaf.vcf.gz

# in R
library(vcfR)
library(vegan)
library(tidyverse)
library(dplyr)

# load in vcf of neutral snps
blgu_vcf = read.vcfR("pruned_nomaf.vcf.gz")
blgu_geno <- extract.gt(blgu_vcf, element = "GT")
blgu_geno <- apply(blgu_geno, 2, as.character)

# convert vcf matrix to more useable form
G1f <- matrix(NA, nrow = nrow(blgu_geno), ncol = ncol(blgu_geno))
G1f[blgu_geno %in% c("0/0", "0|0")] <- 0
G1f[blgu_geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G1f[blgu_geno %in% c("1/1", "1|1")] <- 2

colnames(G1f) <- colnames(blgu_geno)

meta = read_csv("../analysis/locs.csv") %>%
  filter(Sample %in% colnames(blgu_geno))%>%
  base::split(f = as.factor(.$pop.code)) %>%
  map(pull,Sample)

# again with pops removed
keep_pops <- c("COA","COI","GDI","ICE","MLI","NAI","NFL","PLI","QIK","SF")
meta <- meta[keep_pops]

geno_avgs = imap_dfc(meta,function(Sample,pop.code){
    out = G1f[,Sample] %>%
    rowMeans(na.rm = TRUE)/2
    tibble(!!pop.code := out)
    #rowwise %>%
    #summarize(!!pop := mean(everything(),na.rm = TRUE))
})

# remove NAs by row not column
geno_avgs <- geno_avgs %>% filter(if_all(everything(), ~ !is.na(.)))
#### I WANT SNPs COLLASPED VARIATION ONTO POPULATIONS IS EVERY POPULATION SHOULD HAVE so need to tranpose
geno_avgs <- t(geno_avgs)

pca <- rda(geno_avgs, scale=T)
#Column scaling = divide each SNP (column) by its standard deviation.This is equivalent to PCA on a correlation matrix, not a covariance matrix.You do this when SNPs have different variances or when alleles are on different scales.
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")
#Compute explained variance
explained_variance <- (pca$CA$eig / sum(pca$CA$eig)) * 100

PCs <- scores(pca, choices=c(1:2), display="sites", scaling=0)
PopStruct <- PCs
PopStruct <- data.frame(Population = rownames(geno_avgs), PCs)

# rearrnge rows so in same order as individual PCA
PopStruct <- PopStruct[c(9, 1:8, 10:nrow(PopStruct)),]

b <- ggplot(PopStruct, aes(PC1, PC2, colour = Population)) +
  geom_point(size = 3, position = position_jitter(width = 0.02)) +

  scale_colour_manual(
    values = c("#1B9E77", "#7570B3", "#E7298A", "#66A61E",
               "#A6761D", "#666666", "#1F78B4", "#CAB2D6", "#33A02C",
               "#E31A1C", "#A6CEE3"),
    labels = c(
      "Akpaqarvik",
      "Hudson Bay",
      "Country Isl.",
      "Gulf of Maine",
      "Iceland",
      "Middle Lawn Isl.",
      "Labrador Shelf",
      "Conception Bay",
      "Qikiqtarjuaq",
      "Sud Fulgoya"
    ),
    name = "Population"
  ) +

  labs(
    x = "PC1 (14.05%)",
    y = "PC2 (13.25%)"
  ) +

  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 14)
  )
ggsave("allelefreq_pca_final.png",
       plot = b,
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)

# do Kaiser-Guttman criterion for retaining pcs for RDA and variance partitioning
ev <- pca$CA$eig
n <- length(ev)
# Save as PNG

png("kaiser.png", width = 800, height = 600)

# Create the plot
a <- barplot(ev,
             main = "",
             col = "grey",
             las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright",
       legend = "Average eigenvalue",
       lwd = 2,
       col = "red3",
       bty = "n")

dev.off()  # Close the file device

ggsave("kaiser.png",
       plot = a,
       width = 8,
       height = 6,
       units = "in",
       dpi = 300)
# first four PCs are over the average

## save PopStruct

PopStruct <- scores(pca, choices=c(1:4), display="sites", scaling=0)
PopStruct <- data.frame(Population = rownames(geno_avgs), PCs)


saveRDS(PopStruct, file = "PopStruct.rds")
