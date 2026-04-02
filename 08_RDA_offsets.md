Now it's time to run offsets! I started with the RDA offset method, again outlined by Capblanc and Forestor (2021), tutorial found here: https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd

I ran the RDA offset with all SNPs. Lind and Lotterhos (2024) found that there is little difference in offsets generated using all SNPs from across the genome and putativly adaptive SNPs. Local adaptation is most likely harboured in genome-wide variation; traits related to thermal tolerance for instance tend to be polygenic with many SNPs having small effects (Bay et al. 2017). 

Since I did not have samples from many regions across guillemots range, I didn't want to extrapolate the relationship beyond my known colonies. So I did offets per colony, and only across the range for SSP585. 

To calculate genomic offsets using the RDA model, I adjusted the genomic_offset function by Lind and Lotterhos (2021) https://github.com/DrK-Lo/MVP-offsets/blob/main/01_src/MVP_climate_outlier_RDA_offset.R. I removed the first part of the function which scales and centers the raster of environmental data for the species range and converts the raster to points since I already did this and have my data saved as an rds

Below are the functions required to run just the genomic offset - taken from github linked above - which I placed in get_offset_rda_lind.r
```
get_pc_data = function(env_pres){
    #--------------------------------------------------------#
    # Retrieve PCA loadings for individual data pooled data.
    #
    # Parameters
    # ----------
    # env_pres
    #     - data.frame; rownames are individual or popIDs
    #
    # Notes
    # -----
    # - returns only the first two axes
    #     - assumed by train_outlier_rda
    #--------------------------------------------------------#

    if (ind_or_pooled == 'ind'){
        # retrieve PCA - created from MVP-NonClinalAF/src/c-AnalyzeSimOutput.R
#         pcfile = paste0(
#             slimdir,
#             sprintf("/%s_pca.RDS", seed)
#         )
#         pc_object = readRDS(pcfile)
#         pc_object@projDir = paste0(slimdir, '/')  # it is really dumb the LEA needs to hardcode paths,
#                                                         # just incorporate into one file duh!
        subsetfile = paste0(slimdir, '/', seed, '_Rout_ind_subset.txt')
        subset = read.table(subsetfile)
        stopifnot(all(rownames(env_pres) == subset[ , 'indID']))
        pcdata = subset[ , c('PC1', 'PC2')]

    } else {
        # retrieve PCA - created in MVP-Offsets/01_src/MVP_pooled_pca.R
        pcfile = paste0(
            outerdir,
            sprintf("/pca/pca_output/%s_Rout_Gmat_sample_maf-gt-p01_GFready_pooled_all_pca.RDS", seed)
        )
        pc_object = readRDS(pcfile)
        pcdata = pc_object$projections[ , c(1, 2)]
        colnames(pcdata) = c('PC1', 'PC2')
    }

    rownames(pcdata) <- rownames(env_pres)  # rownames are individual or pop IDs

    return(pcdata)    
}

train_outlier_rda = function(outlier_snps, rda, env_pres, ntraits=2){
    #----------------------------------------------------------------------------------#
    # Train new RDA on outliers inferred from script input args.
    # 
    # Parameters
    # ----------
    # outlier_snps
    #     - output from subset_snps
    # rda
    #     - object loaded from rda_file input arg to this script
    # ntraits
    #     - the number of traits (environments) under selection in the simulation seed
    # 
    # Notes
    # -----
    # - assumes any structure correction is done with first two PCs
    #     - see get_pc_data()
    #----------------------------------------------------------------------------------#
    cat(sprintf('\n\nTraining outlier RDA ...'))

    # decide if I need to estimate a new RDA
    if (use_RDA_outliers == FALSE & ntraits==2){
        # just use the original RDA object - this could be structure-corrected or not
        outlier_rda = rda

    } else {  # estimate a new RDA object with `outlier_snps`

        if (grepl('_RDA_structcorr.RDS', basename(rda_file))){  # this is structure corrected
            # get PCs
            pcdata = get_pc_data(env_pres)

            # replicate Variables object from capblancq & forester
            Variables = data.frame(pcdata, env_pres)

            # get RDA
            if (ntraits == 1){  # use only the adaptive trait (temp)
                outlier_rda = rda(outlier_snps ~ temp + Condition(PC1 + PC2), Variables)

            } else if (ntraits == 2) {  # use both traits (sal + temp)
                outlier_rda = rda(outlier_snps ~ sal + temp + Condition(PC1 + PC2), Variables)
                
            } else if (ntraits == 4){
                outlier_rda = rda(outlier_snps ~ sal + temp + ISO + PSsd + Condition(PC1 + PC2), Variables)
                                  
            } else if (ntraits == 5){
                outlier_rda = rda(outlier_snps ~ sal + temp + ISO + PSsd + TSsd + Condition(PC1 + PC2), Variables)

            } else {
                stop(sprintf('ntraits unexpected : %s', ntraits))
            }

        } else {  # this is not structure corrected
            # replicate Variables object from capblancq & forester
            Variables = data.frame(env_pres)

            # get RDA
            if (ntraits == 1){  # use only the adaptive trait (temp)
                outlier_rda = rda(outlier_snps ~ temp, Variables)

            } else if (ntraits == 2) {  # use both traits (sal + temp)
                outlier_rda = rda(outlier_snps ~ sal + temp, Variables)
                
            } else if (ntraits == 4){
                outlier_rda = rda(outlier_snps ~ sal + temp + ISO + PSsd, Variables)

            } else if (ntraits == 5){
                outlier_rda = rda(outlier_snps ~ sal + temp + ISO + PSsd + TSsd, Variables)

            } else {
                stop(sprintf('ntraits unexpected : %s', ntraits))
            }
        }
    }

    return(outlier_rda)
}

ensure_dataframe <- function(df, cols) {

  if (length(cols) == 0) {
    stop("There aren't any cols passed to `ensure_dataframe`")
  }

  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in data frame: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df <- df[, cols, drop = FALSE]

  stopifnot(is.data.frame(df))

  return(df)
}

genomic_offset = function(RDA, env_pres, env_fut, method="loadings"){
    #------------------------------------------------------------------------#
    # Predict genomic offset using RDA.
    #
    # Notes
    # -----
    # - this was modified from capblancq and forester to allow MVP data
    #     - mostly this means using non-raster transformations or
    #       building in code that can accept either one or two envs (which
    #       cause the RDA to have either one or two axes, respectively)
    #     - additionally, `K` in original function was replaced so that
    #       at most 2 RDA axes were used (ie K=2). In our case, when
    #       we used only one trait (temp) there was only one resulting
    #       RDA axis, and therefore replace `K` with `ncol(env_pres)`.
    #       When using two traits (temp + sal), K=2, and therefore
    #       two RDA axes are used for offset estimation.
    # - in the original script, 'predict' could be passed to kwarg `method`
    #     - this code was removed because I didn't see it used in their 
    #       example RDA_landscape_genomics.html
    #
    # Parameters
    # ----------
    # RDA
    #     - RDA object
    # env_pres
    #     - data.frame; current environmental values
    #       rows for individuals or pops, columns for envs
    # env_fut
    #     - data.frame; future environmental values
    #       rows for individuals or pops, columns for envs
    # K
    #     - int; number of RDA axes to use in offset calculation
    #     - set to K=n_env because of convention used in MVP sims/sim processing
    #------------------------------------------------------------------------#
    
    K = ncol(env_pres)

    # Predicting pixels genetic component based on the loadings of the variables
    if(method == "loadings"){
        # Projection for each RDA axis
        Proj_pres = list()
        Proj_fut = list()
        Proj_offset = list()
        for(i in 1:K){
            # Current climates
            ras_pres = as.vector(
                apply(
                    ensure_dataframe(var_env_proj_pres, colnames(env_pres)),
                    1,
                    function(x) sum(x * RDA$CCA$biplot[colnames(env_pres) , i])
                )
            )
            Proj_pres[[i]] <- ras_pres
            names(Proj_pres)[i] <- paste0("RDA", as.character(i))

            # Future climates
            ras_fut <- as.vector(
                apply(
        #             var_env_proj_fut[ , names(RDA$CCA$biplot[ , i])],
                    ensure_dataframe(var_env_proj_fut, colnames(env_pres)),
                    1,
        #             function(x) sum(x * RDA$CCA$biplot[ , i])
                    function(x) sum(x * RDA$CCA$biplot[colnames(env_pres) , i])
                )
            )
            Proj_fut[[i]] <- ras_fut
            names(Proj_fut)[i] <- paste0("RDA", as.character(i))

            # Single axis genetic offset 
            Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
            names(Proj_offset)[i] <- paste0("RDA", as.character(i))
        }
    } else {
        stop('method must be "loading"')
    }

    # Weights based on axis eigen values
    weights <- RDA$CCA$eig/sum(RDA$CCA$eig)

    # Weighing the current and future adaptive indices based on the eigen values of the associated axes
    Proj_offset_pres_ <- do.call(
        cbind,
        lapply(1:K, function(x) {Proj_pres[[x]]})
    )

    Proj_offset_pres <- as.data.frame(
        do.call(
            cbind,
            lapply(1:K, function(x) {Proj_offset_pres_[ , x] * weights[x]})
        )
    )

    Proj_offset_fut_ <- do.call(
        cbind,
        lapply(1:K, function(x) {Proj_fut[[x]]})
    )

    Proj_offset_fut <- as.data.frame(
        do.call(
            cbind,
            lapply(1:K, function(x) {Proj_offset_fut_[ , x] * weights[x]})
        )
    )

    # Predict a global genetic offset, incorporating the K first axes weighted by their eigen values
    Proj_offset_global <- unlist(
        lapply(
            1:nrow(Proj_offset_pres),
            function(x) {dist(
                rbind(Proj_offset_pres[x, ],
                      Proj_offset_fut[x, ]),
                method = "euclidean"
            )}
        )
    )

    # Return projections for current and future climates for each RDA axis,
        # prediction of genetic offset for each RDA axis and a global genetic offset 
    ret = list(Proj_pres = Proj_pres,
               Proj_fut = Proj_fut,
               Proj_offset = Proj_offset,
               Proj_offset_global = Proj_offset_global,
               weights = weights[1:K])

    return(ret)
}
```
### Now use the functions to generate offsets for the colonies

Colonies without population correction
NOTE: I ran this three times, once with o2, siconc, alone and together
```
library(vegan)

source("./get_offset_rda_lind.r")

# load geno_avg_t and environmental variables for each colony - since data is already scaled I cut that part from the function, and current and future data is loaded twice as the same thing
geno_avg_t <- readRDS("./geno_avg_t_FINAL.rds")
all_pred <- readRDS("./all_pred_forRDA.rds")
env_pres <- readRDS("../chelsa_biooracle_current_scaled.rds")
var_env_proj_pres <- readRDS("../chelsa_biooracle_current_scaled.rds")
env_fut <- readRDS("../chelsa_biooracle_ssp126_scaled.rds")
var_env_proj_fut <- readRDS("../chelsa_biooracle_ssp126_scaled.rds")

# make sure environmental data is in the same order as genotype averages
pop_order <- c("COA", "COI", "GDI", "GRI", "ICE", "MLI", "NAI", "NFL",  "PLI", "QIK","SF", "SIG", "SOD")
env_pres$colony <- factor(env_pres$colony, levels = pop_order)
env_pres <- env_pres[order(env_pres$colony), ]

env_fut$colony <- factor(env_fut$colony, levels = pop_order)
env_fut <- env_pres[order(env_fut$colony), ]

var_env_proj_pres$colony <- factor(var_env_proj_pres$colony, levels = pop_order)
var_env_proj_pres <- var_env_proj_pres[order(var_env_proj_pres$colony), ]

var_env_proj_fut$colony <- factor(var_env_proj_fut$colony, levels = pop_order)
var_env_proj_fut <- var_env_proj_fut[order(var_env_proj_fut$colony), ]

# remove three colonies since they are not included in this analysis
env_pres <- env_pres[-c(12,13,14),]
env_fut <- env_fut[-c(12,13,14),]
var_env_proj_pres <- var_env_proj_pres[-c(12,13,14),]
var_env_proj_fut <- var_env_proj_fut[-c(12,13,14),]

# rename future data to same as current data otherwsie model won't recognize
colnames(env_fut)[colnames(env_fut) == "CHELSA_bio01_2071.2100_gfdl.esm4_ssp126_V.2.1"] <- "CHELSA_bio01_1981.2010_V.2.1"
colnames(env_fut)[colnames(env_fut) == "siconc_max_8"] <- "siconc_max_1"
colnames(var_env_proj_fut)[colnames(var_env_proj_fut) == "CHELSA_bio01_2071.2100_gfdl.esm4_ssp126_V.2.1"] <- "CHELSA_bio01_1981.2010_V.2.1"
colnames(var_env_proj_fut)[colnames(var_env_proj_fut) == "siconc_max_8"] <- "siconc_max_1"

RDA <- rda(
  geno_avg_t ~ 
    `CHELSA_bio01_1981.2010_V.2.1` + siconc_max_1,
  data = all_pred
)

rda_offset <- genomic_offset(RDA = RDA, env_pres = env_pres, env_fut = var_env_proj_fut, method="loadings")

saveRDS(rda_offset, file = "offsets_FINAL_si_rda_ssp126_allsnps.rds")

```
Offsets conditioning model with population structure
#NOTE: Ran this three times once with o2, siconc, and together
```
library(vegan)

source("./get_offset_rda_lind.r")

# load geno_avg_t and cropped environmental variables from across whole range
PopStruct <- readRDS("PopStruct.rds")

# load data - since data is already scaled I cut that part from the function, and current and future data is loaded twice as the same thing
geno_avg_t <- readRDS("./geno_avg_t_FINAL.rds")
all_pred <- readRDS("./all_pred_forRDA.rds")
env_pres <- readRDS("../chelsa_biooracle_current_scaled.rds")
var_env_proj_pres <- readRDS("../chelsa_biooracle_current_scaled.rds")
env_fut <- readRDS("../chelsa_biooracle_ssp585_scaled.rds")
var_env_proj_fut <- readRDS("../chelsa_biooracle_ssp585_scaled.rds")

# make sure environmental data is in the same order as genotype averages
pop_order <- c("COA", "COI", "GDI", "GRI", "ICE", "MLI", "NAI", "NFL",  "PLI", "QIK","SF", "SIG", "SOD")
env_pres$colony <- factor(env_pres$colony, levels = pop_order)
env_pres <- env_pres[order(env_pres$colony), ]

env_fut$colony <- factor(env_fut$colony, levels = pop_order)
env_fut <- env_pres[order(env_fut$colony), ]

var_env_proj_pres$colony <- factor(var_env_proj_pres$colony, levels = pop_order)
var_env_proj_pres <- var_env_proj_pres[order(var_env_proj_pres$colony), ]

var_env_proj_fut$colony <- factor(var_env_proj_fut$colony, levels = pop_order)
var_env_proj_fut <- var_env_proj_fut[order(var_env_proj_fut$colony), ]

# remove that three colonies since they are not included in this analysis
env_pres <- env_pres[-c(12,13,14),]
env_fut <- env_fut[-c(12,13,14),]
var_env_proj_pres <- var_env_proj_pres[-c(12,13,14),]
var_env_proj_fut <- var_env_proj_fut[-c(12,13,14),]

# rename future data to same as current data otherwsie model won't recognize
env_pres <- env_pres[, c("CHELSA_bio01_1981.2010_V.2.1", "o2_mean_1", "siconc_max_1")]
var_env_proj_fut <- var_env_proj_fut[, c("CHELSA_bio01_2071.2100_gfdl.esm4_ssp585_V.2.1", "o2_mean_8", "siconc_max_8")]
colnames(env_fut)[colnames(env_fut) == "CHELSA_bio01_2071.2100_gfdl.esm4_ssp585_V.2.1"] <- "CHELSA_bio01_1981.2010_V.2.1"
colnames(env_fut)[colnames(env_fut) == "siconc_max_8"] <- "siconc_max_1"
colnames(env_fut)[colnames(env_fut) == "o2_mean_8"] <- "o2_mean_1"
colnames(var_env_proj_fut)[colnames(var_env_proj_fut) == "CHELSA_bio01_2071.2100_gfdl.esm4_ssp585_V.2.1"] <- "CHELSA_bio01_1981.2010_V.2.1"
colnames(var_env_proj_fut)[colnames(var_env_proj_fut) == "siconc_max_8"] <- "siconc_max_1"
colnames(var_env_proj_fut)[colnames(var_env_proj_fut) == "o2_mean_8"] <- "o2_mean_1"

# combine environmetal and population structure data
all_pred <- cbind(all_pred, PopStruct[, c("PC1", "PC2")])

# amke RDA model
RDA <- rda(
  geno_avg_t ~ 
    `CHELSA_bio01_1981.2010_V.2.1` + o2_mean_1 + siconc_max_1+ Condition(PC1 + PC2),
  data = all_pred)

# offsets!
rda_offset <- genomic_offset(RDA = RDA, env_pres = env_pres, env_fut = var_env_proj_fut, method="loadings")

saveRDS(rda_offset, file = "conditioned_offsets_FINAL_si_o2_rda_ssp585_allsnps.rds")
```
### Plot the offsets!

By range for SSP5-8.5. Have to join with coordinates from current environmental data in order to plot
```
# r/4.4.0
library(ggplot2)
library(rnaturalearth)
library(sf)

env_pres <- readRDS("../current_13_clean_cropped.rds")
all <- readRDS("offsets_range_13_si_o2_rda_ssp585.rds")
env_pres_plot <- env_pres
env_pres_plot$Proj_offset_global <- all$Proj_offset_global

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# plot!
a <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.2) +
  geom_tile(
    data = env_pres_plot,
    aes(x = x, y = y, fill = Proj_offset_global)
  ) +
  coord_sf(
    xlim = range(env_pres_plot$x, na.rm = TRUE),
    ylim = range(env_pres_plot$y, na.rm = TRUE),
    expand = FALSE
  ) +
  scale_fill_viridis_c(name = "Genomic offset", na.value = "transparent") +
  theme_bw() +
    labs(x = "Longitude", y="Latitude")

# save high quality plots
ggsave("range_offsets_si_o2_RDA_585.png", plot = a, width = 8, height = 6, units = "in", dpi = 600)
```
For colonies - Have to join with colony coordinates to plot
```
r/4.4.0
library(ggplot2)
library(maps)

clim.points <- read.csv("../guillemot-colony-locations.csv")

# make sure colony points are in the same order as offsets 
pop_order <- c("COA", "COI", "GDI","ICE", "MLI", "NAI", "NFL",  "PLI", "QIK","SF", "SIG", "SOD", "GRI")
clim.points$pop_code <- factor(clim.points$pop_code, levels = pop_order)
clim.points <- clim.points[order(clim.points$pop_code), ]
coords <- clim.points[, c("lat", "long")] # extract coordinates
coords <- coords[-c(11,12,13,14),] # remove colonies not included

rda_offset <- readRDS("conditioned_offsets_FINAL_o2_rda_ssp585_allsnps.rds")
coords$offset <- rda_offset$Proj_offset_global
#coords$offset <- rda_offset

world <- map_data("world")

a <- ggplot() +
  geom_polygon(data = world,
               aes(long, lat, group = group),
               fill = "gray90",
               colour = "gray50",
               linewidth = 0.2) +
  
  geom_point(data = coords,
             aes(x = long, y = lat, fill = offset),
             shape = 21,
             size = 4,
             colour = "black") +
  
  scale_fill_viridis_c(name = "Offset Z-Score") +
coord_quickmap(
    xlim = c(-100, 60),   # longitude range
    ylim = c(42, 80)      # latitude range
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic()
  
ggsave("conditioned_RDA_o2_allsnps_585.png", plot = a, width = 8, height = 6, units = "in", dpi = 600)
```
