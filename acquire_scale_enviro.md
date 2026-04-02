# Download and Scale CHELSEA data

### download global data for each variable of interest
```
mkdir current
cd current
for VAR in {1..19}; do wget https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio${VAR}_1981-2010_V.2.1.tif; done

mkdir future
cd future
for VAR in {1..19}; do wget https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2071-2100/GFDL-ESM4/ssp585/bio/CHELSA_bio${VAR}_2071-2100_gfdl-esm4_ssp585_V.2.1.tif; done
```
### adjust file names
```
cd chelsa/current

mv CHELSA_bio1_1981-2010_V.2.1.tif CHELSA_bio01_1981-2010_V.2.1.tif
mv CHELSA_bio2_1981-2010_V.2.1.tif CHELSA_bio02_1981-2010_V.2.1.tif
mv CHELSA_bio3_1981-2010_V.2.1.tif CHELSA_bio03_1981-2010_V.2.1.tif
mv CHELSA_bio4_1981-2010_V.2.1.tif CHELSA_bio04_1981-2010_V.2.1.tif
mv CHELSA_bio5_1981-2010_V.2.1.tif CHELSA_bio05_1981-2010_V.2.1.tif
mv CHELSA_bio6_1981-2010_V.2.1.tif CHELSA_bio06_1981-2010_V.2.1.tif
mv CHELSA_bio7_1981-2010_V.2.1.tif CHELSA_bio07_1981-2010_V.2.1.tif
mv CHELSA_bio8_1981-2010_V.2.1.tif CHELSA_bio08_1981-2010_V.2.1.tif
mv CHELSA_bio9_1981-2010_V.2.1.tif CHELSA_bio09_1981-2010_V.2.1.tif

cd ../future

mv CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio01_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio2_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio02_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio3_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio03_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio4_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio04_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio5_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio05_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio6_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio06_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio7_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio07_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio8_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio08_2071-2100_gfdl-esm4_ssp585_V.2.1.tif
mv CHELSA_bio9_2071-2100_gfdl-esm4_ssp585_V.2.1.tif CHELSA_bio09_2071-2100_gfdl-esm4_ssp585_V.2.1.tif

```
### calculate means and sd for scaling

```
names(future_stack) = names(current_stack)

# all_means = list()
# all_sds = list()
# all_stdevs = list()
# all_ns = list()

i = slurm
print(slurm)
print("at layer stack")
layer_stack = c(current_stack[[i]],future_stack[[i]])
print("at global min")
global_min = min(global(layer_stack,"min",na.rm=TRUE))
print("at global n")

global_n = sum(global(layer_stack, "notNA"))
print(global_n)

print("at global mean")
global_mean = mean(as.vector(unlist(global(layer_stack, "mean", na.rm=TRUE))))
print("at numerator raster")
numerator_raster = (layer_stack-global_mean)^2
print("at globals section")
global_numerator = sum(global(numerator_raster, "sum", na.rm = TRUE))
global_std = sqrt(global_numerator/global_n)
global_sd = sqrt(global_numerator/(global_n-1))
print("at append section")
all_means = append(all_means,global_mean)
all_sds = append(all_sds,global_sd)
all_stdevs = append(all_stdevs,global_std)
all_ns = append(all_ns,global_n)

saveRDS(global_n, paste0(for_print,"current_ssp126_ns.rds"))
saveRDS(global_mean,paste0(for_print,"current_ssp126_means.rds"))
saveRDS(global_sd,paste0(for_print,"current_ssp126_sds.rds"))
saveRDS(global_std,paste0(for_print,"current_ssp126_stdevs.rds"))

```
## extract scaler values and put in csv
```
bioclim_id = seq(1:19)
mean_vec = c()
sd_vec = c()
stdev_vec = c()
ns_vec = c()

for(i in 1:19){
  # set in to char
  idx = as.character(i)

  # get stat for each layer
  mean_foradd = readRDS(paste0(idx,"ssp370_means.rds"))
  ns_foradd = readRDS(paste0(idx,"ssp370_ns.rds"))
  stdev_foradd = readRDS(paste0(idx,"ssp370_stdevs.rds"))
  sd_foradd = readRDS(paste0(idx,"ssp370_sds.rds"))

  # merge into single vector
  ns_vec = c(ns_vec,ns_foradd)
  stdev_vec = c(stdev_vec,stdev_foradd)
  mean_vec = c(mean_vec,mean_foradd)
  sd_vec = c(sd_vec ,sd_foradd)
}

# combine all vectors into df
chelsascaling_df = data.frame(bioclim_id = bioclim_id, means=mean_vec, sd = sd_vec, stdev = stdev_vec, n = ns_vec)

# save
write.csv(chelsa_scaling_df,"./chelsa/chelsa_ssp370_future.csv", row.names=F)
```
### extract data for each colony locations, scale data using scaler values. Combine terrestrial data with marine from bio-Oracle
```
library(terra)
library(tidyverse)

# get paths for chelsea variables (1-19)
names_chel_cur = dir("../climate_offset_scaling/chelsa/current/")
names_chel_cur = paste0("../climate_offset_scaling/chelsa/current/",names_chel_cur)

names_chel_fut = dir("../climate_offset_scaling/chelsa/future/")
names_chel_fut = paste0("../climate_offset_scaling/chelsa/future/",names_chel_fut)

# get paths for ocean variables (same order as scaling, alpha order)

current_biooracle = c("../climate_offset_scaling/bioOracle/current/chlorophyll_max_2000.tif", "../climate_offset_scaling/bioOracle/current/oceantemp_mean_2000.tif", "../climate_offset_scaling/bioOracle/current/oxygen_mean_2000.tif", "../climate_offset_scaling/bioOracle/current/phytoplankton_max_2000.tif", "../climate_offset_scaling/bioOracle/current/salinity_mean_2000.tif", "../climate_offset_scaling/bioOracle/current/seaice_conc_max_2000.tif", "../climate_offset_scaling/bioOracle/current/seaicethick_max_2000.tif","../climate_offset_scaling/bioOracle/current/seawaterspeed_range_2000.tif")
future_biooracle = c("../climate_offset_scaling/bioOracle/future/chlorophyll_max_2100.tif", "../climate_offset_scaling/bioOracle/future/oceantemp_mean_2100.tif", "../climate_offset_scaling/bioOracle/future/oxygen_mean_2100.tif", "../climate_offset_scaling/bioOracle/future/phytoplankton_max_2100.tif", "../climate_offset_scaling/bioOracle/future/salinity_mean_2100.tif", "../climate_offset_scaling/bioOracle/future/seaice_conc_max_2100.tif", "../climate_offset_scaling/bioOracle/future/seaicethick_max_2100.tif","../climate_offset_scaling/bioOracle/future/seawaterspeed_range_2100.tif")

current_stack_land = rast(names_chel_cur)
current_stack_ocean = rast(current_biooracle)

future_stack_land = rast(names_chel_fut)
future_stack_ocean = rast(future_biooracle)

# use colony locations to get 10 km bounding circles
locs = read.csv("./guillemot-colony-locations.csv",header=T)
locs2 = locs[,c(7,6)]
circles_land = vect(as.matrix(locs2),atts = data.frame(locs$pop_code),crs=crs(current_stack_land[[1]])) %>% buffer(10000)
# save these circles for future use
# saveRDS(circles_land, "guillemot_colony_land_circle_shape.rds")

# get 60 km bounding circles for oceanic variables
circles_ocean = vect(as.matrix(locs2),atts = data.frame(locs$pop_code),crs=crs(current_stack_ocean[[1]])) %>% buffer(60000)
# save these circles for future use
# saveRDS(circles_ocean, "guillemot_colony_ocean_circle_shape.rds")

# extract current and future land data
land_dat_from_circles_cur = terra::extract(current_stack_land,circles_land)
land_dat_from_circles_fut = terra::extract(future_stack_land,circles_land)

# extract current and future ocean data
ocean_dat_from_circles_cur = terra::extract(current_stack_ocean,circles_ocean)
ocean_dat_from_circles_fut = terra::extract(future_stack_ocean,circles_ocean)

# summarize land data
chel_land_current = land_dat_from_circles_cur %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
chel_land_current$colony = locs$pop_code
chel_land_future = land_dat_from_circles_fut  %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
chel_land_future$colony = locs$pop_code

# summarize ocean data
ocean_current = ocean_dat_from_circles_cur %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
ocean_current$colony = locs$pop_code
ocean_future = ocean_dat_from_circles_fut  %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
ocean_future$colony = locs$pop_code

ocean_land_current = inner_join(chel_land_current[,c(-1)],ocean_current[,c(-1)],by="colony")

ocean_land_future = inner_join(chel_land_future[,c(-1)],ocean_future[,c(-1)],by="colony")

# load scaling constants for chelsa (current_future version)

land_scalers = read_delim("./chelsa/chelsa_ssp545_future.csv")

ocean_scalers = read_delim("./bioOracle/biooracle_scalers_ssp545_future.csv")

land_ocean_means = c(land_scalers$means,ocean_scalers$means)
land_ocean_stdev = c(land_scalers$stdev,ocean_scalers$stdev)

land_ocean_scalers = data.frame(means = land_ocean_means, stdev=land_ocean_stdev )

ocean_land_current_scaled = data.frame(scale(as.matrix(ocean_land_current[,c(-20)]), center = as.numeric(land_ocean_means), scale = as.numeric(land_ocean_stdev)))
ocean_land_current_scaled$colony = locs$pop_code

ocean_land_future_scaled = data.frame(scale(as.matrix(ocean_land_future[,c(-20)]), center = as.numeric(land_ocean_means), scale = as.numeric(land_ocean_stdev)))
ocean_land_future_scaled$colony = locs$pop_code

# saveRDS(ocean_land_current_scaled,"chelsa_biooracle_current_scaled.rds")

saveRDS(ocean_land_future_scaled,"REDO_chelsa_biooracle_ssp585_scaled.rds")

# save scaled chelsa and biooracle separately
ocean_current_scaled = data.frame(scale(as.matrix(ocean_current[,c(-1,-10)]), center = as.numeric(ocean_scalers$means), scale = as.numeric(ocean_scalers$stdev)))
ocean_current_scaled$colony = locs$pop_code
ocean_future_scaled = data.frame(scale(as.matrix(ocean_future[,c(-1,-10)]), center = as.numeric(ocean_scalers$means), scale = as.numeric(ocean_scalers$stdev)))
ocean_future_scaled$colony = locs$pop_code

# saveRDS(ocean_current_scaled,"biooracle_current_scaled.rds")
saveRDS(ocean_future_scaled,"biooracle_ssp585_scaled.rds")

chelsa_current_scaled_2tp = data.frame(scale(as.matrix(chel_land_current[,c(-1,-21)]), center = as.numeric(land_scalers$means), scale = as.numeric(land_scalers$stdev)))
chelsa_current_scaled_2tp$colony = locs$pop_code

chelsa_future_scaled_2tp = data.frame(scale(as.matrix(chel_land_future[,c(-1,-21)]), center = as.numeric(land_scalers$means), scale = as.numeric(land_scalers$stdev)))
chelsa_future_scaled_2tp$colony = locs$pop_code

# saveRDS(chelsa_current_scaled_2tp,"chelsa_current_scaled.rds")
saveRDS(chelsa_future_scaled_2tp,"chelsa_ssp585_scaled.rds")

# extents don't match perfectly across data sets, use resample to fix this
current_stack_ocean_resamp = resample(current_stack_ocean,current_stack_land, method="near")

future_stack_ocean_resamp = resample(future_stack_ocean,future_stack_land, method="near")
current_stack_all = c(current_stack_land, current_stack_ocean_resamp)

future_stack_all = c(future_stack_land, future_stack_ocean_resamp)


