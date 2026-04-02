# Download and Scale CHELSEA data

### Download global data for each variable of interest
```
mkdir current
cd current
for VAR in {1..19}; do wget https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio${VAR}_1981-2010_V.2.1.tif; done

mkdir future
cd future
for VAR in {1..19}; do wget https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2071-2100/GFDL-ESM4/ssp126/bio/CHELSA_bio${VAR}_2071-2100_gfdl-esm4_ssp126_V.2.1.tif; done
```
### Adjust file names
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
### Calculate means and standard deviations for scaling

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

saveRDS(global_n, paste0(for_print,"current_future_ssp126_ns.rds"))
saveRDS(global_mean,paste0(for_print,"current_future_ssp126_means.rds"))
saveRDS(global_sd,paste0(for_print,"current_future_ssp126_sds.rds"))
saveRDS(global_std,paste0(for_print,"current_future_ssp126_stdevs.rds"))

```
## Extract scaler values and place in a csv
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
  mean_foradd = readRDS(paste0(idx,"current_future_ssp126_means.rds"))
  ns_foradd = readRDS(paste0(idx,"current_future_ssp126_means.rds"))
  stdev_foradd = readRDS(paste0(idx,"current_future_ssp126_means.rds"))
  sd_foradd = readRDS(paste0(idx,"current_future_ssp126_means.rds"))

  # merge into single vector
  ns_vec = c(ns_vec,ns_foradd)
  stdev_vec = c(stdev_vec,stdev_foradd)
  mean_vec = c(mean_vec,mean_foradd)
  sd_vec = c(sd_vec ,sd_foradd)
}

# combine all vectors into df
chelsascaling_df = data.frame(bioclim_id = bioclim_id, means=mean_vec, sd = sd_vec, stdev = stdev_vec, n = ns_vec)

# save
write.csv(chelsa_scaling_df,"./chelsa/chelsa_current_future_ssp126.csv", row.names=F)
```
### Extract data for each colony locations, scale data using scaler values. 
```
library(terra)
library(tidyverse)

# get paths for chelsea variables (1-19)
names_chel_cur = dir("../climate_offset_scaling/chelsa/current/")
names_chel_cur = paste0("../climate_offset_scaling/chelsa/current/",names_chel_cur)

names_chel_fut = dir("../climate_offset_scaling/chelsa/future/")
names_chel_fut = paste0("../climate_offset_scaling/chelsa/future/",names_chel_fut)

current_stack_land = rast(names_chel_cur)

future_stack_land = rast(names_chel_fut)

# use colony locations to get 10 km bounding circles
locs = read.csv("./guillemot-colony-locations.csv",header=T)
locs2 = locs[,c(7,6)]
circles_land = vect(as.matrix(locs2),atts = data.frame(locs$pop_code),crs=crs(current_stack_land[[1]])) %>% buffer(10000)
# save these circles for future use
# saveRDS(circles_land, "guillemot_colony_land_circle_shape.rds")

# extract current and future land data
land_dat_from_circles_cur = terra::extract(current_stack_land,circles_land)
land_dat_from_circles_fut = terra::extract(future_stack_land,circles_land)

# summarize land data
chel_land_current = land_dat_from_circles_cur %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
chel_land_current$colony = locs$pop_code
chel_land_future = land_dat_from_circles_fut  %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
chel_land_future$colony = locs$pop_code

# load scaling constants for chelsa (current_future version)

land_scalers = read_delim("./chelsa/chelsa_current_future_ssp126.csv")

chelsa_current_scaled_2tp = data.frame(scale(as.matrix(chel_land_current[,c(-1,-21)]), center = as.numeric(land_scalers$means), scale = as.numeric(land_scalers$stdev)))
chelsa_current_scaled_2tp$colony = locs$pop_code

chelsa_future_scaled_2tp = data.frame(scale(as.matrix(chel_land_future[,c(-1,-21)]), center = as.numeric(land_scalers$means), scale = as.numeric(land_scalers$stdev)))
chelsa_future_scaled_2tp$colony = locs$pop_code

saveRDS(chelsa_current_scaled_2tp,"chelsa_current_scaled.rds")
saveRDS(chelsa_future_scaled_2tp,"chelsa_ssp126_scaled.rds")



