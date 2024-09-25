
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, usdm, raster, spocc, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Raster data
prec <- rast('./tif/wc_5m/prec_wrld.tif')
tmin <- rast('./tif/wc_5m/tmin_wrld.tif')
tmax <- rast('./tif/wc_5m/tmax_wrld.tif')
bioc <- rast('./tif/wc_5m/bioc_wrld.tif')

srad <- dir_ls('//catalogue/workspace-cluster9/DATA/ET_SolRad', regexp = 'et_solrad')
srad <- mixedsort(srad)
srad <- rast(srad)

## Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
wrld <- vect(wrld)

# To calculate the average temperature ------------------------------------
tavg <- (tmax + tmin) / 2
names(tavg) <- glue('tavg_{1:12}')
srad <- terra::resample(srad, tmax, method = 'bilinear')

# To calculate the ETP  ---------------------------------------------------
etps <- 0.0013 * 0.408 * srad * (tavg + 17) * (tmax - tmin - 0.0123 * prec) ^ 0.76
etps <- etps * c(31,29,31,30,31,30,31,31,30,31,30,31)
for(i in 1:12){etps[[i]][which.lyr(is.na(etps[[i]]))] <- 0}
etps <- terra::crop(etps, wrld) %>% terra::mask(., wrld)
names(etps) <- glue('etps_{c(paste0(0, 1:9), 10:12)}')

# To write the ETP  -------------------------------------------------------
terra::writeRaster(x = etps, filename = './tif/wc_5m/etp_harg.tif', overwrite = TRUE)

# To extract by mask  -----------------------------------------------------
prec <- terra::crop(prec, wrld) %>% terra::mask(., wrld)
tmin <- terra::crop(tmin, wrld) %>% terra::mask(., wrld)
tmax <- terra::crop(tmax, wrld) %>% terra::mask(., wrld)
tavg <- terra::crop(tavg, wrld) %>% terra::mask(., wrld)

# To create a stack  ------------------------------------------------------
stck <- c(prec, tmin, tavg, tmax, etps)
tble <- terra::as.data.frame(stck, xy = T)
nrow(tble)
tble <- drop_na(tble)
nrow(tble)
stck <- terra::rast(tble, type = 'xyz')

# Drop the layers for each variable ---------------------------------------
prec <- stck[[grep('prec', names(stck), value = F)]]
tmin <- stck[[grep('tmin', names(stck), value = F)]]
tmax <- stck[[grep('tmax', names(stck), value = F)]]
tavg <- stck[[grep('tavg', names(stck), value = F)]]
etps <- stck[[grep('etps', names(stck), value = F)]]

etps <- map(.x = 1:nlyr(etps), .f = function(i){
  r <- etps[[i]]
  r[r < 0] <- 0
  return(r)
}) %>% 
  reduce(., c)


# To make the bio ETPs ----------------------------------------------------
source('bioclimatic funcitons.R')

# Convert to matrix and create the main matrix -----------------------------
etps.mt <- as.matrix(etps)
prec.mt <- as.matrix(prec)
tavg.mt <- as.matrix(tavg)

## Make a compilate matrix  ------------------------------------------------
etpr <- cbind(etps.mt, prec.mt, tavg.mt)

## To craate the etp bioclimatic variables ---------------------------------
etbi <- apply(etpr, 1, etpvars)
etbi <- t(etbi)

zero <- prec[[1]] * 0 + 1
names(zero) <- 'zero'

## Matrix to raster 
etrs <- purrr::map(.x = 1:ncol(etbi), .f = function(k){
  cat(k, '\n')
  lyer <- zero
  values(lyer) <- etbi[,k]
  cat('Done!\n')
  return(lyer)
})
etrs <- reduce(etrs, c)
names(etrs) <- glue('bioc_{21:29}')

terra::writeRaster(x = etrs, filename = './tif/wc_5m/bioc-etps_wrld.tif')

# Join both rasters into only one -----------------------------------------
names(bioc) <- glue('bioc_{1:19}')
bioc <- terra::crop(bioc, wrld)
bioc <- terra::mask(bioc, wrld)
bioc.all <- c(bioc, etrs)
terra::writeRaster(x = bioc.all, filename = './tif/wc_5m/biocc-all_wrld.tif', overwrite = TRUE)

## Points dataset 
pnts <- read_csv('./tbl/points/points background.csv')
spcs <- unique(pnts$nombre)
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))
pnts <- dplyr::select(pnts, pb:Latitude)

# To extract the climate for the points -----------------------------------

## Worldclim 
pnts.vles <- terra::extract(bioc.all, pnts[,c('Longitude', 'Latitude')])
pnts.vles <- cbind(pnts, pnts.vles)
pnts.vles <- as_tibble(pnts.vles)
pnts.vles <- dplyr::select(pnts.vles, -ID)

write.csv(pnts.vles, './tbl/values/points_worldclim.csv', row.names = FALSE)

# Cleaning the points -----------------------------------------------------
make.clea <- function(spce){
  
  # spce <- 'Theobroma cacao'
  
  ## To filter the specie
  cat('To process: ', spce, '/n')
  pnt <- filter(pnts, nombre == spce & pb == 1)
  bck <- filter(pnts, nombre == spce & pb == 0)
  
  ## To clean the coordinates
  cln <- clean_coordinates(x = as.data.frame(pnt), lon = 'Longitude', lat = 'Latitude', species = 'nombre', tests = c('capitals', 'centroids', 'equal', 'zeros', 'institutions', 'seas'))
  cla <- cln[cln$.summary,]
  cla <- dplyr::select(cla, pb:Latitude)
  cla <- as_tibble(cla)
  
  ## Join with the background 
  cla <- rbind(cla, bck)
  
  ## Return
  cat('Done!/n')
  return(cla)
  
}
pnts.clea <- make.clea(spce = 'Theobroma cacao')
write.csv(pnts.clea, './tbl/values/cocoa_worldclim_clean.csv', row.names = FALSE)

