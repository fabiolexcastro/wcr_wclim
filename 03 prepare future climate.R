
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, usdm, raster, gtools, spocc, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))

## Raster data 
root <- '//ALLIANCEDFS.ALLIANCE.CGIAR.ORG/data_cluster19/GLOBAL/climate/Agroclimas/data/ipcc_6ar_wcl_downscaled/ssp_370/2050s'
dirs <- as.character(dir_ls(root, type = 'directory'))

## Solar radiation
srad <- dir_ls('//catalogue/workspace-cluster9/DATA/ET_SolRad', regexp = 'et_solrad')
srad <- mixedsort(srad)
srad <- rast(srad)

# To calculate bioclimatic variables  -------------------------------------

calc.bioc <- function(dire){
  
  dire <- dirs[1]
  
  ## To list the files
  cat('To process: ', dire, '\n')
  fles <- as.character(dir_ls(dir_ls(dire)))
  
  ## To read as a raster 
  prec <- rast(grep('prec', fles, value = T)) 
  tmin <- rast(grep('tmin', fles, value = T))
  tmax <- rast(grep('tmax', fles, value = T))
  
  ## To change the names
  names(prec) <- glue('prec_{1:12}')
  names(tmin) <- glue('tmin_{1:12}')
  names(tmax) <- glue('tmax_{1:12}')
  
  ## To calculate tavg
  tavg <- (tmax + tmin) / 2
  names(tavg) <- glue('tavg_{1:12}')
  
  ## Bilinear
  # srad <- terra::resample(srad, tmax, method = 'bilinear')
  
  # To calculate the ETP
  etps <- 0.0013 * 0.408 * srad * (tavg + 17) * (tmax - tmin - 0.0123 * prec) ^ 0.76
  etps <- etps * c(31,29,31,30,31,30,31,31,30,31,30,31)
  for(i in 1:12){etps[[i]][which.lyr(is.na(etps[[i]]))] <- 0}
  etps <- terra::crop(etps, wrld) %>% terra::mask(., wrld)
  names(etps) <- glue('etps_{c(paste0(0, 1:9), 10:12)}')
  
  
  
  
  
}
