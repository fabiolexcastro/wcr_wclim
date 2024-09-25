

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, usdm, raster, spocc, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

source('//catalogue/workspace-cluster9/2024/TREECROPSGHANA/02 MAKE MODEL/FunctionsRFclustering.R')

# Load data ---------------------------------------------------------------
lbls <- tibble(specie = c('Theobroma cacao', 'Anacardium occidentale', 'Cocos nucifera', 'Elaeis guineensis', 'Hevea brasiliensis', 'Vitellaria paradoxa', 'Mangifera indica'), common = c('Cocoa', 'Cashew', 'Coconut', 'Oil palm', 'Rubber', 'Shea', 'Mango'))
pnts <- read_csv('./tbl/values/cocoa_worldclim_clean.csv')
bioc <- rast('./tif/wc_5m/biocc-all_wrld.tif')
names(bioc)

# To extract the values ---------------------------------------------------
pnts.vles <- as_tibble(cbind(pnts, terra::extract(bioc, pnts[,c(3, 4)])))
pnts.vles <- dplyr::select(pnts.vles, -ID)
write.csv(pnts.vles, './tbl/values/cocoa_worldclim_clean_swd.csv')

# Make VIF ----------------------------------------------------------------
make.vifs <- function(tble, spce){
  
  # tble <- pnts.vles
  # spce <- 'Theobroma cacao'
  
  ## To filter the specie
  cat('To process: ', spce, '\n')
  pnt <- tble
  occ <- filter(pnt, pb == 1)
  bck <- filter(pnt, pb == 0)
  
  ## To make the VIF 
  vif <- usdm::vifstep(x = as.data.frame(occ[,5:32]), th = 10)
  vrs <- vif@results$Variables
  
  ## To select the variables 
  rsl <- dplyr::select(pnt, pb:Latitude, all_of(vrs))
  
  ## Finish
  cat('Done!\n')
  return(rsl)
  
}
pnts.vifs <- make.vifs(tble = pnts.vles, spce = 'Theobroma cacao')

# Make clustering ---------------------------------------------------------
make.clst.occr <- function(tble, spce){
  
  ## To start the process
  cat('To star the process: ', spce, '\n')
  nme <- filter(lbls, specie == spce) %>% pull(2)
  pnt <- filter(tble, nombre == spce & pb == 1)
  
  ## No Forest / No Trees
  no.forest <- 25
  no.trees <- 100
  nVars <- 8
  
  ## Clustering presences 
  occ <- pnt
  occ.mtx <- occ[,5:ncol(occ)]
  occ.dst <- RFdist(occ.mtx, mtry1 = nVars, no.trees, no.forest, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  occ.ncl <- 5
  occ.lrf <- pamNew(occ.dst$cl1, occ.ncl)
  occ.cls <- hclust(as.dist(occ.dst$cl1), method = 'ward.D2')
  occ.cld <- cbind(pb = as.factor(occ.lrf), occ[2:ncol(occ)])
  occ.clp <- cbind(occ, cluster = occ.lrf) %>% na.omit() %>% as_tibble()
  
  ## To save the results
  dir <- glue('./rData/{nme}/run_'); dir_create(dir)
  save(occ.mtx, file = glue('{dir}/datRF.rData'))
  save(occ.cls, file = glue('{dir}/clusterdata.rData'))
  save(occ, occ.clp, occ.ncl, occ.lrf, file = glue('{dir}/clustereddata.rData'))
  save(occ.cld, file = glue('{dir}/occ_cld.rData'))
  cat('Done!\n')
  
}