# Run GOLUM intercomparison analyses
# BBL November 21, 2017

library(dplyr)
library(readr)
library(tidyr)
library(raster)
library(ggplot2)
theme_set(theme_bw())
library(ncdf4)


# Compute the mean values for all pools by latitude

WHICH_CDO <- "/usr/local/bin/cdo"
GOLUM_FOLDER <- "input/golum/"

tibble(file = list.files(GOLUM_FOLDER, pattern = "*ann.nc", full.names = TRUE)) %>% 
  separate(file, sep = "_", remove = FALSE, into = c("analysis", "group", "run", "forcing", "LUC", "years", "time")) %>% 
  separate(years, into = c("minyear", "maxyear"), convert = TRUE, remove = FALSE) ->
  filelist

results <- list()
hr_data <- list()
for(f in seq_len(nrow(filelist))) {
  cat("Computing new variables and zonal means for", filelist$file[f], "\n")
  
  # compute new variables: 
  #  ROOTC, ROOTN, ROOTP (as fine + coarse)
  #  SOILC, SOILN, SOILP (summing all pools)
  #  TOTVEGC, TOTVEGN, TOTVEGP, WOODC (no change)
  #  TOTLITC, TOTLITN, TOTLITP (no change)
  #  NPP, SMINN_TO_PLANT, SMINP_TO_PLANT (no change)
  
  # CDO compute instructions
  tf <- tempfile()
  computations <- c("ROOTC = FROOTC + DEADCROOTC + LIVECROOTC",
                    "ROOTN = FROOTN + DEADCROOTN + LIVECROOTN",
                    "ROOTP = FROOTP + DEADCROOTP + LIVECROOTP",
                    "SOILC = SOIL1C + SOIL2C + SOIL3C + SOIL4C",
                    "SOILN = SOIL1N + SOIL2N + SOIL3N + SOIL4N",
                    "SOILP = SOIL1P + SOIL2P + SOIL3P + SOIL4P",
                    "SOIL1P = SOIL1P",
                    "SOIL2P = SOIL2P",
                    "SOIL3P = SOIL3P",
                    "LITR1C = LITR1C",
                    "LITR2C = LITR2C",
                    "LITR3C = LITR3C",
                    "LITR1N = LITR1N",
                    "LITR2N = LITR2N",
                    "LITR3N = LITR3N",
                    "LITR1P = LITR1P",
                    "LITR2P = LITR2P",
                    "LITR3P = LITR3P",
                    "TOTVEGC = TOTVEGC",
                    "TOTVEGN = TOTVEGN",
                    "TOTVEGP = TOTVEGP",
                    "WOODC = WOODC",
                    "TOTLITC = TOTLITC",
                    "TOTLITN = TOTLITN",
                    "TOTLITP = TOTLITP",
                    "NPP = NPP * 60 * 60 * 24 * 365",
                    "SMINN_TO_PLANT = SMINN_TO_PLANT * 60 * 60 * 24 * 365",
                    "SMINP_TO_PLANT = SMINP_TO_PLANT * 60 * 60 * 24 * 365")
  
  # LBNL outputs don't have a fourth soil pool - handle this
  outputs <- system2(WHICH_CDO, c("showname", filelist$file[f]), stdout = TRUE)
  if(!grepl("SOIL4", outputs)) {
    cat("Removing 4th soil pool from computations...\n")
    computations <- gsub("SOIL4[CNP]$", "0", computations)  # change to "0"
  } 
  
  # run cdo zonmean chained with variable computation on each file
  ofile <- sub("_ann.nc", "_ann_zonmean.nc", filelist$file[f], fixed = TRUE)
  unlink(ofile)
  
  result <- system2(WHICH_CDO, args = c("zonmean", 
                                        paste0("-expr,'", paste(computations, collapse = ";"), "'"), 
                                        filelist$file[f], ofile))
  if(file.exists(ofile)) cat("Wrote", ofile, "\n") 
  if(result) warning("Non-zero result code:", result, "\n")
  
  # Compute summed area of latitudinal bands
  # December 19, 2017: need to first correct for landfrac though
  cat("Computing area in each ...\n")
  ofile <- sub("_ann.nc", "_ann_area.nc", filelist$file[f], fixed = TRUE)
  unlink(ofile)
  
  result <- system2(WHICH_CDO, args = c("zonsum", 
                                        "-expr,'area=area*landfrac'", 
                                        filelist$file[f], ofile))
  if(file.exists(ofile)) cat("Wrote", ofile, "\n") 
  if(result) warning("Non-zero result code:", result, "\n")
  
}
