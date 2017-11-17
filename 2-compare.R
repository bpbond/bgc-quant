# Run SRDB and FLUXNET BGC intercomparison analyses
# BBL November 15, 2017

library(dplyr)
library(readr)
library(tidyr)
library(raster)
library(ggplot2)
theme_set(theme_bw())

# General purpose function to extract data from netcdf
extract_data <- function(rasterstack, varname, lon, lat, midyear, nyears, 
                         file_startyear, file_layers, months_per_layer = 1) {
  
  
  #  printlog("Starting extraction for varname:", varname)
  
  x <- rep(NA_real_, length(lon))
  
  # Find nearest neighbors for all lon/lat pairs
  for(i in seq_along(lon)) {
    sp <- SpatialPoints(cbind(lon[i], lat[i]))
    
    if(is.null(file_startyear)) {
      start_layer <- nlayers <- 1  # no time dimension
    } else {
      midyear_layer <- (midyear[i] - file_startyear + 1) * (12 / months_per_layer)
      start_layer <- ceiling(midyear_layer - nyears[i] / 2 * (12 / months_per_layer))
      nlayers <- nyears[i] * (12 / months_per_layer)
    }
    
    cat("Extracting", i, lon[i], lat[i], midyear[i], nyears[i], "- layers", start_layer, 
        "to", start_layer + nlayers - 1, "...\n")
    
    # Weirdly, raster::extract does not throw an error if we pass it a negative start layer
    # It just rolls merrily along, returning from the beginning of the file
    if(start_layer < 0 | start_layer + nlayers - 1 > file_layers) {
      #      printlog(varname, "layer out of bounds #", i)
      next
    }
    
    # Extract the information for this point, over as many layers as needed, then average
    rasterstack %>%
      raster::extract(sp, layer = start_layer, nl = nlayers, varname = varname) %>%
      mean(na.rm = TRUE) ->
      x[i]
  }
  
  x
}


get_srdb <- function(srdb_file, bgc_folder) {
  cat("Reading", srdb_file, "\n")
  read_csv(srdb_file, col_types = "dcicicccccdddddccddccccccccddcdddddcddcddddididdddddddddddcccccddddddddcddddddcdcddddddddddddddddddddddc") %>% 
    filter(Manipulation == "None", Meas_method %in% c("IRGA", "Gas chromatography"), !is.na(Latitude), !is.na(Longitude), !is.na(Study_midyear)) %>% 
    dplyr::select(Ecosystem_type, Study_midyear, YearsOfData, Latitude, Longitude,
                  Rs_annual, Rh_annual, Q10_0_10, Q10_5_15, Q10_10_20, Q10_0_20) ->
    srdb
  
#  srdb <- srdb[1:50,]
  
  # Annual values - compare to Rs_annual and Rh_annual
  tibble(file = list.files(bgc_folder, pattern = "*.nc", full.names = TRUE)) %>% 
    separate(file, sep = "_", remove = FALSE, into = c("analysis", "group", "run", "forcing", "LUC", "years", "time")) %>% 
    separate(years, into = c("minyear", "maxyear"), convert = TRUE) ->
    filelist
  
  results <- list()
  for(f in seq_len(nrow(filelist))) {
    cat("Reading", filelist$file[f], "\n")
    
    # For annual files, we compare Rs_annual and Rh_annual
    # SRDB longitudes are -180 to 180, E3SM appears to be 0 to 180
    sp <- SpatialPoints(cbind(srdb$Longitude + 180, srdb$Latitude))
    
    if(filelist$time[f] == "ann.nc") {
      br <- brick(filelist$file[f], varname = "SOILC_HR")
      rh <- raster::extract(br, sp)
      rh <- apply(rh, 1, mean, na.rm = TRUE)   # compute mean over years
      
      # rs <- extract_data(brick(filelist$file[f], varname = "SR"),  varname = "SR",
      #             lon = srdb$Longitude, lat = srdb$Latitude, 
      #             midyear = srdb$Study_midyear, nyears = srdb$YearsOfData, 
      #             file_startyear = filelist$minyear[f], file_layers = filelist$maxyear[f] - filelist$minyear[f] + 1, months_per_layer = 12)
      
      srdb %>% 
        mutate(model_hr = rh * 60 * 60 * 24 * 365,
               group = filelist$group[f], forcing = filelist$forcing[f], LUC = filelist$LUC[f]) ->  # convert to gC/m2/yr
        results[[f]]
    }
    
    # For monthly files, we compute a model Q10 and compare to observations
    
  }
  bind_rows(results)  
}

x <- get_srdb(srdb_file = "input/srdb/srdb-data.csv", bgc_folder = "input/srdb/")

# 1. Observations versus model for HR
p <- ggplot(x, aes(Rh_annual, model_hr, color = group)) + 
  geom_point() + geom_abline() + geom_smooth(method = "lm") + 
  xlim(c(0, 1500)) + 
  xlab("Observed HR (gC/m2/yr)") + ylab("Modeled HR (gC/m2/yr)")
print(p)
ggsave("srdb-1-hr.png")

m <- lm(model_hr ~ Rh_annual * group, data=x)
print(summary(m))

# 2. Global distribution of HR versus Hashimoto data?

# 3. Model SR Q10 versus observations


