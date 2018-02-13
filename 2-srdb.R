# Run SRDB BGC intercomparison analyses
# BBL November 15, 2017

library(dplyr)
library(readr)
library(tidyr)
library(raster)
library(ggplot2)
theme_set(theme_bw())
library(ncdf4)


get_srdb <- function(srdb_file, bgc_folder) {
  cat("Reading", srdb_file, "\n")
  read_csv(srdb_file, col_types = "dcicicccccdddddccddccccccccddcdddddcddcddddididdddddddddddcccccddddddddcddddddcdcddddddddddddddddddddddc") %>% 
    filter(Manipulation == "None", Meas_method %in% c("IRGA", "Gas chromatography"), !is.na(Latitude), !is.na(Longitude), !is.na(Study_midyear)) %>% 
    dplyr::select(Ecosystem_type, Study_midyear, YearsOfData, Latitude, Longitude,
                  Rs_annual, Rh_annual, Q10_0_10, Q10_5_15, Q10_10_20, Q10_0_20) ->
    srdb
  
  # SRDB longitudes are -180 to 180, E3SM appears to be 0 to 180
  sp <- SpatialPoints(cbind(srdb$Longitude + 180, srdb$Latitude))
  
  srdb_q10 <- filter(srdb, !is.na(Q10_5_15))
  sp_q10 <- SpatialPoints(cbind(srdb_q10$Longitude + 180, srdb_q10$Latitude))
  
  # Annual values - compare to Rs_annual and Rh_annual
  tibble(file = list.files(bgc_folder, pattern = "*.nc$", full.names = TRUE)) %>% 
    separate(file, sep = "_", remove = FALSE, into = c("analysis", "group", "run", "forcing", "LUC", "years", "time")) %>% 
    separate(years, into = c("minyear", "maxyear"), convert = TRUE, remove = FALSE) ->
    filelist
  
  results <- list()
  hr_data <- list()
  q10_data <- list()
  globalflux <- list()
  for(f in seq_len(nrow(filelist))) {
    cat("Reading", filelist$file[f], "\n")
    
    # Read rh data from file for SRDB lon/lat points
    brs <- brick(filelist$file[f], varname = "SOILC_HR")
    rhs <- raster::extract(brs, sp)
    brl <- brick(filelist$file[f], varname = "LITTERC_HR")
    rhl <- raster::extract(brl, sp)
    rh <- rhs + rhl  # total RH is soil + litter
    
    if(filelist$time[f] == "ann.nc") {
      # For annual files, we compare Rs_annual and Rh_annual
      rh <- apply(rh, 1, mean, na.rm = TRUE)   # compute mean over years
      
      srdb %>% 
        mutate(model_hr = rh * 60 * 60 * 24 * 365,   # gC/m2/s to gC/m2/yr
               group = filelist$group[f], 
               forcing = filelist$forcing[f],
               years = filelist$years[1],
               LUC = filelist$LUC[f]) ->
        results[[f]]
      
      # Save entire data to look at distribution versus Hashimoto
      nc <- nc_open(filelist$file[f])
      hrd <- ncvar_get(nc, "SOILC_HR") + ncvar_get(nc, "LITTERC_HR")   # took me hours to figure this out  :(
      srd <- ncvar_get(nc, "SR") + ncvar_get(nc, "LITTERC_HR")
      hr_data[[f]] <- tibble(group = filelist$group[f], 
                             forcing = filelist$forcing[f],
                             LUC = filelist$LUC[f],
                             hr = na.omit(as.vector(hrd)) * 60 * 60 * 24 * 365) #,
#                             sr = na.omit(as.vector(srd)) * 60 * 60 * 24 * 365,
#                             hr_to_rs = na.omit(as.vector(hrd / srd)))
      
      # Compute global flux (as soil+litter) and save that too
      area <- ncvar_get(nc, "area")
      landfrac <- ncvar_get(nc, "landfrac")
      time <- ncvar_get(nc, "time")
      # Convert HR flux from gC/m2/s to Pg/yr
      flux <- apply(hrd, 3, function(x) sum(x * 1000 * 1000 * area * landfrac, na.rm = TRUE) * 60 * 60 * 24 * 365 / 1e15 )
      globalflux[[f]] <- tibble(group = filelist$group[f], 
                                forcing = filelist$forcing[f],
                                LUC = filelist$LUC[f],
                                Year = as.integer(floor(time / 365 + 1850)),
                                global_hr = flux)
      
      nc_close(nc)
    } else {
      # For monthly files, we compute a model Q10 and compare to observations
      
      # Read rs data from file for SRDB lon/lat points
      br_rs <- brick(filelist$file[f], varname = "SR")
      rs <- raster::extract(br_rs, sp_q10)
      br_lit <- brick(filelist$file[f], varname = "LITTERC_HR")
      rs_lit <- raster::extract(br_lit, sp_q10)
      br_tmp <- brick(filelist$file[f], varname = "TSOI_10CM")
      tmp <- raster::extract(br_tmp, sp_q10)
      cat("...", nrow(srdb_q10), "points\n")
      
      q10_data[[f]] <- tibble(model_q10 = rep(NA_real_, nrow(srdb_q10)),
                              q10 = rep(NA_real_, nrow(srdb_q10)),
                              group = filelist$group[f], 
                              forcing = filelist$forcing[f],
                              years = filelist$years[1],
                              LUC = filelist$LUC[f])
      
      # We have extracted temperature and HR data; now compute implied Q10 over that range
      for(i in seq_len(dim(rs)[1])) {
        tmp_515 <- tmp[i,] >= 273.1+5 & tmp[i,] <= 273.1+15   # only doing Q10_5_15
        if(sum(tmp_515, na.rm = TRUE) > 3) {
          df <- data.frame(tmp = tmp[tmp_515], rs = rs[tmp_515] + rs_lit[tmp_515])
          m <- lm(rs ~ tmp, data = df)
          pred <- predict(m, newdata = data.frame(tmp = c(273.1+5, 273.1+15))) # predict for 5 & 15 C
          
          q10_data[[f]]$model_q10[i] <- pred[2] / pred[1]
          q10_data[[f]]$q10[i] <- srdb_q10$Q10_5_15[i]
        }
      }
    }
    
    
  }
  list(results = bind_rows(results), 
       hr_data = bind_rows(hr_data), 
       globalflux = bind_rows(globalflux),
       q10_data = bind_rows(q10_data))
}


x <- get_srdb(srdb_file = "input/srdb/srdb-data.csv", bgc_folder = "input/srdb/")
w = 7  # figure width
h = 4  # figure height

# 1. Observations versus model for HR
results <- x$results
p1 <- ggplot(results, aes(Rh_annual, model_hr, color = group)) + 
  geom_point(na.rm=TRUE) + geom_abline() + geom_smooth(method = "lm", na.rm=TRUE) + 
  xlim(c(0, 1500)) + coord_equal() + facet_grid(LUC ~ forcing) +
  scale_color_discrete(labels = c("ECA", "CTC")) + 
  xlab("Observed soil HR (gC/m2/yr)") + ylab(paste0("Modeled HR (", results$years[1], ", gC/m2/yr)"))
print(p1)
ggsave("output/srdb-1-hr.png", width = w, height = h)

cat("Fitting model...")
m <- lm(model_hr ~ Rh_annual * group + forcing * group, data = results)
print(summary(m))

# 2. Global distribution of HR versus Hashimoto data
cat("Global distribution of HR...\n")
hr_data <- x$hr_data

nc <- nc_open("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc")
hashimoto_hr <- ncvar_get(nc, "co2")
nc_close(nc)
nc <- nc_open("~/Data/Hashimoto/RS_yr_Hashimoto2015.nc")
hashimoto_sr <- ncvar_get(nc, "co2")
nc_close(nc)

hashimoto <- tibble(hr = na.omit(as.vector(hashimoto_hr)))

p2 <- ggplot(x$hr_data, aes(x = hr, color = group)) + geom_density(na.rm = TRUE) +
  scale_color_discrete(labels = c("ECA", "CTC")) +
  facet_grid(forcing ~ LUC) +
  geom_density(data = hashimoto, color = "black", na.rm = TRUE) + 
  xlim(c(10, 1000)) +
  ggtitle("Distribution of HR values (black = observations)") + xlab("Soil HR (gC/m2/yr)")
print(p2)
ggsave("output/srdb-2-hr.png", width = w, height = h)


# 3. Model SR Q10 versus observations
results <- x$q10_data
p3 <- ggplot(results, aes(q10, model_q10, color = group)) + 
  geom_point(na.rm=TRUE) + geom_abline() + #geom_smooth(method = "lm", na.rm=TRUE) + 
  facet_grid(LUC ~ forcing) +
  scale_color_discrete(labels = c("ECA", "CTC")) +
  xlab("Observed Q10") + ylab(paste0("Modeled Q10 (", results$years[1], ")"))
print(p3)
ggsave("output/srdb-3-q10.png", width = w, height = h)


# 4. Global flux
# Downloaded August 25, 2017 from http://cse.ffpri.affrc.go.jp/shojih/data/index.html
# Processed using `cdo fldmean` to produce annual values
nc <- nc_open("~/Data/Hashimoto/RH_yr_Hashimoto2015_global.nc")
hashimoto <- ncvar_get(nc, "co2", start = c(1, 1, 1, 1), count = c(-1, -1, 1, -1))
global_area <- 1.293606e+14 # total land area (m2) in Hashimoto file
hashimoto <- tibble(Year = 1901 + ncvar_get(nc, "time"), hashimoto_rh = hashimoto * global_area / 1e15)          
nc_close(nc)
x$globalflux %>% 
  left_join(hashimoto, by = "Year") %>% 
  filter(Year <= 2010) -> 
  globalflux
p4 <- ggplot(globalflux, aes(hashimoto_rh, global_hr, color = group)) + 
  geom_point() + facet_grid(LUC ~ forcing) + geom_abline() + geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("ECA", "CTC")) +
  xlab("Hashimoto (observed) annual HR (Pg C)") + ylab("Modeled annual HR (Pg C)")
print(p4)
ggsave("output/srdb-4-global.png", width = w, height = h)
p4time <- ggplot(globalflux, aes(Year, global_hr, color = group)) + 
  geom_line() + facet_grid(LUC ~ forcing) + geom_line(aes(y = hashimoto_rh), color = "black") +
  scale_color_discrete(labels = c("ECA", "CTC")) +
  ylab("Modeled annual HR (Pg C)")
print(p4time)
ggsave("output/srdb-4-global-time.png", width = w, height = h)
ggsave("output/srdb-4-global-time.pdf", width = w, height = h)

cat("All done.", date())
