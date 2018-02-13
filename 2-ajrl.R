# Run Atul/Ruby 18 December analyses
# BBL December 18, 2017

library(dplyr)
library(readr)
library(tidyr)
library(raster)
library(ggplot2)
theme_set(theme_bw())
library(ncdf4)

WHICH_CDO <- "/usr/local/bin/cdo"

AJRL_FOLDER <- "input/ajrl-ss/"

tibble(file = list.files(AJRL_FOLDER, pattern = "ann.nc$", full.names = TRUE)) %>% 
  separate(file, sep = "_", remove = FALSE, into = c("analysis", "group", "run", "forcing", "LUC", "years", "time")) %>% 
  separate(years, into = c("minyear", "maxyear"), convert = TRUE, remove = FALSE) ->
  filelist

globalflux <- list()
for(f in seq_len(nrow(filelist))) {
  cat("Reading", filelist$file[f], "\n")
  nc <- nc_open(filelist$file[f])
  
  gpp <- ncvar_get(nc, "GPP")
  nee <- ncvar_get(nc, "NEE")
  npp <- ncvar_get(nc, "NPP")
  land_uptake <- ncvar_get(nc, "LAND_UPTAKE")
  area <- ncvar_get(nc, "area")
  landfrac <- ncvar_get(nc, "landfrac")
  time <- ncvar_get(nc, "time")
  nc_close(nc)
  
  # Convert fluxes from gC/m2/s to PgC/yr
  func <- function(x) sum(x * 1000 * 1000 * area * landfrac, na.rm = TRUE) * 60 * 60 * 24 * 365 / 1e15
  gpp <- apply(gpp, 3, func)
  nee <- apply(nee, 3, func)
  npp <- apply(npp, 3, func)
  land_uptake <- apply(land_uptake, 3, func)
  globalflux[[f]] <- tibble(group = filelist$group[f], 
                            forcing = filelist$forcing[f],
                            LUC = filelist$LUC[f],
                            Year = as.integer(floor(time / 365 + 1850)),
                            gpp = gpp, nee = nee, npp = npp, land_uptake = land_uptake)
  
}
bind_rows(globalflux) %>% 
  gather(flux, value, gpp, nee, npp, land_uptake) ->
  globalflux

w = 7  # figure width
h = 4  # figure height

p1 <- ggplot(globalflux, aes(Year, value, color = group, linetype = flux)) +
  geom_line() + facet_grid(LUC ~ forcing) +
  scale_color_discrete(labels = c("ECA", "CTC")) + 
  xlim(c(1851, 1860)) +
  ylab("Flux (Pg C)")
print(p1)
ggsave("output/ajrl-2-fluxes.png", width = w, height = h)


# Atul request: 
# Gridded plots for  (i) Carbon allocation to wood, (ii) Carbon allocation to leaf,
# and (iii) the ratio of (i) and (ii).
# Calculate that ratio
filelist <- list.files("input/ajrl-ss/", pattern = "*ann.nc$", full.names = TRUE)
for(f in filelist) {
  ofile <- gsub("_ann.nc", "_ann_ratio.nc", f, fixed = TRUE)
  result <- system2(WHICH_CDO, c("expr,'wood_to_leaf_alloc=WOODC_ALLOC/LEAFC_ALLOC'", f, ofile))
  if(result) warning("Non-zero result code ", result)
}

# Atul request: NEE over the entire global run
tibble(file = list.files("input/ajrl-nee/", pattern = "*ann.nc$", full.names = TRUE)) %>% 
  separate(file, sep = "_", remove = FALSE, into = c("analysis", "group", "run", "forcing", "LUC", "years", "delete")) %>% 
  dplyr::select(-delete) %>% 
  separate(years, into = c("minyear", "maxyear"), convert = TRUE, remove = FALSE) ->
  filelist

# For each group, run, forcing, and LUC concatenate files and compute global mean NEE

tf <- "temp.nc"
results <- list()
for(f in seq_len(nrow(filelist))) {
  cat(f, filelist$file[f], "\n")
  
  result <- system2(WHICH_CDO, c("fldsum", filelist$file[f], tf))
  nc <- nc_open(tf)
  nee <- ncvar_get(nc, "NEE")
  time <- ncvar_get(nc, "time")
  nc_close(nc)
  unlink(tf)
  
  results[[f]] <- tibble(Year = time / 365 + 1850,
                         nee = nee * 1.28e8 * 60 * 60 * 24 * 365 / 1e15, # gC/m2/s -> PgC/yr
                         group = filelist$group[f],
                         forcing = filelist$forcing[f],
                         LUC = filelist$LUC[f])
}

nee_results <- bind_rows(results)

p_nee <- ggplot(nee_results, aes(Year, nee, color = group)) + geom_line() + facet_grid(LUC ~ forcing) +
  scale_color_discrete(labels = c("ECA", "CTC"))
print(p_nee)
ggsave("output/ajrl-2-nee.png", width = w, height = h)
ggsave("output/ajrl-2-nee.pdf", width = w, height = h)


# Ying-Ping: NEE with LUC - NEE without
nee_results %>% 
  spread(LUC, nee) %>% 
  mutate(nee_diff = LUC - NoLUC) ->
  nee_results

p_nee_luc <- ggplot(nee_results, aes(Year, nee_diff, color = group)) + geom_line() + 
  facet_grid( ~ forcing) +
  ylab("NEE w/ LUC - NEE w/o LUC (Pg C)") +
  scale_color_discrete(labels = c("ECA", "CTC"))
print(p_nee_luc)
ggsave("output/ajrl-2-nee-luc.png", width = w, height = h)
ggsave("output/ajrl-2-nee-luc.pdf", width = w, height = h)
