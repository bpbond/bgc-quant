# Run FLUXNET BGC intercomparison analyses
# BBL November 15, 2017

library(dplyr)
library(readr)
library(tidyr)
library(raster)
library(ggplot2)
theme_set(theme_bw())
library(ncdf4)


FLUXNET_FOLDER_ORNL <- "input/fluxnet/ornl/"

FLUXNET_FOLDER_LBNL <- "input/fluxnet/lbnl/"


# Get the FLUXNET data assembled by the 0-prep_fluxnet script 

fluxnet_monthly <- read_csv("fluxnet_data_monthly.csv", col_types = "cdddddccdddcdd") %>% 
  mutate(YEAR = as.integer(substr(TIMESTAMP, 1, 4)), 
         MONTH = as.integer(substr(TIMESTAMP, 5, 6)))
fluxnet_annual <- read_csv("fluxnet_data_annual.csv")

# Get the ORNL results. Each site is in a separate file

tibble(file = list.files(FLUXNET_FOLDER_ORNL, pattern = "*.nc$", full.names = TRUE)) %>% 
  separate(file, sep = "_", remove = FALSE, into = c("analysis", "group", "run", "forcing", "LUC", "years", "time")) %>% 
  separate(years, into = c("minyear", "maxyear"), convert = TRUE, remove = FALSE) %>% 
  separate(run, into = c("FLUXNET", "Site_ID"), sep = ",") %>% dplyr::select(-FLUXNET) ->
  filelist

results_annual_ornl <- list()
results_monthly_ornl <- list()
for(f in seq_len(nrow(filelist))) {
  cat(f, filelist$file[f], "\n")
  nc <- nc_open(filelist$file[f])  
  
  cat("Extracting netcdf data...\n")
  tibble(group = filelist$group[f],
         Site = filelist$Site_ID[f],
         model_gpp = ncvar_get(nc, "GPP") * 60 * 60 * 24, # to gC/m2/day
         model_nee = ncvar_get(nc, "NEE") * 60 * 60 * 24, # to gC/m2/day
         model_immob = ncvar_get(nc, "ACTUAL_IMMOB") * 60 * 60 * 24, # to gN/m2/day
         model_sminn_to_plant = ncvar_get(nc, "SMINN_TO_PLANT") * 60 * 60 * 24, # to gN/m2/day
         model_frootc = ncvar_get(nc, "FROOTC"),
         model_lh = ncvar_get(nc, "EFLX_LH_TOT")) ->
    x
  time <- ncvar_get(nc, "time")
  
  if(filelist$time[f] == "ann.nc") {
    x %>% 
      mutate(Year = as.integer(time / 365 + 1850),
             model_gpp = model_gpp * 365, # to gC/m2/yr
             model_nee = model_nee * 365) ->
      results_annual_ornl[[f]]
  } else if(filelist$time[f] == "mon.nc") {
    x %>% 
      mutate(Year = as.integer(floor(time / 365 + 1850)),
             Month = rep(1:12, length.out = nrow(x))) ->
      results_monthly_ornl[[f]]
    
  } else {
    stop("Unknown 'time' field!")
  }
  nc_close(nc)
}

# Get the LBNL results. Sites are all in a single file

lbnl_sitelist <- read_delim(file.path(FLUXNET_FOLDER_LBNL, "42_fluxnet_sites_info.txt"), skip=1, delim="\t") %>% 
  dplyr::select(Sitenum = `site_#`, Site = site_code)

# Annual data
nc <- nc_open(file.path(FLUXNET_FOLDER_LBNL, "FLUXNET_LBNL_FLUXNET,all_Local_NoLUC_2000-2015_ann.nc"))  
gpp <- ncvar_get(nc, "GPP")
nee <- ncvar_get(nc, "NEE")
lh <- ncvar_get(nc, "EFLX_LH_TOT")
immob <- ncvar_get(nc, "ACTUAL_IMMOB")
sminn_to_plant <- ncvar_get(nc, "SMINN_TO_PLANT")
frootc <- ncvar_get(nc, "FROOTC")
time <- ncvar_get(nc, "time")
# These data come in as 42 (sites) x 17 (years)
results_ann_lbnl <- list()
for(i in seq_len(dim(gpp)[1])) {
  results_ann_lbnl[[i]] <- tibble(group = "LBNL",
                                  Sitenum = i,
                                  Year = as.integer(time / 365 + 1850),
                                  model_gpp = gpp[i,] * 60 * 60 * 24 * 365,
                                  model_nee = nee[i,] * 60 * 60 * 24 * 365,
                                  model_lh = lh[i,],
                                  model_immob = immob[i,],
                                  model_sminn_to_plant = sminn_to_plant[i,],
                                  model_frootc = frootc[i,])
}
nc_close(nc)

# Monthly data
nc <- nc_open(file.path(FLUXNET_FOLDER_LBNL, "FLUXNET_LBNL_FLUXNET,all_Local_NoLUC_2000-2015_mon.nc"))  
gpp <- ncvar_get(nc, "GPP")
nee <- ncvar_get(nc, "NEE")
lh <- ncvar_get(nc, "EFLX_LH_TOT")
immob <- ncvar_get(nc, "ACTUAL_IMMOB")
sminn_to_plant <- ncvar_get(nc, "SMINN_TO_PLANT")
frootc <- ncvar_get(nc, "FROOTC")
time <- ncvar_get(nc, "time")
# These data come in as 42 (sites) x 192 (12 months * 16 years)
results_mon_lbnl <- list()
for(i in seq_len(dim(gpp)[1])) {
  results_mon_lbnl[[i]] <- tibble(group = "LBNL",
                                  Sitenum = i,
                                  Year = as.integer(floor(time / 365 + 1850)),
                                  Month = rep(1:12, length.out = length(time)),            
                                  model_gpp = gpp[i,] * 60 * 60 * 24,
                                  model_nee = nee[i,] * 60 * 60 * 24,
                                  model_lh = lh[i,],
                                  model_immob = immob[i,],
                                  model_sminn_to_plant = sminn_to_plant[i,],
                                  model_frootc = frootc[i,])
}
nc_close(nc)


# Combine results

bind_rows(results_ann_lbnl) %>% 
  left_join(lbnl_sitelist) %>%    # site number -> FLUXNET site code
  dplyr::select(-Sitenum) %>% 
  bind_rows(results_annual_ornl) %>% 
  left_join(fluxnet_annual, by = c("Site" = "SITE_ID", "Year" = "TIMESTAMP")) %>% 
  filter(model_gpp > 0, GPP_NT_VUT_REF > 0,
         # LBNL ran a crop site--remove
         IGBP != "CRO",
         # ORNL didn't run this site--remove
         Site != "DE-Gri") -> #, Site != "AU-Tum") ->
  results_annual

bind_rows(results_mon_lbnl) %>% 
  left_join(lbnl_sitelist) %>%    # site number -> FLUXNET site code
  dplyr::select(-Sitenum) %>% 
  bind_rows(results_monthly_ornl) %>% 
  left_join(fluxnet_monthly, by = c("Site" = "SITE_ID", "Year" = "YEAR", "Month" = "MONTH")) %>% 
  filter(model_gpp > 0, GPP_NT_VUT_REF > 0, 
         # LBNL ran a crop site--remove
         IGBP != "CRO",
         # ORNL didn't run this site--remove
         Site != "DE-Gri") -> #, Site != "AU-Tum") ->
  results_monthly


# Make graphs etc.

w = 7  # figure width
h = 3  # figure height


# ANNUAL RESULTS

p_gpp <- ggplot(results_annual, aes(GPP_NT_VUT_REF, model_gpp, color = group)) + 
  geom_point() + geom_abline() + geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("ECA", "CTC")) +
  xlab("FLUXNET GPP (GPP_NT_VUT_REF, gC/m2/yr)") + ylab("Model GPP (gC/m2/yr)")
print(p_gpp)
ggsave("output/fluxnet_gpp_ann.png", width = w, height = h)
print(summary(lm(model_gpp ~ GPP_NT_VUT_REF * group, data = results_annual)))

print(p_gpp + facet_wrap(~IGBP))
ggsave("output/fluxnet_gpp_ann_igbp.png")


p_nee <- ggplot(results_annual, aes(NEE_VUT_REF, model_nee, color = group)) + 
  geom_point() + geom_abline() + geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("ECA", "CTC")) +
  xlab("FLUXNET NEE (NEE_VUT_REF, gC/m2/yr)") + ylab("Model NEE (gC/m2/yr)")
print(p_nee)
ggsave("output/fluxnet_nee_ann.png", width = w, height = h)
print(summary(lm(model_nee ~ NEE_VUT_REF * group, data = results_annual)))

print(p_nee + facet_wrap(~IGBP))
ggsave("output/fluxnet_nee_ann_igbp.png")


p_lh <- ggplot(results_annual, aes(LE_F_MDS, model_lh, color = group)) + 
  geom_point() + geom_abline() + geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("ECA", "CTC")) +
  xlab("FLUXNET LE (LE_F_MDS, W/m2)") + ylab("Model LH (W/2)")
print(p_lh)
ggsave("output/fluxnet_lh_ann.png", width = w, height = h)
print(summary(lm(model_lh ~ LE_F_MDS * group, data = results_annual)))

print(p_lh + facet_wrap(~IGBP))
ggsave("output/fluxnet_lh_ann_igbp.png")

# Annual time plots
# Make a template plot for repeated use below
results_annual$y_obs <- results_annual$GPP_NT_VUT_REF
results_annual$y_model <- results_annual$model_gpp
p_ann_time <- ggplot(results_annual, aes(Year, y_model, color = group)) + 
  geom_line() + geom_line(aes(y = y_obs), color = "black") + facet_wrap(~Site) +
  scale_color_discrete(labels = c("ECA", "CTC"))

print(p_ann_time + ylab("GPP (gC/m2/yr)"))
ggsave("output/fluxnet_gpp_ann_time.png")

results_annual$y_obs <- results_annual$NEE_VUT_REF
results_annual$y_model <- results_annual$model_nee
print(p_ann_time %+% results_annual + ylab("NEE (gC/m2/yr)"))
ggsave("output/fluxnet_nee_ann_time.png")

results_annual$y_obs <- results_annual$LE_F_MDS
results_annual$y_model <- results_annual$model_lh
print(p_ann_time %+% results_annual + ylab("LH (W/m2)"))
ggsave("output/fluxnet_lh_ann_time.png")


# MONTHLY RESULTS

# Make a template plot for repeated use below
results_monthly$x <- results_monthly$GPP_NT_VUT_REF
results_monthly$y <- results_monthly$model_gpp
p_mon <- ggplot(results_monthly, aes(x, y, color = group)) + 
  geom_point(alpha = 0.25) + geom_abline() + geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("ECA", "CTC"))
print(p_mon + labs(x = "FLUXNET monthly GPP (GPP_NT_VUT_REF, gC/m2/day)", y = "Model monthly GPP (gC/m2/day)"))
ggsave("output/fluxnet_gpp_mon.png")
print(last_plot() + facet_wrap(~IGBP))
ggsave("output/fluxnet_gpp_mon_igbp.png")
print(last_plot() + facet_wrap(~Month))
ggsave("output/fluxnet_gpp_mon_month.png")
print(last_plot() + facet_wrap(~Site))
ggsave("output/fluxnet_gpp_mon_site.png")

results_monthly$x <- results_monthly$NEE_VUT_REF
results_monthly$y <- results_monthly$model_nee
print(p_mon %+% results_monthly + labs(x = "FLUXNET NEE (GPP_NT_VUT_REF, gC/m2/day)", y = "Model NEE (gC/m2/day)"))
ggsave("output/fluxnet_nee_mon.png")
print(last_plot() + facet_wrap(~IGBP))
ggsave("output/fluxnet_nee_igbp.png")
print(last_plot() + facet_wrap(~Month))
ggsave("output/fluxnet_nee_mon_month.png")
print(last_plot() + facet_wrap(~Site))
ggsave("output/fluxnet_nee_mon_site.png")

results_monthly$x <- results_monthly$LE_F_MDS
results_monthly$y <- results_monthly$model_lh
print(p_mon %+% results_monthly + labs(x = "FLUXNET LH (LE_F_MDS, W/m2)", y = "Model LH (W/m2)"))
ggsave("output/fluxnet_lh_mon.png")
print(last_plot() + facet_wrap(~IGBP))
ggsave("output/fluxnet_lh_igbp.png")
print(last_plot() + facet_wrap(~Month))
ggsave("output/fluxnet_lh_mon_month.png")
print(last_plot() + facet_wrap(~Site))
ggsave("output/fluxnet_lh_mon_site.png")


# Monthly time plots
# Make a template plot for repeated use below
results_monthly$y_obs <- results_monthly$GPP_NT_VUT_REF
results_monthly$y_model <- results_monthly$model_gpp
p_mon_time <- ggplot(results_monthly, aes(Month, y_model, color = group, group = paste(Year, group))) + 
  geom_line(alpha = 0.5) + geom_line(aes(y = y_obs), color = "black", alpha = 0.5) + facet_wrap(~Site) +
  scale_color_discrete(labels = c("ECA", "CTC"))

print(p_mon_time + ylab("GPP (gC/m2/day)"))
ggsave("output/fluxnet_gpp_mon_time.png")

results_monthly$y_obs <- results_monthly$NEE_VUT_REF
results_monthly$y_model <- results_monthly$model_nee
print(p_mon_time %+% results_monthly + ylab("NEE (gC/m2/day)"))
ggsave("output/fluxnet_nee_mon_time.png")

results_monthly$y_obs <- results_monthly$LE_F_MDS
results_monthly$y_model <- results_monthly$model_lh
print(p_mon_time %+% results_monthly + ylab("LH (W/m2)"))
ggsave("output/fluxnet_lh_mon_time.png")

# Compute R2 and RSE values
myfunc <- function(x) {
  m1 <- lm(model_gpp ~ GPP_NT_VUT_REF, data = x) %>% summary
  m2 <- lm(model_nee ~ NEE_VUT_REF, data = x) %>% summary
  m3 <- lm(model_lh ~ LE_F_MDS, data = x) %>% summary
  tibble(r2_gpp = m1$r.squared,
         r2_nee = m2$r.squared,
         r2_lh = m3$r.squared,
         rse_gpp = m1$sigma,
         rse_nee = m2$sigma,
         rse_lh = m3$sigma)
}

results_annual %>% 
  group_by(group) %>% 
  do(myfunc(.)) %>% 
  mutate(time = "annual") ->
  ra_summary

results_monthly %>% 
  group_by(group) %>% 
  do(myfunc(.)) %>% 
  mutate(time = "monthly") %>% 
  bind_rows(ra_summary) %>% 
  gather(label, value, -group, -time) %>% 
  separate(label, into = c("metric", "flux")) ->
  results_summary

results_summary %>% 
  ungroup %>% 
  mutate(timegroup = paste(time, group)) %>% 
  dplyr::select(-time, -group) %>% 
  spread(timegroup, value) %>% 
  readr::write_csv("fluxnet_metrics.csv")

p_metrics <- ggplot(results_summary, aes(flux, value, fill = group)) +
  geom_col(position="dodge") + facet_grid(metric ~ time, scales = "free_y") +
  scale_fill_discrete(labels = c("ECA", "CTC"))
print(p_metrics)
ggsave("output/fluxnet_metrics.png", width = w, height = h)


# Atul's idea - email December 19, 2017
# "I am really keen to understand the difference between  nutrient
# relative demand (RD) approach (i.e., ORNL  approach) and trait-based
# Equilibrium Chemistry Approximation (ECA) approach.  I agree with
# Bill’s argument, but then Peter is saying that Zhu et al. paper  is
# not representing  the relative demand (RD) approach correctly.  Zhu et
# al paper presents the results for 15N Uptake by Microbes / 15N Uptake
# by plant as function of root biomass density (kg/m2), but we don’t
# have two model results for 15N. However, I have come-up with  an idea
# to evaluate the two model results. Instead of 15N, we can use a ratio
# of microbial immobilization flux and plant N uptake flux for
# grasslands and plot this ratio as a function of root density. If
# Bill’s interpretation of RD approach is correct, then the two model
# results should show same or similar pattern as shown in figure 2 of
# Zhu et al. paper. We can test this for steady state conditions for
# two models."
p_aj <- ggplot(filter(results_annual, IGBP == "GRA"), aes(model_frootc, model_immob / model_sminn_to_plant, color = group)) + 
  geom_point() + facet_wrap(~IGBP) + 
  xlab("Fine roots (gC/m2)") + ylab("Ratio of microbial N immobilization to plant N uptake") +
  scale_color_discrete("Model", labels = c("ECA", "CTC"))

zhu2 <- png::readPNG("input/Zhu Figure 2.png")

print(p_aj + annotation_raster(zhu2, xmin = 125, xmax = 385, ymin = 2, ymax = 4.75))
ggsave("output/fluxnet_rd.png")
ggsave("output/fluxnet_rd.pdf")
