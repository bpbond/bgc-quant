# Run FERT BGC intercomparison analyses
# BBL December 19, 2017

library(dplyr)
library(readr)
library(tidyr)
library(raster)
library(ggplot2)
theme_set(theme_bw())
library(ncdf4)


FERT_FOLDER <- "input/fert/"

tibble(filename = list.files(FERT_FOLDER, pattern = "ann.nc$", full.names = TRUE)) %>% 
  separate(filename, into = c("junk1", "group", "run", "forcing", "LUC", "years", "junk2"), sep = "_", remove = FALSE) %>% 
  dplyr::select(-junk1, -junk2) %>% 
  separate(run, into = c("analysis", "run")) ->
  filelist

results <- list()
for(f in seq_len(nrow(filelist))) {
  cat(filelist$filename[f], "\n")
  nc <- nc_open(filelist$filename[f])
  time <- ncvar_get(nc, "time")
  
  # Two groups differ in their start year
  if(filelist$group[f] == "ORNL") {
    startyear = 1951
  } else {
    startyear = 1850
  }
  time = startyear + time / 360
  
  npp <- as_tibble(t(ncvar_get(nc, "NPP") * 60 * 60 * 24 * 365)) %>% 
    mutate(Year = time) %>% gather(site, npp, -Year) %>% 
    mutate(group = filelist$group[f], run = filelist$run[f])
  fpg <- as_tibble(t(ncvar_get(nc, "FPG"))) %>% gather(site, value)
  npp$fpg <- fpg$value
  fpg_p <- as_tibble(t(ncvar_get(nc, "FPG_P"))) %>% gather(site, value)
  npp$fpg_p <- fpg_p$value
  vegn <- as_tibble(t(ncvar_get(nc, "TOTVEGN"))) %>% gather(site, value)
  npp$totvegn <- vegn$value
  vegp <- as_tibble(t(ncvar_get(nc, "TOTVEGP"))) %>% gather(site, value)
  npp$totvegp <- vegp$value
  
  results[[f]] <- npp
  nc_close(nc)
}

bind_rows(results) %>%
  gather(outputvar, value, npp, fpg, fpg_p, totvegn, totvegp) %>% 
  spread(run, value) %>%
  mutate(diff = Experiment - Control, 
         diff_percent = (Experiment - Control) / Control,
         site = as.numeric(gsub("V", "", site))) -> 
  results

cat("Reading Qing's data...\n")
readxl::read_excel("input/fert/N&p fertilization experiment.xlsx") %>% 
  dplyr::select(site_index, PFT = `plant functional type`, lat, lon) %>% 
  right_join(results, by = c("site_index" = "site")) %>% 
  mutate(site = paste(PFT, site_index)) ->
  results

w <- 7
h <- 3

p <- ggplot(filter(results, outputvar == "npp"), aes(Year, diff, group = paste(site, group), color = group)) + 
  geom_line(alpha = 0.75) + 
  xlim(c(1950, 2010)) + coord_cartesian(ylim = c(-500, 500)) +
  facet_wrap(~PFT) + ylab("NPP diff (Experiment - Control) gC/m2/yr") +
  scale_color_discrete(labels = c("ECA", "CTC"))
print(p)
ggsave("output/fert-npp.png", width = w, height = h)
ggsave("output/fert-npp.pdf", width = w, height = h)

library(scales)
p_percent <- ggplot(filter(results, outputvar == "npp"), aes(Year, diff_percent, group = paste(site, group), color = group)) + 
  geom_line(alpha = 0.75) + 
  xlim(c(1950, 2010)) + scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap(~PFT) + ylab("NPP change (Experiment vs. Control) %") +
  scale_color_discrete(labels = c("ECA", "CTC"))
print(p_percent)
ggsave("output/fert-npp_percent.png", width = w, height = h)
ggsave("output/fert-npp_percent.pdf", width = w, height = h)


print(p %+% filter(results, outputvar == "fpg") + coord_cartesian(ylim = c(-0.2, 0.2)) + ylab("FPG diff (Experiment - Control)"))
ggsave("output/fert-fpg.png", width = w, height = h)

print(p %+% filter(results, outputvar == "fpg_p") + ylim(c(-0.1, 0.2)) + ylab("FPG_P diff (Experiment - Control)"))
ggsave("output/fert-fpg_p.png", width = w, height = h)

print(p %+% filter(results, outputvar == "totvegn") + coord_cartesian(ylim = c(-10, 30)) + ylab("TOTVEGN diff (Experiment - Control), gN/m2"))
ggsave("output/fert-totvegn.png", width = w, height = h)

print(p %+% filter(results, outputvar == "totvegp") + coord_cartesian(ylim = c(-2, 8)) + ylab("TOTVEGP diff (Experiment - Control), gP/m2"))
ggsave("output/fert-totvegp.png", width = w, height = h)

# Hawaii chronosequence - requested by Ying-Ping
p <- ggplot(filter(results, lon < -155, lat > 19), aes(Year, diff, group = paste(site, group), color = group)) + 
  geom_line() + 
  xlim(c(1950, 2010)) + ggtitle("Hawaii sites") +
  facet_wrap(~outputvar, scales = "free_y") + ylab("Difference (Experiment - Control)") +
  scale_color_discrete(labels = c("ECA", "CTC"))
print(p)
ggsave("output/fert-hawaii.png")
ggsave("output/fert-hawaii.pdf")
