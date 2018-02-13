# Prep the FLUXNET data
# Ben Bond-Lamberty December 2017
#
# Two tasks here: extract the annual GPP (along with QC, temp/precip,
# NEE, etc) from zipped FLUXNET TIER1 files; and extract the hourly
# NEE, filter for nighttime only, and compute an annual nighttime NEE.

# Downloaded 30 Jan 2017 from http://fluxnet.fluxdata.org (ftp.fluxdata.org/.fluxnet_downloads_86523/)
FLUXNET_DATA <- "~/Data/FLUXNET2015/"

library(dplyr)
library(tidyr)
library(readr)

# ==============================================================================
# Main 

cat("Welcome to 0-prep_fluxnet.R\n")

# Extract the *_SUBSET_MM_* and *_SUBSET_YY_* files for GPP_NT_VUT_REF, NEE_VUT_REF, and LE_F_MDS
td <- tempdir()
td <- "~/Desktop/test/"
files <- list.files(FLUXNET_DATA, pattern = "zip$", full.names = TRUE)
stopifnot(length(files) > 0)
d_annual <- list()
d_monthly <- list()
for(f in files) {
  cat("Unzipping", basename(f), "\n")
  zf <- utils::unzip(f, list = TRUE)
  annual_file <- utils::unzip(f, files = zf$Name[grep("SUBSET_YY", zf$Name)], exdir = td)
  monthly_file <- utils::unzip(f, files = zf$Name[grep("SUBSET_MM", zf$Name)], exdir = td)
  
  # Read in the extracted annual file 
  stopifnot(length(annual_file) == 1)
  cat("Reading", basename(annual_file), "\n")
  readr::read_csv(annual_file, na = "-9999") %>%
    dplyr::select(TIMESTAMP, TA_F, P_F, NEE_VUT_REF, GPP_NT_VUT_REF, LE_F_MDS) %>%
    mutate(filename = annual_file) ->
    d_annual[[f]]
  file.remove(annual_file)
  
  # Read in the extracted hourly files and compute annual nighttime Reco
    cat("Reading", basename(monthly_file))
    readr::read_csv(monthly_file, na = "-9999") %>%
      dplyr::select(TIMESTAMP, TA_F, P_F, NEE_VUT_REF, GPP_NT_VUT_REF, LE_F_MDS) %>%
      mutate(filename = monthly_file) ->
      d_monthly[[f]]
    file.remove(monthly_file)
}


# Combine with site data (in particular lon/lat information)
cat("--------------------------")
sitedata <- read_csv("input/fluxnet/fluxdata_sites.csv", col_types = "ccdddcdi")

cat("Combining flux data and merging with site data...")
bind_rows(d_annual) %>%
  separate(filename, into = c("FLX", "SITE_ID"), extra = "drop", sep = "_") %>%
  dplyr::select(-FLX) %>%
  left_join(sitedata, by = "SITE_ID") %>% 
  write_csv("fluxnet_data_annual.csv")

bind_rows(d_monthly) %>%
  separate(filename, into = c("FLX", "SITE_ID"), extra = "drop", sep = "_") %>%
  dplyr::select(-FLX) %>%
  left_join(sitedata, by = "SITE_ID") %>% 
  write_csv("fluxnet_data_monthly.csv")

cat("All done,", date())
