# Extract variables and times needed for BGC intercomparison analyses
# This script depends on a file specifying location of all the teams' output files,
# and another specifying the years and variables required for each analysis
# BBL November 2017

library(dplyr)
library(tidyr)
library(readr)

cat("Welcome to 1-process_rundata.R", date(), "\n")

runs_list <- read_table2("runs_list.txt", comment = "#")
analyses_list <- read_table2("analyses_list.txt", comment = "#") %>%
  separate(Years_to_use, into = c("minyear", "maxyear"), remove = FALSE, convert = TRUE)

# Run loop: step through the teams' different runs
for(r in seq_len(nrow(runs_list))) {
  cat("=========================\n")
  cat("Team:", runs_list$Team[r], "Run:", runs_list$Run[r], "\n")
  
  # Analysis loop: step through the various analyses the panel is doing
  for(a in which(analyses_list$Run == runs_list$Run[r])) {
    cat("-------------------------\n")
    cat("Analysis:", analyses_list$Analysis[a], "Run:", analyses_list$Run[a], "\n")

    ofile <- paste0(paste(analyses_list$Analysis[a],
                          runs_list$Team[r],
                          runs_list$Run[r],
                          runs_list$Forcing[r],
                          runs_list$LUC[r],
                          analyses_list$Years_to_use[a],
                          "mon", sep = "_"),
                    ".nc")

	if(file.exists(ofile)) {
		cat("File exists--skipping\n")
		next
	}
     
    YEARMON_PATTERN <- "[0-9]{4}-[0-9]{2}.nc$"
    ifiles <- list.files(runs_list$Directory[r], pattern = YEARMON_PATTERN, full.names = TRUE)
    cat("Files:", length(ifiles), "\n")
    if(length(ifiles) == 0) {
    	stop("Looking in ", runs_list$Directory[r], " and not seeing any files!")
    }
    
    # Select just the files corresponding to years we want
    cat("Limiting to", analyses_list$Years_to_use[a], "\n")
    yrpos <- regexpr(ifiles, pattern = YEARMON_PATTERN)
    yrs <- as.integer(substr(ifiles, yrpos, yrpos+3))
    ifiles <- ifiles[yrs %in% analyses_list$minyear[a]:analyses_list$maxyear[a]]
    cat("Files:", length(ifiles), "\n")
    
    cat("Variables to extract:", analyses_list$Vars_to_use[a], "\n")
    
    # Use CDO (for speed) to extract variables we need, from files we need
    # First we make a monthly file, and then an annual file
    cat("Extracting and concatenating data...\n")
    WHICH_CDO <- "/usr/common/software/cdo/1.9.0/bin/cdo"
    result <- system2(WHICH_CDO, 
                      args = c(paste0("select,name=", analyses_list$Vars_to_use[a]),
                               ifiles,
                               ofile))
    if(file.exists(ofile)) cat("Wrote", ofile, "\n") 
    if(result) warning("Non-zero result code:", result, "\n")
    
    # Make an annually-averaged file while we're at it
    if(!result) {
		cat("Computing annual average file...\n")
		tf1 <- tempfile()
		system2(WHICH_CDO, args = c("yearmonmean", ofile, tf1))

		# It turns out that yearmonmean screws up the 'area' and 'landfrac' variables
   		# (they don't have a time dimensions, get multiplied by 28! wtf?)
   		# Copy them over from the first of the monthly files to a tempfile
   		cat("Copying area and landfrac into annual file...\n")
    	tf2 <- tempfile()
    	system2(WHICH_CDO, args = c("select,name=area,landfrac", ifiles[1], tf2))
    	
		# ...and then merge in
		ofile_ann <- sub("mon", "ann", ofile)
		result <- system2(WHICH_CDO, args = c("merge", tf1, tf2, ofile_ann))
		if(file.exists(ofile_ann)) cat("Wrote", ofile_ann, "\n") 
		if(result) warning("Non-zero result code:", result, "\n")
    }
    
  }
  
}

cat("All done.", date(), "\n")
