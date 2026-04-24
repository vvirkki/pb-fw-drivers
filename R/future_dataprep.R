# ------------------------------------------------------------------------------
#
# future_dataprep
#
# Supplementary future projection analysis script that calls for required
# libraries, initialises analysis functions, and prepares and calls configs
# for parallel data pre-processing.
#
# Analysis script for publication:
#
# "Regionally divergent drivers behind transgressions of the freshwater change
# planetary boundary"
#
# Vili Virkki, Lauren Seaby Andersen, Sofie te Wierik, Dieter Gerten,
# Miina Porkka
#
# Published in Nature Communications:
# ADD LINK
#
# Output data available in a Zenodo repository:
# https://doi.org/10.5281/zenodo.19663530
#
# Code availability:
# https://github.com/vvirkki/pb-fw-drivers
#
# Corresponding author & script author
# Vili Virkki (vili.virkki@uef.fi)
#
# Date:
# April 24th, 2026
#
# ------------------------------------------------------------------------------

# Load scripts -----------------------------------------------------------------

source("R/libraries.R")

# data preparation scripts
source("R/00_configure_run.R")
source("R/01_prepare_monthly_data.R")
source("R/03_set_local_baseline_range.R")
source("R/04_detect_local_deviations.R")
source("R/05_build_ice_areas.R")
source("R/06_get_area_aggregates.R")

# save sessionInfo()
sink(paste0("sessionInfo ", Sys.time(), ".txt"))
sessionInfo()
sf::sf_extSoftVersion()
sink()

# Prepare configs for parallel running -----------------------------------------

cmd_args <- commandArgs(trailingOnly = TRUE)
run_variable <- cmd_args[1]
ncores <- cmd_args[2]

# template config
cfg <- configure_run(paste0("configs/", run_variable, " main isimip3b"))

# isimip 3b scenarios
scenarios <- c("historical_histsoc_default",
               "ssp126_2015soc_default",
               "ssp370_2015soc_default",
               "ssp585_2015soc_default")

# isimip 3b climate forcings
gcms <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")

parallel_configs <- list()

for (s in 1:length(scenarios)) {

  if (run_variable == "dis") {

    ghms <- c("h08", "miroc-integ-land", "watergap2-2e")

  } else if (run_variable == "rootmoist") {

    ghms <-  c("miroc-integ-land")

  }

  for (i in 1:length(ghms)) {
    for (j in 1:length(gcms)) {

      ensemble_member <- paste0(run_variable, " ", ghms[i], " ", gcms[j])
      new_config <- cfg
      new_config$variable <- run_variable
      new_config$impactmodel <- ghms[i]
      new_config$forcing <- gcms[j]
      new_config$period_label <- c("picontrol_1850soc_default", scenarios[s])

      if (scenarios[s] == "historical_histsoc_default") {
        new_config$period_begin <- c(1850, 1850)
        new_config$period_end <- c(2014, 2014)
      } else if (grepl("ssp", scenarios[s])) {
        new_config$period_begin <- c(1850, 2015)
        new_config$period_end <- c(2014, 2100)
      } else {
        stop("isimip 3b scenario not given correctly")
      }

      # open log file
      if (!dir.exists("logs")) { dir.create("logs") }
      new_config$log_file <- paste0("logs/log ", Sys.time(),
                                    " ", ensemble_member, ".txt")
      write(format(Sys.time(), "%a %b %d %X %Y"), new_config$log_file, append = TRUE)
      write("opening log...", new_config$log_file, append = TRUE)

      parallel_configs[[length(parallel_configs) + 1]] <- new_config
      message(paste0("prepared config: ", ensemble_member, " ", scenarios[s]))

    }

  }
}

# avoid messing up with globalenv in parallel processing
unlink(cfg$log_file)
remove(cfg, ensemble_member, run_variable, gcms, ghms, scenarios, i, j, s, new_config)

# Prepare main analysis data with prescribed configs ---------------------------

registerDoParallel(cores = ncores)
foreach(i = 1:length(parallel_configs),
        .errorhandling = "pass",
        .verbose = TRUE) %dopar% {

          prepare_monthly_data(parallel_configs[[i]])
          set_local_baseline_range(parallel_configs[[i]])
          detect_local_deviations(parallel_configs[[i]])
          build_ice_areas(parallel_configs[[i]])
          parallel_configs[[i]]$areal_division <- ""
          get_area_aggregates(parallel_configs[[i]])

        }
stopImplicitCluster()
