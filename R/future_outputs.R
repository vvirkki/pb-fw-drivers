# ------------------------------------------------------------------------------
#
# future_outputs
#
# Supplementary intermediate output script that creates future projection
# outputs that are further post-processed into figures used in the manuscript.
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

# data export & plotting scripts
source("R/00_configure_run.R")
source("R/output01_land_area_with_local_deviations.R")

# save sessionInfo()
sink(paste0("sessionInfo ", Sys.time(), ".txt"))
sessionInfo()
sf::sf_extSoftVersion()
sink()

.find_ensmems <- function(cfg) {

  dirs <- list.dirs(paste0("data/", cfg$variable, "/", cfg$period_label[2],
                           "/", cfg$variable, "_monthly"))

  imdirs <- dirs %>%
    lapply(str_split, "/") %>%
    lapply(unlist) %>%
    lapply(tail, 1) %>%
    unlist()

  imdirs <- dirs[which(imdirs %in% (str_replace_all(cfg$impactmodel, "-", "_")))]

  emdirs <- c()
  for (i in 1:length(imdirs)) {
    emdirs <- c(emdirs, list.dirs(imdirs[i])[-1])
  }

  ensmems_lst <- emdirs %>%
    lapply(str_split, "/") %>%
    lapply(unlist) %>%
    lapply(tail, 2) %>%
    lapply(str_replace_all, "_", "-")

  return (ensmems_lst)

}

# Create outputs ---------------------------------------------------------------

if (!dir.exists("output")) { dir.create("output") }

# ISIMIP 3b --------------------------------------------------------------------

run_variable <- commandArgs(trailingOnly = TRUE)

# template config
cfg <- configure_run(paste0("configs/", run_variable, " main isimip3b"))

scenarios <- c("historical_histsoc_default",
               "ssp126_2015soc_default",
               "ssp370_2015soc_default",
               "ssp585_2015soc_default")

areal_divisions <- c("")

if (run_variable == "dis") {

  cfg$impactmodel <- c("h08", "miroc-integ-land", "watergap2-2e")
  cfg$forcing <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr",
                   "mri-esm2-0", "ukesm1-0-ll")

} else if (run_variable == "rootmoist") {

  cfg$impactmodel <- c("miroc-integ-land")
  cfg$forcing <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr",
                   "mri-esm2-0", "ukesm1-0-ll")

}

# outputs global ---------------------------------------------------------------

global_output_out <- paste0("output/", cfg$variable, "_future_ensemble")
ensmem_output_out <- paste0("output/", cfg$variable, "_future_single_ensemble_members")

if (!dir.exists(global_output_out)) { dir.create(global_output_out) }
if (!dir.exists(ensmem_output_out)) { dir.create(ensmem_output_out) }

storecfg <- cfg # for individual ensemble member looping
for (s in 1:length(scenarios)) {

  scenario_output_global <- paste0(global_output_out, "/", scenarios[s], "/")
  scenario_output_ensmem <- paste0(ensmem_output_out, "/", scenarios[s], "/")

  if (!dir.exists(scenario_output_global)) { dir.create(scenario_output_global) }
  if (!dir.exists(scenario_output_ensmem)) { dir.create(scenario_output_ensmem) }

  cfg$period_label <- c("picontrol_1850soc_default", scenarios[s])
  cfg$ensemble_selection <- ""

  if (scenarios[s] == "historical_histsoc_default") {
    cfg$period_begin <- c(1850, 1850)
    cfg$period_end <- c(2014, 2014)
  } else if (grepl("ssp", scenarios[s])) {
    cfg$period_begin <- c(1850, 2015)
    cfg$period_end <- c(2014, 2100)
  } else {
    stop("isimip 3b scenario not given correctly")
  }

  # area aggregations for each scenario and areal division
  for (a in 1:length(areal_divisions)) {

    message(paste0(scenarios[s], " ", areal_divisions[a]))
    cfg$areal_division <- areal_divisions[a]
    cfg$draw_ts_var <- "landShareInClass"

    if (cfg$areal_division == "") {
      cfg$normalise_ts <- FALSE
      fldr_created <- output_land_area_with_local_deviations(cfg)
    } else {
      stop("areal divisions should not be done for future")
      cfg$normalise_ts <- TRUE
      fldr_created <- output_land_area_with_local_deviations(cfg)
    }
    file.rename(fldr_created,
                str_replace(fldr_created, "output/", scenario_output_global))

  }

  # outputs ensemble members ---------------------------------------------------
  # run ensemble members exhaustively; all that are available after data prep

  message(scenarios[s])
  ensmems <- .find_ensmems(cfg)
  for (em in 1:length(ensmems)) {

    em_out <- paste0(ensmem_output_out, "/", scenarios[s], "/",
                     paste(ensmems[[em]], collapse = "_"), "/")
    if (!dir.exists(em_out)) { dir.create(em_out) }

    # message(paste(ensmems[[em]][1], ensmems[[em]][2], collapse = " "))
    cfg$impactmodel <- ensmems[[em]][1]
    cfg$forcing <- ensmems[[em]][2]
    cfg$areal_division <- ""
    cfg$draw_ts_var <- "landShareInClass"
    cfg$normalise_ts <- FALSE
    cfg$ensemble_selection <- ""

    if (scenarios[s] == "historical_histsoc_default") {
      cfg$period_begin <- c(1850, 1850)
      cfg$period_end <- c(2014, 2014)
    } else if (grepl("ssp", scenarios[s])) {
      cfg$period_begin <- c(1850, 2015)
      cfg$period_end <- c(2014, 2100)
    } else {
      stop("isimip 3b scenario not given correctly")
    }

    fldr_created <- output_land_area_with_local_deviations(cfg)
    file.rename(fldr_created,
                str_replace(fldr_created, "output/", em_out))

    cfg <- storecfg
    cfg$period_label <- c("picontrol_1850soc_default", scenarios[s])

  }

}
